// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/ParticleId.hpp"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <ostream>
#include <ranges>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/Print.h>
#include <boost/container/flat_map.hpp>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsExamples {

TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  m_outputJets.initialize(m_cfg.outputJets);

  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);

  if (m_cfg.doJetLabeling && !m_cfg.inputHepMC3Event.has_value()) {
    throw std::invalid_argument("Input HepMC3 event is not configured ");
  }

  m_inputHepMC3Event.initialize(m_cfg.inputHepMC3Event.value());
}

namespace {

JetLabel jetLabelFromHadronType(ParticleId::HadronType hType) {
  using enum ActsExamples::ParticleId::HadronType;
  switch (hType) {
    case BBbarMeson:
    case BottomMeson:
    case BottomBaryon:
      return JetLabel::BJet;
    case CCbarMeson:
    case CharmedMeson:
    case CharmedBaryon:
      return JetLabel::CJet;
    case StrangeMeson:
    case StrangeBaryon:
    case LightMeson:
    case LightBaryon:
      return JetLabel::LightJet;
    default:
      return JetLabel::Unknown;
  }
}
}  // namespace

ProcessCode ActsExamples::TruthJetAlgorithm::initialize() {
  std::ofstream outfile;
  outfile.open("particles.csv");
  outfile << "event,pt,eta,phi,pdg" << std::endl;

  outfile.flush();
  outfile.close();

  outfile.open("jets.csv");
  outfile << "event,pt,eta,phi,label" << std::endl;

  outfile.flush();
  outfile.close();

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  TrackJetContainer outputJets;

  const SimParticleContainer& truthParticlesRaw = m_inputTruthParticles(ctx);
  std::vector<const SimParticle*> truthParticles;
  truthParticles.reserve(truthParticlesRaw.size());
  std::ranges::transform(truthParticlesRaw, std::back_inserter(truthParticles),
                         [](const auto& particle) { return &particle; });

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition =
      fastjet::JetDefinition(fastjet::antikt_algorithm, m_cfg.jetClusteringR);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  static std::mutex mtxPseudoJets;
  {
    std::lock_guard lock(mtxPseudoJets);
    std::ofstream outfile;
    outfile.open("particles.csv",
                 std::ios_base::app);  // append instead of overwrite

    for (unsigned int i = 0; i < truthParticles.size(); i++) {
      const auto* particle = truthParticles.at(i);
      fastjet::PseudoJet pseudoJet(
          particle->momentum().x(), particle->momentum().y(),
          particle->momentum().z(), particle->energy());

      outfile << ctx.eventNumber << "," << pseudoJet.pt() << ","
              << pseudoJet.eta() << "," << pseudoJet.phi() << ","
              << static_cast<int>(particle->pdg());
      outfile << std::endl;

      pseudoJet.set_user_index(i);
      inputPseudoJets.push_back(pseudoJet);
    }

    outfile.flush();
    outfile.close();
  }

  ACTS_DEBUG("Number of input jet input particles: " << inputPseudoJets.size());

  std::vector<fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clusterSeq;
  {
    Acts::ScopedTimer timer("Jet clustering", logger(), Acts::Logging::DEBUG);
    // Run the jet clustering, only once
    clusterSeq =
        fastjet::ClusterSequence(inputPseudoJets, defaultJetDefinition);

    // Get the jets above a certain pt threshold
    jets = sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
    ACTS_DEBUG("Number of clustered jets: " << jets.size());
  }

  std::vector<std::pair<JetLabel, std::shared_ptr<const HepMC3::GenParticle>>>
      hadrons;

  if (m_cfg.doJetLabeling) {
    ACTS_DEBUG("Jet labeling is enabled");
    const auto& genEvent = *m_inputHepMC3Event(ctx);

    auto hadronView =
        genEvent.particles() | std::views::filter([](const auto& particle) {
          return ParticleId::isHadron(particle->pdg_id());
        }) |
        std::views::transform([](const auto& particle) {
          auto type = ActsExamples::ParticleId::hadronType(particle->pdg_id());
          auto label = jetLabelFromHadronType(type);
          return std::pair{label, particle};
        }) |
        std::views::filter([](const auto& hadron) {
          return hadron.first > JetLabel::Unknown;
        });

    std::ranges::copy(hadronView, std::back_inserter(hadrons));

    // deduplicate hadrons
    std::ranges::sort(hadrons, [](const auto& a, const auto& b) {
      return a.second->id() < b.second->id();
    });
    auto unique = std::ranges::unique(hadrons);
    hadrons.erase(unique.begin(), unique.end());
  }

  // Prepare jets for the storage - conversion of jets to custom track jet
  // class (and later add here the jet classification)

  constexpr static auto deltaR = [](const fastjet::PseudoJet& a,
                                    const fastjet::PseudoJet& b) {
    double dphi = abs(a.phi() - b.phi());
    if (dphi > std::numbers::pi) {
      dphi = std::numbers::pi * 2 - dphi;
    }
    double drap = a.rap() - b.rap();
    return std::sqrt(dphi * dphi + drap * drap);
  };

  auto classifyJet = [&](const fastjet::PseudoJet& jet) {
    auto hadronsInJetView =
        hadrons | std::views::filter([&jet, this](const auto& hadron) {
          const auto& momentum = hadron.second->momentum();
          fastjet::PseudoJet hadronJet(momentum.px(), momentum.py(),
                                       momentum.pz(), momentum.e());
          return deltaR(jet, hadronJet) < m_cfg.jetLabelingDeltaR;
        }) |
        std::views::transform([](const auto& hadron) {
          return std::pair{
              hadron.second,
              jetLabelFromHadronType(ActsExamples::ParticleId::hadronType(
                  hadron.second->pdg_id()))};
        });

    std::vector<std::pair<std::shared_ptr<const HepMC3::GenParticle>, JetLabel>>
        hadronsInJet;
    std::ranges::copy(hadronsInJetView, std::back_inserter(hadronsInJet));

    ACTS_VERBOSE("-> hadrons in jet: " << hadronsInJet.size());
    for (const auto& hadron : hadronsInJet) {
      ACTS_VERBOSE(
          "  - " << hadron.first->pdg_id() << " "
                 << Acts::findName(hadron.first->pdg_id()).value_or("UNKNOWN")
                 << " label=" << hadron.second);
    }

    auto maxHadronIt = std::ranges::max_element(
        hadronsInJet, [](const auto& a, const auto& b) { return a < b; },
        [](const auto& a) {
          const auto& [hadron, label] = a;
          return label;
        });

    if (maxHadronIt == hadronsInJet.end()) {
      // Now hadronic "jet"
      return JetLabel::Unknown;
    }

    const auto& [maxHadron, maxHadronLabel] = *maxHadronIt;

    ACTS_VERBOSE("-> max hadron type="
                 << Acts::findName(maxHadron->pdg_id()).value_or("UNKNOWN")
                 << " label=" << maxHadronLabel);

    return maxHadronLabel;
  };

  boost::container::flat_map<JetLabel, std::size_t> jetLabelCounts;

  Acts::AveragingScopedTimer timer("Jet classification", logger(),
                                   Acts::Logging::DEBUG);

  static std::mutex mtxJets;
  {
    std::lock_guard lock(mtxJets);
    std::ofstream outfile;
    outfile.open("jets.csv", std::ios_base::app);

    for (unsigned int i = 0; i < jets.size(); i++) {
      // Get information on the jet constituents
      const auto& jet = jets[i];
      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
      std::vector<int> constituentIndices;
      constituentIndices.reserve(jetConstituents.size());

      // Get the jet classification label later here! For now, we use "unknown"

      Acts::Vector4 jetFourMomentum(jets[i].px(), jets[i].py(), jets[i].pz(),
                                    jets[i].e());
      ACTS_VERBOSE("Found jet "
                   << i << " with 4-momentum: " << jetFourMomentum(0) << ", "
                   << jetFourMomentum(1) << ", " << jetFourMomentum(2) << ", "
                   << jetFourMomentum(3) << " and " << constituentIndices.size()
                   << " constituents.");

      JetLabel label = JetLabel::Unknown;
      if (m_cfg.doJetLabeling) {
        timer.sample();
        label = classifyJet(jet);
      }

      outfile << ctx.eventNumber << "," << jet.pt() << "," << jet.eta() << ","
              << jet.phi() << "," << static_cast<int>(label);
      outfile << std::endl;

      // Initialize the (track) jet with 4-momentum and jet label
      ActsExamples::TrackJet storedJet(jetFourMomentum, label);

      // Add the jet constituents to the (track)jet
      for (unsigned int j = 0; j < jetConstituents.size(); j++) {
        // Get the index of the constituent in the original input pseudo jets
        constituentIndices.push_back(jetConstituents[j].user_index());
      }

      storedJet.setConstituents(constituentIndices);

      outputJets.push_back(storedJet);

      jetLabelCounts[label] += 1;

      ACTS_VERBOSE("-> jet label: " << label);
      ACTS_VERBOSE("-> jet constituents: ");

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        for (const auto& constituent : constituentIndices) {
          const auto& particle = truthParticles.at(constituent);
          ACTS_VERBOSE("- " << particle);
        }
      }
    }

    outfile.flush();
    outfile.close();
  }

  ACTS_DEBUG("-> jet label counts: ");
  for (const auto& [label, count] : jetLabelCounts) {
    ACTS_DEBUG("  - " << label << ": " << count);
  }

  m_outputJets(ctx, std::move(outputJets));
  return ProcessCode::SUCCESS;
}
};  // namespace ActsExamples
