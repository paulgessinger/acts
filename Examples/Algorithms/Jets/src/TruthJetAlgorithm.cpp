// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

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

using namespace Acts::UnitLiterals;

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

JetLabel jetLabelFromHadronType(Acts::HadronType hType) {
  using enum Acts::HadronType;
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
  if (m_cfg.debugCsvOutput) {
    std::ofstream outfile;
    outfile.open("particles.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("jets.csv");
    outfile << "event,pt,eta,phi,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("hadrons.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();
  }

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  TrackJetContainer outputJets;
  Acts::ScopedTimer globalTimer("TruthJetAlgorithm", logger(),
                                Acts::Logging::DEBUG);

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
    std::ofstream outfile;
    Acts::ScopedTimer timer("Input particle building", logger(),
                            Acts::Logging::DEBUG);

    std::unique_lock lock(mtxPseudoJets, std::defer_lock);
    if (m_cfg.debugCsvOutput) {
      lock.lock();
      outfile.open("particles.csv",
                   std::ios_base::app);  // append instead of overwrite
    }

    for (unsigned int i = 0; i < truthParticles.size(); i++) {
      const auto* particle = truthParticles.at(i);

      const auto* gp = particle->genParticle();
      // Convention is that idx=0 is hard-scatter, check if we need to skip it
      if (m_cfg.clusterHardScatterParticlesOnly && gp != nullptr &&
          HepMC3Util::eventGeneratorIndex(*gp) != 0) {
        continue;
      }

      fastjet::PseudoJet pseudoJet(
          particle->momentum().x(), particle->momentum().y(),
          particle->momentum().z(), particle->energy());

      if (m_cfg.debugCsvOutput) {
        outfile << ctx.eventNumber << "," << pseudoJet.pt() << ","
                << pseudoJet.eta() << "," << pseudoJet.phi() << ","
                << static_cast<int>(particle->pdg()) << ","
                << static_cast<int>(jetLabelFromHadronType(
                       Acts::ParticleId::hadronType(particle->pdg())));
        outfile << std::endl;
      }

      pseudoJet.set_user_index(i);
      inputPseudoJets.push_back(pseudoJet);
    }

    if (m_cfg.debugCsvOutput) {
      outfile.flush();
      outfile.close();
    }
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

    if (m_cfg.jetEtaRange.first.has_value() ||
        m_cfg.jetEtaRange.second.has_value()) {
      double minEta = m_cfg.jetEtaRange.first.value_or(
          std::numeric_limits<double>::lowest());
      double maxEta =
          m_cfg.jetEtaRange.second.value_or(std::numeric_limits<double>::max());
      std::erase_if(jets, [minEta, maxEta](const auto& jet) {
        return jet.eta() < minEta || jet.eta() > maxEta;
      });
    }

    ACTS_DEBUG("Number of clustered jets: " << jets.size());
  }

  std::vector<std::pair<JetLabel, std::shared_ptr<const HepMC3::GenParticle>>>
      hadrons;

  if (m_cfg.doJetLabeling) {
    Acts::ScopedTimer timer("Hadron finding", logger(), Acts::Logging::DEBUG);
    ACTS_DEBUG("Jet labeling is enabled");
    const auto& genEvent = *m_inputHepMC3Event(ctx);

    auto hadronView =
        genEvent.particles() | std::views::filter([this](const auto& particle) {
          if (m_cfg.jetLabelingHardScatterHadronsOnly) {
            if (HepMC3Util::eventGeneratorIndex(*particle) != 0) {
              // Convention is that idx=0 is hard-scatter
              return false;
            }
          }

          if (!Acts::ParticleId::isHadron(particle->pdg_id())) {
            return false;
          }

          if (particle->status() != HepMC3Util::kDecayedParticleStatus &&
              particle->status() != HepMC3Util::kUndecayedParticleStatus) {
            return false;
          }

          // Apply pt cut only on B or C hadrons
          auto label = jetLabelFromHadronType(
              Acts::ParticleId::hadronType(particle->pdg_id()));
          using enum JetLabel;

          if (label == BJet || label == CJet) {
            if (particle->momentum().pt() < m_cfg.jetLabelingHadronPtMin) {
              return false;
            }
          }

          return true;
        }) |
        std::views::transform([](const auto& particle) {
          auto type = Acts::ParticleId::hadronType(particle->pdg_id());
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

    if (m_cfg.debugCsvOutput) {
      static std::mutex mtxHadrons;
      std::lock_guard lock(mtxHadrons);
      std::ofstream outfile;
      outfile.open("hadrons.csv", std::ios_base::app);
      for (const auto& hadron : hadrons) {
        outfile << ctx.eventNumber << "," << hadron.second->momentum().pt()
                << "," << hadron.second->momentum().eta() << ","
                << hadron.second->momentum().phi() << ","
                << static_cast<int>(hadron.second->pdg_id()) << ","
                << static_cast<int>(hadron.first);
        outfile << std::endl;
      }
    }
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
          return std::pair{hadron.second,
                           jetLabelFromHadronType(Acts::ParticleId::hadronType(
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

  static std::mutex mtxJets;
  {
    Acts::AveragingScopedTimer timer("Jet classification", logger(),
                                     Acts::Logging::DEBUG);

    std::ofstream outfile;
    std::unique_lock lock(mtxJets, std::defer_lock);
    if (m_cfg.debugCsvOutput) {
      lock.lock();
      outfile.open("jets.csv", std::ios_base::app);
    }

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
        ACTS_DEBUG("Classifying jet " << i);
        auto sample = timer.sample();
        label = classifyJet(jet);
      }

      if (m_cfg.debugCsvOutput) {
        outfile << ctx.eventNumber << "," << jet.pt() << "," << jet.eta() << ","
                << jet.phi() << "," << static_cast<int>(label);
        outfile << std::endl;
      }

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

    if (m_cfg.debugCsvOutput) {
      outfile.flush();
      outfile.close();
    }
  }

  m_numJets += outputJets.size();
  if (m_cfg.doOverlapRemoval) {
    overlapRemoval(truthParticlesRaw, outputJets);
  }

  ACTS_DEBUG("-> jet label counts: ");
  for (const auto& [label, count] : jetLabelCounts) {
    ACTS_DEBUG("  - " << label << ": " << count);
  }

  m_numLightJets += jetLabelCounts[JetLabel::LightJet];
  m_numCJets += jetLabelCounts[JetLabel::CJet];
  m_numBJets += jetLabelCounts[JetLabel::BJet];
  m_numJetsAfterOverlapRemoval += outputJets.size();

  m_outputJets(ctx, std::move(outputJets));
  return ProcessCode::SUCCESS;
}

void TruthJetAlgorithm::overlapRemoval(
    const SimParticleContainer& truthParticles, TrackJetContainer& jets) const {
  ACTS_DEBUG("Running overlap removal for jets against isolated truth leptons");
  Acts::ScopedTimer timer("Overlap removal", logger(), Acts::Logging::DEBUG);

  std::vector<const SimParticle*> isolatedLeptons;
  // Empirical size estimate
  isolatedLeptons.reserve(
      static_cast<std::size_t>(5 * std::log(truthParticles.size())));

  std::vector<double> phis;
  std::vector<double> etas;

  phis.reserve(truthParticles.size());
  etas.reserve(truthParticles.size());

  {
    Acts::ScopedTimer t{"Calculating phis and etas", logger(),
                        Acts::Logging::DEBUG};

    for (const auto& particle : truthParticles) {
      phis.push_back(Acts::VectorHelpers::phi(particle.direction()));
      etas.push_back(Acts::VectorHelpers::eta(particle.direction()));
    }
  }

  for (const auto& [i, particle] : Acts::enumerate(truthParticles)) {
    bool accept = Acts::ParticleId::isMuon(particle.pdg()) ||
                  Acts::ParticleId::isElectron(particle.pdg()) ||
                  Acts::ParticleId::isTau(particle.pdg()) ||
                  Acts::ParticleId::isPhoton(particle.pdg());

    if (!accept) {
      continue;
    }

    // For this lepton, sum up all total momenta inside a cone
    double totalMomentum = 0.;
    for (const auto& [j, otherParticle] : Acts::enumerate(truthParticles)) {
      // exclude self
      if (particle.particleId() == otherParticle.particleId()) {
        continue;
      }

      double deltaR =
          Acts::VectorHelpers::deltaR(phis[i], etas[i], phis[j], etas[j]);

      if (deltaR < m_cfg.overlapRemovalIsolationDeltaR) {
        // Add the momentum of the other particle to the total momentum
        totalMomentum += otherParticle.absoluteMomentum();
      }
    }

    double isolation = totalMomentum / particle.absoluteMomentum();
    if (isolation < m_cfg.overlapRemovalIsolation) {
      isolatedLeptons.push_back(&particle);
    }
  }

  ACTS_DEBUG("Number of isolated leptons: " << isolatedLeptons.size());
  ACTS_DEBUG("Number of jets before overlap removal: " << jets.size());

  // erase if jet is within deltaR of any isolated lepton
  std::erase_if(jets, [&](const auto& jet) -> bool {
    for (const auto& lepton : isolatedLeptons) {
      double deltaR =
          Acts::VectorHelpers::deltaR(jet.getDirection(), lepton->direction());
      if (deltaR < m_cfg.overlapRemovalDeltaR) {
        return true;
      }
    }
    return false;
  });

  ACTS_DEBUG("Number of jets after overlap removal: " << jets.size());
}

ProcessCode TruthJetAlgorithm::finalize() {
  ACTS_INFO("TruthJetAlgorithm::finalize:");

  ACTS_INFO("- Total jets         : " << m_numJets);
  ACTS_INFO("- Total jets after OR: " << m_numJetsAfterOverlapRemoval);
  ACTS_INFO("- Light jets         : " << m_numLightJets);
  ACTS_INFO("- C jets             : " << m_numCJets);
  ACTS_INFO("- B jets             : " << m_numBJets);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
