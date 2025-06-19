// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootJetWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsFatras/EventData/Barcode.hpp"


#include <cmath>
#include <ios>
#include <limits>
#include <numbers>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>


const Acts::MagneticFieldContext magFieldContext;

Acts::MagneticFieldProvider::Cache magFieldCache() {
  return Acts::NullBField{}.makeCache(magFieldContext);
}

namespace ActsExamples {

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

RootJetWriter::RootJetWriter(const RootJetWriter::Config& config,
                             Acts::Logging::Level level)
    : WriterT(config.inputTrackJets, "RootJetWriter", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
  }
  if (m_cfg.inputTrackJets.empty()) {
    throw std::invalid_argument("Missing track jets input collection");
  }
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Missing vertices input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Initialize the input handles
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackJets.initialize(m_cfg.inputTrackJets);
  m_inputVertices.initialize(m_cfg.inputVertices);
  m_inputTrackParticleMatching.initialize(
      m_cfg.inputTrackParticleMatching);
  m_inputParticles.initialize(m_cfg.inputParticles);

  // Tools
  Acts::EigenStepper<> stepper(m_cfg.field);
  m_propagator = std::make_shared<Propagator>(stepper);
  Acts::ImpactPointEstimator::Config ipEstCfg(m_cfg.field, m_propagator);
  m_ipEst = std::make_shared<Acts::ImpactPointEstimator>(ipEstCfg);

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  m_outputTree->Branch("EventNumber", &m_eventNr);

  m_outputTree->Branch("recovertex_x", &m_recovtx_x);
  m_outputTree->Branch("recovertex_y", &m_recovtx_y);
  m_outputTree->Branch("recovertex_z", &m_recovtx_z);
  m_outputTree->Branch("recovertex_t", &m_recovtx_t);
  m_outputTree->Branch("recovertex_sumPt2", &m_recovtx_sumPt2);
  m_outputTree->Branch("recovertex_isHS", &m_recovtx_isHS);
  m_outputTree->Branch("recovertex_isPU", &m_recovtx_isPU);
  m_outputTree->Branch("recovertex_isSec", &m_recovtx_isSec);

  m_outputTree->Branch("secvertex_x", &m_secvtx_x);
  m_outputTree->Branch("secvertex_y", &m_secvtx_y);
  m_outputTree->Branch("secvertex_z", &m_secvtx_z);
  m_outputTree->Branch("secvertex_t", &m_secvtx_t);

  // The jets
  m_outputTree->Branch("jet_pt", &m_jet_pt);
  m_outputTree->Branch("jet_eta", &m_jet_eta);
  m_outputTree->Branch("jet_phi", &m_jet_phi);
  m_outputTree->Branch("jet_ncomponents", &m_jet_ncomponents);
  m_outputTree->Branch("jet_components", &m_jet_components);
  m_outputTree->Branch("jet_tracks_idx", &m_jet_tracks_idx); // for each jet, the indices of the tracks it contains
  m_outputTree->Branch("jet_ntracks", &m_jet_ntracks);
  m_outputTree->Branch("jet_isPU", &m_jet_isPU);
  m_outputTree->Branch("jet_isHS", &m_jet_isHS);
  m_outputTree->Branch("jet_label", &m_jet_label);

  // Tracks in jets
  m_outputTree->Branch("track_matched_secvtx_idx", &m_matched_secvtx_idx); // for each track (that is matched to a jet), the index of the vertex it belongs to
  m_outputTree->Branch("track_matched_jet_idx", &m_matched_jet_idx); // for each track, the index of the jet it belongs to
  m_outputTree->Branch("track_prob", &m_tracks_prob);
  m_outputTree->Branch("track_d0", &m_trk_d0);
  m_outputTree->Branch("track_z0", &m_trk_z0);
  m_outputTree->Branch("track_signedd0", &m_trk_signed_d0);
  m_outputTree->Branch("track_signedd0sig", &m_trk_signed_d0sig);
  m_outputTree->Branch("track_signedz0sinTheta", &m_trk_signed_z0sinTheta);
  m_outputTree->Branch("track_signedz0sinThetasig",
                       &m_trk_signed_z0sinThetasig);
  m_outputTree->Branch("track_eta", &m_trk_eta);
  m_outputTree->Branch("track_theta", &m_trk_theta);
  m_outputTree->Branch("track_phi", &m_trk_phi);
  m_outputTree->Branch("track_pt", &m_trk_pt);
  m_outputTree->Branch("track_qOverP", &m_trk_qOverP);
  m_outputTree->Branch("track_t", &m_trk_t);

  m_outputTree->Branch("track_var_d0",     &m_trk_var_d0);
  m_outputTree->Branch("track_var_z0",     &m_trk_var_z0);
  m_outputTree->Branch("track_var_phi",    &m_trk_var_phi);
  m_outputTree->Branch("track_var_theta",  &m_trk_var_theta);
  m_outputTree->Branch("track_var_qOverP", &m_trk_var_qOverP);

  m_outputTree->Branch("track_cov_d0z0"       ,&m_trk_cov_d0z0);
  m_outputTree->Branch("track_cov_d0phi"      ,&m_trk_cov_d0phi);
  m_outputTree->Branch("track_cov_d0theta"    ,&m_trk_cov_d0theta);
  m_outputTree->Branch("track_cov_d0qOverP"   ,&m_trk_cov_d0qOverP);
  m_outputTree->Branch("track_cov_z0phi"      ,&m_trk_cov_z0phi);
  m_outputTree->Branch("track_cov_z0theta"    ,&m_trk_cov_z0theta);
  m_outputTree->Branch("track_cov_z0qOverP"   ,&m_trk_cov_z0qOverP);
  m_outputTree->Branch("track_cov_phitheta"   ,&m_trk_cov_phitheta);
  m_outputTree->Branch("track_cov_phiqOverP"  ,&m_trk_cov_phiqOverP);
  m_outputTree->Branch("track_cov_tehtaqOverP",&m_trk_cov_thetaqOverP);
}

RootJetWriter::~RootJetWriter() {
 // m_outputFile->Close();
}

ProcessCode RootJetWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote states of trajectories to tree '"
            << m_cfg.treeName << "' in '" << m_cfg.treeName << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootJetWriter::writeT(const AlgorithmContext& ctx,
                                  const TrackJetContainer& tp) {
  // Exclusive access to all the writing procedure
  std::lock_guard<std::mutex> lock(m_writeMutex);

  Clear();

  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS; // make it an error condition
  }

  // Read input collections

  const auto& tracks = m_inputTracks(ctx);
  const auto& trackJets = m_inputTrackJets(ctx);
  const auto& vertices = m_inputVertices(ctx);
  const auto& trackParticleMatching =
      m_inputTrackParticleMatching(ctx);
  const auto& particles = m_inputParticles(ctx);

  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputTracks << " "
                                        << tracks.size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputTrackJets << " "
                                        << trackJets.size());
  ACTS_VERBOSE("RootWriter::Number of " << m_cfg.inputVertices << " "
                                        << m_inputVertices(ctx).size());

  m_eventNr = ctx.eventNumber;

  // Read the vertices and find the HS vertex
  // Find the vertex with the highest sumPt2
  double maxSumPt2 = -1.;
  int HS_idx = -1;

  for (int v = 0; v < vertices.size(); v++) {
    double sumPt2 = calcSumPt2(vertices.at(v));
    if (sumPt2 > maxSumPt2) {
      maxSumPt2 = sumPt2;
      HS_idx = v;
    }
  }

  ACTS_VERBOSE("Vertex " << HS_idx << " with sumPt2: " << maxSumPt2 << " has "
                         << vertices.at(HS_idx).tracks().size()
                         << " tracks associated");

  Acts::Vertex hs_vtx = vertices.at(HS_idx);
  Acts::Vector3 vertexPosition =
      vertices.at(HS_idx).position();  // Hard scatter vertex position
  Acts::Vector4 stddev;
  stddev[Acts::ePos0] = 10 * Acts::UnitConstants::um;
  stddev[Acts::ePos1] = 10 * Acts::UnitConstants::um;
  stddev[Acts::ePos2] = 75 * Acts::UnitConstants::um;
  stddev[Acts::eTime] = 1 * Acts::UnitConstants::ns;
  Acts::SquareMatrix4 vertexCov = stddev.cwiseProduct(stddev).asDiagonal();

  Acts::ImpactPointEstimator::State state{magFieldCache()};
  std::vector<Acts::Vector4> secondaryVertices;

  // Loop over the tracks
  for (std::size_t itrk = 0; itrk < tracks.size(); itrk++) {
    double signed_d0 = -9999;
    double signed_z0SinTheta = -9999;
    double signed_d0_err = 1.;
    double signed_z0SinTheta_err = 1.;

    int matched_jet_idx = -999;
    int isecvtx = -999; // flag to check if we found a secondary vertex

    double covd0 = -999;
    double covz0 = -999;
    double covphi = -999;
    double covtheta = -999;
    double covqOverP = -999;
    double covd0z0 = -999;
    double covd0phi = -999;
    double covd0theta = -999;
    double covd0qOverP = -999;
    double covz0phi = -999;
    double covz0theta = -999;
    double covz0qOverP = -999;
    double covphitheta = -999;
    double covphiqOverP = -999;
    double covthetaqOverP = -999;

    const auto trk = tracks.at(itrk);
    const auto params = trk.parameters();

    auto boundParams = trk.createParametersAtReference();

    // Check if this track belongs to a jet and compute the IPs
    for (std::size_t ijet = 0; ijet < trackJets.size(); ++ijet) {
      std::vector<int> jtrks = trackJets[ijet].getTracks();
      

      if (std::find(jtrks.begin(), jtrks.end(), itrk) != jtrks.end()) {
        ACTS_DEBUG("Track " << itrk << " is in jet " << ijet);
        matched_jet_idx = ijet;

        Acts::Vector3 jetDir{trackJets[ijet].getFourMomentum()[0],
                             trackJets[ijet].getFourMomentum()[1],
                             trackJets[ijet].getFourMomentum()[2]};

        jetDir = jetDir.normalized();  // normalize the jet direction

        auto ipAndSigma =
            m_ipEst->getImpactParameters(boundParams, hs_vtx, gctx_, mctx_, 1);
        auto vszs = m_ipEst->getLifetimeSignOfTrack(boundParams, hs_vtx, jetDir,
                                                    gctx_, mctx_);

        if (!ipAndSigma.ok() || !vszs.ok()) {
          continue;
        }

        signed_d0 = std::fabs((*ipAndSigma).d0) * (*vszs).first;
        auto z0 = std::fabs((*ipAndSigma).z0);
        auto theta = boundParams.theta();
        signed_z0SinTheta = z0 * std::sin(theta) * (*vszs).second;
        signed_d0_err = (*ipAndSigma).sigmaD0;
        signed_z0SinTheta_err =
            (*ipAndSigma).sigmaZ0;  // not correct yet - have to take into
                                    // account the sin(theta) factor

        // Get the associated truth particle for this track
        // if matched, take this truth particle's initial ? position 
        // make a vector for 2ndvtxpos (out of the loop) am I close enough to the position of the truth particle?

        auto match = trackParticleMatching.find(itrk);
        if(match != trackParticleMatching.end() && match->second.particle.has_value()) {
          auto barcode = match->second.particle.value();
          auto findParticle = particles.find(barcode);

          if(findParticle != particles.end()) {
            auto initParticle = *findParticle;
            Acts::Vector4 initParticlePosition = initParticle.fourPosition();

              auto it = std::ranges::find_if(
                  secondaryVertices.begin(), secondaryVertices.end(),
                  [&initParticlePosition](const Acts::Vector4& vtx) {
                    return (vtx - initParticlePosition).norm() < 1e-3;  // tolerance // can also make it configurable
                  });
              if(it == secondaryVertices.end()) {
                secondaryVertices.push_back(initParticlePosition);
                isecvtx = secondaryVertices.size() - 1; // index of the new secondary vertex
              }
              else {
              isecvtx = std::distance(secondaryVertices.begin(), it);
              }

            } // if found
          } // if match
      }  // track in jet
    }  // loop on jets

    float trk_theta = params[Acts::eBoundTheta];
    float trk_eta = std::atanh(std::cos(trk_theta));
    float trk_qop = params[Acts::eBoundQOverP];
    float trk_p = std::abs(1.0 / trk_qop);
    float trk_pt = trk_p * std::sin(trk_theta);

    m_tracks_prob.push_back(1.);  // todo
    m_trk_d0.push_back(params[Acts::eBoundLoc0]);
    m_trk_z0.push_back(params[Acts::eBoundLoc1]);

    m_trk_signed_d0.push_back(signed_d0);
    m_trk_signed_d0sig.push_back(signed_d0 / signed_d0_err);
    m_trk_signed_z0sinTheta.push_back(signed_z0SinTheta);
    m_trk_signed_z0sinThetasig.push_back(signed_z0SinTheta /
                                         signed_z0SinTheta_err);

    m_matched_jet_idx.push_back(matched_jet_idx); // for each track, the index of the jet it belongs to
    m_matched_secvtx_idx.push_back(isecvtx); // fill number of secondary vertices for each track that are matched to a jet HINT: this might also include the hard scatter vertex

    if (isecvtx > -999) {
      m_secvtx_x.push_back(secondaryVertices[isecvtx][Acts::ePos0]);
      m_secvtx_y.push_back(secondaryVertices[isecvtx][Acts::ePos1]);
      m_secvtx_z.push_back(secondaryVertices[isecvtx][Acts::ePos2]);
      m_secvtx_t.push_back(secondaryVertices[isecvtx][Acts::eTime]); }
    else {
      m_secvtx_x.push_back(-9999.);
      m_secvtx_y.push_back(-9999.);
      m_secvtx_z.push_back(-9999.);
      m_secvtx_t.push_back(-9999.);
    }

    m_trk_eta.push_back(trk_eta);
    m_trk_theta.push_back(trk_theta);
    m_trk_phi.push_back(params[Acts::eBoundPhi]);
    m_trk_pt.push_back(trk_pt);
    m_trk_qOverP.push_back(trk_p);
    m_trk_t.push_back(params[Acts::eBoundTime]);

    // //This give un-initialized warnings

    auto cov = trk.covariance();

    covd0          = cov(Acts::eBoundLoc0,Acts::eBoundLoc0);
    covz0          = cov(Acts::eBoundLoc1,Acts::eBoundLoc1);
    covphi         = cov(Acts::eBoundPhi,Acts::eBoundPhi);
    covtheta       = cov(Acts::eBoundTheta,Acts::eBoundTheta);
    covqOverP      = cov(Acts::eBoundQOverP,Acts::eBoundQOverP);
    covd0z0        = cov(Acts::eBoundLoc0,Acts::eBoundLoc1);
    covd0phi       = cov(Acts::eBoundLoc0,Acts::eBoundPhi);
    covd0theta     = cov(Acts::eBoundLoc0,Acts::eBoundTheta);
    covd0qOverP    = cov(Acts::eBoundLoc0,Acts::eBoundQOverP);
    covz0phi       = cov(Acts::eBoundLoc1,Acts::eBoundPhi);
    covz0theta     = cov(Acts::eBoundLoc1,Acts::eBoundTheta);
    covz0qOverP    = cov(Acts::eBoundLoc1,Acts::eBoundQOverP);
    covphitheta    = cov(Acts::eBoundPhi,Acts::eBoundTheta);
    covphiqOverP   = cov(Acts::eBoundPhi,Acts::eBoundQOverP);
    covthetaqOverP = cov(Acts::eBoundTheta,Acts::eBoundQOverP);

    m_trk_var_d0.push_back(covd0);
    m_trk_var_z0.push_back(covz0);
    m_trk_var_phi.push_back(covphi);
    m_trk_var_theta.push_back(covtheta);
    m_trk_var_qOverP.push_back(covqOverP);
    m_trk_cov_d0z0.push_back(covd0z0);
    m_trk_cov_d0phi.push_back(covd0phi);
    m_trk_cov_d0theta.push_back(covd0theta);
    m_trk_cov_d0qOverP.push_back(covd0qOverP);
    m_trk_cov_z0phi.push_back(covz0phi);
    m_trk_cov_z0theta.push_back(covz0theta);
    m_trk_cov_z0qOverP.push_back(covz0qOverP);
    m_trk_cov_phitheta.push_back(covphitheta);
    m_trk_cov_phiqOverP.push_back(covphiqOverP);
    m_trk_cov_thetaqOverP.push_back(covthetaqOverP);

  

  }  // loop over tracks

  for (std::size_t ijets = 0; ijets < trackJets.size(); ++ijets) {
    Acts::Vector4 jet_4mom = trackJets[ijets].getFourMomentum();
    Acts::Vector3 jet_3mom{jet_4mom[0], jet_4mom[1], jet_4mom[2]};
    float jet_theta = theta(jet_3mom);
    m_jet_pt.push_back(perp(jet_4mom));
    m_jet_eta.push_back(std::atanh(std::cos(jet_theta)));
    m_jet_phi.push_back(phi(jet_4mom));

    m_jet_ncomponents.push_back(trackJets[ijets].getConstituents().size());
    m_jet_components.push_back(trackJets[ijets].getConstituents());
    m_jet_tracks_idx.push_back(trackJets[ijets].getTracks());
    m_jet_ntracks.push_back(trackJets[ijets].getTracks().size());
    m_jet_label.push_back(static_cast<int>(trackJets[ijets].getLabel()));

    m_jet_isPU.push_back(0);  // Need to check....
    m_jet_isHS.push_back(1);

  }  // jets

  // Fill the reconstructed vertex information
  for (std::size_t ivtx = 0; ivtx < vertices.size(); ivtx++) {
    const auto& vtx = vertices.at(ivtx);
    m_recovtx_x.push_back(vtx.position()[Acts::ePos0]);
    m_recovtx_y.push_back(vtx.position()[Acts::ePos1]);
    m_recovtx_z.push_back(vtx.position()[Acts::ePos2]);
    m_recovtx_t.push_back(vtx.time());
    m_recovtx_sumPt2.push_back(calcSumPt2(vtx));

    if (&vtx == &m_inputVertices(ctx).at(HS_idx)) {
      m_recovtx_isHS.push_back(1);
      m_recovtx_isPU.push_back(0);
    } else {
      m_recovtx_isHS.push_back(0);
      m_recovtx_isPU.push_back(1);
    }
    
  }  // vertices

  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

void ActsExamples::RootJetWriter::Clear() {
  // Vertices
  m_recovtx_x.clear();
  m_recovtx_y.clear();
  m_recovtx_z.clear();
  m_recovtx_t.clear();
  m_recovtx_sumPt2.clear();
  m_recovtx_isHS.clear();
  m_recovtx_isPU.clear();
  m_recovtx_isSec.clear();
  m_matched_secvtx_idx.clear(); // for each track (that is matched to a jet), the index of the vertex it belongs to

  // clear secondary vertex information
  m_secvtx_x.clear();
  m_secvtx_y.clear();
  m_secvtx_z.clear();
  m_secvtx_t.clear();

  // Jets
  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_ncomponents.clear();
  m_jet_components.clear();
  m_jet_tracks_idx.clear();
  m_matched_jet_idx.clear(); // for each track, the index of the jet it belongs to
  m_jet_ntracks.clear();
  m_jet_isPU.clear();
  m_jet_isHS.clear();
  m_jet_label.clear();

  // Tracks
  m_tracks_prob.clear();
  m_trk_d0.clear();
  m_trk_z0.clear();
  m_trk_signed_d0.clear();
  m_trk_signed_d0sig.clear();
  m_trk_signed_z0sinTheta.clear();
  m_trk_signed_z0sinThetasig.clear();
  m_trk_eta.clear();
  m_trk_theta.clear();
  m_trk_phi.clear();
  m_trk_pt.clear();
  m_trk_qOverP.clear();
  m_trk_t.clear();
    m_trk_var_d0.clear();
    m_trk_var_z0.clear();
    m_trk_var_phi.clear();
    m_trk_var_theta.clear();
    m_trk_var_qOverP.clear();
    m_trk_cov_d0z0.clear();
    m_trk_cov_d0phi.clear();
    m_trk_cov_d0theta.clear();
    m_trk_cov_d0qOverP.clear();
    m_trk_cov_z0phi.clear();
    m_trk_cov_z0theta.clear();
    m_trk_cov_z0qOverP.clear();
    m_trk_cov_phitheta.clear();
    m_trk_cov_phiqOverP.clear();
    m_trk_cov_thetaqOverP.clear();
}

double ActsExamples::RootJetWriter::calcSumPt2(const Acts::Vertex& vtx) {
  double sumPt2 = 0;
  const double minTrkWeight = 0.0;
  for (const auto& trk : vtx.tracks()) {
    if (trk.trackWeight > minTrkWeight) {
      double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                      ->transverseMomentum();
      sumPt2 += pt * pt;
    }
  }
  return sumPt2;
}
