// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "BeamSpotConditionsData/BeamSpotData.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "StoreGate/ReadCondHandle.h"
#include "StoreGate/ReadHandle.h"
#include "TrkParameters/TrackParameters.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/Vertex.h"

class AmvfTrackExport final : public AthReentrantAlgorithm {
 public:
  using AthReentrantAlgorithm::AthReentrantAlgorithm;
  StatusCode initialize() override {
    ATH_CHECK(m_tracks.initialize());
    ATH_CHECK(m_event.initialize());
    ATH_CHECK(m_beam.initialize());
    ATH_CHECK(m_selector.retrieve());
    std::filesystem::create_directories(m_output.value());
    return StatusCode::SUCCESS;
  }
  StatusCode execute(const EventContext& ctx) const override {
    SG::ReadHandle<xAOD::TrackParticleContainer> tracks(m_tracks, ctx);
    SG::ReadHandle<xAOD::EventInfo> event(m_event, ctx);
    SG::ReadCondHandle<InDet::BeamSpotData> beam(m_beam, ctx);
    ATH_CHECK(tracks.isValid() && event.isValid() && beam.isValid());
    const auto& beamVertex = beam->beamVtx();
    xAOD::Vertex constraint;
    constraint.makePrivateStore();
    constraint.setPosition(beamVertex.position());
    constraint.setCovariancePosition(beamVertex.covariancePosition());
    if (!m_useBeam) {
      constraint.setPosition(Amg::Vector3D::Zero());
      constraint.setCovariancePosition(AmgSymMatrix(3)::Zero());
    }
    std::ostringstream stem;
    stem << m_output.value() << "/event" << std::setw(9) << std::setfill('0')
         << ctx.evt();
    std::ofstream csv(stem.str() + "-tracks.csv");
    std::ofstream meta(stem.str() + "-metadata.txt");
    ATH_CHECK(csv.good() && meta.good());
    csv << std::setprecision(17);
    csv << "trackId,d0,z0,phi,theta,qop,var_d0,var_z0,var_phi,var_theta,var_"
           "qop";
    constexpr std::array names{"d0", "z0", "phi", "theta", "qop"};
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        if (i != j)
          csv << ",cov_" << names[i] << names[j];
      }
    }
    csv << '\n';
    std::size_t selected = 0;
    std::size_t missingCovariance = 0;
    Amg::Transform3D reference = Amg::Transform3D::Identity();
    bool haveReference = false;
    for (std::size_t id = 0; id < tracks->size(); ++id) {
      const auto& track = *tracks->at(id);
      if (!static_cast<bool>(m_selector->accept(track, &constraint)))
        continue;
      const auto& perigee = track.perigeeParameters();
      // Athena uses the first selected track's surface for the whole
      // collection.
      if (!haveReference) {
        reference = perigee.associatedSurface().transform();
        haveReference = true;
      }
      if (!perigee.covariance()) {
        ++missingCovariance;
        continue;
      }
      ATH_CHECK(reference.matrix().isApprox(
          perigee.associatedSurface().transform().matrix(), 0.));
      auto params = perigee.parameters().eval();
      auto cov = perigee.covariance()->eval();
      // Match Athena's conversion from MeV-based q/p to ACTS GeV-based q/p.
      params[4] *= 1000.;
      cov.row(4) *= 1000.;
      cov.col(4) *= 1000.;
      csv << id;
      for (int i = 0; i < 5; ++i)
        csv << ',' << params[i];
      for (int i = 0; i < 5; ++i)
        csv << ',' << cov(i, i);
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
          if (i != j)
            csv << ',' << cov(i, j);
        }
      }
      csv << '\n';
      ++selected;
    }
    // Simple named records, row-major matrices, mm/GeV units. Kept separate
    // from the standard ACTS TrackParameterData CSV schema.
    meta << std::setprecision(17) << "schema 1\nrun " << event->runNumber()
         << "\nevent " << event->eventNumber() << "\nlumi "
         << event->lumiBlock() << "\nmu "
         << event->averageInteractionsPerCrossing() << "\nactual_mu "
         << event->actualInteractionsPerCrossing() << "\ninput_tracks "
         << tracks->size() << "\nselected_tracks " << selected
         << "\nmissing_covariance " << missingCovariance
         << "\nuse_beam_constraint " << m_useBeam.value() << "\nbeam_position";
    for (int i = 0; i < 3; ++i)
      meta << ' ' << beamVertex.position()[i];
    meta << "\nbeam_covariance";
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        meta << ' ' << beamVertex.covariancePosition()(i, j);
    }
    meta << "\nperigee_transform";
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j)
        meta << ' ' << reference.matrix()(i, j);
    }
    meta << '\n';
    ATH_CHECK(csv.good() && meta.good());
    ATH_MSG_INFO("AMVF export event=" << event->eventNumber() << " mu="
                                      << event->averageInteractionsPerCrossing()
                                      << " tracks=" << tracks->size()
                                      << " selected=" << selected);
    return StatusCode::SUCCESS;
  }

 private:
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_tracks{
      this, "Tracks", "InDetTrackParticles"};
  SG::ReadHandleKey<xAOD::EventInfo> m_event{this, "EventInfo", "EventInfo"};
  SG::ReadCondHandleKey<InDet::BeamSpotData> m_beam{this, "BeamSpotKey",
                                                    "BeamSpotData"};
  ToolHandle<InDet::IInDetTrackSelectionTool> m_selector{this, "TrackSelector",
                                                         ""};
  Gaudi::Property<std::string> m_output{this, "OutputDirectory", "amvf-tracks"};
  Gaudi::Property<bool> m_useBeam{this, "UseBeamConstraint", true};
};
DECLARE_COMPONENT(AmvfTrackExport)
