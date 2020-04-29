// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Utilities/BinAdjustment.hpp"
#include "Acts/Utilities/BinUtility.hpp"

Acts::SurfaceMaterialMapper::SurfaceMaterialMapper(
    const Config& cfg, StraightLinePropagator propagator,
    std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_propagator(std::move(propagator)),
      m_logger(std::move(slogger)) {}

Acts::SurfaceMaterialMapper::State Acts::SurfaceMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const TrackingGeometry& tGeometry) const {
  // Parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();

  // The Surface material mapping state
  State mState(gctx, mctx);
  resolveMaterialSurfaces(mState, *world);

  ACTS_DEBUG(mState.accumulatedMaterial.size()
             << " Surfaces with PROXIES collected ... ");
  for (auto& smg : mState.accumulatedMaterial) {
    ACTS_VERBOSE(" -> Surface in with id " << smg.first);
  }
  return mState;
}

void Acts::SurfaceMaterialMapper::resolveMaterialSurfaces(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- boundary surfaces ...");
  // Check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces()) {
    checkAndInsert(mState, bSurface->surfaceRepresentation());
  }

  ACTS_VERBOSE("- confined layers ...");
  // Check the confined layers
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // Take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // Check the representing surface
        checkAndInsert(mState, cLayer->surfaceRepresentation());
        // Get the approach surfaces if present
        if (cLayer->approachDescriptor() != nullptr) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces()) {
            if (aSurface != nullptr) {
              checkAndInsert(mState, *aSurface);
            }
          }
        }
        // Get the sensitive surface is present
        if (cLayer->surfaceArray() != nullptr) {
          // Sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->surfaces()) {
            if (sSurface != nullptr) {
              checkAndInsert(mState, *sSurface);
            }
          }
        }
      }
    }
  }
  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      resolveMaterialSurfaces(mState, *sVolume);
    }
  }
}

void Acts::SurfaceMaterialMapper::checkAndInsert(State& mState,
                                                 const Surface& surface) const {
  auto surfaceMaterial = surface.surfaceMaterial();
  // Check if the surface has a proxy
  if (surfaceMaterial != nullptr) {
    auto geoID = surface.geoID();
    size_t volumeID = geoID.volume();
    ACTS_DEBUG("Material surface found with volumeID " << volumeID);
    ACTS_DEBUG("       - surfaceID is " << geoID);

    // We need a dynamic_cast to either a surface material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoSurfaceMaterial*>(surfaceMaterial);

    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, surface);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
      mState.accumulatedMaterial[geoID] =
          AccumulatedSurfaceMaterial(buAdjusted);
      return;
    }

    // Second attempt: binned material
    auto bmp = dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
    bu = (bmp != nullptr) ? (&bmp->binUtility()) : nullptr;
    // Creaete a binned type of material
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning is " << *bu);
      mState.accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial(*bu);
    } else {
      // Create a homogeneous type of material
      ACTS_DEBUG("       - this is homogeneous material.");
      mState.accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial();
    }
  }
}

void Acts::SurfaceMaterialMapper::finalizeMaps(State& mState) const {
  // iterate over the map to call the total average
  for (auto& accMaterial : mState.accumulatedMaterial) {
    ACTS_DEBUG("Finalizing map for Surface " << accMaterial.first);
    mState.surfaceMaterial[accMaterial.first] =
        accMaterial.second.totalAverage();
  }
}

void Acts::SurfaceMaterialMapper::mapMaterialTrack(
    State& mState, RecordedMaterialTrack& mTrack) const {
  // Neutral curvilinear parameters
  NeutralCurvilinearParameters start(std::nullopt, mTrack.first.first,
                                     mTrack.first.second, 0.);

  // Prepare Action list and abort list
  using DebugOutput = DebugOutputActor;
  using MaterialSurfaceCollector = SurfaceCollector<MaterialSurface>;
  using MaterialVolumeCollector = VolumeCollector<MaterialVolume>;
  using ActionList = ActionList<MaterialSurfaceCollector,
                                MaterialVolumeCollector, DebugOutput>;
  using AbortList = AbortList<EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options(mState.geoContext,
                                                   mState.magFieldContext);
  options.debug = m_cfg.mapperDebugOutput;

  // Now collect the material layers by using the straight line propagator
  const auto& result = m_propagator.propagate(start, options).value();
  auto mcResult = result.get<MaterialSurfaceCollector::result_type>();
  auto mvcResult = result.get<MaterialVolumeCollector::result_type>();
  // Massive screen output
  if (m_cfg.mapperDebugOutput) {
    auto debugOutput = result.get<DebugOutput::result_type>();
    ACTS_VERBOSE("Debug propagation output.");
    ACTS_VERBOSE(debugOutput.debugString);
  }

  auto mappingSurfaces = mcResult.collected;
  auto mappingVolumes = mvcResult.collected;

  // Retrieve the recorded material from the recorded material track
  auto& rMaterial = mTrack.second.materialInteractions;
  std::map<GeometryID, unsigned int> assignedMaterial;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material steps to map.")

  // These should be mapped onto the mapping surfaces found
  ACTS_VERBOSE("Found     " << mappingSurfaces.size()
                            << " mapping surfaces for this track.");
  ACTS_VERBOSE("Mapping surfaces are :")
  for (auto& mSurface : mappingSurfaces) {
    ACTS_VERBOSE(" - Surface : " << mSurface.surface->geoID()
                                 << " at position = (" << mSurface.position.x()
                                 << ", " << mSurface.position.y() << ", "
                                 << mSurface.position.z() << ")");
    assignedMaterial[mSurface.surface->geoID()] = 0;
  }

  // Run the mapping process, i.e. take the recorded material and map it
  // onto the mapping surfaces:
  // - material steps and surfaces are assumed to be ordered along the
  // mapping ray
  //- do not record the material inside a volume with material
  auto rmIter = rMaterial.begin();
  auto sfIter = mappingSurfaces.begin();
  auto volIter = mappingVolumes.begin();
  bool encounterVolume = false;

  // Use those to minimize the lookup
  GeometryID lastID = GeometryID();
  GeometryID currentID = GeometryID();
  Vector3D currentPos(0., 0., 0);
  double currentPathCorrection = 0.;
  auto currentAccMaterial = mState.accumulatedMaterial.end();

  // To remember the bins of this event
  using MapBin = std::pair<AccumulatedSurfaceMaterial*, std::array<size_t, 3>>;
  std::multimap<AccumulatedSurfaceMaterial*, std::array<size_t, 3>>
      touchedMapBins;

  // Assign the recorded ones, break if you hit an end
  while (rmIter != rMaterial.end() && sfIter != mappingSurfaces.end()) {
    if (volIter != mappingVolumes.end() && encounterVolume == true &&
        !volIter->volume->inside(rmIter->position)) {
      encounterVolume = false;
      // Switch to next material volume
      ++volIter;
    }
    /// check if we are inside a material volume
    if (volIter != mappingVolumes.end() &&
        volIter->volume->inside(rmIter->position)) {
      encounterVolume = true;
      ++rmIter;
      continue;
    }
    if (sfIter != mappingSurfaces.end() - 1 &&
        (rmIter->position - sfIter->position).norm() >
            (rmIter->position - (sfIter + 1)->position).norm()) {
      // Switch to next assignment surface
      ++sfIter;
    }
    // get the current Surface ID
    currentID = sfIter->surface->geoID();
    // We have work to do: the assignemnt surface has changed
    if (not(currentID == lastID)) {
      // Let's (re-)assess the information
      lastID = currentID;
      currentPos = (sfIter)->position;
      currentPathCorrection = sfIter->surface->pathCorrection(
          mState.geoContext, currentPos, sfIter->direction);
      currentAccMaterial = mState.accumulatedMaterial.find(currentID);
    }
    // Now assign the material for the accumulation process
    auto tBin = currentAccMaterial->second.accumulate(
        currentPos, rmIter->materialProperties, currentPathCorrection);
    touchedMapBins.insert(MapBin(&(currentAccMaterial->second), tBin));
    ++assignedMaterial[currentID];
    // Update the material interaction with the associated surface and direction
    rmIter->direction = mTrack.first.second.normalized();
    rmIter->surface = sfIter->surface;
    // Switch to next material
    ++rmIter;
  }

  ACTS_VERBOSE("Surfaces have following number of assigned hits :")
  for (auto& [key, value] : assignedMaterial) {
    ACTS_VERBOSE(" + Surface : " << key << " has " << value << " hits.");
  }

  // After mapping this track, average the touched bins
  for (auto tmapBin : touchedMapBins) {
    std::vector<std::array<size_t, 3>> trackBins = {tmapBin.second};
    tmapBin.first->trackAverage(trackBins);
  }

  // After mapping this track, average the untouched but intersected bins
  if (m_cfg.emptyBinCorrection) {
    // Use the assignedMaterial map to account for empty hits, i.e.
    // the material surface has been intersected by the mapping ray
    // but no material step was assigned to this surface
    for (auto& mSurface : mappingSurfaces) {
      auto mgID = mSurface.surface->geoID();
      // Count an empty hit only if the surface does not appear in the
      // list of assigned surfaces
      if (assignedMaterial[mgID] == 0) {
        auto missedMaterial = mState.accumulatedMaterial.find(mgID);
        missedMaterial->second.trackAverage(mSurface.position, true);
      }
    }
  }
}
