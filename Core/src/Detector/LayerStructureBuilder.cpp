// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/LayerStructureBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Detector/detail/SupportHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <cmath>
#include <cstddef>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

namespace {

/// Check autorange for a given binning
///
/// @param pBinning the proto binning
/// @param extent the extent
/// @param fullPhi indicates whether the full phi range is used
///
void adaptBinningRage(std::vector<Acts::Experimental::ProtoBinning>& pBinning,
                      const Acts::Extent& extent, bool fullPhiBinning) {
  for (auto& pb : pBinning) {
    // Starting values
    Acts::ActsScalar vmin = pb.edges.front();
    Acts::ActsScalar vmax = pb.edges.back();
    // Get the number of bins
    std::size_t nBins = pb.bins();
    // Check if extent overwrites that
    if (extent.constrains(pb.binValue)) {
      const auto& range = extent.range(pb.binValue);
      // Patch the edges values from the range
      vmin = range.min();
      vmax = range.max();
    } else if (pb.binValue == Acts::binPhi && fullPhiBinning) {
      vmin = -M_PI;
      vmax = M_PI;
      pb.boundaryType = Acts::detail::AxisBoundaryType::Closed;
    }
    // Possibly update the edges
    if (pb.axisType == Acts::detail::AxisType::Equidistant) {
      Acts::ActsScalar binWidth = (vmax - vmin) / nBins;
      // Fill the edges
      pb.edges = {vmin};
      pb.edges.resize(nBins + 1);
      for (std::size_t ib = 0; ib <= nBins; ++ib) {
        pb.edges[ib] = vmin + ib * binWidth;
      }
    } else {
      pb.edges.front() = vmin;
      pb.edges.back() = vmax;
    }
  }
}

/// Helper for 1-dimensional generators
///
/// @tparam aType is the axis boundary type: closed or bound
///
/// @param gctx the geometry context
/// @param lSurfaces the surfaces of the layer
/// @param assignToAll the indices assigned to all
/// @param binning the binning struct
///
/// @return a configured surface candidate updators
template <Acts::detail::AxisBoundaryType aType>
Acts::Experimental::SurfaceCandidatesUpdater createUpdater(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Surface>> lSurfaces,
    std::vector<std::size_t> assignToAll,
    const Acts::Experimental::ProtoBinning& binning) {
  // The surface candidate updator & a generator for polyhedrons
  Acts::Experimental::SurfaceCandidatesUpdater sfCandidates;
  Acts::Experimental::detail::PolyhedronReferenceGenerator rGenerator;
  // Indexed Surface generator for this case
  Acts::Experimental::detail::IndexedSurfacesGenerator<
      decltype(lSurfaces), Acts::Experimental::IndexedSurfacesImpl>
      isg{std::move(lSurfaces),
          std::move(assignToAll),
          {binning.binValue},
          {binning.expansion}};
  if (binning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant
    Acts::Experimental::detail::GridAxisGenerators::Eq<aType> aGenerator{
        {binning.edges.front(), binning.edges.back()}, binning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else {
    // Variable
    Acts::Experimental::detail::GridAxisGenerators::Var<aType> aGenerator{
        binning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  }
  return sfCandidates;
}

/// Helper for 2-dimensional generators
///
/// @tparam aType is the axis boundary type, axis a: closed or bound
/// @tparam bType is the axis boundary type, axis b: closed or bound
///
/// @param gctx the geometry context
/// @param lSurfaces the surfaces of the layer
/// @param assignToAll the indices assigned to all
/// @param aBinning the binning struct of axis a
/// @param bBinning the binning struct of axis b
///
/// @return a configured surface candidate updators
template <Acts::detail::AxisBoundaryType aType,
          Acts::detail::AxisBoundaryType bType>
Acts::Experimental::SurfaceCandidatesUpdater createUpdater(
    const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Surface>>& lSurfaces,
    const std::vector<std::size_t>& assignToAll,
    const Acts::Experimental::ProtoBinning& aBinning,
    const Acts::Experimental::ProtoBinning& bBinning) {
  // The surface candidate updator & a generator for polyhedrons
  Acts::Experimental::SurfaceCandidatesUpdater sfCandidates;
  Acts::Experimental::detail::PolyhedronReferenceGenerator rGenerator;
  // Indexed Surface generator for this case
  Acts::Experimental::detail::IndexedSurfacesGenerator<
      decltype(lSurfaces), Acts::Experimental::IndexedSurfacesImpl>
      isg{lSurfaces,
          assignToAll,
          {aBinning.binValue, bBinning.binValue},
          {aBinning.expansion, bBinning.expansion}};
  // Run through the cases
  if (aBinning.axisType == Acts::detail::AxisType::Equidistant &&
      bBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant-Equidistant
    Acts::Experimental::detail::GridAxisGenerators::EqEq<aType, bType>
        aGenerator{{aBinning.edges.front(), aBinning.edges.back()},
                   aBinning.bins(),
                   {bBinning.edges.front(), bBinning.edges.back()},
                   bBinning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else if (bBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Variable-Equidistant
    Acts::Experimental::detail::GridAxisGenerators::VarEq<aType, bType>
        aGenerator{aBinning.edges,
                   {bBinning.edges.front(), bBinning.edges.back()},
                   bBinning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else if (aBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant-Variable
    Acts::Experimental::detail::GridAxisGenerators::EqVar<aType, bType>
        aGenerator{{aBinning.edges.front(), aBinning.edges.back()},
                   aBinning.bins(),
                   bBinning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else {
    // Variable-Variable
    Acts::Experimental::detail::GridAxisGenerators::VarVar<aType, bType>
        aGenerator{aBinning.edges, bBinning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  }
  // Return the candidates
  return sfCandidates;
}

}  // namespace

Acts::Experimental::LayerStructureBuilder::LayerStructureBuilder(
    const Acts::Experimental::LayerStructureBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IInternalStructureBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.surfacesProvider == nullptr) {
    throw std::invalid_argument(
        "LayerStructureBuilder: surfaces provider is nullptr.");
  }
}

Acts::Experimental::InternalStructure
Acts::Experimental::LayerStructureBuilder::construct(
    const Acts::GeometryContext& gctx) const {
  // Trivialities first: internal volumes
  std::vector<std::shared_ptr<DetectorVolume>> internalVolumes = {};
  DetectorVolumeUpdater internalVolumeUpdater = tryNoVolumes();

  // Print the auxiliary information
  if (!m_cfg.auxiliary.empty()) {
    ACTS_DEBUG(m_cfg.auxiliary);
  }

  // Retrieve the layer surfaces
  SurfaceCandidatesUpdater internalCandidatesUpdater =
      tryAllPortalsAndSurfaces();
  auto internalSurfaces = m_cfg.surfacesProvider->surfaces(gctx);
  ACTS_DEBUG("Building internal layer structure from "
             << internalSurfaces.size() << " provided surfaces.");

  // Check whether support structure is scheduled to be built, and if so
  // collect those that should be assigned to all bins
  std::vector<std::size_t> assignToAll = {};
  if (!m_cfg.supports.empty()) {
    ACTS_DEBUG("Adding " << m_cfg.supports.size() << " support structures.")
    // The surface candidate updator
    for (const auto& support : m_cfg.supports) {
      // Check if the supportsurface has already been built
      if (support.surface != nullptr) {
        ACTS_VERBOSE("- Use provided support surface directly.");
        if (support.assignToAll) {
          assignToAll.push_back(internalSurfaces.size());
          ACTS_VERBOSE("  Support surface is assigned to all bins.");
        }
        internalSurfaces.push_back(support.surface);
        continue;
      }

      // Throw an exception is misconfigured
      if (support.type == Surface::SurfaceType::Other) {
        throw std::invalid_argument(
            "LayerStructureBuilder: support surface type not specified.");
      }
      ACTS_VERBOSE("- Build support of type '"
                   << Acts::Surface::s_surfaceTypeNames[support.type] << "'.");
      if (support.splits > 1u) {
        ACTS_VERBOSE("  Support surface is modelled with " << support.splits
                                                           << " planes.");
      }
      // To correctly attach the support structures, estimate the extent
      Extent internalExtent;
      if (m_cfg.extent.has_value()) {
        internalExtent = m_cfg.extent.value();
      } else {
        // Estimate the extent from the surfaces
        for (const auto& s : internalSurfaces) {
          auto sPolyhedron = s->polyhedronRepresentation(gctx, m_cfg.nSegments);
          internalExtent.extend(sPolyhedron.extent(), support.constraints);
        }
      }
      // Use the support bulder helper to add support surfaces
      detail::SupportHelper::addSupport(
          internalSurfaces, assignToAll, internalExtent, support.type,
          support.values, support.transform, support.splits);
    }
  }

  if (internalSurfaces.size() >= m_cfg.nMinimalSurfaces) {
    // Copy as we might patch it with the surface extent
    auto binnings = m_cfg.binnings;

    if (binnings.empty()) {
      ACTS_DEBUG(
          "No surface binning provided, navigation will be 'tryAll' "
          "(potentially slow).");
    } else if (binnings.size() == 1u) {
      // Check if autorange for binning applies
      if (m_cfg.extent.has_value()) {
        ACTS_DEBUG("- adapting the proto binning range to the surface extent.");
        adaptBinningRage(binnings, m_cfg.extent.value(), m_cfg.fullPhiBinning);
      }
      ACTS_DEBUG("- 1-dimensional surface binning detected.");
      // Capture the binning
      auto binning = binnings[0u];
      if (binning.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
        ACTS_VERBOSE("-- closed binning option.");
        internalCandidatesUpdater =
            createUpdater<Acts::detail::AxisBoundaryType::Closed>(
                gctx, internalSurfaces, assignToAll, binning);
      } else {
        ACTS_VERBOSE("-- bound binning option.");
        internalCandidatesUpdater =
            createUpdater<Acts::detail::AxisBoundaryType::Bound>(
                gctx, internalSurfaces, assignToAll, binning);
      }
    } else if (binnings.size() == 2u) {
      // Check if autorange for binning applies
      if (m_cfg.extent.has_value()) {
        ACTS_DEBUG(
            "- adapting the proto binning range(s) to the surface extent.");
        adaptBinningRage(binnings, m_cfg.extent.value(), m_cfg.fullPhiBinning);
      }
      // Sort the binning for conventions
      std::sort(binnings.begin(), binnings.end(),
                [](const ProtoBinning& a, const ProtoBinning& b) {
                  return a.binValue < b.binValue;
                });

      ACTS_DEBUG("- 2-dimensional surface binning detected.");
      // Capture the binnings
      const auto& binning0 = binnings[0u];
      const auto& binning1 = binnings[1u];

      if (binning0.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
        ACTS_VERBOSE("-- closed/bound binning option.");
        internalCandidatesUpdater =
            createUpdater<Acts::detail::AxisBoundaryType::Closed,
                          Acts::detail::AxisBoundaryType::Bound>(
                gctx, internalSurfaces, assignToAll, binning0, binning1);
      } else if (binning1.boundaryType ==
                 Acts::detail::AxisBoundaryType::Closed) {
        ACTS_VERBOSE("-- bound/closed binning option.");
        internalCandidatesUpdater =
            createUpdater<Acts::detail::AxisBoundaryType::Bound,
                          Acts::detail::AxisBoundaryType::Closed>(
                gctx, internalSurfaces, assignToAll, binning0, binning1);
      } else {
        ACTS_VERBOSE("-- bound/bound binning option.");
        internalCandidatesUpdater =
            createUpdater<Acts::detail::AxisBoundaryType::Bound,
                          Acts::detail::AxisBoundaryType::Bound>(
                gctx, internalSurfaces, assignToAll, binning0, binning1);
      }
    }
  } else {
    ACTS_DEBUG("Only " << internalSurfaces.size() << " surfaces provided, "
                       << "navigation will be 'tryAll'");
    ACTS_DEBUG("Per configuration " << m_cfg.nMinimalSurfaces
                                    << " surfaces are "
                                    << "required to use the surface binning.");
  }

  // Check if everything went ok
  if (!internalCandidatesUpdater.connected()) {
    throw std::runtime_error(
        "LayerStructureBuilder: could not connect surface candidate updator.");
  }

  // Return the internal structure
  return InternalStructure{internalSurfaces, internalVolumes,
                           std::move(internalCandidatesUpdater),
                           std::move(internalVolumeUpdater)};
}
