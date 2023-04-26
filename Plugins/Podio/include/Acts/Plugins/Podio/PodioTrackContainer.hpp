// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/detail/DynamicColumn.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsPodioEdm/Track.h"
#include "ActsPodioEdm/TrackCollection.h"

#include <mutex>
#include <stdexcept>

namespace Acts {

class MutablePodioTrackContainer;
class ConstPodioTrackContainer;

template <>
struct IsReadOnlyTrackContainer<MutablePodioTrackContainer> : std::false_type {
};

class MutablePodioTrackContainer {
 public:
  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
  static constexpr auto MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

  using Parameters =
      typename detail_lt::Types<eBoundSize, false>::CoefficientsMap;
  using Covariance =
      typename detail_lt::Types<eBoundSize, false>::CovarianceMap;

  using ConstParameters =
      typename detail_lt::Types<eBoundSize, true>::CoefficientsMap;
  using ConstCovariance =
      typename detail_lt::Types<eBoundSize, true>::CovarianceMap;

 public:
  MutablePodioTrackContainer(const PodioUtil::ConversionHelper& helper,
                             ActsPodioEdm::TrackCollection& collection)
      : m_collection{&collection}, m_helper{helper} {
    m_surfaces.reserve(m_collection->size());
    for (ActsPodioEdm::Track track : *m_collection) {
      m_surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          m_helper, track.getReferenceSurface()));
    }
  }

  MutablePodioTrackContainer(const MutablePodioTrackContainer& other);
  MutablePodioTrackContainer(MutablePodioTrackContainer&& other) = default;

  MutablePodioTrackContainer(const ConstPodioTrackContainer& other);

  // BEGIN INTERFACE HELPER

 private:
  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType itrack) {
    using namespace Acts::HashedStringLiteral;
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }

    using namespace Acts::HashedStringLiteral;
    auto track = instance.m_collection->at(itrack);
    auto& data = track.data();
    switch (key) {
      case "tipIndex"_hash:
        return &data.tipIndex;
      case "params"_hash:
        return data.parameters.data();
      case "cov"_hash:
        return data.covariance.data();
      case "referenceSurface"_hash:
        if constexpr (EnsureConst) {
          return instance.getSurface(itrack);
        } else {
          return instance.getOrCreateSurface(itrack);
        }
      case "nMeasurements"_hash:
        return &data.nMeasurements;
      case "nHoles"_hash:
        return &data.nHoles;
      case "chi2"_hash:
        return &data.chi2;
      case "ndf"_hash:
        return &data.ndf;
      case "nOutliers"_hash:
        return &data.nOutliers;
      case "nSharedHits"_hash:
        return &data.nSharedHits;
      default:
        // auto it = instance.m_dynamic.find(key);
        // if (it == instance.m_dynamic.end()) {
        throw std::runtime_error("Unable to handle this component");
    }

    // std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
    // detail::DynamicColumnBase*>
    // col = it->second.get();
    // assert(col && "Dynamic column is null");
    // return col->get(itrack);
    // }
  }

  std::shared_ptr<const Surface> getOrCreateSurface(IndexType itrack) {
    std::shared_ptr<const Surface>& ptr = m_surfaces.at(itrack);
    if (!ptr) {
      ActsPodioEdm::Track track = m_collection->at(itrack);
      ptr = PodioUtil::convertSurfaceFromPodio(m_helper,
                                               track.getReferenceSurface());
    }
    return ptr;
  }

  std::shared_ptr<const Surface> getSurface(IndexType itrack) const {
    return m_surfaces.at(itrack);
  }

 public:
  std::any component_impl(HashedString key, IndexType itrack) {
    return component_impl<false>(*this, key, itrack);
  }

  std::any component_impl(HashedString key, IndexType itrack) const {
    return component_impl<true>(*this, key, itrack);
  }

  constexpr bool hasColumn_impl(HashedString /*key*/) const { return false; }

  std::size_t size_impl() const { return m_collection->size(); }
  // END INTERFACE HELPER

  const Surface& referenceSurface_impl(IndexType itrack) const {
    return *m_surfaces.at(itrack);
  }

  void setReferenceSurface_impl(IndexType itrack,
                                std::shared_ptr<const Surface> surface) {
    auto track = m_collection->at(itrack);
    track.setReferenceSurface(
        PodioUtil::convertSurfaceToPodio(m_helper, *surface));
    m_surfaces.at(itrack) = std::move(surface);
  }

 public:
  // BEGIN INTERFACE

  IndexType addTrack_impl() {
    auto track = m_collection->create();
    m_surfaces.emplace_back();
    return m_collection->size() - 1;
  };

  void removeTrack_impl(IndexType itrack);

  template <typename T>
  constexpr void addColumn_impl(const std::string& /*key*/) {
    throw std::runtime_error{"addColumn not implemented"};
  }

  Parameters parameters(IndexType itrack) {
    return Parameters{m_collection->at(itrack).data().parameters.data()};
  }

  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{m_collection->at(itrack).data().parameters.data()};
  }

  Covariance covariance(IndexType itrack) {
    return Covariance{m_collection->at(itrack).data().covariance.data()};
  }

  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{m_collection->at(itrack).data().covariance.data()};
  }

  // @TODO What's the equivalent of this?
  // void copyDynamicFrom_impl(IndexType dstIdx, const PodioTrackContainer& src,
  // IndexType srcIdx);

  // void ensureDynamicColumns_impl(const PodioTrackContainer& other);

  void reserve(IndexType /*size*/) {}

  // END INTERFACE

 private:
  ActsPodioEdm::TrackCollection* m_collection;
  std::reference_wrapper<const PodioUtil::ConversionHelper> m_helper;
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};
}  // namespace Acts
