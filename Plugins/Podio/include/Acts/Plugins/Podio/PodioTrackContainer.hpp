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
#include "ActsPodioEdm/TrackCollection.h"

#include <stdexcept>

namespace Acts {

class PodioTrackContainer;
class ConstPodioTrackContainer;

template <>
struct IsReadOnlyTrackContainer<PodioTrackContainer> : std::false_type {};

class PodioTrackContainer {
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
  PodioTrackContainer(ActsPodioEdm::TrackCollection& collection)
      : m_collection{&collection} {}

  PodioTrackContainer(const PodioTrackContainer& other);

  PodioTrackContainer(PodioTrackContainer&& other) = default;

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
      // case "referenceSurface"_hash:
      // return &instance.m_referenceSurfaces[itrack];
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

 public:
  std::any component_impl(HashedString key, IndexType itrack) {
    return component_impl<false>(*this, key, itrack);
  }

  std::any component_impl(HashedString key, IndexType itrack) const {
    return component_impl<true>(*this, key, itrack);
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    return false;
    // using namespace Acts::HashedStringLiteral;
    // switch (key) {
    // default:
    // return m_dynamic.find(key) != m_dynamic.end();
    // }
  }

  std::size_t size_impl() const { return m_collection->size(); }
  // END INTERFACE HELPER

  // std::vector<IndexType> m_tipIndex;
  // std::vector<typename detail_lt::Types<eBoundSize>::Coefficients> m_params;
  // std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_cov;
  // std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  // std::vector<unsigned int> m_nMeasurements;
  // std::vector<unsigned int> m_nHoles;
  // std::vector<float> m_chi2;
  // std::vector<unsigned int> m_ndf;
  // std::vector<unsigned int> m_nOutliers;
  // std::vector<unsigned int> m_nSharedHits;

  // std::unordered_map<HashedString,
  // std::unique_ptr<detail::DynamicColumnBase>> m_dynamic;

  PodioTrackContainer(const ConstPodioTrackContainer& other);

 public:
  // BEGIN INTERFACE

  IndexType addTrack_impl() {
    auto track = m_collection->create();
    return m_collection->size() - 1;
  };

  void removeTrack_impl(IndexType itrack);

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {
    throw std::runtime_error{"Not implemented"};
    // m_dynamic.insert(
    // {hashString(key), std::make_unique<detail::DynamicColumn<T>>()});
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

  void copyDynamicFrom_impl(IndexType dstIdx, const PodioTrackContainer& src,
                            IndexType srcIdx);

  void ensureDynamicColumns_impl(const PodioTrackContainer& other);

  void reserve(IndexType size) {}

  // END INTERFACE

 private:
  ActsPodioEdm::TrackCollection* m_collection;
};
}  // namespace Acts
