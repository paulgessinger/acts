// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryBackendConcept.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Plugins/Podio/PodioTrackContainer.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>

namespace Acts {

class MutablePodioTrackStateContainer;
class ConstPodioTrackStateContainer;

template <>
struct IsReadOnlyMultiTrajectory<MutablePodioTrackStateContainer>
    : std::false_type {};

template <>
struct IsReadOnlyMultiTrajectory<ConstPodioTrackStateContainer>
    : std::true_type {};

class PodioTrackStateContainerBase {
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

 protected:
  template <typename T>
  static constexpr bool has_impl(T& instance, HashedString key,
                                 IndexType istate) {
    return false;
    // using namespace Acts::HashedStringLiteral;
    // switch (key) {
    // case "predicted"_hash:
    // return instance.m_index[istate].ipredicted != kInvalid;
    // case "filtered"_hash:
    // return instance.m_index[istate].ifiltered != kInvalid;
    // case "smoothed"_hash:
    // return instance.m_index[istate].ismoothed != kInvalid;
    // case "calibrated"_hash:
    // return instance.m_measOffset[istate] != kInvalid;
    // case "calibratedCov"_hash:
    // return instance.m_measCovOffset[istate] != kInvalid;
    // case "jacobian"_hash:
    // return instance.m_index[istate].ijacobian != kInvalid;
    // case "projector"_hash:
    // return instance.m_index[istate].iprojector != kInvalid;
    // case "uncalibratedSourceLink"_hash:
    // return instance.m_sourceLinks[instance.m_index[istate].iuncalibrated]
    // .has_value();
    // case "previous"_hash:
    // case "measdim"_hash:
    // case "referenceSurface"_hash:
    // case "chi2"_hash:
    // case "pathLength"_hash:
    // case "typeFlags"_hash:
    // return true;
    // default:
    // return instance.m_dynamic.find(key) != instance.m_dynamic.end();
  }

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType istate) {
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }
    // using namespace Acts::HashedStringLiteral;
    // switch (key) {
    // case "previous"_hash:
    // return &instance.m_index[istate].iprevious;
    // case "predicted"_hash:
    // return &instance.m_index[istate].ipredicted;
    // case "filtered"_hash:
    // return &instance.m_index[istate].ifiltered;
    // case "smoothed"_hash:
    // return &instance.m_index[istate].ismoothed;
    // case "calibrated"_hash:
    // return &instance.m_measOffset[istate];
    // case "calibratedCov"_hash:
    // return &instance.m_measCovOffset[istate];
    // case "jacobian"_hash:
    // return &instance.m_index[istate].ijacobian;
    // case "projector"_hash:
    // return &instance.m_projectors[instance.m_index[istate].iprojector];
    // case "measdim"_hash:
    // return &instance.m_index[istate].measdim;
    // case "chi2"_hash:
    // return &instance.m_index[istate].chi2;
    // case "pathLength"_hash:
    // return &instance.m_index[istate].pathLength;
    // case "typeFlags"_hash:
    // return &instance.m_index[istate].typeFlags;
    // default:
    // auto it = instance.m_dynamic.find(key);
    // if (it == instance.m_dynamic.end()) {
    // throw std::runtime_error("Unable to handle this component");
    // }
    // std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
    // detail::DynamicColumnBase*>
    // col = it->second.get();
    // assert(col && "Dynamic column is null");
    // return col->get(istate);
    // }
  }

  template <typename T>
  static constexpr bool hasColumn_impl(T& instance, HashedString key) {
    return false;
    // using namespace Acts::HashedStringLiteral;
    // switch (key) {
    // case "predicted"_hash:
    // case "filtered"_hash:
    // case "smoothed"_hash:
    // case "calibrated"_hash:
    // case "calibratedCov"_hash:
    // case "jacobian"_hash:
    // case "projector"_hash:
    // case "previous"_hash:
    // case "uncalibratedSourceLink"_hash:
    // case "referenceSurface"_hash:
    // case "measdim"_hash:
    // case "chi2"_hash:
    // case "pathLength"_hash:
    // case "typeFlags"_hash:
    // return true;
    // default:
    // return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    // }
  }

 public:
  IndexType calibratedSize_impl(IndexType istate) const { return 0; }

  SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return SourceLink{0, 5};
  }

  const Surface* referenceSurface_impl(IndexType istate) const {
    return nullptr;
  }
};

class ConstPodioTrackStateContainer final
    : public PodioTrackStateContainerBase,
      public MultiTrajectory<ConstPodioTrackStateContainer> {
 public:
  ConstParameters parameters_impl(IndexType parIdx) const {
    return ConstParameters{nullptr};
  }

  ConstCovariance covariance_impl(IndexType parIdx) const {
    return ConstCovariance{nullptr};
  }

  ConstCovariance jacobian_impl(IndexType parIdx) const {
    return ConstCovariance{nullptr};
  }

  template <size_t measdim>
  ConstTrackStateProxy::Measurement<measdim> measurement_impl(
      IndexType offset) const {
    return ConstTrackStateProxy::Measurement<measdim>{nullptr};
  }

  template <size_t measdim>
  ConstTrackStateProxy::MeasurementCovariance<measdim>
  measurementCovariance_impl(IndexType offset) const {
    return ConstTrackStateProxy::MeasurementCovariance<measdim>{nullptr};
  }

  IndexType size_impl() const { return 0; }

  std::any component_impl(HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::component_impl<true>(*this, key,
                                                              istate);
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    return PodioTrackStateContainerBase::hasColumn_impl(*this, key);
  }

  constexpr bool has_impl(HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::has_impl(*this, key, istate);
  }
};

ACTS_STATIC_CHECK_CONCEPT(ConstMultiTrajectoryBackend,
                          ConstPodioTrackStateContainer);

}  // namespace Acts
