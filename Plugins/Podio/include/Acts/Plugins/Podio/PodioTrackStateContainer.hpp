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
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Plugins/Podio/PodioTrackContainer.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"

#include <any>

namespace Acts {

class MutablePodioTrackStateContainer;
class ConstPodioTrackStateContainer;

class PodioTrackStateContainerBase {
 public:
  using IndexType = MultiTrajectoryTraits::IndexType;

  using Parameters =
      typename detail_lt::Types<eBoundSize, false>::CoefficientsMap;
  using Covariance =
      typename detail_lt::Types<eBoundSize, false>::CovarianceMap;

  using ConstParameters =
      typename detail_lt::Types<eBoundSize, true>::CoefficientsMap;
  using ConstCovariance =
      typename detail_lt::Types<eBoundSize, true>::CovarianceMap;

 protected:
  // PodioTrackStateContainerBase(const PodioUtil::ConversionHelper& helper)
  // : m_helper{helper} {}

  template <typename T>
  static constexpr bool has_impl(T& instance, HashedString key,
                                 IndexType istate) {
    constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
    using namespace Acts::HashedStringLiteral;
    auto trackState = instance.m_collection.at(istate);
    const auto& data = trackState.getData();
    switch (key) {
      case "predicted"_hash:
        return data.ipredicted != kInvalid;
      case "filtered"_hash:
        return data.ifiltered != kInvalid;
      case "smoothed"_hash:
        return data.ismoothed != kInvalid;
      case "calibrated"_hash:
        return data.measdim != 0;
      case "calibratedCov"_hash:
        return data.measdim != 0;
      case "jacobian"_hash:
        return data.ijacobian != kInvalid;
      case "projector"_hash:
        return data.hasProjector;
      case "uncalibratedSourceLink"_hash:
        return data.uncalibratedIdentifier != PodioUtil::kNoIdentifier;
      case "previous"_hash:
      case "measdim"_hash:
      case "referenceSurface"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        throw std::runtime_error("Unable to handle this component");
        // return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    }

    return false;
  }

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType istate) {
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }
    using namespace Acts::HashedStringLiteral;
    auto trackState = instance.m_collection.at(istate);
    std::conditional_t<EnsureConst, const ActsPodioEdm::TrackStateInfo*,
                       ActsPodioEdm::TrackStateInfo*>
        dataPtr;
    if constexpr (EnsureConst) {
      dataPtr = &trackState.getData();
    } else {
      dataPtr = &trackState.data();
    }
    auto& data = *dataPtr;
    switch (key) {
      case "previous"_hash:
        return &data.previous;
      case "predicted"_hash:
        return &data.ipredicted;
      case "filtered"_hash:
        return &data.ifiltered;
      case "smoothed"_hash:
        return &data.ismoothed;
      case "projector"_hash:
        return &data.projector;
      case "measdim"_hash:
        return &data.measdim;
      case "chi2"_hash:
        return &data.chi2;
      case "pathLength"_hash:
        return &data.pathLength;
      case "typeFlags"_hash:
        return &data.typeFlags;
      default:
        throw std::runtime_error("Unable to handle this component");
        // auto it = instance.m_dynamic.find(key);
        // if (it == instance.m_dynamic.end()) {
        // throw std::runtime_error("Unable to handle this component");
        // }
        // std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
        // detail::DynamicColumnBase*>
        // col = it->second.get();
        // assert(col && "Dynamic column is null");
        // return col->get(istate);
    }
  }

  template <typename T>
  static constexpr bool hasColumn_impl(T& instance, HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      // case "predicted"_hash:
      // case "filtered"_hash:
      // case "smoothed"_hash:
      // case "calibrated"_hash:
      // case "calibratedCov"_hash:
      // case "jacobian"_hash:
      // case "projector"_hash:
      // case "previous"_hash:
      case "uncalibratedSourceLink"_hash:
      // case "referenceSurface"_hash:
      // case "measdim"_hash:
      // case "chi2"_hash:
      // case "pathLength"_hash:
      // case "typeFlags"_hash:
      // return true;
      default:
        throw std::runtime_error("Unable to handle this component");
        // return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    }
    return false;
  }
};

template <>
struct IsReadOnlyMultiTrajectory<ConstPodioTrackStateContainer>
    : std::true_type {};

class ConstPodioTrackStateContainer final
    : public PodioTrackStateContainerBase,
      public MultiTrajectory<ConstPodioTrackStateContainer> {
 public:
  ConstPodioTrackStateContainer(const PodioUtil::ConversionHelper& helper)
      : m_helper{helper} {
    populateSurfaceBuffer();
  }

  ConstPodioTrackStateContainer(const MutablePodioTrackStateContainer& other);
  ConstPodioTrackStateContainer(MutablePodioTrackStateContainer&& other);

  ConstParameters parameters_impl(IndexType istate) const {
    return ConstParameters{m_params.at(istate).getData().values.data()};
  }

  ConstCovariance covariance_impl(IndexType istate) const {
    return ConstCovariance{m_params.at(istate).getData().covariance.data()};
  }

  ConstCovariance jacobian_impl(IndexType istate) const {
    IndexType ijacobian = m_collection.at(istate).getData().ijacobian;
    return ConstCovariance{m_jacs.at(ijacobian).getData().values.data()};
  }

  template <size_t measdim>
  ConstTrackStateProxy::Measurement<measdim> measurement_impl(
      IndexType index) const {
    return ConstTrackStateProxy::Measurement<measdim>{
        m_collection.at(index).getData().measurement.data()};
  }

  template <size_t measdim>
  ConstTrackStateProxy::MeasurementCovariance<measdim>
  measurementCovariance_impl(IndexType index) const {
    return ConstTrackStateProxy::MeasurementCovariance<measdim>{
        m_collection.at(index).getData().measurementCovariance.data()};
  }

  IndexType size_impl() const { return m_collection.size(); }

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

  MultiTrajectoryTraits::IndexType calibratedSize_impl(IndexType istate) const {
    return m_collection.at(istate).getData().measdim;
  }

  SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return m_helper.get().identifierToSourceLink(
        m_collection.at(istate).getData().uncalibratedIdentifier);
  }

  const Surface* referenceSurface_impl(IndexType istate) const {
    return m_surfaces.at(istate).get();
  }

 private:
  void populateSurfaceBuffer() noexcept {
    m_surfaces.reserve(m_collection->size());
    for (ActsPodioEdm::TrackState trackState : m_collection) {
      m_surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          m_helper, trackState.getReferenceSurface()));
    }
  }

  friend class PodioTrackStateContainerBase;

  std::reference_wrapper<const PodioUtil::ConversionHelper> m_helper;
  ActsPodioEdm::TrackStateCollection m_collection;
  ActsPodioEdm::BoundParametersCollection m_params;
  ActsPodioEdm::JacobianCollection m_jacs;
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};

static_assert(IsReadOnlyMultiTrajectory<ConstPodioTrackStateContainer>::value,
              "MutablePodioTrackStateContainer should not be read-only");

ACTS_STATIC_CHECK_CONCEPT(ConstMultiTrajectoryBackend,
                          ConstPodioTrackStateContainer);

template <>
struct IsReadOnlyMultiTrajectory<MutablePodioTrackStateContainer>
    : std::false_type {};

class MutablePodioTrackStateContainer final
    : public PodioTrackStateContainerBase,
      public MultiTrajectory<MutablePodioTrackStateContainer> {
 public:
  MutablePodioTrackStateContainer(PodioUtil::ConversionHelper& helper)
      : m_helper{helper} {
    m_surfaces.reserve(m_collection->size());
    for (ActsPodioEdm::TrackState trackState : m_collection) {
      m_surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          m_helper, trackState.getReferenceSurface()));
    }
  }

  ConstParameters parameters_impl(IndexType istate) const {
    return ConstParameters{m_params.at(istate).getData().values.data()};
  }

  Parameters parameters_impl(IndexType istate) {
    return Parameters{m_params.at(istate).data().values.data()};
  }

  ConstCovariance covariance_impl(IndexType istate) const {
    return ConstCovariance{m_params.at(istate).getData().covariance.data()};
  }

  Covariance covariance_impl(IndexType istate) {
    return Covariance{m_params.at(istate).data().covariance.data()};
  }

  ConstCovariance jacobian_impl(IndexType istate) const {
    IndexType ijacobian = m_collection.at(istate).getData().ijacobian;
    return ConstCovariance{m_jacs.at(ijacobian).getData().values.data()};
  }

  Covariance jacobian_impl(IndexType istate) {
    IndexType ijacobian = m_collection.at(istate).getData().ijacobian;
    return Covariance{m_jacs.at(ijacobian).data().values.data()};
  }

  template <size_t measdim>
  ConstTrackStateProxy::Measurement<measdim> measurement_impl(
      IndexType index) const {
    return ConstTrackStateProxy::Measurement<measdim>{
        m_collection.at(index).getData().measurement.data()};
  }

  template <size_t measdim>
  TrackStateProxy::Measurement<measdim> measurement_impl(IndexType index) {
    return TrackStateProxy::Measurement<measdim>{
        m_collection->at(index).data().measurement.data()};
  }

  template <size_t measdim>
  ConstTrackStateProxy::MeasurementCovariance<measdim>
  measurementCovariance_impl(IndexType index) const {
    return ConstTrackStateProxy::MeasurementCovariance<measdim>{
        m_collection.at(index).getData().measurementCovariance.data()};
  }

  template <size_t measdim>
  TrackStateProxy::MeasurementCovariance<measdim> measurementCovariance_impl(
      IndexType index) {
    return TrackStateProxy::MeasurementCovariance<measdim>{
        m_collection->at(index).data().measurementCovariance.data()};
  }

  IndexType size_impl() const { return m_collection.size(); }

  std::any component_impl(HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::component_impl<true>(*this, key,
                                                              istate);
  }

  std::any component_impl(HashedString key, IndexType istate) {
    return PodioTrackStateContainerBase::component_impl<false>(*this, key,
                                                               istate);
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    return PodioTrackStateContainerBase::hasColumn_impl(*this, key);
  }

  constexpr bool has_impl(HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::has_impl(*this, key, istate);
  }

  IndexType addTrackState_impl(
      TrackStatePropMask mask = TrackStatePropMask::All,
      TrackIndexType iprevious = kTrackIndexInvalid) {
    auto trackState = m_collection->create();
    auto& data = trackState.data();
    data.previous = iprevious;
    data.ipredicted = kInvalid;
    data.ifiltered = kInvalid;
    data.ismoothed = kInvalid;
    data.ijacobian = kInvalid;
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
      m_params->create();
      data.ipredicted = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
      m_params->create();
      data.ifiltered = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
      m_params->create();
      data.ismoothed = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
      m_jacs->create();
      data.ijacobian = m_jacs->size() - 1;
    }
    data.measdim = 0;
    data.hasProjector = false;
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
      data.hasProjector = true;
    }
    m_surfaces.emplace_back();

    data.uncalibratedIdentifier = PodioUtil::kNoIdentifier;
    assert(m_collection->size() == m_surfaces.size() &&
           "Inconsistent surface buffer");
    return m_collection->size() - 1;
  }

  void shareFrom_impl(TrackIndexType iself, TrackIndexType iother,
                      TrackStatePropMask shareSource,
                      TrackStatePropMask shareTarget) {
    auto& self = m_collection->at(iself).data();
    auto& other = m_collection->at(iother).data();

    assert(ACTS_CHECK_BIT(getTrackState(iother).getMask(), shareSource) &&
           "Source has incompatible allocation");

    using PM = TrackStatePropMask;

    IndexType sourceIndex{kInvalid};
    switch (shareSource) {
      case PM::Predicted:
        sourceIndex = other.ipredicted;
        break;
      case PM::Filtered:
        sourceIndex = other.ifiltered;
        break;
      case PM::Smoothed:
        sourceIndex = other.ismoothed;
        break;
      case PM::Jacobian:
        sourceIndex = other.ijacobian;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }

    assert(sourceIndex != kInvalid);

    switch (shareTarget) {
      case PM::Predicted:
        assert(shareSource != PM::Jacobian);
        self.ipredicted = sourceIndex;
        break;
      case PM::Filtered:
        assert(shareSource != PM::Jacobian);
        self.ifiltered = sourceIndex;
        break;
      case PM::Smoothed:
        assert(shareSource != PM::Jacobian);
        self.ismoothed = sourceIndex;
        break;
      case PM::Jacobian:
        assert(shareSource == PM::Jacobian);
        self.ijacobian = sourceIndex;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }
  }

  void unset_impl(TrackStatePropMask target, TrackIndexType istate) {
    auto& data = m_collection->at(istate).data();
    switch (target) {
      case TrackStatePropMask::Predicted:
        data.ipredicted = kInvalid;
        break;
      case TrackStatePropMask::Filtered:
        data.ifiltered = kInvalid;
        break;
      case TrackStatePropMask::Smoothed:
        data.ismoothed = kInvalid;
        break;
      case TrackStatePropMask::Jacobian:
        data.ijacobian = kInvalid;
        break;
      case TrackStatePropMask::Calibrated:
        data.measdim = 0;
        break;
      default:
        throw std::domain_error{"Unable to unset this component"};
    }
  }

  void clear_impl() {
    m_collection->clear();
    m_params->clear();
    m_surfaces.clear();
  }

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {}

  void allocateCalibrated_impl(IndexType istate, size_t measdim) {
    assert(measdim > 0 && "Zero measdim not supported");
    auto& data = m_collection->at(istate).data();
    data.measdim = measdim;
  }

  void setUncalibratedSourceLink_impl(IndexType istate, SourceLink sourceLink) {
    PodioUtil::Identifier id =
        m_helper.get().sourceLinkToIdentifier(std::move(sourceLink));
    m_collection->at(istate).data().uncalibratedIdentifier = id;
  }

  void setReferenceSurface_impl(IndexType istate,
                                std::shared_ptr<const Surface> surface) {
    auto trackState = m_collection->at(istate);
    trackState.setReferenceSurface(
        PodioUtil::convertSurfaceToPodio(m_helper, *surface));
    m_surfaces.at(istate) = std::move(surface);
  }

  MultiTrajectoryTraits::IndexType calibratedSize_impl(IndexType istate) const {
    return m_collection.at(istate).getData().measdim;
  }

  SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return m_helper.get().identifierToSourceLink(
        m_collection.at(istate).getData().uncalibratedIdentifier);
  }

  const Surface* referenceSurface_impl(IndexType istate) const {
    return m_surfaces.at(istate).get();
  }

 private:
  friend class PodioTrackStateContainerBase;
  friend class ConstPodioTrackStateContainer;

  std::reference_wrapper<PodioUtil::ConversionHelper> m_helper;
  ActsPodioEdm::TrackStateCollection m_collection;
  ActsPodioEdm::BoundParametersCollection m_params;
  ActsPodioEdm::JacobianCollection m_jacs;
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};

static_assert(
    !IsReadOnlyMultiTrajectory<MutablePodioTrackStateContainer>::value,
    "MutablePodioTrackStateContainer should not be read-only");

static_assert(!MutablePodioTrackStateContainer::ReadOnly,
              "MutablePodioTrackStateContainer should not be read-only");

ACTS_STATIC_CHECK_CONCEPT(MutableMultiTrajectoryBackend,
                          MutablePodioTrackStateContainer);

ConstPodioTrackStateContainer::ConstPodioTrackStateContainer(
    MutablePodioTrackStateContainer&& other)
    : m_helper{other.m_helper},
      m_collection{std::move(other.m_collection)},
      m_params{std::move(other.m_params)},
      m_jacs{std::move(other.m_jacs)},
      m_surfaces{std::move(other.m_surfaces)} {}

ConstPodioTrackStateContainer::ConstPodioTrackStateContainer(
    const MutablePodioTrackStateContainer& other)
    : m_helper{other.m_helper},
      m_surfaces{other.m_surfaces.begin(), other.m_surfaces.end()} {
  for (auto src : other.m_collection) {
    auto dst = m_collection.create();
    dst = src.clone();
  }
  for (auto src : other.m_params) {
    auto dst = m_params.create();
    dst = src.clone();
  }
  for (auto src : other.m_jacs) {
    auto dst = m_jacs.create();
    dst = src.clone();
  }
}

}  // namespace Acts
