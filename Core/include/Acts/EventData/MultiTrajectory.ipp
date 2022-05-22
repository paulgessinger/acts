// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

namespace Acts {
namespace detail_lt {
template <typename D, size_t M, bool ReadOnly>
inline TrackStateProxy<D, M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory<D>, ReadOnly>& trajectory, size_t istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <typename D, size_t M, bool ReadOnly>
TrackStatePropMask TrackStateProxy<D, M, ReadOnly>::getMask() const {
  using PM = TrackStatePropMask;

  PM mask = PM::None;
  if (hasPredicted()) {
    mask |= PM::Predicted;
  }
  if (hasFiltered()) {
    mask |= PM::Filtered;
  }
  if (hasSmoothed()) {
    mask |= PM::Smoothed;
  }
  if (hasJacobian()) {
    mask |= PM::Jacobian;
  }
  if (hasCalibrated()) {
    mask |= PM::Calibrated;
  }
  return mask;
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::parameters() const -> Parameters {
  IndexType idx;
  if (hasSmoothed()) {
    return smoothed();
  } else if (hasFiltered()) {
    return filtered();
  } else {
    return predicted();
  }

  return m_traj->parameters(idx);
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::covariance() const -> Covariance {
  IndexType idx;
  if (hasSmoothed()) {
    return smoothed();
  } else if (hasFiltered()) {
    return filtered();
  } else {
    return predicted();
  }
  return m_traj->covariance(idx);
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::predicted() const -> Parameters {
  assert(has<hashString("predicted")>());
  return m_traj->parameters(component<IndexType, hashString("predicted")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::predictedCovariance() const
    -> Covariance {
  assert(has<hashString("predicted")>());
  return m_traj->covariance(component<IndexType, hashString("predicted")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::filtered() const -> Parameters {
  assert(has<hashString("filtered")>());
  return m_traj->parameters(component<IndexType, hashString("filtered")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::filteredCovariance() const
    -> Covariance {
  assert(has<hashString("filtered")>());
  return m_traj->covariance(component<IndexType, hashString("filtered")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::smoothed() const -> Parameters {
  assert(has<hashString("smoothed")>());
  return m_traj->parameters(component<IndexType, hashString("smoothed")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::smoothedCovariance() const
    -> Covariance {
  assert(has<hashString("smoothed")>());
  return m_traj->covariance(component<IndexType, hashString("smoothed")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::jacobian() const -> Covariance {
  assert(has<hashString("jacobian")>());
  return m_traj->jacobian(component<IndexType, hashString("jacobian")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::projector() const -> Projector {
  assert(has<hashString("projector")>());
  return bitsetToMatrix<Projector>(
      component<ProjectorBitset, hashString("projector")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::uncalibrated() const
    -> const SourceLink& {
  assert(has<hashString("sourceLink")>());
  using T = const SourceLink*;
  const T& sl = component<const SourceLink*, hashString("sourceLink")>();
  assert(sl != nullptr);
  return *sl;
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::calibrated() const -> Measurement {
  assert(has<hashString("calibrated")>());
  return m_traj->measurement(component<IndexType, hashString("calibrated")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::calibratedSourceLink() const
    -> const SourceLink& {
  assert(has<hashString("calibratedSourceLink")>());
  using T = const SourceLink*;
  const T& sl =
      component<const SourceLink*, hashString("calibratedSourceLink")>();
  assert(sl != nullptr);
  return *sl;
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::calibratedCovariance() const
    -> MeasurementCovariance {
  assert(has<hashString("calibrated")>());
  return m_traj->measurementCovariance(
      component<IndexType, hashString("calibrated")>());
}

}  // namespace detail_lt

template <typename D>
template <typename F>
void MultiTrajectory<D>::visitBackwards(size_t iendpoint, F&& callable) const {
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    auto ts = getTrackState(iendpoint);
    if constexpr (std::is_same_v<std::invoke_result_t<F, ConstTrackStateProxy>,
                                 bool>) {
      bool proceed = callable(ts);
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (!ts.hasPrevious() || !proceed) {
        break;
      }
    } else {
      callable(ts);
      // this point has no parent and ends the trajectory
      if (!ts.hasPrevious()) {
        break;
      }
    }
    iendpoint = ts.previous();
  }
}

template <typename D>
template <typename F>
void MultiTrajectory<D>::applyBackwards(size_t iendpoint, F&& callable) {
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    auto ts = getTrackState(iendpoint);
    if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                 bool>) {
      bool proceed = callable(ts);
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (!ts.hasPrevious() || !proceed) {
        break;
      }
    } else {
      callable(ts);
      // this point has no parent and ends the trajectory
      if (!ts.hasPrevious()) {
        break;
      }
    }
    iendpoint = ts.previous();
  }
}
}  // namespace Acts
