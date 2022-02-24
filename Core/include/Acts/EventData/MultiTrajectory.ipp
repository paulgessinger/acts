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
template <size_t M, bool ReadOnly>
inline TrackStateProxy<M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory, ReadOnly>& trajectory, size_t istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::parameters() const -> Parameters {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
    return Parameters(m_traj->m_smth.data.col(idx).data());
  } else if (hasFiltered()) {
    idx = data().ifiltered;
    return Parameters(m_traj->m_filt.data.col(idx).data());
  } else {
    idx = data().ipredicted;
    return Parameters(m_traj->m_pred.data.col(idx).data());
  }
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::covariance() const -> Covariance {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
    return Covariance(m_traj->m_smthCov.data.col(idx).data());
  } else if (hasFiltered()) {
    idx = data().ifiltered;
    return Covariance(m_traj->m_filtCov.data.col(idx).data());
  } else {
    idx = data().ipredicted;
    return Covariance(m_traj->m_predCov.data.col(idx).data());
  }
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::predicted() const -> Parameters {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Predicted));
  return Parameters(m_traj->m_pred.col(data().ipredicted).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::predictedCovariance() const
    -> Covariance {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Predicted));
  return Covariance(m_traj->m_predCov.col(data().ipredicted).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::filtered() const -> Parameters {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Filtered));
  return Parameters(m_traj->m_filt.col(data().ifiltered).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::filteredCovariance() const
    -> Covariance {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Filtered));
  return Covariance(m_traj->m_filtCov.col(data().ifiltered).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::smoothed() const -> Parameters {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Smoothed));
  return Parameters(m_traj->m_smth.col(data().ismoothed).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::smoothedCovariance() const
    -> Covariance {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Smoothed));
  return Covariance(m_traj->m_smthCov.col(data().ismoothed).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::jacobian() const -> Covariance {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Jacobian));
  return Covariance(m_traj->m_jac.col(data().ijacobian).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::projector() const -> Projector {
  return bitsetToMatrix<Projector>(m_traj->m_projectors[data().iprojector]);
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::uncalibrated() const
    -> const SourceLink& {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Uncalibrated));
  assert(m_traj->m_sourceLinks[data().iuncalibrated] != nullptr);
  return *m_traj->m_sourceLinks[data().iuncalibrated];
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibrated() const -> Measurement {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Calibrated));
  return Measurement(m_traj->m_meas.col(data().icalibrated).data());
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibratedSourceLink() const
    -> const SourceLink& {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Calibrated));
  assert(m_traj->m_sourceLinks[data().icalibratedsourcelink] != nullptr);
  return *m_traj->m_sourceLinks[data().icalibratedsourcelink];
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibratedCovariance() const
    -> MeasurementCovariance {
  assert(ACTS_CHECK_BIT(data().mask, TrackStatePropMask::Calibrated));
  return MeasurementCovariance(
      m_traj->m_measCov.col(data().icalibrated).data());
}

}  // namespace detail_lt

inline size_t MultiTrajectory::addTrackState(TrackStatePropMask mask,
                                             size_t iprevious) {
  m_index.emplace_back();
  detail_lt::IndexData& p = m_index.back();
  p.mask = mask;

  size_t index = m_index.size() - 1;

  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
  }

  // always set, but can be null
  m_referenceSurfaces.emplace_back(nullptr);
  p.irefsurface = m_referenceSurfaces.size() - 1;

  m_pred.addCol();
  m_predCov.addCol();
  p.ipredicted = m_pred.size() - 1;

  m_filt.addCol();
  m_filtCov.addCol();
  p.ifiltered = m_filt.size() - 1;

  m_smth.addCol();
  m_smthCov.addCol();
  p.ismoothed = m_smth.size() - 1;

  m_jac.addCol();
  p.ijacobian = m_jac.size() - 1;

  m_sourceLinks.emplace_back();
  p.iuncalibrated = m_sourceLinks.size() - 1;

  m_meas.addCol();
  m_measCov.addCol();
  p.icalibrated = m_meas.size() - 1;

  m_sourceLinks.emplace_back();
  p.icalibratedsourcelink = m_sourceLinks.size() - 1;

  m_projectors.emplace_back();
  p.iprojector = m_projectors.size() - 1;

  return index;
}

template <typename F>
void MultiTrajectory::visitBackwards(size_t iendpoint, F&& callable) const {
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, ConstTrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid ||
          !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

template <typename F>
void MultiTrajectory::applyBackwards(size_t iendpoint, F&& callable) {
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid ||
          !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}
}  // namespace Acts
