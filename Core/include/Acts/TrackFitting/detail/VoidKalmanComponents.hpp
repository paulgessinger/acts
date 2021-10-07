// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

// /// @brief void Measurement calibrator and converter
// struct VoidKalmanComponents {
//   /// @brief Public call mimicking a calibrator
//   ///
//   /// @tparam measurement_t Type of the measurement
//   /// @tparam parameter_t Type of the parameters for calibration
//   ///
//   /// @param measurement Measurement to be moved through
//   /// @param parameters Parameters to be used for calibration
//   ///
//   /// @return void-calibrated measurement
//   template <typename measurement_t, typename parameters_t>
//   Result<measurement_t> operator()(measurement_t measurement,
//                                    const parameters_t& parameters) const {
//     (void)parameters;
//     return measurement;
//   }
// };

struct VoidKalmanCalibrator {
  template <typename source_link_t>
  void operator()(const GeometryContext& gctx,
                  typename MultiTrajectory<source_link_t>::TrackStateProxy
                      trackState) const {
    throw std::runtime_error{"VoidKalmanCalibrator should not ever execute"};
  }
};

template <typename source_link_t>
void voidKalmanCalibrator(
    const GeometryContext& gctx,
    typename MultiTrajectory<source_link_t>::TrackStateProxy trackState) {
  throw std::runtime_error{"VoidKalmanCalibrator should not ever execute"};
}

/// @brief void Kalman updater
struct VoidKalmanUpdater {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam track_state_t Type of the track state to be used
  /// @tparam predicted_state_t Type of the (bound) predicted state
  ///
  /// @param trackState The track state
  /// @param predicted The predicted parameters
  ///
  /// @return The copied predicted parameters
  template <typename source_link_t>
  Result<void> operator()(
      const GeometryContext&,
      typename MultiTrajectory<source_link_t>::TrackStateProxy trackState,
      NavigationDirection direction,
      LoggerWrapper logger = getDummyLogger()) const {
    trackState.filtered() = trackState.predicted();
    trackState.filteredCovariance() = trackState.predictedCovariance();
    return Result<void>::success();
  }
};

template <typename source_link_t>
Result<void> voidKalmanUpdater(
    const GeometryContext&,
    typename MultiTrajectory<source_link_t>::TrackStateProxy trackState,
    NavigationDirection direction, LoggerWrapper logger = getDummyLogger()) {
  trackState.filtered() = trackState.predicted();
  trackState.filteredCovariance() = trackState.predictedCovariance();
  return Result<void>::success();
}

/// @brief void Kalman smoother
struct VoidKalmanSmoother {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam track_states_t Type of the track states
  ///
  /// @param trackStates The track states to be smoothed
  ///
  /// @return The resulting
  template <typename source_link_t>
  Result<void> operator()(const GeometryContext& gctx,
                          MultiTrajectory<source_link_t>& trajectory,
                          size_t entryIndex,
                          LoggerWrapper logger = getDummyLogger()) const {
    // trackState.filtered() = trackState.predicted();
    // trackState.filteredCovariance() = trackState.predictedCovariance();
    return Result<void>::success();
  }
};

template <typename source_link_t>
Result<void> voidKalmanSmoother(const GeometryContext& gctx,
                                MultiTrajectory<source_link_t>& trajectory,
                                size_t entryIndex,
                                LoggerWrapper logger = getDummyLogger()) {
  return Result<void>::success();
}

/// @brief void outlier finder
struct VoidOutlierFinder {
  /// @brief Public call mimicking an outlier finder
  ///
  /// @tparam track_state_t Type of the track state
  ///
  /// @param trackState The trackState to investigate
  ///
  /// @return Whether it's outlier or not
  template <typename source_link_t, size_t kMeasurementSizeMax>
  constexpr bool operator()(
      detail_lt::TrackStateProxy<source_link_t, kMeasurementSizeMax, false>
          trackState) const {
    (void)trackState;
    return false;
  }
};

template <typename source_link_t, size_t kMeasurementSizeMax>
bool voidOutlierFinder(
    detail_lt::TrackStateProxy<source_link_t, kMeasurementSizeMax, false>
        trackState) {
  (void)trackState;
  return false;
}

/// @brief void smoothing logic
struct VoidReverseFilteringLogic {
  /// @brief Public call mimicking an outlier finder
  ///
  /// @tparam track_state_t Type of the track state
  ///
  /// @param trackState The trackState of the last measurement
  ///
  /// @return Whether to run filtering in reversed direction as smoothing or not
  template <typename source_link_t, size_t kMeasurementSizeMax>
  constexpr bool operator()(
      detail_lt::TrackStateProxy<source_link_t, kMeasurementSizeMax, false>
          trackState) const {
    (void)trackState;
    return false;
  }
};

template <typename source_link_t, size_t kMeasurementSizeMax>
constexpr bool voidReverseFilteringLogic(
    detail_lt::TrackStateProxy<source_link_t, kMeasurementSizeMax, false>
        trackState) {
  (void)trackState;
  return false;
}

}  // namespace Acts
