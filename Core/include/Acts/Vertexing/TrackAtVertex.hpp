// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

#include <any>
#include <cstring>
#include <functional>
#include <memory>
#include <typeindex>

namespace Acts {

/// Type-erased wrapper around an input track pointer.
struct InputTrack {
  /// Construct from input track pointer
  /// @param inputTrack Pointer to the input track
  template <typename input_track_t>
  explicit InputTrack(const input_track_t* inputTrack)
      : m_type{typeid(inputTrack)}, m_ptr{inputTrack} {}

  InputTrack() = delete;
  /// Copy constructor
  InputTrack(const InputTrack&) = default;
  /// Move constructor
  InputTrack(InputTrack&&) = default;
  /// Copy assignment
  /// @return Reference to this object
  InputTrack& operator=(const InputTrack&) = default;
  /// Move assignment
  /// @return Reference to this object
  InputTrack& operator=(InputTrack&&) = default;

  /// Equality comparison with another InputTrack
  /// @param other The other InputTrack to compare with
  /// @return True if equal
  bool operator==(const InputTrack& other) const {
    return m_ptr == other.m_ptr;
  }

  /// Equality comparison with a raw pointer
  /// @param other The pointer to compare with
  /// @return True if equal
  template <typename input_track_t>
  bool operator==(const input_track_t* other) const {
    return m_ptr == other;
  }

  /// Comparison operator for ordering
  /// @param other The other InputTrack to compare with
  /// @return True if this is less than other
  bool operator<(const InputTrack& other) const { return m_ptr < other.m_ptr; }

  /// Cast to the original pointer type
  /// @return Pointer to the original track object
  template <typename T>
  const T* as() const {
    using ptr_t = const T*;
    if (m_type != typeid(ptr_t)) {
      throw std::bad_any_cast();
    }
    return static_cast<ptr_t>(m_ptr);
  }

  friend std::ostream& operator<<(std::ostream& os, const InputTrack& track) {
    os << track.m_ptr;
    return os;
  }

  friend struct std::hash<Acts::InputTrack>;

  /// Extract track parameters from an InputTrack
  /// @param track The input track
  /// @return The bound track parameters
  static BoundTrackParameters extractParameters(const InputTrack& track) {
    return *track.as<BoundTrackParameters>();
  }

  /// Function type for extracting track parameters
  using Extractor = Acts::Delegate<BoundTrackParameters(const InputTrack&)>;

 private:
  std::type_index m_type;
  const void* m_ptr;
};

/// @class TrackAtVertex
///
/// @brief Defines a track at vertex object
struct TrackAtVertex {
  /// Deleted default constructor
  TrackAtVertex() = delete;

  /// @brief Parameterized constructor
  ///
  /// @param chi2perTrack Chi2 of track
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(double chi2perTrack, const BoundTrackParameters& paramsAtVertex,
                InputTrack originalTrack)
      : fittedParams(paramsAtVertex),
        originalParams(originalTrack),
        chi2Track(chi2perTrack) {}

  /// @brief Constructor with default chi2
  ///
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const BoundTrackParameters& paramsAtVertex,
                InputTrack originalTrack)
      : fittedParams(paramsAtVertex), originalParams(originalTrack) {}

  /// Fitted perigee
  BoundTrackParameters fittedParams;

  /// Original input parameters
  InputTrack originalParams;

  /// Chi2 of track
  double chi2Track = 0;

  /// Number degrees of freedom
  /// Note: Can be different from integer value
  /// since annealing can result in effective
  /// non-interger values
  double ndf = 0;

  /// Value of the compatibility of the track to the actual vertex, based
  /// on the estimation of the 3d distance between the track and the vertex
  double vertexCompatibility = 0;

  /// Weight of track in fit
  double trackWeight = 1;

  /// The linearized state of the track at vertex
  LinearizedTrack linearizedState;

  /// Is already linearized
  bool isLinearized = false;

  /// Return the inverse covariance block for the current linearized track.
  /// Mutating the covariance or switching fit dimension invalidates the cache.
  /// Cache storage is allocated on first use; copies share it until mutation.
  /// @tparam nParams Number of fitted parameters (5 spatial or 6 spacetime)
  template <unsigned int nParams>
  SquareMatrix<nParams> parameterWeight() {
    static_assert(nParams == 5 || nParams == 6);
    const auto& covariance = linearizedState.covarianceAtPCA;
    if (!m_weightCache || m_weightCache->dimension != nParams ||
        std::memcmp(covariance.data(), m_weightCache->covariance.data(),
                    eBoundSize * eBoundSize * sizeof(double)) != 0) {
      if (!m_weightCache || !m_weightCache.unique()) {
        m_weightCache = std::make_shared<WeightCache>();
      }
      const SquareMatrix<nParams> weight =
          covariance.template block<nParams, nParams>(0, 0).inverse();
      m_weightCache->covariance = covariance;
      m_weightCache->weight.template block<nParams, nParams>(0, 0) = weight;
      m_weightCache->dimension = nParams;
    }
    return m_weightCache->weight.template block<nParams, nParams>(0, 0);
  }

 private:
  struct WeightCache {
    BoundMatrix covariance = BoundMatrix::Zero();
    BoundMatrix weight = BoundMatrix::Zero();
    unsigned int dimension = 0;
  };
  std::shared_ptr<WeightCache> m_weightCache;
};

}  // namespace Acts

template <>
struct std::hash<Acts::InputTrack> {
  std::size_t operator()(const Acts::InputTrack& track) const noexcept {
    return std::hash<const void*>{}(track.m_ptr);
  }
};
