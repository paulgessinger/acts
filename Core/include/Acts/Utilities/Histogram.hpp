// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <boost/histogram.hpp>
#include <boost/histogram/accumulators/weighted_sum.hpp>
#include <boost/histogram/axis/category.hpp>
#include <boost/histogram/axis/integer.hpp>
#include <boost/histogram/axis/regular.hpp>
#include <boost/histogram/axis/variable.hpp>

namespace Acts {

// Convenience alias to avoid boost:: qualification when filling with weights
using boost::histogram::weight;

/// A histogram class wrapping boost::histogram with ACTS-specific features
///
/// This class provides a thin wrapper around boost::histogram with dynamic
/// axes, adding histogram-level metadata (title, label) and following ACTS
/// conventions. It uses weight storage for variance tracking and supports
/// various axis types through a variant.
///
/// The class uses composition rather than inheritance to wrap boost::histogram,
/// allowing for better encapsulation and ACTS-specific extensions while
/// leveraging boost::histogram's type erasure capabilities.
///
/// Example usage:
/// @code
/// using namespace Acts;
/// using namespace boost::histogram;
///
/// // Create a 2D histogram
/// Histogram::AxesType axes;
/// axes.emplace_back(axis::regular<>(50, 0.0, 10.0));
/// axes.emplace_back(axis::category<std::string>({"A", "B", "C"}));
///
/// Histogram hist(std::move(axes), "My Histogram", "hist_label");
///
/// // Fill with weights
/// hist(3.5, "A", weight(2.0));
/// hist(7.2, "B");  // weight defaults to 1.0
///
/// // Access statistics
/// std::cout << "Sum: " << hist.sum() << "\n";
/// std::cout << "Variance: " << hist.variance() << "\n";
///
/// // Iterate over bins
/// for (auto&& bin : hist.indexed()) {
///   std::cout << "Bin value: " << bin->value()
///             << " ± " << std::sqrt(bin->variance()) << "\n";
/// }
/// @endcode
class Histogram {
 public:
  // ==================== Type Aliases ====================

  /// Axis variant supporting multiple axis types for dynamic histograms
  using AxisVariant = boost::histogram::axis::variant<
      boost::histogram::axis::regular<>,            // Equidistant bins
      boost::histogram::axis::variable<>,           // Variable-width bins
      boost::histogram::axis::integer<>,            // Integer axis
      boost::histogram::axis::category<std::string>,  // String categories
      boost::histogram::axis::category<>            // Integer categories
      >;

  /// Dynamic axes container type
  using AxesType = std::vector<AxisVariant>;

  /// Storage type with weight and variance tracking
  using StorageType = boost::histogram::storage_adaptor<
      std::vector<boost::histogram::accumulators::weighted_sum<>>>;

  /// Underlying histogram type with dynamic axes and weight storage
  using HistogramType = boost::histogram::histogram<AxesType, StorageType>;

  // ==================== Constructors ====================

  /// Default constructor creates an empty histogram
  Histogram() = default;

  /// Construct from axes with optional title and label
  /// @param axes Vector of axis variants defining the histogram structure
  /// @param title Optional histogram title for documentation/output
  /// @param label Optional histogram label for identification
  Histogram(AxesType axes, std::string title = "", std::string label = "")
      : m_hist(std::move(axes), StorageType{}),
        m_title(std::move(title)),
        m_label(std::move(label)) {}

  /// Construct directly from a boost histogram
  /// @param hist The boost histogram to wrap
  /// @param title Optional histogram title
  /// @param label Optional histogram label
  explicit Histogram(HistogramType hist, std::string title = "",
                     std::string label = "")
      : m_hist(std::move(hist)),
        m_title(std::move(title)),
        m_label(std::move(label)) {}

  // ==================== Accessors ====================

  /// Access the underlying boost histogram (const)
  /// @return Const reference to the internal histogram
  const HistogramType& histogram() const { return m_hist; }

  /// Access the underlying boost histogram (mutable)
  /// @return Mutable reference to the internal histogram
  HistogramType& histogram() { return m_hist; }

  /// Get histogram title
  /// @return The histogram title
  const std::string& title() const { return m_title; }

  /// Set histogram title
  /// @param title New title for the histogram
  void setTitle(std::string title) { m_title = std::move(title); }

  /// Get histogram label
  /// @return The histogram label
  const std::string& label() const { return m_label; }

  /// Set histogram label
  /// @param label New label for the histogram
  void setLabel(std::string label) { m_label = std::move(label); }

  /// Get number of dimensions
  /// @return Number of axes in the histogram
  std::size_t rank() const { return m_hist.rank(); }

  /// Get total number of bins (including under/overflow)
  /// @return Total number of bins
  std::size_t size() const { return m_hist.size(); }

  // ==================== Fill Operations ====================

  /// Fill histogram with optional weight using parameter pack
  ///
  /// This enables the ergonomic syntax: h(x, y, weight(w))
  /// Arguments are forwarded directly to the underlying boost::histogram
  ///
  /// @param args Arguments to pass to boost::histogram (coordinates and
  /// optional weight)
  /// @return Reference to the bin that was filled
  template <typename... Args>
  decltype(auto) operator()(Args&&... args) {
    return m_hist(std::forward<Args>(args)...);
  }

  /// Access histogram bin at specific indices (mutable)
  /// @param indices Bin indices to access
  /// @return Reference to the bin at the given indices
  template <typename... Indices>
  decltype(auto) at(Indices&&... indices) {
    return m_hist.at(std::forward<Indices>(indices)...);
  }

  /// Access histogram bin at specific indices (const)
  /// @param indices Bin indices to access
  /// @return Const reference to the bin at the given indices
  template <typename... Indices>
  decltype(auto) at(Indices&&... indices) const {
    return m_hist.at(std::forward<Indices>(indices)...);
  }

  // ==================== Iteration ====================

  /// Get indexed iterator range for mutable access
  /// @return Indexed range that can be used in range-based for loops
  auto indexed() { return boost::histogram::indexed(m_hist); }

  /// Get indexed iterator range for const access
  /// @return Const indexed range that can be used in range-based for loops
  auto indexed() const { return boost::histogram::indexed(m_hist); }

  // ==================== Output ====================

  /// Write histogram to output stream
  /// @param os Output stream
  /// @param indent Number of spaces to indent (for nested output)
  void toStream(std::ostream& os, std::size_t indent = 0) const;

  /// Stream output operator
  /// @param os Output stream
  /// @param hist Histogram to output
  /// @return Reference to the output stream
  friend std::ostream& operator<<(std::ostream& os, const Histogram& hist) {
    hist.toStream(os);
    return os;
  }

  // ==================== Statistics Operations ====================

  /// Get sum of all weights
  /// @return Sum of weights across all bins
  double sum() const {
    double total = 0.0;
    for (auto&& bin : boost::histogram::indexed(m_hist)) {
      total += bin->value();
    }
    return total;
  }

  /// Get sum of weight variances
  /// @return Sum of variances across all bins
  double variance() const {
    double total = 0.0;
    for (auto&& bin : boost::histogram::indexed(m_hist)) {
      total += bin->variance();
    }
    return total;
  }

  /// Reset all bin contents to zero
  void reset() { m_hist.reset(); }

 private:
  HistogramType m_hist;  ///< Underlying boost histogram
  std::string m_title;   ///< Histogram title
  std::string m_label;   ///< Histogram label
};

}  // namespace Acts
