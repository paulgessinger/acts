// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>

#include <boost/histogram/axis.hpp>

namespace Acts {

void Histogram::toStream(std::ostream& os, std::size_t indent) const {
  std::string prefix(indent, ' ');

  // Print title if set
  if (!m_title.empty()) {
    os << prefix << "Histogram: " << m_title << "\n";
  }

  // Print label if set
  if (!m_label.empty()) {
    os << prefix << "Label: " << m_label << "\n";
  }

  // Print basic histogram information
  os << prefix << "Dimensions: " << m_hist.rank() << "\n";
  os << prefix << "Total bins: " << m_hist.size() << "\n";
  os << prefix << "Sum: " << sum() << "\n";

  // Print axis information
  for (std::size_t i = 0; i < m_hist.rank(); ++i) {
    os << prefix << "Axis " << i << ": ";

    // Use visit to handle the variant axis types
    boost::histogram::axis::visit(
        [&](const auto& ax) {
          using AxisType = std::decay_t<decltype(ax)>;

          os << ax.size() << " bins";

          // For numeric axes, print the range
          if constexpr (requires { ax.bin(0); }) {
            // Check if the axis has bin() method that returns an interval
            if constexpr (requires {
                            ax.bin(0).lower();
                            ax.bin(0).upper();
                          }) {
              // Regular or variable axis with interval bins
              os << " [" << ax.bin(0).lower() << ", "
                 << ax.bin(ax.size() - 1).upper() << "]";
            } else if constexpr (std::is_same_v<
                                     AxisType,
                                     boost::histogram::axis::integer<>>) {
              // Integer axis - print first and last values
              if (ax.size() > 0) {
                os << " [" << ax.value(0) << ", " << ax.value(ax.size()) << ")";
              }
            }
          } else if constexpr (requires { ax.value(0); }) {
            // Category axis - mention it's categorical
            if constexpr (std::is_same_v<AxisType, boost::histogram::axis::
                                                        category<std::string>>) {
              os << " (string categories)";
            } else {
              os << " (integer categories)";
            }
          }
        },
        m_hist.axis(i));

    os << "\n";
  }
}

}  // namespace Acts
