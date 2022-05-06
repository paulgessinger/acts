// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackStatePropMask.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <ostream>

namespace Acts {

std::ostream& operator<<(std::ostream& os, TrackStatePropMask mask) {
  using PM = TrackStatePropMask;
  auto check = [mask](auto attrib) -> const char* {
    if (ACTS_CHECK_BIT(mask, attrib)) {
      return "x";
    }

    return " ";
  };
  os << "TrackStatePropMask(";
  if (mask == PM::None) {
    os << "None";
  } else {
    os << "\n  [" << check(PM::Predicted) << "] predicted";
    os << "\n  [" << check(PM::Filtered) << "] filtered";
    os << "\n  [" << check(PM::Smoothed) << "] smoothed";
    os << "\n  [" << check(PM::Jacobian) << "] jacobian";
    os << "\n  [" << check(PM::Uncalibrated) << "] uncalibrated";
    os << "\n  [" << check(PM::Calibrated) << "] calibrated";
    os << "\n";
  }
  os << ")";
  return os;
}

}  // namespace Acts
