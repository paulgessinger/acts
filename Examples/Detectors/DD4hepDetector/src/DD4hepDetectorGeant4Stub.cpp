// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

namespace ActsExamples {

std::shared_ptr<Geant4DetectorConstructionFactory>
DD4hepDetector::buildGeant4DetectorConstructionFactory() {
  throw std::runtime_error("Not implemented");
}

}  // namespace ActsExamples