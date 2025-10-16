// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

/// Fuzzer entry point for testing material interaction calculations.
///
/// This fuzzes computeEnergyLossLandauSigmaQOverP to test for:
/// - Memory safety issues (via AddressSanitizer)
/// - Undefined behavior (via UndefinedBehaviorSanitizer)
/// - Crashes, assertions, or other runtime errors
///
/// The fuzzer generates pseudo-random but deterministic inputs to explore
/// edge cases in the material interaction calculations.
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Need at least 32 bytes to generate meaningful test data
  // 5 floats * 4 bytes + 3 floats * 4 bytes = 32 bytes
  if (size < 32) {
    return 0;
  }

  // Parse fuzzer input into parameters
  float params[8];
  std::memcpy(params, data, sizeof(params));

  // Material parameters
  float X0 = std::abs(params[0]);           // Radiation length
  float L0 = std::abs(params[1]);           // Nuclear interaction length
  float Ar = std::abs(params[2]);           // Relative atomic mass
  float Z = std::abs(params[3]);            // Atomic number
  float rho = std::abs(params[4]);          // Mass density
  // Clamp to reasonable ranges to avoid extreme values
  X0 = std::clamp(X0, 0.1f, 1000.0f);       // mm
  L0 = std::clamp(L0, 0.1f, 1000.0f);       // mm
  Ar = std::clamp(Ar, 1.0f, 300.0f);        // g/mol
  Z = std::clamp(Z, 1.0f, 100.0f);          // atomic number
  rho = std::clamp(rho, 0.001f, 30.0f);     // g/cm³

  // Particle parameters
  float m = std::abs(params[5]);            // Particle mass
  float qOverP = params[6];                 // Charge over momentum
  float absQ = std::abs(params[7]);         // Absolute charge

  // Clamp particle parameters to physically reasonable ranges
  m = std::clamp(m, 0.1f, 1000.0f);         // MeV (electron to muon range)
  absQ = std::clamp(absQ, 0.1f, 10.0f);     // Elementary charges

  // Prevent division by zero and extreme values
  if (std::abs(qOverP) < 1e-6f) {
    qOverP = (qOverP >= 0) ? 1e-6f : -1e-6f;
  }
  qOverP = std::clamp(qOverP, -1.0f, 1.0f); // 1/GeV

  // Thickness from remaining data if available
  float thickness = 1.0f;  // mm
  if (size >= 36) {
    std::memcpy(&thickness, data + 32, sizeof(float));
    thickness = std::abs(thickness);
    thickness = std::clamp(thickness, 0.001f, 100.0f);  // mm
  }

  // Create material and slab
  try {
    Acts::Material material =
        Acts::Material::fromMassDensity(X0, L0, Ar, Z, rho);
    Acts::MaterialSlab slab(material, thickness);

    // Test the target function - computeEnergyLossLandauSigmaQOverP
    // The return value is checked for NaN/Inf which would indicate issues
    float result = Acts::computeEnergyLossLandauSigmaQOverP(slab, m, qOverP, absQ);

    // Check for invalid results (NaN or Inf might indicate problems)
    // Note: We don't crash on these, but they might indicate issues
    if (std::isnan(result) || std::isinf(result)) {
      // Just continue - this might be valid for extreme inputs
    }

    // Also test related functions to increase coverage
    volatile float resultLandau = Acts::computeEnergyLossLandau(slab, m, qOverP, absQ);
    volatile float resultBethe = Acts::computeEnergyLossBethe(slab, m, qOverP, absQ);
    volatile float resultSigma = Acts::computeEnergyLossLandauSigma(slab, m, qOverP, absQ);

    // Use volatile to prevent optimization
    (void)resultLandau;
    (void)resultBethe;
    (void)resultSigma;

  } catch (const std::exception& e) {
    // Expected exceptions are fine - we're looking for crashes and UB
    // The sanitizers will catch real issues like use-after-free, etc.
  } catch (...) {
    // Catch any other exceptions
  }

  return 0;
}
