// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include <boost/histogram.hpp>
#include <boost/histogram/axis.hpp>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_SUITE(HistogramSuite)

// Test default construction
BOOST_AUTO_TEST_CASE(DefaultConstruction) {
  Histogram hist;
  BOOST_CHECK_EQUAL(hist.rank(), 0);
  BOOST_CHECK_EQUAL(hist.size(), 0);
  BOOST_CHECK(hist.title().empty());
  BOOST_CHECK(hist.label().empty());
  BOOST_CHECK_EQUAL(hist.sum(), 0.0);
  BOOST_CHECK_EQUAL(hist.variance(), 0.0);
}

// Test construction with axes
BOOST_AUTO_TEST_CASE(ConstructionWithAxes) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));
  axes.emplace_back(axis::regular<>(20, -5.0, 5.0));

  Histogram hist(std::move(axes), "Test Histogram", "test_label");

  BOOST_CHECK_EQUAL(hist.rank(), 2);
  BOOST_CHECK_EQUAL(hist.title(), "Test Histogram");
  BOOST_CHECK_EQUAL(hist.label(), "test_label");
  BOOST_CHECK_GT(hist.size(), 0);
  BOOST_CHECK_EQUAL(hist.sum(), 0.0);
}

// Test construction without title/label
BOOST_AUTO_TEST_CASE(ConstructionWithAxesNoMetadata) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(5, 0.0, 5.0));

  Histogram hist(std::move(axes));

  BOOST_CHECK_EQUAL(hist.rank(), 1);
  BOOST_CHECK(hist.title().empty());
  BOOST_CHECK(hist.label().empty());
}

// Test filling with operator()
BOOST_AUTO_TEST_CASE(FillingBasic) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  // Fill some values
  hist(2.5);
  hist(2.5);
  hist(7.3);

  BOOST_CHECK_EQUAL(hist.sum(), 3.0);
}

// Test filling with weights
BOOST_AUTO_TEST_CASE(FillingWithWeights) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  // Fill with weights
  hist(2.5, weight(2.0));
  hist(2.5, weight(1.5));
  hist(7.3, weight(3.0));

  BOOST_CHECK_EQUAL(hist.sum(), 6.5);

  // Variance should be sum of weight^2
  // 2.0^2 + 1.5^2 + 3.0^2 = 4.0 + 2.25 + 9.0 = 15.25
  BOOST_CHECK_CLOSE(hist.variance(), 15.25, 1e-10);
}

// Test 2D histogram filling
BOOST_AUTO_TEST_CASE(Filling2D) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  hist(2.5, 3.5);
  hist(2.5, 3.5, weight(2.0));
  hist(7.3, 8.2);

  BOOST_CHECK_EQUAL(hist.sum(), 4.0);  // 1 + 2 + 1 = 4
}

// Test regular axis type
BOOST_AUTO_TEST_CASE(RegularAxis) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes), "Regular Axis Test");

  hist(0.5);
  hist(5.5);
  hist(9.5);

  BOOST_CHECK_EQUAL(hist.sum(), 3.0);
  BOOST_CHECK_EQUAL(hist.rank(), 1);
}

// Test variable axis type
BOOST_AUTO_TEST_CASE(VariableAxis) {
  using namespace boost::histogram;

  std::vector<double> edges = {0.0, 1.0, 3.0, 6.0, 10.0};

  Histogram::AxesType axes;
  axes.emplace_back(axis::variable<>(edges));

  Histogram hist(std::move(axes), "Variable Axis Test");

  hist(0.5);
  hist(2.0);
  hist(7.0);

  BOOST_CHECK_EQUAL(hist.sum(), 3.0);
  BOOST_CHECK_EQUAL(hist.rank(), 1);
}

// Test integer axis type
BOOST_AUTO_TEST_CASE(IntegerAxis) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::integer<>(0, 10));  // [0, 10)

  Histogram hist(std::move(axes), "Integer Axis Test");

  hist(0);
  hist(5);
  hist(9);

  BOOST_CHECK_EQUAL(hist.sum(), 3.0);
  BOOST_CHECK_EQUAL(hist.rank(), 1);
}

// Test category axis with strings
BOOST_AUTO_TEST_CASE(StringCategoryAxis) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::category<std::string>({"A", "B", "C"}));

  Histogram hist(std::move(axes), "String Category Test");

  hist("A");
  hist("A");
  hist("B");
  hist("C", weight(2.0));

  BOOST_CHECK_EQUAL(hist.sum(), 5.0);  // 1 + 1 + 1 + 2 = 5
  BOOST_CHECK_EQUAL(hist.rank(), 1);
}

// Test category axis with integers
BOOST_AUTO_TEST_CASE(IntegerCategoryAxis) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::category<>({0, 1, 2, 3}));

  Histogram hist(std::move(axes), "Integer Category Test");

  hist(0);
  hist(1);
  hist(1);
  hist(3);

  BOOST_CHECK_EQUAL(hist.sum(), 4.0);
  BOOST_CHECK_EQUAL(hist.rank(), 1);
}

// Test mixed axis types
BOOST_AUTO_TEST_CASE(MixedAxisTypes) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));
  axes.emplace_back(axis::integer<>(0, 5));
  axes.emplace_back(axis::category<std::string>({"X", "Y", "Z"}));

  Histogram hist(std::move(axes), "Mixed Axes Test");

  hist(2.5, 2, "X");
  hist(5.0, 3, "Y");
  hist(8.0, 1, "Z", weight(1.5));

  BOOST_CHECK_EQUAL(hist.sum(), 3.5);  // 1 + 1 + 1.5 = 3.5
  BOOST_CHECK_EQUAL(hist.rank(), 3);
}

// Test title and label accessors
BOOST_AUTO_TEST_CASE(MetadataAccessors) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes), "Initial Title", "initial_label");

  BOOST_CHECK_EQUAL(hist.title(), "Initial Title");
  BOOST_CHECK_EQUAL(hist.label(), "initial_label");

  hist.setTitle("New Title");
  hist.setLabel("new_label");

  BOOST_CHECK_EQUAL(hist.title(), "New Title");
  BOOST_CHECK_EQUAL(hist.label(), "new_label");
}

// Test sum calculation
BOOST_AUTO_TEST_CASE(SumCalculation) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  BOOST_CHECK_EQUAL(hist.sum(), 0.0);

  hist(1.0);
  BOOST_CHECK_EQUAL(hist.sum(), 1.0);

  hist(2.0, weight(2.5));
  BOOST_CHECK_CLOSE(hist.sum(), 3.5, 1e-10);

  hist(3.0);
  hist(4.0);
  BOOST_CHECK_CLOSE(hist.sum(), 5.5, 1e-10);
}

// Test variance calculation
BOOST_AUTO_TEST_CASE(VarianceCalculation) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  BOOST_CHECK_EQUAL(hist.variance(), 0.0);

  hist(1.0);  // weight = 1, variance contribution = 1^2 = 1
  BOOST_CHECK_EQUAL(hist.variance(), 1.0);

  hist(2.0, weight(2.0));  // variance contribution = 2^2 = 4
  BOOST_CHECK_CLOSE(hist.variance(), 5.0, 1e-10);

  hist(3.0, weight(3.0));  // variance contribution = 3^2 = 9
  BOOST_CHECK_CLOSE(hist.variance(), 14.0, 1e-10);
}

// Test reset functionality
BOOST_AUTO_TEST_CASE(ResetFunctionality) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  hist(1.0);
  hist(2.0);
  hist(3.0, weight(2.0));

  BOOST_CHECK_EQUAL(hist.sum(), 4.0);
  BOOST_CHECK_GT(hist.variance(), 0.0);

  hist.reset();

  BOOST_CHECK_EQUAL(hist.sum(), 0.0);
  BOOST_CHECK_EQUAL(hist.variance(), 0.0);
  BOOST_CHECK_EQUAL(hist.rank(), 1);  // Structure unchanged
}

// Test indexed iteration
BOOST_AUTO_TEST_CASE(IndexedIteration) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(5, 0.0, 5.0));

  Histogram hist(std::move(axes));

  // Fill bins 0, 1, 2
  hist(0.5);
  hist(1.5, weight(2.0));
  hist(2.5, weight(3.0));

  // Count non-zero bins
  std::size_t nonZeroCount = 0;
  double totalSum = 0.0;

  for (auto&& bin : hist.indexed()) {
    if (bin->value() > 0.0) {
      nonZeroCount++;
      totalSum += bin->value();
    }
  }

  BOOST_CHECK_EQUAL(nonZeroCount, 3);
  BOOST_CHECK_CLOSE(totalSum, 6.0, 1e-10);
}

// Test at() method for bin access
BOOST_AUTO_TEST_CASE(BinAccessWithAt) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  hist(2.5, weight(5.0));

  // Find which bin contains 2.5
  // For regular axis from 0 to 10 with 10 bins, bin 2 covers [2, 3)
  auto bin = hist.at(2);
  BOOST_CHECK_CLOSE(bin.value(), 5.0, 1e-10);
  BOOST_CHECK_CLOSE(bin.variance(), 25.0, 1e-10);  // 5^2
}

// Test toStream output
BOOST_AUTO_TEST_CASE(ToStreamOutput) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));
  axes.emplace_back(axis::category<std::string>({"A", "B"}));

  Histogram hist(std::move(axes), "Test Histogram", "test_hist");

  hist(2.5, "A");
  hist(7.5, "B", weight(2.0));

  std::ostringstream oss;
  hist.toStream(oss);

  std::string output = oss.str();

  // Check that output contains expected information
  BOOST_CHECK(output.find("Test Histogram") != std::string::npos);
  BOOST_CHECK(output.find("test_hist") != std::string::npos);
  BOOST_CHECK(output.find("Dimensions: 2") != std::string::npos);
  BOOST_CHECK(output.find("Axis 0") != std::string::npos);
  BOOST_CHECK(output.find("Axis 1") != std::string::npos);
}

// Test stream operator
BOOST_AUTO_TEST_CASE(StreamOperator) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(5, 0.0, 5.0));

  Histogram hist(std::move(axes), "Stream Test");

  hist(1.0);

  std::ostringstream oss;
  oss << hist;

  std::string output = oss.str();
  BOOST_CHECK(output.find("Stream Test") != std::string::npos);
  BOOST_CHECK(output.find("Dimensions: 1") != std::string::npos);
}

// Test toStream with indent
BOOST_AUTO_TEST_CASE(ToStreamWithIndent) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(5, 0.0, 5.0));

  Histogram hist(std::move(axes), "Indented Output");

  std::ostringstream oss;
  hist.toStream(oss, 4);

  std::string output = oss.str();

  // Check that lines are indented
  BOOST_CHECK(output.find("    Histogram:") != std::string::npos ||
              output.find("    Dimensions:") != std::string::npos);
}

// Test histogram() accessor
BOOST_AUTO_TEST_CASE(HistogramAccessor) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  // Access underlying histogram
  auto& underlyingHist = hist.histogram();
  BOOST_CHECK_EQUAL(underlyingHist.rank(), 1);

  // Fill through accessor
  underlyingHist(5.0, weight(3.0));
  BOOST_CHECK_EQUAL(hist.sum(), 3.0);

  // Const access
  const Histogram& constHist = hist;
  const auto& constUnderlyingHist = constHist.histogram();
  BOOST_CHECK_EQUAL(constUnderlyingHist.rank(), 1);
}

// Test empty histogram (edge case)
BOOST_AUTO_TEST_CASE(EmptyHistogram) {
  Histogram hist;

  BOOST_CHECK_EQUAL(hist.rank(), 0);
  BOOST_CHECK_EQUAL(hist.size(), 0);
  BOOST_CHECK_EQUAL(hist.sum(), 0.0);
  BOOST_CHECK_EQUAL(hist.variance(), 0.0);

  // Should be able to call toStream without crashing
  std::ostringstream oss;
  hist.toStream(oss);
  BOOST_CHECK(!oss.str().empty());
}

// Test single bin histogram
BOOST_AUTO_TEST_CASE(SingleBinHistogram) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(1, 0.0, 10.0));

  Histogram hist(std::move(axes));

  hist(5.0);
  hist(7.0, weight(2.0));

  BOOST_CHECK_EQUAL(hist.rank(), 1);
  BOOST_CHECK_EQUAL(hist.sum(), 3.0);
}

// Test overflow/underflow handling (boost::histogram handles this automatically)
BOOST_AUTO_TEST_CASE(OverflowUnderflow) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  // Fill values outside the range
  hist(-1.0);   // Underflow
  hist(15.0);   // Overflow
  hist(5.0);    // Normal

  // Only the in-range value is counted by sum() (which uses indexed() that
  // excludes overflow bins by default)
  BOOST_CHECK_EQUAL(hist.sum(), 1.0);

  // But the histogram does store all values including overflow/underflow
  BOOST_CHECK_EQUAL(hist.histogram().size(), 12);  // 10 regular + 2 overflow bins
}

// Test large number of fills
BOOST_AUTO_TEST_CASE(LargeNumberOfFills) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(100, 0.0, 100.0));

  Histogram hist(std::move(axes));

  for (int i = 0; i < 10000; ++i) {
    hist(static_cast<double>(i % 100));
  }

  BOOST_CHECK_EQUAL(hist.sum(), 10000.0);
}

// Test variance with uniform weights
BOOST_AUTO_TEST_CASE(VarianceUniformWeights) {
  using namespace boost::histogram;

  Histogram::AxesType axes;
  axes.emplace_back(axis::regular<>(10, 0.0, 10.0));

  Histogram hist(std::move(axes));

  // Fill with uniform weight
  for (int i = 0; i < 100; ++i) {
    hist(5.0);  // All in same bin, weight = 1
  }

  BOOST_CHECK_EQUAL(hist.sum(), 100.0);
  BOOST_CHECK_EQUAL(hist.variance(), 100.0);  // 100 * 1^2 = 100
}

BOOST_AUTO_TEST_SUITE_END()  // HistogramSuite

BOOST_AUTO_TEST_SUITE_END()  // UtilitiesSuite

}  // namespace Acts::Test
