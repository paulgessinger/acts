// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <exception>

#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>

#include "spdlog/common.h"
#include "spdlog/pattern_formatter.h"

namespace Acts {

namespace Logging {

#if defined(ACTS_ENABLE_LOG_FAILURE_THRESHOLD) and \
    not defined(ACTS_LOG_FAILURE_THRESHOLD)
namespace {
Level& getFailureThresholdMutable() {
  static Level _level = []() {
    Level level = Level::MAX;

    const char* envvar = std::getenv("ACTS_LOG_FAILURE_THRESHOLD");
    if (envvar == nullptr) {
      return level;
    }
    std::string slevel = envvar;
    if (slevel == "VERBOSE") {
      level = std::min(level, Level::VERBOSE);
    } else if (slevel == "DEBUG") {
      level = std::min(level, Level::DEBUG);
    } else if (slevel == "INFO") {
      level = std::min(level, Level::INFO);
    } else if (slevel == "WARNING") {
      level = std::min(level, Level::WARNING);
    } else if (slevel == "ERROR") {
      level = std::min(level, Level::ERROR);
    } else if (slevel == "FATAL") {
      level = std::min(level, Level::FATAL);
    } else {
      std::cerr << "ACTS_LOG_FAILURE_THRESHOLD environment variable is set to "
                   "unknown value: "
                << slevel << std::endl;
    }
    return level;
  }();

  return _level;
}
}  // namespace

Level getFailureThreshold() {
  return getFailureThresholdMutable();
}

void setFailureThreshold(Level level) {
  getFailureThresholdMutable() = level;
}

#else

void setFailureThreshold(Level /*lvl*/) {
  throw std::logic_error{
      "Compile-time log failure threshold defined (ACTS_LOG_FAILURE_THRESHOLD "
      "is set or ACTS_ENABLE_LOG_FAILURE_THRESHOLD is OFF), unable to "
      "override. See "
      "https://acts.readthedocs.io/en/latest/core/"
      "logging.html#logging-thresholds"};
}

#endif

}  // namespace Logging

namespace {

class CustomLevelFormatter : public spdlog::custom_flag_formatter {
 public:
  void format(const spdlog::details::log_msg& msg, const std::tm&,
              spdlog::memory_buf_t& dest) override {
    const static std::array<std::string_view, spdlog::level::n_levels> levels =
        {"VERBOSE  ", "DEBUG    ", "INFO     ", "WARNING  ",
         "ERROR    ", "FATAL    ", "MAX      "};

    const std::string_view& levelString = levels[msg.level];
    dest.append(levelString.data(), levelString.data() + levelString.size());
  }

  std::unique_ptr<spdlog::custom_flag_formatter> clone() const override {
    return spdlog::details::make_unique<CustomLevelFormatter>();
  }
};

}  // namespace

Logger::Logger(Acts::Logging::Level level, const std::string& name,
               std::ostream& os) {
  auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(os);
  m_logger = std::make_shared<spdlog::logger>(name, std::move(sink));
  m_logger->set_level(toSpdlog(level));

  auto formatter = std::make_unique<spdlog::pattern_formatter>();
  formatter->add_flag<CustomLevelFormatter>('Z').set_pattern("%-29!n %Z %v");
  m_logger->set_formatter(std::move(formatter));
}

Logger::Logger(std::shared_ptr<spdlog::logger> logger)
    : m_logger{std::move(logger)} {}

Logger::~Logger() = default;

bool Logger::doPrint(const Logging::Level& lvl) const {
  return m_logger->should_log(toSpdlog(lvl));
}

void Logger::log(Logging::Level lvl, const std::string& input) const {
  m_logger->log(toSpdlog(lvl), input);
}

Logging::Level Logger::level() const {
  return fromSpdlog(m_logger->level());
}

const std::string& Logger::name() const {
  return m_logger->name();
}

std::unique_ptr<Logger> Logger::clone(
    const std::optional<std::string>& _name,
    const std::optional<Logging::Level>& _level) const {
  auto newLogger = m_logger->clone(_name.value_or(name()));
  if (_level.has_value()) {
    newLogger->set_level(toSpdlog(_level.value()));
  }
  return std::make_unique<Logger>(std::move(newLogger));
}

std::unique_ptr<Logger> Logger::clone(Logging::Level _level) const {
  return clone(std::nullopt, _level);
}

std::unique_ptr<Logger> Logger::cloneWithSuffix(
    const std::string& suffix, std::optional<Logging::Level> _level) const {
  return clone(name() + suffix, _level);
}

std::unique_ptr<const Logger> getDefaultLogger(const std::string& name,
                                               const Logging::Level& lvl,
                                               std::ostream* log_stream) {
  using namespace Logging;
  return std::make_unique<const Logger>(lvl, name, *log_stream);
}

const Logger& getDummyLogger() {
  static const std::unique_ptr<const Logger> logger =
      std::make_unique<const Logger>(Logging::Level::MAX, "Dummy");
  return *logger;
}

spdlog::level::level_enum Logger::toSpdlog(Logging::Level level) {
  switch (level) {
    case Logging::VERBOSE:
      return spdlog::level::level_enum::trace;
    case Logging::DEBUG:
      return spdlog::level::level_enum::debug;
    case Logging::INFO:
      return spdlog::level::level_enum::info;
    case Logging::WARNING:
      return spdlog::level::level_enum::warn;
    case Logging::ERROR:
      return spdlog::level::level_enum::err;
    case Logging::FATAL:
      return spdlog::level::level_enum::critical;
    case Logging::MAX:
      return spdlog::level::level_enum::off;
    default:
      std::terminate();
  }
}

Logging::Level Logger::fromSpdlog(spdlog::level::level_enum level) {
  switch (level) {
    case spdlog::level::level_enum::trace:
      return Logging::VERBOSE;
    case spdlog::level::level_enum::debug:
      return Logging::DEBUG;
    case spdlog::level::level_enum::info:
      return Logging::INFO;
    case spdlog::level::level_enum::warn:
      return Logging::WARNING;
    case spdlog::level::level_enum::err:
      return Logging::ERROR;
    case spdlog::level::level_enum::critical:
      return Logging::FATAL;
    case spdlog::level::level_enum::off:
      return Logging::MAX;
    default:
      std::terminate();
  }
}

}  // namespace Acts
