// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <cassert>
#include <ctime>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

namespace spdlog {
class logger;
}

// clang-format off
/// @brief macro to use a local Acts::Logger object
/// @ingroup Logging
///
/// @param log_object logger instance of type
//         <tt>std::unique_ptr<const Acts::Logger></tt>
///
/// @pre In the current scope, the symbol @c logger is not yet defined.
/// @post The ownership of the given @c log_object is transferred and
///       @c log_object should not be used directly any more.
///
/// This macro allows to use a locally defined logging object with the ACTS_*
/// logging macros. The envisaged usage is the following:
///
/// @code{.cpp}
/// void myFunction() {
///    std::unique_ptr<const Acts::Logger> myLogger
///        = /* .. your initialization .. */;
///    ACTS_LOCAL_LOGGER(std::move(myLogger));
///
///    ACTS_VERBOSE("hello world!");
/// }
/// @endcode
#define ACTS_LOCAL_LOGGER(log_object)                                          \
  std::unique_ptr __log_object_owned = log_object; \
  const Acts::Logger& logger = *__log_object_owned \

// Debug level agnostic implementation of the ACTS_XYZ logging macros
#define ACTS_LOG(level, x)                                                     \
  if (logger().doPrint(level)) {                                               \
    std::ostringstream os;                                                     \
    os << x;                                                                   \
    logger().log(level, os.str());                                             \
  }

/// @brief macro for verbose debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::VERBOSE.
#define ACTS_VERBOSE(x)  ACTS_LOG(Acts::Logging::VERBOSE, x)

/// @brief macro for debug debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::DEBUG.
#define ACTS_DEBUG(x)  ACTS_LOG(Acts::Logging::DEBUG, x)

/// @brief macro for info debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::INFO.
#define ACTS_INFO(x)  ACTS_LOG(Acts::Logging::INFO, x)

/// @brief macro for warning debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::WARNING.
#define ACTS_WARNING(x)  ACTS_LOG(Acts::Logging::WARNING, x)

/// @brief macro for error debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::ERROR.
#define ACTS_ERROR(x)  ACTS_LOG(Acts::Logging::ERROR, x)

/// @brief macro for fatal debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::FATAL.
#define ACTS_FATAL(x)  ACTS_LOG(Acts::Logging::FATAL, x)
// clang-format on

namespace Acts {

/// @brief debug output related helper classes and functions
/// @ingroup Logging
namespace Logging {
/// @brief constants steering the debug output
///
/// All messages with a debug level equal or higher than the currently set
/// debug output level will be printed.
enum Level {
  VERBOSE = 0,  ///< VERBOSE level
  DEBUG,        ///< DEBUG level
  INFO,         ///< INFO level
  WARNING,      ///< WARNING level
  ERROR,        ///< ERROR level
  FATAL,        ///< FATAL level
  MAX           ///< Must be kept above the maximum supported debug level
};

inline std::string_view levelName(Level level) {
  switch (level) {
    case Level::VERBOSE:
      return "VERBOSE";
    case Level::DEBUG:
      return "DEBUG";
    case Level::INFO:
      return "INFO";
    case Level::WARNING:
      return "WARNING";
    case Level::ERROR:
      return "ERROR";
    case Level::FATAL:
      return "FATAL";
    case Level::MAX:
      return "MAX";
    default:
      throw std::invalid_argument{"Unknown level"};
  }
}

#ifdef DOXYGEN
/// @brief Get debug level above which an exception will be thrown after logging
///
/// All messages with a debug level equal or higher than the return value of
/// this function will cause an exception to be thrown after log emission.
///
/// @note Depending on preprocessor settings @c ACTS_ENABLE_LOG_FAILURE_THRESHOLD
///       and @c ACTS_LOG_FAILURE_THRESHOLD, this operations is either constexpr
///       or a runtime operation.
Level getFailureThreshold();

#else

#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
#ifdef ACTS_LOG_FAILURE_THRESHOLD
// We have a fixed compile time log failure threshold
constexpr Level getFailureThreshold() {
  return Level::ACTS_LOG_FAILURE_THRESHOLD;
}
#else
Level getFailureThreshold();
#endif
#else
constexpr Level getFailureThreshold() {
  // Default "NO" failure threshold
  return Level::MAX;
}
#endif

#endif

/// @brief Set debug level above which an exception will be thrown after logging
///
/// All messages with a debug level equal or higher than @p level will
/// cause an exception to be thrown after log emission.
///
/// @warning The runtime log failure threshold is **global state**, therefore
///          this function is  **not threadsafe**. The intention is that this
///          level is set once, before multi-threaded execution begins, and then
///          not modified before the end of the job.
/// @note This function is only available if @c ACTS_LOG_FAILURE_THRESHOLD is
///       unset, i.e. no compile-time threshold is used. Otherwise an
///       exception is thrown.
void setFailureThreshold(Level level);

/// Custom exception class so threshold failures can be caught
class ThresholdFailure : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

}  // namespace Logging

/// @brief class for printing debug output
///
/// This class provides the user interface for printing debug messages with
/// different levels of severity.
///
/// @ingroup Logging
class Logger {
 public:
  /// @brief construct from output print and filter policy
  ///
  /// @param [in] pPrint  policy for printing debug messages
  /// @param [in] pFilter policy for filtering debug messages
  Logger(Acts::Logging::Level level, const std::string& name,
         std::ostream& os = std::cout);

  Logger(std::shared_ptr<spdlog::logger> logger);

  ~Logger();

  /// @brief decide whether a message with a given debug level has to be printed
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true if debug message should be printed, otherwise @c false
  bool doPrint(const Logging::Level& lvl) const;

  /// @brief log a debug message
  ///
  /// @param [in] lvl debug level of debug message
  /// @param [in] input text of debug message
  void log(Logging::Level lvl, const std::string& input) const;

  /// Return the level of the filter policy of this logger
  /// @return the level
  Logging::Level level() const;

  /// Return the name of the print policy of this logger
  /// @return the name
  const std::string& name() const;

  /// Make a copy of this logger, optionally changing the name or the level
  /// @param _name the optional new name
  /// @param _level the optional new level
  std::unique_ptr<Logger> clone(
      const std::optional<std::string>& _name = std::nullopt,
      const std::optional<Logging::Level>& _level = std::nullopt) const;

  /// Make a copy of the logger, with a new level. Convenience function for
  /// if you only want to change the level but not the name.
  /// @param _level the new level
  /// @return the new logger
  std::unique_ptr<Logger> clone(Logging::Level _level) const;

  /// Make a copy of the logger, with a suffix added to the end of it's
  /// name. You can also optionally supply a new level
  /// @param suffix the suffix to add to the end of the name
  /// @param _level the optional new level
  std::unique_ptr<Logger> cloneWithSuffix(
      const std::string& suffix,
      std::optional<Logging::Level> _level = std::nullopt) const;

  /// Helper function so a logger reference can be used as is with the logging
  /// macros
  const Logger& operator()() const { return *this; }

 private:
  std::shared_ptr<spdlog::logger> m_logger;
};

/// @brief get default debug output logger
///
/// @param [in] name       name of the logger instance
/// @param [in] lvl        debug threshold level
/// @param [in] log_stream output stream used for printing debug messages
///
/// This function returns a pointer to a Logger instance with the following
/// decorations enabled:
/// - time stamps
/// - name of logging instance
/// - debug level
///
/// @return pointer to logging instance
std::unique_ptr<const Logger> getDefaultLogger(
    const std::string& name, const Logging::Level& lvl,
    std::ostream* log_stream = &std::cout);

const Logger& getDummyLogger();

}  // namespace Acts
