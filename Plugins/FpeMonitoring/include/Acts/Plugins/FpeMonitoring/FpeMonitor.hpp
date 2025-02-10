// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <atomic>
#include <csignal>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <stack>
#include <vector>

#include <boost/container/static_vector.hpp>
#include <boost/stacktrace/stacktrace_fwd.hpp>

namespace Acts {

enum class FpeType : std::uint32_t {
  INTDIV = FPE_INTDIV,
  INTOVF = FPE_INTOVF,
  FLTDIV = FPE_FLTDIV,
  FLTOVF = FPE_FLTOVF,
  FLTUND = FPE_FLTUND,
  FLTRES = FPE_FLTRES,
  FLTINV = FPE_FLTINV,
  FLTSUB = FPE_FLTSUB,
};

std::ostream &operator<<(std::ostream &os, FpeType type);

class FpeMonitor {
 public:
  static constexpr std::size_t kBufferSize = 65536;
  struct Buffer {
    explicit Buffer()
        : m_data{std::make_unique<std::array<std::byte, kBufferSize>>()} {}
    Buffer(const Buffer &) = delete;
    Buffer(Buffer &&other) noexcept {
      m_data = std::move(other.m_data);
      m_offset = other.m_offset;
      other.m_offset = 0;
    }

    std::pair<void *, std::size_t> next() {
      return {m_data.get() + m_offset, size() - m_offset};
    }

    void pushOffset(std::size_t offset) {
      assert(m_offset + offset < size());
      m_offset = offset;
    }

    void reset() { m_offset = 0; }

    std::size_t size() const { return m_data->size(); }
    std::size_t offset() const { return m_offset; }

    std::byte *data() { return m_data->data(); }

   private:
    std::unique_ptr<std::array<std::byte, kBufferSize>> m_data;
    std::size_t m_offset{};
  };

  struct Result {
    struct FpeInfo {
      std::size_t count;
      FpeType type;
      std::shared_ptr<const boost::stacktrace::stacktrace> st;

      FpeInfo(std::size_t countIn, FpeType typeIn,
              std::shared_ptr<const boost::stacktrace::stacktrace> stIn);
      ~FpeInfo();
    };

    Result merged(const Result &with) const;
    void merge(const Result &with);

    bool encountered(FpeType type) const;
    unsigned int count(FpeType type) const;

    const std::vector<FpeInfo> &stackTraces() const;
    unsigned int numStackTraces() const;

    void deduplicate();

    bool contains(FpeType type, const boost::stacktrace::stacktrace &st) const;

    void summary(
        std::ostream &os,
        std::size_t depth = std::numeric_limits<std::size_t>::max()) const;

    Result() = default;

    bool hasStackTraces() const { return !m_stackTraces.empty(); }

    void add(Acts::FpeType type, void *stackPtr, std::size_t bufferSize);

   private:
    std::vector<FpeInfo> m_stackTraces;
    std::array<unsigned int, 32> m_counts{};

    friend FpeMonitor;
  };

  FpeMonitor();
  explicit FpeMonitor(int excepts);
  FpeMonitor(FpeMonitor &&other) = default;
  ~FpeMonitor();

  Result &result();

  void consumeRecorded();

  void rearm();

  static std::string stackTraceToString(const boost::stacktrace::stacktrace &st,
                                        std::size_t depth);
  static std::string getSourceLocation(const boost::stacktrace::frame &frame);

  static bool canSymbolize();

 private:
  void enable();
  void disable();

  static void ensureSignalHandlerInstalled();
  static void signalHandler(int signal, siginfo_t *si, void *ctx);

  struct GlobalState {
    std::atomic_bool isSignalHandlerInstalled{false};
    std::mutex mutex{};
  };

  static std::stack<FpeMonitor *> &stack();
  static GlobalState &globalState();

  int m_excepts = 0;

  Result m_result;

  Buffer m_buffer{};

  boost::container::static_vector<std::tuple<FpeType, void *, std::size_t>, 128>
      m_recorded;
};

}  // namespace Acts
