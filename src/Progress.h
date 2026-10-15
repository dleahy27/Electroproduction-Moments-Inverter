#pragma once

// Cosmetic progress for ROOT forked workers. A normal parent-owned counter
// would be copied by fork, so one lock-free atomic lives in shared mmap storage.
// This is the only process-shared state; fitted data still use worker ROOT files.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <sstream>
#include <string>
#include <thread>
#include <utility>

#include <sys/mman.h>
#include <unistd.h>

namespace emi::detail {

// Non-owning worker handle: Advance counts completed trials (including failed
// fits). It must not outlive the monitor in the parent. A null handle is a no-op.
// Relaxed ordering is sufficient: this count does not publish any fit data.
struct ProgressCounter {
  std::atomic<std::uint64_t>* completed = nullptr;
  void Advance() const {
    if (completed)
      completed->fetch_add(1, std::memory_order_relaxed);
  }
};

// TProcessExecutor children inherit this MAP_SHARED counter.  Only the parent
// renders output, so parallel workers never interleave terminal messages.
// The monitor must surround the entire pool.Map call, including child joins.
// mmap failure disables cosmetic reporting without aborting a physics fit.
// Only lock-free atomics are safe here: fallback locks need not be shared
// between forked processes.
class ProgressMonitor {
  static_assert(std::atomic<std::uint64_t>::is_always_lock_free,
                "Progress requires a process-shared lock-free counter");

public:
  ProgressMonitor(bool enabled, std::uint64_t total, unsigned workers,
                  std::string label)
      : enabled_(enabled), total_(total), workers_(workers), label_(std::move(label)) {
    if (!enabled_)
      return;
    storage_ = static_cast<std::atomic<std::uint64_t>*>(
        mmap(nullptr, sizeof(*storage_), PROT_READ | PROT_WRITE,
             MAP_SHARED | MAP_ANONYMOUS, -1, 0));
    if (storage_ == MAP_FAILED) {
      storage_ = nullptr;
      enabled_ = false;
      return;
    }
    new (storage_) std::atomic<std::uint64_t>(0);
    interactive_ = isatty(STDERR_FILENO);
    started_ = std::chrono::steady_clock::now();
    thread_ = std::thread([this] { Run(); });
  }
  ProgressMonitor(const ProgressMonitor&) = delete;
  ~ProgressMonitor() {
    if (!enabled_)
      return;
    stop_.store(true, std::memory_order_relaxed);
    thread_.join();
    Render(true);
    storage_->~atomic();
    munmap(storage_, sizeof(*storage_));
  }
  ProgressCounter Counter() const {
    return {storage_};
  }

private:
  void Run() {
    while (!stop_.load(std::memory_order_relaxed)) {
      Render(false);
      std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
  }
  static std::string Time(double seconds) {
    const auto value = static_cast<std::uint64_t>(seconds < 0 ? 0 : seconds);
    std::ostringstream out;
    if (value >= 3600)
      out << value / 3600 << 'h' << std::setw(2) << std::setfill('0')
          << (value / 60) % 60 << 'm';
    else
      out << value / 60 << 'm' << std::setw(2) << std::setfill('0') << value % 60
          << 's';
    return out.str();
  }
  void Render(bool final) const {
    const auto completed = storage_->load(std::memory_order_relaxed);
    // Batch logs receive a new line only when the count changes. Interactive
    // terminals redraw in place every 250 ms to refresh elapsed time and ETA.
    if (!interactive_ && completed == lastRendered_)
      return;
    lastRendered_ = completed;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(now - started_).count();
    const double fraction =
        total_ ? std::min(1.0, static_cast<double>(completed) / total_) : 1.0;
    const int filled = static_cast<int>(fraction * 30.0);
    std::ostringstream out;
    out << label_ << " [" << std::string(filled, '=') << std::string(30 - filled, ' ')
        << "] " << completed << '/' << total_ << " ("
        << static_cast<int>(fraction * 100.0) << "%)"
        << " workers=" << workers_ << " elapsed=" << Time(elapsed) << " ETA=";
    if (completed == 0)
      out << "--";
    else
      out << Time(elapsed / completed * (total_ > completed ? total_ - completed : 0));
    // Return to start of line and erase the previous contents.
    if (interactive_)
      std::cerr << "\r\033[2K";
    std::cerr << out.str();
    if (final || !interactive_)
      std::cerr << '\n';
    std::cerr << std::flush;
  }
  bool interactive_ = false;
  mutable std::uint64_t lastRendered_ = std::numeric_limits<std::uint64_t>::max();
  bool enabled_ = false;
  std::uint64_t total_ = 0;
  unsigned workers_ = 0;
  std::string label_;
  std::atomic<std::uint64_t>* storage_ = nullptr;
  std::chrono::steady_clock::time_point started_{};
  mutable std::atomic<bool> stop_{false};
  std::thread thread_;
};

} // namespace emi::detail
