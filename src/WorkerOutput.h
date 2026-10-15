#pragma once

// File mechanics shared by the fit, mass-fit, and bootstrap runners.
// Each fork owns a ROOT file; the parent merges only after every worker returns.
// These helpers deliberately leave trial scheduling and physics in each runner.
#include "TFileMerger.h"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace emi::detail {

inline std::filesystem::path WorkerPartFile(const std::filesystem::path& output,
                                            unsigned worker) {
  return output.parent_path() / (output.stem().string() + ".part_" +
                                 std::to_string(worker) + output.extension().string());
}

// Copy the same trees and metadata as ROOT's ordinary file merger. Mass-fit
// parts are already local and skip staging; the other runners keep their
// existing ROOT staging behavior. Leave parts available until the caller's
// final processing (such as writing mass predictions) also succeeds.
inline void MergeWorkerFiles(const std::filesystem::path& output,
                             const std::vector<std::string>& parts,
                             bool stageFiles = true) {
  TFileMerger merger(stageFiles, false);
  merger.OutputFile(output.c_str(), "RECREATE");
  for (const auto& part : parts)
    merger.AddFile(part.c_str());
  if (!merger.Merge())
    throw std::runtime_error("Could not merge worker output into " + output.string());
}

} // namespace emi::detail
