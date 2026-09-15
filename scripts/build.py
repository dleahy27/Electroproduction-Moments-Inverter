#!/usr/bin/env python3

"""Configure and build EMI."""

import subprocess
from pathlib import Path


# Settings: edit these values before running the script.
BUILD_TYPE = "Release"
PARALLEL_JOBS = 0  # Use 0 to let CMake choose the number of jobs.


# Locate the repository and its build directory.
PROJECT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = PROJECT_DIR / "build"

# Configure, then compile, the project.
subprocess.run(
    [
        "cmake", "-S", str(PROJECT_DIR), "-B", str(BUILD_DIR),
        f"-DCMAKE_BUILD_TYPE={BUILD_TYPE}",
    ],
    check=True,
)
build_command = ["cmake", "--build", str(BUILD_DIR), "--parallel"]
if PARALLEL_JOBS:
    build_command.append(str(PARALLEL_JOBS))
subprocess.run(build_command, check=True)

print(f"Built {BUILD_DIR / 'emi'}")
