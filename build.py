#!/usr/bin/env python3

"""Configure, build, and optionally test EMI with CMake.

The helper is intentionally a thin wrapper: every command it runs is printed
by CMake and can also be copied from the README and run manually.
"""

import argparse
import subprocess
from pathlib import Path


def arguments():
    """Parse convenience options without hiding any CMake configuration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-type", default="Release")
    parser.add_argument("--jobs", type=int, default=0,
                        help="parallel jobs; zero lets CMake choose")
    parser.add_argument("--root-dir", type=Path,
                        help="directory containing ROOTConfig.cmake")
    parser.add_argument("--test", action="store_true",
                        help="run CTest after a successful build")
    return parser.parse_args()


options = arguments()
project_dir = Path(__file__).resolve().parent
build_dir = project_dir / "build"

# An explicit ROOT_DIR is useful on systems with several ROOT installations.
configure = [
    "cmake", "-S", str(project_dir), "-B", str(build_dir),
    f"-DCMAKE_BUILD_TYPE={options.build_type}",
]
if options.root_dir:
    configure.append(f"-DROOT_DIR={options.root_dir.resolve()}")
subprocess.run(configure, check=True)

build_command = ["cmake", "--build", str(build_dir), "--parallel"]
if options.jobs:
    build_command.append(str(options.jobs))
subprocess.run(build_command, check=True)

if options.test:
    subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "--output-on-failure"],
        check=True,
    )

print(f"Built {build_dir / 'emi'}")
