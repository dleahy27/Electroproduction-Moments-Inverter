#!/usr/bin/env python3

"""Generate PhotoTest truth points and fit each one with the single-k model."""

import argparse
import csv
import shlex
import subprocess
from pathlib import Path


parser = argparse.ArgumentParser(
    description="Run the 20-by-5 PhotoTest mass and k-minus-scale scan."
)
parser.add_argument("--mass-bins", type=int, default=20)
parser.add_argument("--mass-min", type=float, default=0.5)
parser.add_argument("--mass-max", type=float, default=12.0)
parser.add_argument("--scale-count", type=int, default=5)
parser.add_argument("--scale-min", type=float, default=0.0)
parser.add_argument("--scale-max", type=float, default=0.3)
parser.add_argument("--starts", type=int, default=1000)
parser.add_argument("--workers", type=int, default=8)
parser.add_argument("--seed", type=int, default=12345)
parser.add_argument("--settings", type=Path, default=None)
parser.add_argument("--name", default="photo_test_mass_scan")
parser.add_argument("--no-background", action="store_true")
parser.add_argument(
    "--overwrite", action="store_true",
    help="regenerate and refit points whose ROOT files already exist",
)
parser.add_argument(
    "--dry-run", action="store_true",
    help="print commands and the manifest without running EMI",
)
args = parser.parse_args()

if args.mass_bins < 2 or args.scale_count < 2:
    parser.error("--mass-bins and --scale-count must both be at least two")
if args.mass_min <= 0.0 or args.mass_max <= args.mass_min:
    parser.error("mass limits must satisfy 0 < mass-min < mass-max")
if args.scale_min < 0.0 or args.scale_max < args.scale_min:
    parser.error("scale limits must satisfy 0 <= scale-min <= scale-max")
if args.starts < 1 or args.workers < 0:
    parser.error("--starts must be positive and --workers must be nonnegative")
allowed_name_characters = (
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
)
if not args.name or any(
    character not in allowed_name_characters for character in args.name
):
    parser.error("--name may contain only letters, digits, underscores, and hyphens")

project_dir = Path(__file__).resolve().parents[1]
emi = project_dir / "build" / "emi"
settings = (
    args.settings.resolve()
    if args.settings
    else project_dir / "app" / "UserSettings.h"
)
if not emi.is_file():
    raise SystemExit(f"Missing {emi}; build the project first")
if not settings.is_file():
    raise SystemExit(f"Missing settings file {settings}")

truth_dir = project_dir / "InputFiles" / "Generated" / args.name
fit_dir = project_dir / "OutputFiles" / args.name
truth_dir.mkdir(parents=True, exist_ok=True)
fit_dir.mkdir(parents=True, exist_ok=True)

mass_ratio = args.mass_max / args.mass_min
masses = [
    args.mass_min * mass_ratio ** (index / (args.mass_bins - 1))
    for index in range(args.mass_bins)
]
scales = [
    args.scale_min
    + (args.scale_max - args.scale_min) * index / (args.scale_count - 1)
    for index in range(args.scale_count)
]

manifest_rows = []
total = len(masses) * len(scales)
point = 0
for scale_index, scale in enumerate(scales):
    for mass_index, mass in enumerate(masses):
        point += 1
        stem = f"scale_{scale_index:02d}_mass_{mass_index:02d}"
        truth_file = truth_dir / f"{stem}_truth.root"
        fit_file = fit_dir / f"{stem}_fit.root"
        run_seed = args.seed + scale_index * len(masses) + mass_index

        manifest_rows.append({
            "scale_index": scale_index,
            "mass_index": mass_index,
            "mass_GeV": f"{mass:.17g}",
            "k_minus_scale": f"{scale:.17g}",
            "seed": run_seed,
            "truth_file": str(truth_file.relative_to(project_dir)),
            "fit_file": str(fit_file.relative_to(project_dir)),
        })

        print(
            f"[{point:3d}/{total}] mass={mass:.7g} GeV, "
            f"k-minus scale={scale:.5g}",
            flush=True,
        )

        truth_needs_generation = args.overwrite or not truth_file.exists()
        if truth_needs_generation:
            generate_command = [
                str(emi), "generate-fixed",
                "--settings", str(settings),
                "--output", str(truth_file),
                "--mass", f"{mass:.17g}",
                "--k-minus-scale", f"{scale:.17g}",
                "--quiet",
            ]
            if args.no_background:
                generate_command.append("--no-background")
            print("  " + shlex.join(generate_command), flush=True)
            if not args.dry_run:
                subprocess.run(generate_command, cwd=project_dir, check=True)
        else:
            print(f"  using existing {truth_file.name}", flush=True)

        fit_needs_run = (
            args.overwrite or truth_needs_generation or not fit_file.exists()
        )
        if fit_needs_run:
            fit_command = [
                str(emi), "fit",
                "--settings", str(settings),
                "--input", str(truth_file),
                "--tree", "genMoments",
                "--output", str(fit_file),
                "--photo",
                "--polarization", "none",
                "--starts", str(args.starts),
                "--workers", str(args.workers),
                "--seed", str(run_seed),
                "--no-hesse",
                "--quiet",
            ]
            print("  " + shlex.join(fit_command), flush=True)
            if not args.dry_run:
                subprocess.run(fit_command, cwd=project_dir, check=True)
        else:
            print(f"  using existing {fit_file.name}", flush=True)

manifest = fit_dir / "scan_points.csv"
with manifest.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=manifest_rows[0].keys())
    writer.writeheader()
    writer.writerows(manifest_rows)

print(f"\nWrote scan manifest to {manifest}")
print(
    "Analyse it with: python3 AnalysisScripts/analyse-photo-test-mass-scan.py"
)
