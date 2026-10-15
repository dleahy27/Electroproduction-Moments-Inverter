"""Compare fitted unnormalised intensity at each k asymmetry with s=1."""

import matplotlib.pyplot as plt
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis
from tutorials.paths import PROJECT_DIR
from .scan_data import load_scan, unnormalised, intensity_bounds, comparison_errors

CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/photo_test_mass_scan"
MANIFEST_FILE = CASE_DIR / "output/root/scan_points.csv"
OUTPUT_DIR = CASE_DIR / "output/plots/asymmetry_intensity_comparisons"
SHOW_PLOTS = False
FIGURE_DPI = 100
REFERENCE_ASYMMETRY = 1.0


def finish_axis(axis, title, ylabel):
    """Apply shared scientific styling to an intensity-comparison panel."""
    axis.set(xlabel=r"Invariant mass $M$ (GeV)", ylabel=ylabel, title=title)
    style_axis(axis, grid=True)


def main():
    data = load_scan(MANIFEST_FILE, require_unnormalised=True)
    masses, asymmetries = data.masses, data.asymmetries
    matches = np.flatnonzero(np.isclose(asymmetries, REFERENCE_ASYMMETRY, atol=1.0e-12))
    if len(matches) != 1:
        raise RuntimeError("The manifest must contain exactly one s=1")
    reference_index = int(matches[0])
    display = unnormalised(data)
    fitted_amplitudes = display.fitted_amplitudes
    lower, upper = intensity_bounds(display)
    if not np.all(np.isfinite(fitted_amplitudes)):
        raise RuntimeError(
            "The intensity comparison requires all mass/asymmetry points"
        )
    summed_intensities = np.sum(np.abs(fitted_amplitudes) ** 2, axis=2)
    reference_intensity = summed_intensities[reference_index]
    figure_count = 0
    for asymmetry_index, asymmetry in enumerate(asymmetries):
        if asymmetry_index == reference_index:
            continue
        comparison_intensity = summed_intensities[asymmetry_index]
        intensity_difference = comparison_intensity - reference_intensity
        intensity_ratio = np.divide(
            comparison_intensity,
            reference_intensity,
            out=np.full_like(comparison_intensity, np.nan),
            where=reference_intensity > 0.0,
        )
        difference_errors, ratio_errors = comparison_errors(
            comparison_intensity, reference_intensity,
            lower[asymmetry_index], upper[asymmetry_index],
            lower[reference_index], upper[reference_index],
        )
        total_difference = np.sum(intensity_difference)

        figure, axes = plt.subplots(
            1, 2, figsize=(15, 6), sharex=True, constrained_layout=True
        )
        difference_axis, ratio_axis = axes
        difference_axis.errorbar(
            masses, intensity_difference, yerr=difference_errors,
            marker="o", color="tab:green", capsize=3,
        )
        difference_axis.axhline(0.0, color="0.25", linewidth=0.9)
        finish_axis(
            difference_axis,
            rf"Fitted intensity difference: $s={asymmetry:g}$ minus $s=1$; best 0.1% ranges",
            rf"$\Delta I_{{\rm fit}}^{{s={asymmetry:g}}}(M) = "
            r"\sum_{\rm waves} (|A_{\rm fit}^{s}|^2-|A_{\rm fit}^{s=1}|^2)$",
        )
        difference_axis.text(
            0.03,
            0.97,
            rf"$\sum_M \Delta I_{{\rm fit}}^{{s={asymmetry:g}}}="
            rf"{total_difference:.5g}$",
            transform=difference_axis.transAxes,
            verticalalignment="top",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9},
        )

        ratio_axis.errorbar(masses, intensity_ratio, yerr=ratio_errors,
                            marker="o", color="tab:purple", capsize=3)
        finish_axis(
            ratio_axis,
            rf"Fitted intensity ratio: $s={asymmetry:g}$ divided by $s=1$; best 0.1% ranges",
            rf"$R_{{\rm fit}}^{{s={asymmetry:g}}}(M) = "
            r"\frac{\sum_{\rm waves}|A_{\rm fit}^{s}|^2}"
            r"{\sum_{\rm waves}|A_{\rm fit}^{s=1}|^2}$",
        )
        output_name = f"asymmetry_{asymmetry_index:02d}_vs_s1_intensity_comparison.pdf"
        save_figure(figure, OUTPUT_DIR / output_name, show=SHOW_PLOTS, dpi=FIGURE_DPI)
        figure_count += 1

    print(f"Wrote {figure_count} figures beneath {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
