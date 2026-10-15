"""Compare the unnormalised PhotoTest limits s=0 and s=1."""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis, wave_label, phase_branch
from tutorials.paths import PROJECT_DIR
from .scan_data import load_scan, unnormalised, sample_errors, intensity_bounds, comparison_errors

CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/photo_test_mass_scan"
MANIFEST_FILE = CASE_DIR / "output/root/scan_points.csv"
OUTPUT_DIR = CASE_DIR / "output/plots/ideal_cases_unnormalised"
SHOW_PLOTS = False
FIGURE_DPI = 100
TARGET_ASYMMETRIES = (0.0, 1.0)


def finish_axis(axis, title, ylabel):
    """Apply shared labels and styling to one ideal-case panel."""
    axis.set(xlabel=r"Invariant mass $M$ (GeV)", ylabel=ylabel, title=title)
    style_axis(axis, grid=True)


def amplitude_legend():
    """Return the marker legend for combined ideal-case amplitude figures."""
    return [
        Line2D(
            [],
            [],
            marker="*",
            linestyle="none",
            color="tab:green",
            markerfacecolor="tab:green",
            markersize=11,
            label=r"$s=1$; best 0.1% range",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            color="black",
            markerfacecolor="tab:green",
            label=r"$s=0$; best 0.1% range",
        ),
    ]


def plot_intensity(masses, intensities):
    """Plot the total and S/P/D model intensities for the two ideal cases."""
    figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
    for column, label, colour in (
        (0, "Total", "black"),
        (1, "S wave", "tab:blue"),
        (2, "P waves", "tab:orange"),
        (3, "D waves", "tab:green"),
    ):
        axis.plot(
            masses,
            intensities[1, :, column],
            marker="*",
            markersize=10,
            color=colour,
            label=label,
            zorder=4,
        )
        axis.plot(
            masses,
            intensities[0, :, column],
            marker="o",
            color=colour,
            zorder=3,
        )
    finish_axis(
        axis,
        "Unnormalised PhotoTest invariant-mass intensity",
        "Unnormalised model intensity (arb. units)",
    )
    axis.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="*",
                linestyle="none",
                color="black",
                markersize=10,
                label=r"$s=1$",
            ),
            Line2D([], [], marker="o", linestyle="none", color="black", label=r"$s=0$"),
        ]
        + [
            Line2D([], [], color=colour, label=label)
            for _, label, colour in (
                (0, "Total", "black"),
                (1, "S wave", "tab:blue"),
                (2, "P waves", "tab:orange"),
                (3, "D waves", "tab:green"),
            )
        ],
        frameon=False,
        ncol=2,
    )
    save_figure(
        figure,
        OUTPUT_DIR / "unnormalised_mass_intensity.pdf",
        show=SHOW_PLOTS,
        dpi=FIGURE_DPI,
    )


def plot_moments(masses, moment_pairs, set_moments, reconstructed_moments, samples):
    """Compare generated and fitted moments in separate s=0 and s=1 panels."""
    errors = sample_errors(reconstructed_moments, samples)
    moment_dir = OUTPUT_DIR / "moments"
    for moment_index, (_, moment) in enumerate(moment_pairs):
        figure, axes = plt.subplots(
            1, 2, figsize=(15, 6), sharex=True, sharey=True, constrained_layout=True
        )
        for panel, (axis, asymmetry) in enumerate(zip(axes, TARGET_ASYMMETRIES)):
            axis.scatter(
                masses,
                set_moments[panel, :, moment_index],
                marker="*",
                s=90,
                color="tab:orange",
                edgecolors="black",
                zorder=4,
                label="Set moment",
            )
            axis.errorbar(
                masses,
                reconstructed_moments[panel, :, moment_index],
                yerr=errors[:, panel, :, moment_index],
                fmt="o", ms=7.7, capsize=3,
                color="tab:orange",
                zorder=3,
                label="Reconstructed moment; best 0.1% range",
            )
            axis.axhline(0.0, color="0.45", linewidth=0.8)
            finish_axis(
                axis,
                rf"Unnormalised {moment}, $s={asymmetry:g}$",
                "Unnormalised moment value (arb. units)" if panel == 0 else "",
            )
        axes[0].legend(frameon=False)
        save_figure(
            figure, moment_dir / f"{moment}.pdf", show=SHOW_PLOTS, dpi=FIGURE_DPI
        )


def plot_intensity_difference(masses, fitted_amplitudes, lower, upper):
    """Compare wave-summed fitted intensities by difference and ratio."""
    # Compare the total squared magnitude recovered by the unpolarized minimizer
    # in the two ideal asymmetry limits. The sum is over all reflectivity and
    # orbital-wave amplitudes at each invariant mass.
    summed_fit_intensity_s0 = np.sum(np.abs(fitted_amplitudes[0]) ** 2, axis=1)
    summed_fit_intensity_s1 = np.sum(np.abs(fitted_amplitudes[1]) ** 2, axis=1)
    summed_fit_intensity_difference = summed_fit_intensity_s0 - summed_fit_intensity_s1
    summed_fit_intensity_ratio = np.divide(
        summed_fit_intensity_s0,
        summed_fit_intensity_s1,
        out=np.full_like(summed_fit_intensity_s0, np.nan),
        where=summed_fit_intensity_s1 > 0.0,
    )
    difference_errors, ratio_errors = comparison_errors(
        summed_fit_intensity_s0, summed_fit_intensity_s1,
        lower[0], upper[0], lower[1], upper[1],
    )
    total_fit_intensity_difference = np.sum(summed_fit_intensity_difference)

    figure, axes = plt.subplots(
        1, 2, figsize=(15, 6), sharex=True, constrained_layout=True
    )
    difference_axis, ratio_axis = axes
    difference_axis.errorbar(
        masses, summed_fit_intensity_difference, yerr=difference_errors,
        marker="o", color="tab:green", capsize=3,
    )
    difference_axis.axhline(0.0, color="0.25", linewidth=0.9)
    finish_axis(
        difference_axis,
        r"Fitted intensity difference: $s=0$ minus $s=1$",
        r"$\Delta I_{\rm fit}(M) = \sum_{\rm waves} (|A_{\rm fit}^{s=0}|^2-|A_{\rm fit}^{s=1}|^2)$",
    )
    difference_axis.text(
        0.03,
        0.97,
        rf"$\sum_M \Delta I_{{\rm fit}}={total_fit_intensity_difference:.5g}$",
        transform=difference_axis.transAxes,
        verticalalignment="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9},
    )
    ratio_axis.errorbar(masses, summed_fit_intensity_ratio, yerr=ratio_errors,
                        marker="o", color="tab:purple", capsize=3)
    finish_axis(
        ratio_axis,
        r"Fitted intensity ratio: $s=0$ divided by $s=1$",
        r"$R_{\rm fit}(M) = \frac{\sum_{\rm waves}|A_{\rm fit}^{s=0}|^2}{\sum_{\rm waves}|A_{\rm fit}^{s=1}|^2}$",
    )
    save_figure(
        figure,
        OUTPUT_DIR / "summed_fit_intensity_difference.pdf",
        show=SHOW_PLOTS,
        dpi=FIGURE_DPI,
    )


def plot_amplitudes(masses, waves, fitted_amplitudes, samples):
    """Plot each fitted magnitude and phase with the original case markers."""
    magnitude_errors = sample_errors(np.abs(fitted_amplitudes), samples, np.abs)
    phase_errors = sample_errors(np.angle(fitted_amplitudes), samples, np.angle)
    amplitude_dir = OUTPUT_DIR / "amplitudes"
    for wave_index, wave in enumerate(waves):
        figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
        axis.errorbar(
            masses,
            np.abs(fitted_amplitudes[1, :, wave_index]),
            yerr=magnitude_errors[:, 1, :, wave_index],
            fmt="*", ms=9.5, capsize=3,
            color="tab:green",
            markeredgecolor="black",
            zorder=4,
        )
        axis.errorbar(
            masses,
            np.abs(fitted_amplitudes[0, :, wave_index]),
            yerr=magnitude_errors[:, 0, :, wave_index],
            fmt="o", ms=7.7, capsize=3,
            color="tab:green",
            zorder=3,
        )
        finish_axis(axis, wave_label(wave), "Unnormalised amplitude magnitude")
        axis.legend(handles=amplitude_legend(), loc="upper right", frameon=False)
        save_figure(
            figure,
            amplitude_dir / f"{wave}_ideal_magnitude.pdf",
            show=SHOW_PLOTS,
            dpi=FIGURE_DPI,
        )

        figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
        axis.errorbar(
            masses,
            np.angle(fitted_amplitudes[1, :, wave_index]),
            yerr=phase_errors[:, 1, :, wave_index],
            fmt="*", ms=9.5, capsize=3,
            color="tab:green",
            markeredgecolor="black",
            zorder=4,
        )
        axis.errorbar(
            masses,
            np.angle(fitted_amplitudes[0, :, wave_index]),
            yerr=phase_errors[:, 0, :, wave_index],
            fmt="o", ms=7.7, capsize=3,
            color="tab:green",
            zorder=3,
        )
        finish_axis(axis, wave_label(wave), "Phase (rad)")
        axis.set_ylim(-np.pi, np.pi)
        axis.set_yticks(
            [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
            [r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"],
        )
        axis.legend(handles=amplitude_legend(), loc="upper right", frameon=False)
        save_figure(
            figure,
            amplitude_dir / f"{phase_branch(wave)}_ideal_phase.pdf",
            show=SHOW_PLOTS,
            dpi=FIGURE_DPI,
        )


def main():
    data = load_scan(MANIFEST_FILE, require_unnormalised=True)
    indices = []
    for asymmetry in TARGET_ASYMMETRIES:
        matches = np.flatnonzero(np.isclose(data.asymmetries, asymmetry, atol=1.0e-12))
        if len(matches) != 1:
            raise RuntimeError(f"The manifest must contain exactly one s={asymmetry:g}")
        indices.append(int(matches[0]))
    display = unnormalised(data)
    masses, waves, moment_pairs = data.masses, data.waves, data.moment_pairs
    fitted_amplitudes = display.fitted_amplitudes[indices]
    intensities = data.unnormalised_intensities[indices]
    set_moments = display.truth_moments[indices]
    reconstructed_moments = display.fitted_moments[indices]
    if not np.all(np.isfinite(fitted_amplitudes)):
        raise RuntimeError(
            "The ideal-case comparison requires all requested mass points"
        )
    plot_intensity(masses, intensities)
    plot_moments(masses, moment_pairs, set_moments, reconstructed_moments,
                 display.moment_samples[indices])
    lower, upper = intensity_bounds(display)
    plot_intensity_difference(masses, fitted_amplitudes, lower[indices], upper[indices])
    plot_amplitudes(masses, waves, fitted_amplitudes, display.amplitude_samples[indices])
    print(
        f"Wrote {2 + len(moment_pairs) + 2 * len(waves)} figures beneath {OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()
