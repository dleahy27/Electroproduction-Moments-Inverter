"""Matplotlib plots for one loaded PhotoTest scan.

Data loading is in scan_data.py. Each function below produces one plot family;
the normalized and unnormalized entry points differ only in the display scale.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis, wave_label, phase_branch, range_errors
from tutorials.paths import PROJECT_DIR
from .scan_data import load_scan, unnormalised, parameter_label, sample_errors

CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/photo_test_mass_scan"
MANIFEST_FILE = CASE_DIR / "output/root/scan_points.csv"
SHOW_PLOTS = False
FIGURE_DPI = 100
POINT_OFFSET_FRACTION = 0.06
UNNORMALISED_INTENSITY_ASYMMETRY = 0.0


def offset_masses(masses, number_of_asymmetries):
    """Separate coincident asymmetry markers while retaining mass as x.

    The displacement is purely graphical and is expressed in GeV.  It does not
    modify the mass used to generate or fit any scan point.
    """
    if len(masses) < 2 or number_of_asymmetries < 2:
        return np.tile(masses, (number_of_asymmetries, 1))
    offsets = np.linspace(
        -POINT_OFFSET_FRACTION / 2,
        POINT_OFFSET_FRACTION / 2,
        number_of_asymmetries,
    )
    return masses[None, :] + offsets[:, None]


def asymmetry_handles(asymmetries, colours):
    """Make a compact colour legend shared by the scan figures."""
    return [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            color=colour,
            label=rf"$s={asymmetry:.3g}$",
        )
        for asymmetry, colour in zip(asymmetries, colours)
    ]


def finish_mass_axis(axis, title, ylabel):
    """Apply the common mass-scan labels and styling."""
    axis.set(xlabel=r"Invariant mass $M$ (GeV)", ylabel=ylabel, title=title)
    style_axis(axis, grid=True)


def save_scan_figure(figure, path):
    """Save and close one scan figure."""
    save_figure(figure, path, show=SHOW_PLOTS, dpi=FIGURE_DPI)


def plot_style(data):
    """Use the same colour and marker offsets in every mass-scan plot."""
    colours = plt.colormaps["viridis"](np.linspace(0.08, 0.92, len(data.asymmetries)))
    return (
        offset_masses(data.masses, len(data.asymmetries)),
        colours,
        asymmetry_handles(data.asymmetries, colours),
    )


def plot_moment_comparisons(data, directory, generated, fitted, ylabel, title_prefix, samples):
    """Write one generated-versus-fitted figure for every moment."""
    plot_masses, colours, colour_legend = plot_style(data)
    errors = sample_errors(fitted, samples)
    for index, (_, moment) in enumerate(data.moment_pairs):
        figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
        for row, colour in enumerate(colours):
            axis.scatter(
                plot_masses[row],
                generated[row, :, index],
                marker="*",
                s=46,
                color=colour,
                zorder=4,
                edgecolors="black",
            )
            axis.errorbar(
                plot_masses[row],
                fitted[row, :, index],
                yerr=errors[:, row, :, index],
                fmt="o", ms=7.7, capsize=3,
                color=colour,
                zorder=3,
            )
        axis.axhline(0.0, color="0.45", linewidth=0.8)
        finish_mass_axis(axis, f"{title_prefix}{moment}", ylabel)
        marker_legend = [
            Line2D(
                [],
                [],
                marker="*",
                linestyle="none",
                color="black",
                markersize=10,
                label="Input moment",
            ),
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                color="black",
                label="Fitted moment; best 0.1% range",
            ),
        ]
        axis.legend(handles=colour_legend + marker_legend, frameon=False, ncol=2)
        save_scan_figure(figure, directory / f"{moment}.pdf")


def plot_amplitudes(data, display, directory, ylabel):
    """Plot each wave's magnitude and phase with both injected k sectors."""
    plot_masses, colours, colour_legend = plot_style(data)
    amplitude_dir = directory / "amplitudes"
    magnitude_errors = sample_errors(np.abs(display.fitted_amplitudes), display.amplitude_samples, np.abs)
    phase_errors = sample_errors(np.angle(data.fitted_amplitudes), data.amplitude_samples, np.angle)
    for index, wave in enumerate(data.waves):
        title = wave_label(wave)

        figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
        for row, colour in enumerate(colours):
            axis.scatter(
                plot_masses[row],
                np.abs(display.truth_plus[row, :, index]),
                marker="^",
                s=46,
                color=colour,
                zorder=4,
                edgecolors="black",
            )
            axis.scatter(
                plot_masses[row],
                np.abs(display.truth_minus[row, :, index]),
                marker="v",
                s=46,
                color=colour,
                zorder=4,
                edgecolors="black",
            )
            axis.errorbar(
                plot_masses[row],
                np.abs(display.fitted_amplitudes[row, :, index]),
                yerr=magnitude_errors[:, row, :, index],
                fmt="o", ms=7.7, capsize=3,
                color=colour,
                zorder=3,
            )
        finish_mass_axis(axis, title, ylabel)
        amplitude_markers = [
            Line2D(
                [],
                [],
                marker="^",
                linestyle="none",
                color="black",
                label=r"Generated $k=+1$",
            ),
            Line2D(
                [],
                [],
                marker="v",
                linestyle="none",
                color="black",
                label=r"Generated $k=-1$",
            ),
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                color="black",
                label="Unpolarized fit; best 0.1% range",
            ),
        ]
        axis.legend(handles=colour_legend + amplitude_markers, frameon=False, ncol=2)
        save_scan_figure(figure, amplitude_dir / f"{wave}_magnitude.pdf")

        figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
        for row, colour in enumerate(colours):
            plus_phase = np.angle(data.truth_plus[row, :, index])
            minus_phase = np.angle(data.truth_minus[row, :, index])
            plus_phase[np.abs(data.truth_plus[row, :, index]) < 1.0e-12] = np.nan
            minus_phase[np.abs(data.truth_minus[row, :, index]) < 1.0e-12] = np.nan
            axis.scatter(
                plot_masses[row],
                plus_phase,
                marker="^",
                s=46,
                color=colour,
                zorder=4,
                edgecolors="black",
            )
            axis.scatter(
                plot_masses[row],
                minus_phase,
                marker="v",
                s=46,
                color=colour,
                zorder=4,
                edgecolors="black",
            )
            axis.errorbar(
                plot_masses[row],
                np.angle(data.fitted_amplitudes[row, :, index]),
                yerr=phase_errors[:, row, :, index],
                fmt="o", ms=7.7, capsize=3,
                color=colour,
                zorder=3,
            )
        finish_mass_axis(axis, title, "Phase (rad)")
        axis.set_ylim(-np.pi, np.pi)
        axis.set_yticks(
            [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
            [r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"],
        )
        axis.legend(handles=colour_legend + amplitude_markers, frameon=False, ncol=2)
        save_scan_figure(figure, amplitude_dir / f"{phase_branch(wave)}.pdf")


def plot_correlations(data, directory):
    """Plot recorded Hessian correlations without synthesizing missing entries."""
    correlation_dir = directory / "correlations"
    correlation_figure_count = 0
    labels = [parameter_label(name) for name in data.correlation_parameters]
    for asymmetry_row, asymmetry in enumerate(data.asymmetries):
        for mass_column, mass in enumerate(data.masses):
            matrix = data.correlations[asymmetry_row, mass_column]
            off_diagonal = ~np.eye(len(matrix), dtype=bool)
            if not np.any(np.isfinite(matrix[off_diagonal])):
                continue
            figure, axis = plt.subplots(figsize=(18, 16), constrained_layout=False)
            image = axis.imshow(
                matrix, vmin=-1.0, vmax=1.0, cmap="coolwarm", interpolation="none"
            )
            axis.set_xticks(np.arange(len(labels)), labels, rotation=90)
            axis.set_yticks(np.arange(len(labels)), labels)
            axis.set_title(
                rf"Amplitude correlations: $M={mass:.3f}$ GeV, $s={asymmetry:.3g}$"
            )
            figure.colorbar(image, ax=axis, label="Correlation coefficient")
            save_scan_figure(
                figure,
                correlation_dir
                / f"asymmetry_{asymmetry_row:02d}_mass_{mass_column:02d}_correlation.pdf",
            )
            correlation_figure_count += 1

    return correlation_figure_count


def plot_signed_errors(data, display, directory):
    """Keep the natural/unnatural, linear/squared, per-wave/sum panel layout."""
    for mass_column, mass in enumerate(data.masses):
        quadrature = np.hypot(
            np.abs(display.truth_plus[:, mass_column]),
            np.abs(display.truth_minus[:, mass_column]),
        )
        differences = [
            np.abs(display.fitted_amplitudes[:, mass_column]) - quadrature,
            np.abs(display.fitted_amplitudes[:, mass_column]) ** 2 - quadrature**2,
        ]
        selected_samples = [np.abs(display.amplitude_samples[row, mass_column])
                            for row in range(len(data.asymmetries))]
        figure, axes = plt.subplots(2, 4, figsize=(28, 10), sharex=True)
        axes[0, 1].sharey(axes[0, 0])
        axes[0, 3].sharey(axes[0, 2])
        for power, errors in enumerate(differences):
            for sector, (reflectivity, title) in enumerate(
                (("a", "Natural"), ("b", "Unnatural"))
            ):
                column = 2 * power + sector
                top, bottom = axes[:, column]
                selected = [
                    (index, wave)
                    for index, wave in enumerate(data.waves)
                    if wave.startswith(reflectivity + "_")
                ]
                colours = plt.colormaps["tab10"](
                    np.linspace(0.0, 0.9, max(len(selected), 1))
                )
                for colour, (index, wave) in zip(colours, selected):
                    wave_errors = np.column_stack([
                        range_errors(errors[row, index], samples[:, index] ** (power + 1)
                                     - quadrature[row, index] ** (power + 1))
                        for row, samples in enumerate(selected_samples)
                    ])
                    top.errorbar(
                        data.asymmetries,
                        errors[:, index],
                        yerr=wave_errors, capsize=3,
                        marker="o",
                        color=colour,
                        label=wave_label(wave),
                    )
                indices = [index for index, _ in selected]
                summed = np.sum(errors[:, indices], axis=1)
                sum_errors = np.column_stack([
                    range_errors(summed[row], np.sum(
                        samples[:, indices] ** (power + 1)
                        - quadrature[row, indices] ** (power + 1), axis=1))
                    for row, samples in enumerate(selected_samples)
                ])
                bottom.errorbar(
                    data.asymmetries,
                    summed, yerr=sum_errors, capsize=3,
                    fmt="o",
                    color="tab:blue",
                )
                top.set_title(
                    title
                    + " exchange — "
                    + (
                        "squared-amplitude difference"
                        if power
                        else "amplitude difference"
                    )
                )
                bottom.set_title(
                    title
                    + " exchange — "
                    + (
                        "summed squared-amplitude difference"
                        if power
                        else "summed difference"
                    )
                )
                top.legend(frameon=False, ncol=2)
                if sector == 0:
                    top.set_ylabel(
                        r"$|A_{\rm fit}|^2-|A_{\rm quad}|^2$"
                        if power
                        else r"$|A_{\rm fit}|-\sqrt{|A_{+}|^2+|A_{-}|^2}$"
                    )
                    bottom.set_ylabel(
                        r"Summed $|A_{\rm fit}|^2-|A_{\rm quad}|^2$"
                        if power
                        else "Summed signed error"
                    )
        for axis in axes.flat:
            axis.axhline(0.0, color="black", linestyle="--", linewidth=0.9)
            axis.set(xlabel=r"K-sector asymmetry $s$", xlim=(-1.05, 1.05))
            style_axis(axis, grid=True)
        figure.suptitle(rf"Invariant mass $M={mass:.3f}$ GeV; bars: best 0.1% range")
        save_scan_figure(
            figure, directory / "signed_errors" / f"signed_error_mass_{mass:.3f}.pdf"
        )


def plot_chi2(data, directory):
    """Preserve the mass-bin sum of the original selected objective values."""
    summed_chi2 = np.sum(data.minimum_chi2, axis=1)

    figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
    axis.plot(data.asymmetries, summed_chi2, marker="o", color="tab:blue")
    axis.set(
        xlabel=r"$s$",
        ylabel=r"$\sum_M \chi^2_{\min}(s, M)$",
        title=r"Mass bin summed reduced $\chi^2$ vs $k$-asymmetry factor",
    )
    axis.set_xlim(-1.05, 1.05)
    style_axis(axis, grid=True)
    save_scan_figure(figure, directory / "chi2_vs_asymmetry.pdf")


def plot_intensity(data, directory):
    """Show the total/S/P/D raw intensity at the configured balanced point."""
    # The generator exports this quantity before H_0000 normalization. It is
    # an angularly integrated intensity in arbitrary units, useful for showing
    # resonance content but not a calibrated cross section.
    if UNNORMALISED_INTENSITY_ASYMMETRY not in {
        value: index for index, value in enumerate(data.asymmetries)
    }:
        raise RuntimeError(
            "The manifest has no scan point at "
            f"s={UNNORMALISED_INTENSITY_ASYMMETRY:g}"
        )
    unnormalised_intensity_row = {
        value: index for index, value in enumerate(data.asymmetries)
    }[UNNORMALISED_INTENSITY_ASYMMETRY]
    figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
    for column, label, colour in (
        (0, "Total", "black"),
        (1, "S wave", "tab:blue"),
        (2, "P waves", "tab:orange"),
        (3, "D waves", "tab:green"),
    ):
        axis.plot(
            data.masses,
            data.unnormalised_intensities[unnormalised_intensity_row, :, column],
            marker="o",
            color=colour,
            label=label,
        )
    finish_mass_axis(
        axis,
        "Unnormalised PhotoTest mass intensity, "
        rf"$s={UNNORMALISED_INTENSITY_ASYMMETRY:g}$",
        "Unnormalised model intensity (arb. units)",
    )
    axis.legend(frameon=False)
    save_scan_figure(figure, directory / "unnormalised_mass_intensity.pdf")


def main(plot_unnormalised=False):
    data = load_scan(MANIFEST_FILE, require_unnormalised=plot_unnormalised)
    display = unnormalised(data) if plot_unnormalised else data
    directory = (
        CASE_DIR
        / "output/plots"
        / ("unnormalised" if plot_unnormalised else "normalised")
    )
    font_style = {"font.family": "serif", "mathtext.fontset": "stix", "font.size": 11}
    plot_chi2(data, directory)
    plt.rcParams.update(font_style)
    if plot_unnormalised:
        plot_intensity(data, directory)
    plot_moment_comparisons(
        data,
        directory / "moments",
        display.truth_moments,
        display.fitted_moments,
        (
            "Unnormalised moment value (arb. units)"
            if plot_unnormalised
            else "Moment value"
        ),
        "Unnormalised " if plot_unnormalised else "",
        display.moment_samples,
    )
    plot_amplitudes(
        data,
        display,
        directory,
        (
            "Unnormalised amplitude magnitude"
            if plot_unnormalised
            else "Amplitude magnitude"
        ),
    )
    correlation_count = plot_correlations(data, directory)
    plot_signed_errors(data, display, directory)
    count = (
        len(data.moment_pairs)
        + 2 * len(data.waves)
        + len(data.masses)
        + correlation_count
        + int(plot_unnormalised)
        + 1  # The chi-square summary is a separate figure.
    )
    print(f"Read {data.read_points} of {data.total_points} scan points")
    print(f"Wrote {count} figures beneath {directory}")


if __name__ == "__main__":
    main()
