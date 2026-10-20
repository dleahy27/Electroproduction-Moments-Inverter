"""Plot and verify the f0(980) Flatte + f2(1270) BW closure test."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tutorials.analysis_utils import read_tree, save_figure, style_axis
from tutorials.photo_test_m.analysis import load_closure, moment_indices, coupling_name, pole_names
from tutorials.paths import PROJECT_DIR

CASE = PROJECT_DIR / "tutorials/photo_test_m/flatte"
PION = 0.13957039
KAON = 0.493677


def shapes(mass, parameters):
    """Independent reconstruction in the same denominator convention as C++."""
    mass = np.asarray(mass)
    rho_pi = np.sqrt(1 - (2 * PION / mass)**2 + 0j)
    rho_k = np.sqrt(1 - (2 * KAON / mass)**2 + 0j)
    scalar = 1 / (parameters['f0_980_mass']**2 - mass**2 - 1j * (
        parameters['f0_980_pipi_g2'] * rho_pi + parameters['f0_980_KK_g2'] * rho_k))
    m, width = parameters['f2_1270_mass'], parameters['f2_1270_width']
    tensor = m * width / (m*m - mass**2 - 1j*m*width)
    return scalar, tensor


def intensities(mass, parameters):
    scalar, tensor = shapes(mass, parameters)
    weights = {}
    for resonance in ('f0_980', 'f2_1270'):
        weights[resonance] = 2 * sum(value**2 for name, value in parameters.items()
                                   if name.startswith(resonance + '_') and name.endswith('_mag'))
    return weights['f0_980'] * abs(scalar)**2, weights['f2_1270'] * abs(tensor)**2


def parameter_plots(data, best, plots):
    for suffix, ylabel, filename in (
        ('_mass', 'mass parameter [GeV]', 'masses'),
        ('_width', 'BW width [GeV]', 'widths'),
        ('_g2', r'decay coefficient $g_i^2$ [GeV$^2$]', 'decay_couplings'),
        ('_mag', 'production magnitude', 'magnitudes'),
        ('_phase', 'production phase [rad]', 'phases'),
    ):
        names = [n for n in data['coordinates'] if n.endswith(suffix)]
        figure, axis = plt.subplots(figsize=(max(7, len(names)*0.8), 4.5))
        x = np.arange(len(names))
        axis.scatter(x, [data['truth'][n] for n in names], marker='x', s=55, label='generated')
        axis.scatter(x, [best[n] for n in names], marker='_', s=170, label='best fit')
        labels = (pole_names(names) if suffix in ('_mass', '_width') else
                  coupling_name(names) if suffix in ('_mag', '_phase') else names)
        axis.set(xticks=x, xticklabels=labels, ylabel=ylabel)
        axis.tick_params(axis='x', rotation=45)
        axis.legend()
        style_axis(axis)
        save_figure(figure, plots / (filename + '.pdf'))


def spectrum_plots(data, best, curves, plots):
    truth, prediction = data['moments'], data['predictions']
    mass = curves['mass_GeV']
    intensity = next(n for n in data['spectra'] if moment_indices(n) == (0, 0, 0))
    figure, (top, residual) = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
                                         gridspec_kw={'height_ratios': [3, 1]})
    top.plot(truth['mass_GeV'], truth[intensity], 'o', ms=3, label='generated points')
    top.plot(mass, curves['model__'+intensity], '-', label='best-fit mass model')
    residual.plot(truth['mass_GeV'], prediction['model__'+intensity]-truth[intensity], 'o', ms=3)
    top.set_ylabel(r'$H^0(00)$')
    residual.set(xlabel='mass [GeV]', ylabel='model − input')
    top.legend()
    for axis in (top, residual): style_axis(axis)
    save_figure(figure, plots/'intensity_and_residuals.pdf')

    groups = {}
    for name in data['spectra']:
        indices = moment_indices(name)
        if indices is not None:
            alpha, l, m = indices
            groups.setdefault((l, m), {})[alpha] = name
    for (l, m), names in sorted(groups.items()):
        alphas = (0, 1) if m == 0 else (0, 1, 2, 3)
        figure, axes = plt.subplots(1 if m == 0 else 2, 2, squeeze=False,
                                   figsize=(10, 3.8 if m == 0 else 7), constrained_layout=True)
        for axis, alpha in zip(axes.flat, alphas):
            name = names[alpha]
            axis.plot(truth['mass_GeV'], truth[name], 'o', ms=3, label='generated points')
            axis.plot(mass, curves['model__'+name], '-', lw=1.5, label='best-fit mass model')
            axis.set(xlabel='mass [GeV]', ylabel=rf'$H^{{{alpha}}}({l}{m})$')
            style_axis(axis)
        axes.flat[0].legend(fontsize=8)
        figure.suptitle(f'S-wave Flatté + D-wave BW; L={l}, M={m}')
        save_figure(figure, plots/'moments'/f'moments_L{l}_M{m}.pdf')

    figure, axis = plt.subplots(figsize=(8, 4.5))
    for index, label in enumerate(('S: '+r'$f_0(980)$ Flatté', 'D: '+r'$f_2(1270)$ BW')):
        axis.plot(truth['mass_GeV'], intensities(truth['mass_GeV'], data['truth'])[index],
                  'o', ms=3, color=f'C{index}', label=label+' generated')
        axis.plot(mass, intensities(mass, best)[index], color=f'C{index}', label=label+' fitted')
    axis.axvline(2*KAON, color='grey', ls=':', label=r'$K^+K^-$ threshold')
    axis.set(xlabel='mass [GeV]', ylabel='partial-wave intensity')
    axis.legend(fontsize=8)
    style_axis(axis)
    save_figure(figure, plots/'partial_wave_intensities.pdf')

    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for col, label in enumerate((r'$f_0(980)$ Flatté', r'$f_2(1270)$ BW')):
        target = shapes(truth['mass_GeV'], data['truth'])[col]
        fitted = shapes(mass, best)[col]
        axes[0, col].plot(truth['mass_GeV'], abs(target)**2, 'o', ms=3, label='generated shape')
        axes[0, col].plot(mass, abs(fitted)**2, '-', label='fitted shape')
        axes[1, col].plot(truth['mass_GeV'], np.angle(target), 'o', ms=3)
        axes[1, col].plot(mass, np.angle(fitted), '-')
        axes[0, col].set(title=label, ylabel=r'$|F(m)|^2$')
        axes[1, col].set(ylabel='shape phase [rad]')
        for axis in axes[:, col]:
            axis.set_xlabel('mass [GeV]')
            axis.axvline(2*KAON, color='grey', ls=':')
            style_axis(axis)
    axes[0, 0].legend()
    save_figure(figure, plots/'lineshapes.pdf')
    # H00 is a sum over orthogonal waves; this also checks C++ curve output
    # against the independent Flatte threshold continuation above.
    np.testing.assert_allclose(sum(intensities(mass, best)), curves['model__'+intensity],
                               rtol=1e-10, atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=CASE)
    options = parser.parse_args()
    data = load_closure(options.case)
    best = {n: float(data['results'][n][data['best']]) for n in data['coordinates']}
    columns = ['mass_GeV', *['model__'+n for n in data['spectra']]]
    curves = read_tree(options.case/'output/fit.root', 'massFitCurves', columns)
    # Include original model evaluations too, preserving the exact KK cusp
    # which need not fall on the uniform dense output grid.
    combined = {name: np.concatenate((curves[name], data['predictions'][name]))
                for name in columns}
    _, indices = np.unique(combined['mass_GeV'], return_index=True)
    curves = {name: values[indices] for name, values in combined.items()}
    assert np.isfinite(np.column_stack(list(curves.values()))).all()
    plots = options.case/'output/plots'
    parameter_plots(data, best, plots)
    spectrum_plots(data, best, curves, plots)
    good = np.asarray(data['results']['fit_ok'], dtype=bool)
    chi2 = float(data['results']['chi2'][data['best']])
    report = [f'Starts: {len(good)}; accepted: {good.sum()}; best chi2: {chi2:.12g}']
    for name in data['coordinates']:
        if name.endswith(('_mass', '_width', '_g2')):
            report.append(f'{name}: fitted {best[name]:.10g}; generated {data["truth"][name]:.10g}')
            np.testing.assert_allclose(best[name], data['truth'][name], rtol=2e-4, atol=2e-6)
    for sector in ('a', 'b'):
        name = f'f0_980_{sector}_0_0_phase'
        assert np.all(data['results'][name] == 0), name
    report.append('PASS: recovered shape parameters, both phase references, finite dense curves, independent H00 reconstruction.')
    report.append('Plots show exact generated points and best-fit curves without uncertainty bars.')
    text = '\n'.join(report)+'\n'
    print(text)
    (options.case/'output/verification_report.txt').write_text(text)

if __name__ == '__main__':
    main()
