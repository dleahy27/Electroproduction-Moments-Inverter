"""Check mass-fit ROOT output without repeating the minimizations.

Checks accepted rows for finite physical parameters and ordered masses, all rows
for fixed gauge phases, and available covariances against the full Jacobian map.
Use individual fit files: merged scan trees do not carry each case's metadata.
"""
import argparse
from pathlib import Path

import numpy as np
import ROOT


def require(condition, message):
    if not condition:
        raise ValueError(message)


def metadata(source, name):
    # An empty TObjString is valid (e.g. a single-resonance ordering map),
    # but PyROOT converts it to False. Check the key rather than its contents.
    require(bool(source.GetKey(name)), f"Missing metadata: {name}")
    return source.Get(name)


def verify(path, starts=None, workers=None, require_covariance=False, mass_branch="mass_GeV"):
    source = ROOT.TFile.Open(str(path), "READ")
    require(bool(source) and not source.IsZombie(), f"Cannot open {path}")
    try:
        physical = str(metadata(source, 'mass_fit_physical_coordinate_names').GetString()).splitlines()
        fitted = str(metadata(source, 'mass_fit_coordinate_names').GetString()).splitlines()
        names = list(dict.fromkeys(['fit_ok', 'chi2', 'start', 'worker', 'covariance_ok', *physical, *fitted]))
        data = ROOT.RDataFrame('massFitResults', source).AsNumpy(names)
        count = len(data['chi2'])
        require(count > 0, 'No saved starts')
        require(starts is None or count == starts, f'Expected {starts} starts, got {count}')
        require(np.array_equal(np.sort(data['start']), np.arange(count)), 'Missing or duplicate start indices')
        if workers is not None:
            require(np.array_equal(np.unique(data['worker']), np.arange(workers)), 'Worker indices differ from expectation')
        good = np.asarray(data['fit_ok'], dtype=bool) & np.isfinite(data['chi2'])
        require(good.any(), 'No accepted finite minimum')
        for name in physical:
            values = data[name][good]
            require(np.isfinite(values).all(), f'Nonfinite {name}')
            if name.endswith(('_mass', '_width')):
                require((values > 0).all(), f'Nonpositive {name}')
            elif name.endswith('_g2'):
                require((values >= 0).all(), f'Negative {name}')
        for name in fitted:
            if name.startswith('positive__'):
                np.testing.assert_allclose(data[name.removeprefix('positive__')][good],
                                           np.exp(data[name][good]), rtol=1e-12)
        gauges = metadata(source, 'mass_fit_phase_gauge_count').GetVal()
        for index in range(gauges):
            coordinate = metadata(source, f'mass_fit_phase_gauge_coordinate_{index}').GetVal()
            name = fitted[coordinate]
            require(np.all(data[name] == data[name][0]), f'Reference phase changed: {name}')
        gap_floor = metadata(source, 'mass_fit_minimum_gap_GeV').GetVal()
        for line in str(metadata(source, 'mass_fit_ordering_map').GetString()).splitlines():
            group, resonance, current, previous, coordinate, kind = line.split()
            current, previous = int(current), int(previous)
            if previous < 0:
                continue
            high, low = data[physical[current]][good], data[physical[previous]][good]
            gaps = high - low
            # Subtraction at very large free masses can round below the floor.
            tolerance = 16 * np.finfo(float).eps * (np.abs(high) + np.abs(low))
            require(np.all((gaps > 0) & (gaps + tolerance >= gap_floor)),
                    f'Mass order failed in {group}: {resonance}')
        tree = metadata(source, 'massFitResults')
        covariances = np.flatnonzero(good & np.asarray(data['covariance_ok'], dtype=bool))
        require(not require_covariance or len(covariances) > 0, 'No accepted covariance available')
        for entry in covariances:
            tree.GetEntry(int(entry))
            jacobian = np.asarray(tree.mass_jacobian).reshape(len(physical), len(fitted))
            fit_covariance = np.asarray(tree.fit_covariance).reshape(len(fitted), len(fitted))
            physical_covariance = np.asarray(tree.physical_covariance).reshape(len(physical), len(physical))
            require(np.isfinite(physical_covariance).all(), 'Nonfinite physical covariance')
            np.testing.assert_allclose(physical_covariance, jacobian @ fit_covariance @ jacobian.T,
                                       rtol=1e-8, atol=1e-10)
        for name in ('massFitPredictions', 'massFitCurves'):
            tree = source.Get(name)
            if not tree:
                require(name == 'massFitCurves', 'Missing original-grid predictions')
                continue
            columns = [mass_branch, *[b.GetName() for b in tree.GetListOfBranches()
                                     if b.GetName().startswith('model__')]]
            require(len(columns) > 1 and tree.GetEntries() > 0, f'Empty {name}')
            values = ROOT.RDataFrame(name, source).AsNumpy(columns)
            require(all(np.isfinite(v).all() for v in values.values()), f'Nonfinite {name}')
            if name == 'massFitCurves':
                require(np.all(np.diff(values[mass_branch]) > 0), 'Dense mass grid is not increasing')
        best = float(np.min(data['chi2'][good]))
        np.testing.assert_allclose(best, metadata(source, 'mass_fit_best_chi2').GetVal(), rtol=1e-12)
        return (f'{path}: PASS; starts={count}, accepted={good.sum()}, best chi2={best:.12g}; '
                f'fixed references={gauges}, covariances checked={len(covariances)}')
    finally:
        source.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', type=Path, nargs='+')
    parser.add_argument('--starts', type=int)
    parser.add_argument('--workers', type=int)
    parser.add_argument('--require-covariance', action='store_true')
    parser.add_argument('--mass-branch', default='mass_GeV')
    options = parser.parse_args()
    for path in options.files:
        print(verify(path, options.starts, options.workers,
                     options.require_covariance, options.mass_branch))


if __name__ == '__main__':
    main()
