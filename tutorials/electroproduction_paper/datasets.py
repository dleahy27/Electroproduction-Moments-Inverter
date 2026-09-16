"""Shared metadata for the experimental electroproduction channels.

Keeping these numbers in one module prevents the fit runners and the Hessian
and bootstrap plotting scripts from silently using different binning.  Values
of ``q2`` and ``invariant_mass`` are in GeV squared and GeV, respectively.
The published ``R`` rows contain ``(central value, uncertainty)``.
"""

from dataclasses import dataclass
from pathlib import Path

from tutorials.paths import PROJECT_DIR


INPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/input"


@dataclass(frozen=True)
class Channel:
    """Physics metadata shared by runners and analysis programs."""

    key: str
    title: str
    q2: tuple[float, ...]
    invariant_mass: float
    epsilon: float
    published_r: tuple[tuple[float, float], ...] = ()
    published_r_source: str = ""

    @property
    def input_file(self) -> Path:
        """Return the tracked ROOT moment table for this channel."""
        return INPUT_DIR / f"{self.key}_moments.root"


CHANNELS = {
    "e_rho": Channel(
        "e_rho", r"$e^-\rho^0$", (0.82, 1.19, 1.66, 3.06), 4.8, 0.8,
        ((0.649, 0.188), (0.671, 0.065), (0.794, 0.083), (1.068, 0.087)),
        "Published HERMES",
    ),
    "mu_rho": Channel(
        "mu_rho", r"$\mu^-\rho^0$", (1.14, 1.60, 2.80, 6.02), 9.9, 0.9,
        ((0.724, 0.078), (0.930, 0.081), (1.227, 0.122), (1.200, 0.729)),
        "Published COMPASS",
    ),
    "e_omega": Channel(
        "e_omega", r"$e^-\omega$", (1.28, 2.00, 4.00), 4.8, 0.8,
        ((0.226, 0.066), (0.237, 0.091), (0.283, 0.091)),
        "Published HERMES",
    ),
    "mu_omega": Channel(
        "mu_omega", r"$\mu^-\omega$", (1.16, 1.64, 3.58), 7.6, 0.96,
        ((0.475, 0.085), (0.495, 0.124), (0.739, 0.197)),
        "Published COMPASS",
    ),
    "e_phi": Channel(
        "e_phi", r"$e^-\phi$", (1.2, 1.7, 4.5), 21.89 ** 0.5, 0.8,
    ),
}

# The batch runners reproduce the four channels used in the main comparison.
MAIN_CHANNEL_KEYS = ("e_rho", "mu_rho", "e_omega", "mu_omega")
