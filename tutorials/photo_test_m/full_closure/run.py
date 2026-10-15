"""Generate and fit the full PhotoTest closure from its supplied truth start."""

import subprocess
from tutorials.paths import PROJECT_DIR


def main():
    settings = PROJECT_DIR / "tutorials/photo_test_m/full_closure/settings.h"
    settings.with_name("output").mkdir(parents=True, exist_ok=True)
    for command in ("generate-mass-model", "mass-fit"):
        subprocess.run(
            [str(PROJECT_DIR / "build/emi"), command, "--settings", str(settings)],
            cwd=PROJECT_DIR,
            check=True,
        )


if __name__ == "__main__":
    main()
