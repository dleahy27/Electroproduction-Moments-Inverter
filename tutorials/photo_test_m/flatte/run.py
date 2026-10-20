"""Generate exact S-wave Flatte plus D-wave BW moments and fit them."""
import argparse
import subprocess
from tutorials.paths import PROJECT_DIR

CASE = PROJECT_DIR / "tutorials/photo_test_m/flatte"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int)
    parser.add_argument("--workers", type=int)
    options = parser.parse_args()
    output = CASE / "output"
    output.mkdir(parents=True, exist_ok=True)
    common = [str(PROJECT_DIR / "build/emi"), "--settings", str(CASE / "settings.h")]
    subprocess.run([common[0], "generate-mass-model", *common[1:],
                    "--output", str(output / "moments.root")], cwd=PROJECT_DIR, check=True)
    command = [common[0], "mass-fit", *common[1:], "--input", str(output / "moments.root"),
               "--output", str(output / "fit.root")]
    for flag in ("starts", "workers"):
        value = getattr(options, flag)
        if value is not None: command.extend(["--" + flag, str(value)])
    subprocess.run(command, cwd=PROJECT_DIR, check=True)

if __name__ == "__main__":
    main()
