"""Plot and check the full PhotoTest closure using the shared closure plotter."""

import argparse
from tutorials.paths import PROJECT_DIR
from tutorials.photo_test_m.closure.analyse import analyse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    options = parser.parse_args()
    analyse(PROJECT_DIR / "tutorials/photo_test_m/full_closure", options.check)


if __name__ == "__main__":
    main()
