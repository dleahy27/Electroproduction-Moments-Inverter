"""Runnable EMI studies and the small helpers shared between them.

The package deliberately performs no eager imports.  In particular, importing
``tutorials`` should not initialise PyROOT; individual analysis modules import
ROOT only when they need to read a ROOT file.

Generated ROOT files and plots remain local to each study; this package contains
only reusable Python code and the modules needed to reproduce those outputs.
"""

__all__ = ["analysis_utils", "paths"]
