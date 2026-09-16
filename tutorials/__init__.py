"""Runnable EMI studies and the small helpers shared between them.

The package deliberately performs no eager imports.  In particular, importing
``tutorials`` should not initialise PyROOT; individual analysis modules import
ROOT only when they need to read a ROOT file.
"""

__all__ = ["analysis_utils", "paths"]
