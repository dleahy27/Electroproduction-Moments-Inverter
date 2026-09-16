"""Repository paths shared by lightweight tutorial runners.

``__file__`` points inside the ``tutorials`` package, so its parent's parent is
the checkout containing ``build/emi`` and the tutorial data.  This remains true
when a runner is launched from a different working directory.
"""

from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
