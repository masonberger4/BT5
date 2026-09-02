"""`python -m bt5` -- same entry point as the `bt5` console script."""

from __future__ import annotations

import sys

from bt5.cli import main

if __name__ == "__main__":
    sys.exit(main())
