#!/usr/bin/env python3
"""VisIt renderer for the 2D plane-strain copper foam compression example."""
from copperFoam_visit_common import main

if __name__ == "__main__":
    raise SystemExit(main(default_case_name="copperFoamCompression2D", default_view="xy", default_width=1200, default_height=1600))
