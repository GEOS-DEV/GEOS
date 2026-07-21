#!/usr/bin/env python3
"""Render ParticleColor at four states for the rigid-body jamming example."""
from __future__ import annotations

import sys

import pfw_visit_render


def main() -> int:
    # mpm_example_postprocess passes initial/final variable hints meant for the
    # generic renderer.  This case always renders ParticleColor and intentionally
    # overrides the state list so the VisIt output shows color at four times.
    passthrough = []
    skip_next = False
    for i, arg in enumerate(sys.argv[1:]):
        if skip_next:
            skip_next = False
            continue
        if arg in ("--initial-variable", "--final-variable", "--states"):
            skip_next = True
            continue
        passthrough.append(arg)
    argv = passthrough + [
        "--variable", "particleColor",
        "--strict-variable",
        "--states", "initial,quarter,threequarter,final",
        "--range-mode", "explicit",
        "--color-min", "0",
        "--color-max", "5",
        "--point-size-pixels", "7",
        "--view", "xy",
        "--mesh", "CellRegion1",
        "--colortable", "hot_desaturated",
    ]
    return pfw_visit_render.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
