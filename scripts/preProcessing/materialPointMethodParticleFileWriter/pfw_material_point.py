#!/usr/bin/env python3
"""Compatibility entry point for the isolated material-point driver tool."""
from material_point_driver.pfw_material_point import *  # noqa: F401,F403
from material_point_driver.pfw_material_point import main

if __name__ == "__main__":
    raise SystemExit(main())
