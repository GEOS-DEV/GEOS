# -*- coding: utf-8 -*-
"""Uniform particle cooling, exact adiabatic dT/dt.

Units: meters, kg, seconds, Kelvin, watts. See ../README.md.
Requires thermal PFW integration from geos_mpm_thermal_extension.patch.
"""
#[pfw_dependency] input:../thermal_pfw.py
#[pfw_dependency] input:../thermal_cases.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from thermal_pfw import build
pfw = build('sink')
