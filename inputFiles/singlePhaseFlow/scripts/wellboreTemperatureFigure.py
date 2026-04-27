"""
wellboreTemperatureFigure.py

Reads GEOS results (produced by wellboreTemperatureQueries.py)
and produces a 2x2 subplot figure comparing at 4 time steps:
  - GEOS results (from model-results.txt)
  - FDM non-linear reference
  - Linear analytic solution for infinite domain(Wang & Papamichos 1994)
  - Steady-state log solution

Supports three test cases via xmlFilePrefix:
  - 'thermalCompressible_2d'  (linear / constant properties)
  - 'thermalCompressible_temperatureDependentVolumetricHeatCapacity'
  - 'thermalCompressible_temperatureDependentSinglePhaseThermalConductivity'

Usage:
    python wellboreTemperatureFigure.py --xmlFilePrefix <prefix> [--geosDir /path/to/GEOS] [--savePng]
"""

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from xml.etree import ElementTree

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import wellboreTemperatureSolutions as wts


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def _xmlList(s):
    return s.replace('{', '').replace('}', '').strip().split(',')


def getWellboreGeometryFromXML(xmlPath):
    tree  = ElementTree.parse(xmlPath)
    radii = _xmlList(tree.find('Mesh/InternalWellbore').get('radius'))
    return float(radii[0]), float(radii[-1])


def getParametersFromXML(xmlPath, xmlFilePrefix):
    """Parse constitutive / BC parameters from the base XML.

    Parameters
    ----------
    xmlPath : str
        Path to the base XML file (thermalCompressible_2d_base.xml).
    xmlFilePrefix : str
        Test-case prefix, determines which constitutive models to look for.

    Returns
    -------
    tuple
        (Tin, Tout, lambda0, lambda_gradient, c0, c_gradient, Tref, porosity)
    """
    tree     = ElementTree.parse(xmlPath)
    fsParams = tree.findall('FieldSpecifications/FieldSpecification')

    Pin = Pout = Tin = Tout = None
    for fs in fsParams:
        if fs.get('fieldName') == 'pressure' and fs.get('initialCondition') != '1':
            if fs.get('setNames') == '{ rneg }':
                Pin  = float(fs.get('scale'))
            if fs.get('setNames') == '{ rpos }':
                Pout = float(fs.get('scale'))
    for fs in fsParams:
        if fs.get('fieldName') == 'temperature' and fs.get('initialCondition') != '1':
            if fs.get('setNames') == '{ rneg }':
                Tin  = float(fs.get('scale'))
            if fs.get('setNames') == '{ rpos }':
                Tout = float(fs.get('scale'))

    lambda0 = lambda_gradient = Tref = c0 = c_gradient = None

    if 'ThermalConductivity' in xmlFilePrefix:
        # T-dependent conductivity, constant heat capacity
        for tc in tree.findall('Constitutive/SinglePhaseThermalConductivity'):
            if tc.get('name') == 'thermalCond_nonLinear':
                lambda0         = float(_xmlList(tc.get('defaultThermalConductivityComponents'))[0])
                lambda_gradient = float(_xmlList(tc.get('thermalConductivityGradientComponents'))[0])
                Tref            = float(tc.get('referenceTemperature'))
        for sie in tree.findall('Constitutive/SolidInternalEnergy'):
            if sie.get('name') == 'rockInternalEnergy_linear':
                c0         = float(sie.get('referenceVolumetricHeatCapacity'))
                c_gradient = 0.0
                if Tref is None:
                    Tref = float(sie.get('referenceTemperature'))
    elif 'VolumetricHeatCapacity' in xmlFilePrefix:
        # T-dependent heat capacity, constant conductivity
        for tc in tree.findall('Constitutive/SinglePhaseThermalConductivity'):
            if tc.get('name') == 'thermalCond_linear':
                lambda0         = float(_xmlList(tc.get('defaultThermalConductivityComponents'))[0])
                lambda_gradient = 0.0
        for sie in tree.findall('Constitutive/SolidInternalEnergy'):
            if sie.get('name') == 'rockInternalEnergy_nonLinear':
                c0         = float(sie.get('referenceVolumetricHeatCapacity'))
                c_gradient = float(sie.get('dVolumetricHeatCapacity_dTemperature'))
                Tref       = float(sie.get('referenceTemperature'))
    else:
        # Linear / constant properties
        for tc in tree.findall('Constitutive/SinglePhaseThermalConductivity'):
            if tc.get('name') == 'thermalCond_linear':
                lambda0         = float(_xmlList(tc.get('defaultThermalConductivityComponents'))[0])
                lambda_gradient = 0.0
        for sie in tree.findall('Constitutive/SolidInternalEnergy'):
            if sie.get('name') == 'rockInternalEnergy_linear':
                c0         = float(sie.get('referenceVolumetricHeatCapacity'))
                c_gradient = 0.0
                Tref       = float(sie.get('referenceTemperature'))

    porosity = float(tree.find('Constitutive/PressurePorosity').get('defaultReferencePorosity'))

    return Tin, Tout, lambda0, lambda_gradient, c0, c_gradient, Tref, porosity


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(xmlFilePrefix='', outputDir=None):
    _scripts_dir      = os.path.dirname(os.path.abspath(__file__))
    _default_geos_dir = os.path.normpath(os.path.join(_scripts_dir, '..', '..', '..'))

    parser = argparse.ArgumentParser(
        description='Plot wellbore temperature validation (reads model-results.txt).')
    parser.add_argument('--geosDir', default=_default_geos_dir,
                        help='Path to the GEOS repository root')
    parser.add_argument('--xmlFilePrefix', default='',
                        help='Test-case XML prefix')
    parser.add_argument('--outputDir', default=outputDir if outputDir is not None else '.',
                        help='Path to the directory containing model-results.txt and where the figure is saved')
    parser.add_argument('--savePng', action='store_true',
                        help='Save PNG to outputDir using <xmlFilePrefix>.png')
    args = parser.parse_args()

    if not xmlFilePrefix:
        xmlFilePrefix = args.xmlFilePrefix
    if not xmlFilePrefix:
        raise ValueError('xmlFilePrefix must be provided either as argument or via --xmlFilePrefix')

    xmlBase      = os.path.join(args.geosDir, 'inputFiles', 'singlePhaseFlow',
                                'thermalCompressible_2d_base.xml')
    xmlBenchmark = os.path.join(args.geosDir, 'inputFiles', 'singlePhaseFlow',
                                xmlFilePrefix + '_benchmark.xml')

    Rin, Rout = getWellboreGeometryFromXML(xmlBenchmark)
    (Tin, Tout,
     lambda0, lambda_gradient,
     c0, c_gradient, Tref, porosity) = getParametersFromXML(xmlBase, xmlFilePrefix)

    # Whether properties are temperature-dependent (non-linear)
    isNonLinear = (lambda_gradient != 0.0) or (c_gradient != 0.0)

    # Porosity-weighted reference diffusivity for the linear analytic solution
    effectiveRefCap    = (1.0 - porosity) * c0
    thermalDiffusivity = lambda0 / effectiveRefCap

    # -----------------------------------------------------------------------
    # Load results
    # -----------------------------------------------------------------------
    dataPath = os.path.join(args.outputDir, 'model-results.txt')
    data     = np.loadtxt(dataPath, skiprows=1)
    # columns: time, r, temperature
    unique_times = np.unique(data[:, 0])

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    plt.rc('font', size=16)

    for chart_idx, t_val in enumerate(unique_times):
        mask             = data[:, 0] == t_val
        r_cells          = data[mask, 1]
        temperature_geos = data[mask, 2]

        # FDM non-linear reference (only when λ or c depends on T)
        if isNonLinear:
            T_fdm = wts.solveRadialDiffusion(
                r_cells, Rin, Rout,
                t_val, t_val / 100.0,
                Tin, Tout,
                lambda0, lambda_gradient,
                c0, c_gradient, Tref, porosity)

        # Linear analytic (Wang & Papamichos 1994)
        T_linear = wts.transientTemperature(Tin, Tout, Rin, r_cells, thermalDiffusivity, t_val)

        # Steady-state log profile
        T_ss = wts.steadyState(Tin, Tout, Rin, Rout, r_cells)

        ax = axes.flat[chart_idx]
        ax.plot(r_cells, temperature_geos, 'k+',  label='GEOS')
        if isNonLinear:
            ax.plot(r_cells, T_fdm,        'g-',  label='FDM Non-Linear')
        ax.plot(r_cells, T_linear,         'r-',  label='Analytic Linear, infinite domain')
        ax.plot(r_cells, T_ss,             'b-',  label='Steady State')

        if chart_idx == 1:
            ax.legend(fontsize=11)
        if chart_idx in [2, 3]:
            ax.set_xlabel('Radial distance from well centre (m)')
        if chart_idx in [0, 2]:
            ax.set_ylabel('Temperature (°C)')

        ax.set_ylim(-30, 110)
        ax.set_xlim(0.0, 1.0 if not isNonLinear else 0.5)
        ax.set_title(f't = {t_val:g} s')
        ax.grid(True)

    plt.tight_layout()
    if args.savePng:
        pngName = xmlFilePrefix + '.png'
        pngPath = os.path.join(args.outputDir, pngName)
        plt.savefig(pngPath, dpi=150)
        print(f'Saved {pngPath}')
    plt.show()


if __name__ == '__main__':
    main()
