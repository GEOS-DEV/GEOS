#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only
"""Independent real-solver regressions for hydrostatic capillary correctness.

Requires Python VTK. The baseline is the unmodified/contract-only solver; every
run uses a fresh directory. All generated material and boundary values belong
to these explicit test cases and are not application defaults.
"""
import argparse
import copy
import json
import math
import pathlib
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
import vtk

IDENTITY = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
UPWARD = ((1., 0., 0.), (0., -1., 0.), (0., 0., -1.))
SIDEWAYS = ((0., 0., 1.), (0., 1., 0.), (-1., 0., 0.))


def oblique_rotation():
    axis = (1., 2., 3.)
    length = math.sqrt(sum(v*v for v in axis))
    x, y, z = (v / length for v in axis)
    c, s = math.cos(.731), math.sin(.731)
    t = 1-c
    return ((t*x*x+c, t*x*y-s*z, t*x*z+s*y),
            (t*x*y+s*z, t*y*y+c, t*y*z-s*x),
            (t*x*z-s*y, t*y*z+s*x, t*z*z+c))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def transform(matrix, point):
    return tuple(dot(row, point) for row in matrix)


def vector(values):
    return '{ ' + ', '.join(format(x, '.17g') for x in values) + ' }'


MODEL = '''<Problem>
<Solvers gravityVector="{0,0,-9.81}"><CompositionalMultiphaseFVM name="flow" logLevel="0" discretization="tpfa" temperature="350" useMass="1" targetRegions="{reservoir}">
<NonlinearSolverParameters newtonTol="1e-10" newtonMaxIter="15" lineSearchAction="None"/>
<LinearSolverParameters solverType="direct" directParallel="0"/>
</CompositionalMultiphaseFVM></Solvers>
<Mesh><VTKMesh name="mesh" file="mesh.vtu" useGlobalIds="1"/></Mesh>
<Events maxTime="1"><PeriodicEvent name="output" timeFrequency="1" target="/Outputs/vtk"/></Events>
<NumericalMethods><FiniteVolume><TwoPointFluxApproximation name="tpfa"/></FiniteVolume></NumericalMethods>
<ElementRegions><CellElementRegion name="reservoir" cellBlocks="{*}" materialList="{fluid,rock,relperm}"/></ElementRegions>
<Constitutive>
<InvariantImmiscibleFluid name="fluid" componentNames="{CO2,H2O}" phaseNames="{gas,water}" densities="{500,1000}" componentMolarWeight="{0.044,0.018}" viscosities="{0.001,0.001}"/>
<BrooksCoreyRelativePermeability name="relperm" phaseNames="{gas,water}" phaseMinVolumeFraction="{0,0}" phaseRelPermExponent="{2,2}" phaseRelPermMaxValue="{1,1}"/>
<CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="porosity" permeabilityModelName="permeability"/>
<NullModel name="nullSolid"/><PressurePorosity name="porosity" defaultReferencePorosity="0.1" referencePressure="0" compressibility="0"/>
<ConstantPermeability name="permeability" permeabilityComponents="{1e-13,1e-13,1e-13}"/>
</Constitutive>
<FieldSpecifications><HydrostaticEquilibrium name="equilibrium" objectPath="ElementRegions/reservoir" datumElevation="10" datumPressure="10000000" phaseContacts="{5}" componentNames="{CO2,H2O}" componentFractionVsElevationTableNames="{gasFraction,waterFraction}" temperatureVsElevationTableName="temperature" elevationIncrementInHydrostaticPressureTable="0.1" equilibrationTolerance="1e-6" maxNumberOfEquilibrationIterations="20"/></FieldSpecifications>
<Functions>
<TableFunction name="gasFraction" coordinates="{0,4.9,5,5.1,10}" values="{0,0,1,1,1}"/>
<TableFunction name="waterFraction" coordinates="{0,4.9,5,5.1,10}" values="{1,1,0,0,0}"/>
<TableFunction name="temperature" coordinates="{0,10}" values="{300,350}"/>
</Functions>
<Outputs><VTK name="vtk" plotFileRoot="state"/></Outputs>
</Problem>'''


def write_mesh(path, rotation, translation, nz):
    nx = ny = 4
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    point_ids = vtk.vtkIdTypeArray()
    point_ids.SetName('GLOBAL_ID')
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                point = transform(rotation, (10*i/nx, 10*j/ny, 10*k/nz))
                points.InsertNextPoint(*(x+t for x, t in zip(point, translation)))
                point_ids.InsertNextValue(point_ids.GetNumberOfTuples())
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.GetPointData().SetGlobalIds(point_ids)
    cell_ids = vtk.vtkIdTypeArray()
    cell_ids.SetName('GLOBAL_ID')
    attribute = vtk.vtkIntArray()
    attribute.SetName('attribute')
    def node(i, j, k):
        return i+(nx+1)*(j+(ny+1)*k)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = vtk.vtkHexahedron()
                vertices = ((i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k),
                            (i,j,k+1),(i+1,j,k+1),(i+1,j+1,k+1),(i,j+1,k+1))
                for q, xyz in enumerate(vertices):
                    cell.GetPointIds().SetId(q, node(*xyz))
                grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
                cell_ids.InsertNextValue(cell_ids.GetNumberOfTuples())
                attribute.InsertNextValue(0)
    grid.GetCellData().SetGlobalIds(cell_ids)
    grid.GetCellData().AddArray(attribute)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    assert writer.Write() == 1


def records(directory, rotation=IDENTITY, translation=(0.,0.,0.), final=False):
    collections = sorted((directory/'state').glob('*.vtm'))
    assert collections, f'No actual solver output in {directory}'
    collection = collections[-1] if final else collections[0]
    transpose = tuple(zip(*rotation))
    result = {}
    for dataset in ET.parse(collection).iter('DataSet'):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(collection.parent / dataset.get('file')))
        reader.Update()
        grid = reader.GetOutput()
        fields = grid.GetCellData()
        center = fields.GetArray('elementCenter')
        assert center is not None
        for i in range(grid.GetNumberOfCells()):
            world = center.GetTuple3(i)
            local = transform(transpose, tuple(x-t for x,t in zip(world, translation)))
            key = tuple(round(x, 6) for x in local)
            assert key not in result, ('Duplicated owned cell', key)
            result[key] = {fields.GetArrayName(a): fields.GetArray(a).GetTuple(i)
                           for a in range(fields.GetNumberOfArrays())}
    return result


EVIDENCE = {}


def configure_fixed_capillary(tree, three=False, entry=0., minimum=0., maximum=20000., tolerance=1e-6):
    """Explicit invariant-fluid fixture, independent of later initializer patches."""
    phases = ('gas', 'oil', 'water') if three else ('gas', 'oil')
    names = '{' + ','.join(phases) + '}'
    n = len(phases)
    fluid = tree.find('./Constitutive/InvariantImmiscibleFluid')
    fluid.set('phaseNames', names)
    fluid.set('componentNames', names)
    fluid.set('densities', vector([500., 800., 1000.][:n]))
    fluid.set('componentMolarWeight', vector([.044, .114, .018][:n]))
    fluid.set('viscosities', vector([.001] * n))
    relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
    relperm.set('phaseNames', names)
    for key, values in [('phaseMinVolumeFraction', [0.] * n),
                        ('phaseRelPermExponent', [2.] * n), ('phaseRelPermMaxValue', [1.] * n)]:
        relperm.set(key, vector(values))
    tree.find('./ElementRegions/CellElementRegion').set('materialList', '{fluid,rock,relperm,capillary}')
    tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping', '0')
    eq = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
    eq.set('componentNames', names)
    eq.set('componentFractionVsElevationTableNames', '{' + ','.join('z' + phase for phase in phases) + '}')
    eq.set('phaseContacts', '{3,7}' if three else '{5}')
    eq.set('equilibrationTolerance', str(tolerance))
    eq.set('maxNumberOfEquilibrationIterations', '100')
    functions = tree.find('Functions')
    for table in list(functions):
        if table.get('name') != 'temperature':
            functions.remove(table)
    for phase in phases:
        ET.SubElement(functions, 'TableFunction', name='z' + phase, coordinates='{0,10}', values=vector([1 / n] * 2))
    if three:
        options = dict(wettingIntermediateCapPressureTableName='pcwater', nonWettingIntermediateCapPressureTableName='pcgas')
        ET.SubElement(functions, 'TableFunction', name='pcwater', coordinates='{0,1}', values=vector([entry + maximum, entry]))
        ET.SubElement(functions, 'TableFunction', name='pcgas', coordinates='{0,1}', values=vector([entry, entry + maximum]))
    else:
        options = dict(wettingNonWettingCapPressureTableName='pc')
        ET.SubElement(functions, 'TableFunction', name='pc', coordinates=vector([minimum, 1.]), values=vector([entry, entry + maximum]))
    ET.SubElement(tree.find('Constitutive'), 'TableCapillaryPressure', name='capillary', phaseNames=names, **options)


def configure_table_relperm(tree, hidden_oil=False):
    relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
    relperm.tag = 'TableRelativePermeability'
    for key in ('phaseMinVolumeFraction', 'phaseRelPermExponent', 'phaseRelPermMaxValue'):
        relperm.attrib.pop(key)
    relperm.set('wettingIntermediateRelPermTableNames', '{krw,krow}')
    relperm.set('nonWettingIntermediateRelPermTableNames', '{krg,krog}')
    relperm.set('threePhaseInterpolator', 'BAKER')
    for name in ('krw', 'krow', 'krg', 'krog'):
        coordinates = '{.8,1}' if hidden_oil and name == 'krog' else '{0,1}'
        ET.SubElement(tree.find('Functions'), 'TableFunction', name=name, coordinates=coordinates, values='{0,1}')


def phase_flux_metrics(values):
    """Actual vertical TPFA potential and upwind mass flux, with no mobility cutoff."""
    levels = sorted({key[2] for key in values})
    assert len(levels) > 1, 'A flux check needs adjacent cell levels'
    dz = levels[1] - levels[0]
    potential = flux = 0.
    faces = 0
    for (x, y, z), lower in values.items():
        upper = values.get((x, y, round(z + dz, 6)))
        if upper is None:
            continue
        faces += 1
        for ip in range(len(lower['phaseVolumeFraction'])):
            p0 = lower['pressure'][0] - lower['capillary_phaseCapPressure'][ip]
            p1 = upper['pressure'][0] - upper['capillary_phaseCapPressure'][ip]
            density = .5 * (lower['fluid_phaseMassDensity'][ip] + upper['fluid_phaseMassDensity'][ip])
            residual = p0 - p1 - density * 9.81 * dz
            mobility = (lower if residual >= 0. else upper)['phaseMobility'][ip]
            assert all(math.isfinite(value) for value in (p0, p1, density, residual, mobility))
            assert density > 0. and mobility >= 0.
            if mobility > 0.:
                potential = max(potential, abs(residual))
            flux = max(flux, abs(1e-13 * mobility * residual / dz))
    assert faces > 0, 'No actual adjacent cell pairs were checked'
    return dict(max_mobile_potential_pa=potential, max_phase_mass_flux_per_area=flux, checked_faces=faces)



class HydrostaticCapillary(unittest.TestCase):

    count = 0

    def case(self, name, mode=None, rotation=IDENTITY, translation=(0.,0.,0.),
             gravity=None, nz=4, mutate=None, baseline=False, success=True, mpi=False, run_solver=False, expected_cells=None):
        HydrostaticCapillary.count += 1
        directory = ROOT / f'{HydrostaticCapillary.count:02d}-{name}'
        directory.mkdir()
        write_mesh(directory/'mesh.vtu', rotation, translation, nz)
        tree = ET.fromstring(MODEL)
        actual_gravity = gravity if gravity is not None else transform(rotation, (0.,0.,-9.81))
        tree.find('Solvers').set('gravityVector', vector(actual_gravity))
        magnitude = math.sqrt(dot(actual_gravity, actual_gravity))
        up = tuple(-g/magnitude for g in actual_gravity) if magnitude else (0.,0.,1.)
        offset = dot(up, translation) if mode == 'gravityAligned' else translation[2]
        equilibrium = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
        if mode is None:
            equilibrium.attrib.pop('coordinateSystem', None)
        else:
            equilibrium.set('coordinateSystem', mode)
        equilibrium.set('datumElevation', format(10+offset, '.17g'))
        equilibrium.set('phaseContacts', vector((5+offset,)))
        for table in tree.find('Functions'):
            coords = [float(x) for x in table.get('coordinates').strip('{}').split(',')]
            table.set('coordinates', vector(x+offset for x in coords))
        if run_solver:
            # The deliberately incompressible test fluid needs a pressure gauge.
            # The top-cell value is the analytic hydrostatic pressure, not a
            # pressure chosen to conceal an initialization residual.
            geometry = ET.SubElement(tree, 'Geometry')
            specs = tree.find('FieldSpecifications')
            # Pure, zero-capillary phases have disconnected mobility blocks;
            # each needs its own analytic pressure reference in this test.
            for region, lower, upper, pressure, gas in (
                    ('topReference',7.5,10.,1e7+500*9.81*1.25,1),
                    ('bottomReference',0.,2.5,1e7+500*9.81*5+1000*9.81*3.75,0)):
                ET.SubElement(geometry, 'Box', name=region, xMin=vector((0,0,lower)), xMax=vector((10,10,upper)))
                for name, field, component, value in (('Pressure','pressure','-1',str(pressure)),
                                                      ('Gas','globalCompFraction','0',str(gas)),
                                                      ('Water','globalCompFraction','1',str(1-gas))):
                    ET.SubElement(specs, 'FieldSpecification', name=region+name, setNames='{'+region+'}',
                                  objectPath='ElementRegions/reservoir', fieldName=field, component=component, scale=value)
            ET.SubElement(tree.find('Events'), 'PeriodicEvent', name='solve', forceDt='1', target='/Solvers/flow')
        if mutate:
            mutate(tree)
        ET.ElementTree(tree).write(directory/'input.xml', encoding='utf-8', xml_declaration=True)
        program = BASELINE if baseline else GEOS
        command = [program, '-i', str(directory/'input.xml'), '-o', str(directory)]
        if mpi:
            command = [MPIEXEC, '-np', '2', *command, '-x', '2']
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120, text=True)
        (directory/'run.log').write_text(completed.stdout)
        if success is True or (success is None and completed.returncode == 0):
            self.assertEqual(completed.returncode, 0, completed.stdout[-7000:])
            values = records(directory, rotation, translation)
            self.assertEqual(len(values), 4*4*nz if expected_cells is None else expected_cells)
            return directory, values
        self.assertNotEqual(completed.returncode, 0, completed.stdout[-7000:])
        self.assertFalse(list(directory.glob('*.pvd')))
        return directory, completed.stdout


    def compare_fields(self, first, second, tolerance=0.):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            for field in ('pressure','temperature','globalCompFraction','phaseVolumeFraction',
                          'fluid_phaseMassDensity','fluid_phaseDensity'):
                self.assertIn(field, first[key])
                self.assertEqual(len(first[key][field]), len(second[key][field]))
                for a,b in zip(first[key][field], second[key][field]):
                    self.assertTrue(math.isfinite(a) and math.isfinite(b))
                    self.assertLessEqual(abs(a-b), tolerance, (key,field,a,b))


    def test_oil_capillary_pressure_and_rest_state(self):
        def configure(tree, three=False, solve=False):
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','0')
            fluid = tree.find('./Constitutive/InvariantImmiscibleFluid')
            relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
            phases = '{gas,oil,water}' if three else '{gas,oil}'
            fluid.set('phaseNames', phases)
            fluid.set('densities', '{500,800,1000}' if three else '{500,800}')
            relperm.set('phaseNames', phases)
            functions = tree.find('Functions')
            if three:
                fluid.set('componentNames','{C0,C1,C2}')
                fluid.set('componentMolarWeight','{0.044,0.114,0.018}')
                fluid.set('viscosities','{0.001,0.001,0.001}')
                for key,value in {'phaseMinVolumeFraction':'{0,0,0}', 'phaseRelPermExponent':'{2,2,2}',
                                  'phaseRelPermMaxValue':'{1,1,1}'}.items(): relperm.set(key,value)
                equilibrium = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
                equilibrium.set('componentNames','{C0,C1,C2}')
                equilibrium.set('componentFractionVsElevationTableNames','{gasFraction,oilFraction,waterFraction}')
                equilibrium.set('phaseContacts','{3,7}')
                for table in list(functions):
                    if table.get('name') != 'temperature': functions.remove(table)
                for name,values in (('gasFraction','{0,0,0,0,1,1}'),('oilFraction','{0,0,1,1,0,0}'),('waterFraction','{1,1,0,0,0,0}')):
                    ET.SubElement(functions,'TableFunction',name=name,coordinates='{0,2.9,3,6.9,7,10}',values=values)
            tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock,relperm,capillary}')
            options = {'wettingIntermediateCapPressureTableName':'pcow','nonWettingIntermediateCapPressureTableName':'pcgo'} if three else {'wettingNonWettingCapPressureTableName':'pcgo'}
            ET.SubElement(tree.find('Constitutive'),'TableCapillaryPressure',name='capillary',phaseNames=phases,**options)
            ET.SubElement(functions,'TableFunction',name='pcgo',coordinates='{0,1}',values='{0,10000}' if three else '{0,20000}')
            if three: ET.SubElement(functions,'TableFunction',name='pcow',coordinates='{0,1}',values='{10000,0}')
            if solve:
                # A small physical rock compressibility removes the incompressible
                # pressure nullspace without injecting gauge boundary composition.
                # The fluid remains constant-density, so the analytic hydrostatics
                # and capillary equilibria above are unchanged.
                tree.find('./Constitutive/PressurePorosity').set('compressibility','1e-9')
                tree.find('./Solvers/CompositionalMultiphaseFVM').set('logLevel','1')
                # Let the solver accept an already converged rest state. Forcing
                # a Newton correction of a machine-zero residual can introduce
                # signed roundoff into exactly absent components when chopping
                # is intentionally disabled. Pressure/flux checks below remain
                # independent of the nonlinear solver's convergence decision.
                tree.find('./Solvers/CompositionalMultiphaseFVM/NonlinearSolverParameters').set('newtonMinIter','0')
                ET.SubElement(tree.find('Events'),'PeriodicEvent',name='solve',forceDt='1',target='/Solvers/flow')

        def maximum_phase_flux(values):
            # Actual TPFA phase potential/upwind mobility, using solver-output
            # primary pressure, constitutive Pc, density and mobility fields.
            maximum = 0.
            for (x,y,z), lower in values.items():
                upper = values.get((x,y,round(z+1.,6)))
                if upper is None: continue
                for ip in range(len(lower['phaseVolumeFraction'])):
                    p0 = lower['pressure'][0]-lower['capillary_phaseCapPressure'][ip]
                    p1 = upper['pressure'][0]-upper['capillary_phaseCapPressure'][ip]
                    density = .5*(lower['fluid_phaseMassDensity'][ip]+upper['fluid_phaseMassDensity'][ip])
                    potential = p0-p1-density*9.81
                    upstream = lower if potential >= 0 else upper
                    flux_per_area = 1e-13*upstream['phaseMobility'][ip]*potential
                    maximum = max(maximum,abs(flux_per_area))
            return maximum

        for three in (False, True):
            prefix = 'three-capillary' if three else 'gas-oil-capillary'
            _, old = self.case(prefix+'-old',mode=None,baseline=True,nz=10,mutate=lambda t:configure(t,three))
            self.assertGreater(maximum_phase_flux(old),1e-7)  # Reproduces the old physical defect.
            for mode in (None,):
                directory, current = self.case(prefix+'-corrected',mode=mode,nz=10,mutate=lambda t:configure(t,three,True))
                self.assertLess(maximum_phase_flux(current),1e-11)
                times = [float(item.get('timestep')) for item in ET.parse(directory/'state.pvd').iter('DataSet')]
                self.assertIn(0., times)
                self.assertIn(1., times)
                final = records(directory,final=True)
                self.assertLess(maximum_phase_flux(final),1e-11)
                for key in current:
                    self.assertLess(abs(current[key]['pressure'][0]-final[key]['pressure'][0]),1e-3)
                    for a,b in zip(current[key]['phaseVolumeFraction'],final[key]['phaseVolumeFraction']):
                        self.assertLess(abs(a-b),1e-8)
                for (_,_,z),fields in current.items():
                    sg = min(1.,max(0.,300*9.81*(z-(7. if three else 5.))/(10000. if three else 20000.)))
                    self.assertAlmostEqual(fields['phaseVolumeFraction'][0],sg,delta=1e-8)
                    if z > (7. if three else 5.):
                        expected = 1e7+500*9.81*(10-(7. if three else 5.))-800*9.81*(z-(7. if three else 5.))
                        self.assertAlmostEqual(fields['pressure'][0],expected,delta=1e-4)



    def test_capillary_entry_pressure_and_inconsistent_mobile_endpoints(self):
        def capillary(tree, impossible=False):
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','0')
            tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock,relperm,capillary}')
            ET.SubElement(tree.find('Constitutive'),'TableCapillaryPressure',name='capillary',phaseNames='{gas,water}',wettingNonWettingCapPressureTableName='pc')
            ET.SubElement(tree.find('Functions'),'TableFunction',name='pc',coordinates='{0.1,1}' if impossible else '{0,1}',values='{1000,0}' if impossible else '{21000,1000}')
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('phaseContacts','{5.4}')
        for mode in (None,):
            _, values = self.case('capillary-entry',mode=mode,nz=10,mutate=capillary)
            for (_,_,z),fields in values.items():
                desired_pc = 500*9.81*(z-5.4)
                expected_water = 1.-max(0.,min(1.,(desired_pc-1000.)/20000.))
                self.assertAlmostEqual(fields['phaseVolumeFraction'][1],expected_water,delta=1e-8)
                # At z=5.5 gas is geometrically above its contact, but absent
                # until entry pressure is reached; water fixes primary pressure.
                if z == 5.5:
                    expected = 1e7+500*9.81*(10-5.4)+1000*9.81*(5.4-z)+1000.
                    self.assertAlmostEqual(fields['pressure'][0],expected,delta=1e-4)
            _, log = self.case('capillary-inconsistent-endpoint',mode=mode,nz=10,mutate=lambda t:capillary(t,True),success=False)
            self.assertIn('inconsistent mobile-phase capillary pressure',log)


    def test_tiny_saturation_can_have_substantial_mobility(self):
        def configure(tree):
            configure_fixed_capillary(tree, entry=1000., minimum=5e-9)
            tree.find('./Constitutive/BrooksCoreyRelativePermeability').set('phaseRelPermExponent', '{.1,2}')
        _, log = self.case('tiny-but-mobile-endpoint', nz=10, mutate=configure, success=False)
        self.assertIn('inconsistent mobile-phase capillary pressure', log)
        EVIDENCE['tiny_mobile_endpoint'] = dict(result='explicit rejection', saturation=5e-9,
                                               analytic_relative_permeability=(5e-9) ** .1)


    def test_table_baker_oil_mobility_below_reported_minimum(self):
        def configure(tree):
            configure_fixed_capillary(tree, three=True)
            configure_table_relperm(tree, hidden_oil=True)
            water = tree.find("./Functions/TableFunction[@name='pcwater']")
            water.set('coordinates', '{.5,1}')
            water.set('values', '{1000,0}')
            tree.find("./Functions/TableFunction[@name='pcgas']").set('values', '{1000000,1001000}')
        _, log = self.case('table-baker-hidden-mobile-oil', nz=10, mutate=configure, success=False)
        self.assertIn('inconsistent mobile-phase capillary pressure', log)
        EVIDENCE['table_baker_hidden_oil'] = dict(result='explicit rejection', oil_saturation=.5,
                                                 reported_oil_minimum=.8, actual_oil_relative_permeability=.5)


    def test_realized_residual_endpoint_and_default_chopping(self):
        def residual_roundtrip(tree):
            configure_fixed_capillary(tree, entry=1000., minimum=.01)
            relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
            relperm.set('phaseMinVolumeFraction', '{.01,0}')
            relperm.set('phaseRelPermExponent', '{.1,2}')
        _, log = self.case('residual-one-ulp-mobile', nz=10, mutate=residual_roundtrip, success=False)
        self.assertIn('Realized hydrostatic EOS/capillary state violates', log)
        # This second case intentionally leaves the authoritative chopping default
        # untouched. A successful absence fixture must not hide that setting.
        def default_chopping(tree):
            configure_fixed_capillary(tree, entry=1000.)
            tree.find('./Solvers/CompositionalMultiphaseFVM').attrib.pop('allowLocalCompDensityChopping')
        _, log = self.case('default-chopping-creates-mobile-phase', nz=10, mutate=default_chopping, success=False)
        self.assertIn('Realized hydrostatic EOS/capillary state violates', log)
        self.assertIn('allowLocalCompDensityChopping', log)
        EVIDENCE['realized_endpoints'] = dict(residual_roundtrip='explicit rejection', default_chopping='explicit rejection')


    def test_representable_near_entry_complementarity(self):
        entry = 1471.5 - 5e-7
        tolerance = 1e-8
        def configure(tree):
            configure_fixed_capillary(tree, entry=entry, tolerance=tolerance)
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('datumPressure', '100000')
            tree.find('./Constitutive/BrooksCoreyRelativePermeability').set('phaseRelPermExponent', '{.1,2}')
        _, result = self.case('representable-just-super-entry', nz=10, mutate=configure, success=None)
        if isinstance(result, str):
            self.assertTrue('inconsistent mobile-phase capillary pressure' in result or
                            'inconsistent capillary endpoint complementarity' in result or
                            'Realized hydrostatic EOS/capillary state violates' in result, result[-7000:])
            EVIDENCE['near_entry'] = dict(result='explicit endpoint incompatibility rejection')
            return
        # A future pressure-controlled inverse may solve this valid case. A
        # snapped zero gas saturation is not a successful alternative: its
        # pressure error is 50 times the authored, representable tolerance.
        for (_, _, z), fields in result.items():
            expected = max(0., min(1., (300. * 9.81 * (z - 5.) - entry) / 20000.))
            allowance = tolerance + 128. * math.ulp(1.) * max(1., abs(fields['pressure'][0]))
            self.assertLessEqual(abs(fields['phaseVolumeFraction'][0] - expected) * 20000., allowance)
            if z == 5.5:
                self.assertGreater(fields['phaseVolumeFraction'][0], 0.)
        metrics = phase_flux_metrics(result)
        allowance = tolerance + 128. * math.ulp(1.) * max(abs(f['pressure'][0]) for f in result.values())
        self.assertLessEqual(metrics['max_mobile_potential_pa'], 2. * allowance)
        EVIDENCE['near_entry'] = dict(result='pressure-correct successful state', **metrics)


    def test_admissible_table_and_stone2_relative_permeability(self):
        for family in ('table', 'stone2'):
            def configure(tree, family=family):
                configure_fixed_capillary(tree, three=True, maximum=10000.)
                if family == 'table':
                    configure_table_relperm(tree)
                else:
                    relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
                    relperm.tag = 'BrooksCoreyStone2RelativePermeability'
                    for key in ('phaseRelPermExponent', 'phaseRelPermMaxValue'):
                        relperm.attrib.pop(key)
                    relperm.set('phaseMinVolumeFraction', '{0,.1,0}')
                    for key, value in [('waterOilRelPermExponent', '{2,.5}'), ('gasOilRelPermExponent', '{2,.5}'),
                                       ('waterOilRelPermMaxValue', '{1,1}'), ('gasOilRelPermMaxValue', '{1,1}')]:
                        relperm.set(key, value)
            _, values = self.case(family + '-admissible', nz=10, mutate=configure)
            metrics = phase_flux_metrics(values)
            self.assertLess(metrics['max_mobile_potential_pa'], 2e-6)
            self.assertLess(metrics['max_phase_mass_flux_per_area'], 1e-11)
            EVIDENCE[family + '_admissible'] = metrics


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geos',required=True)
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--mpiexec')
    parser.add_argument('--output',required=True,help='New directory for durable numerical evidence')
    parser.add_argument('--test',action='append',help='Run selected test methods; omit for the full acceptance suite')
    args=parser.parse_args()
    GEOS=str(pathlib.Path(args.geos).absolute())
    BASELINE=str(pathlib.Path(args.baseline).absolute())
    MPIEXEC=args.mpiexec
    ROOT=pathlib.Path(args.output).absolute()
    ROOT.mkdir(parents=True,exist_ok=False)
    suite = unittest.TestSuite(HydrostaticCapillary(name) for name in args.test) if args.test else unittest.defaultTestLoader.loadTestsFromTestCase(HydrostaticCapillary)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    (ROOT / 'numerical-evidence.json').write_text(json.dumps(EVIDENCE, indent=2) + '\n')
    raise SystemExit(not result.wasSuccessful())
