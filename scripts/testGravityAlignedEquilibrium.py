#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only
"""Real-solver regression for opt-in potential-distance hydrostatic initialization.

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
<FieldSpecifications><HydrostaticEquilibrium name="equilibrium" objectPath="ElementRegions/reservoir" coordinateSystem="gravityAligned" datumElevation="10" datumPressure="10000000" phaseContacts="{5}" componentNames="{CO2,H2O}" componentFractionVsElevationTableNames="{gasFraction,waterFraction}" temperatureVsElevationTableName="temperature" elevationIncrementInHydrostaticPressureTable="0.1" equilibrationTolerance="1e-6" maxNumberOfEquilibrationIterations="20"/></FieldSpecifications>
<Functions>
<TableFunction name="gasFraction" coordinates="{0,4.9,5,5.1,10}" values="{0,0,1,1,1}"/>
<TableFunction name="waterFraction" coordinates="{0,4.9,5,5.1,10}" values="{1,1,0,0,0}"/>
<TableFunction name="temperature" coordinates="{0,10}" values="{300,350}"/>
</Functions>
<Outputs><VTK name="vtk" plotFileRoot="state"/></Outputs>
</Problem>'''


def write_mesh(path, rotation, translation, nz, z_coordinates=None):
    nx = ny = 4
    levels = z_coordinates if z_coordinates is not None else [10*k/nz for k in range(nz+1)]
    assert len(levels) == nz+1
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    point_ids = vtk.vtkIdTypeArray()
    point_ids.SetName('GLOBAL_ID')
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                point = transform(rotation, (10*i/nx, 10*j/ny, levels[k]))
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


class GravityAlignedEquilibrium(unittest.TestCase):
    count = 0

    def case(self, name, mode='gravityAligned', rotation=IDENTITY, translation=(0.,0.,0.),
             gravity=None, nz=4, mutate=None, baseline=False, success=True, mpi=False, run_solver=False, expected_cells=None, z_coordinates=None):
        GravityAlignedEquilibrium.count += 1
        directory = ROOT / f'{GravityAlignedEquilibrium.count:02d}-{name}'
        directory.mkdir()
        write_mesh(directory/'mesh.vtu', rotation, translation, nz, z_coordinates)
        tree = ET.fromstring(MODEL)
        actual_gravity = gravity if gravity is not None else transform(rotation, (0.,0.,-9.81))
        tree.find('Solvers').set('gravityVector', vector(actual_gravity))
        magnitude = math.sqrt(dot(actual_gravity, actual_gravity))
        up = tuple(-g/magnitude for g in actual_gravity) if magnitude else (0.,0.,1.)
        offset = dot(up, translation) if mode == 'gravityAligned' else translation[2]
        equilibrium = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
        if mode is None:
            equilibrium.attrib.pop('coordinateSystem')
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
        if success:
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

    def check_analytic(self, values, gravity=9.81, selected=False):
        for (x,y,z), fields in values.items():
            if selected and x > 5:
                self.assertEqual(fields['pressure'][0], 9e6)
                continue
            contact_pressure = 1e7 + 500*gravity*5
            pressure = contact_pressure + 1000*gravity*(5-z) if z < 5 else 1e7+500*gravity*(10-z)
            self.assertAlmostEqual(fields['pressure'][0], pressure, delta=1e-4)
            self.assertAlmostEqual(fields['temperature'][0], 300+5*z, delta=1e-9)
            gas = 0. if z < 5 else 1.
            self.assertAlmostEqual(fields['globalCompFraction'][0], gas, delta=1e-12)
            self.assertAlmostEqual(fields['phaseVolumeFraction'][0], gas, delta=1e-12)
            self.assertAlmostEqual(sum(fields['globalCompFraction']), 1., delta=1e-12)

    def test_declared_contract(self):
        import json
        capabilities = subprocess.run([GEOS,'--capabilities','--format=json'],capture_output=True,text=True,check=True)
        document=json.loads(capabilities.stdout)
        self.assertTrue(document['gravity-aligned-hydrostatic-initialization'])
        self.assertEqual(document['hydrostaticInitializationVersion'],1)
        self.assertFalse(document['arbitrary-plane-phase-initialization'])
        self.assertFalse(document['generated-set-phase-initialization'])
        path=ROOT/'catalog.json'
        result=subprocess.run([GEOS,'--input-catalog',str(path)],capture_output=True,text=True)
        (ROOT/'catalog.log').write_text(result.stdout+result.stderr)
        self.assertEqual(result.returncode,0,result.stdout[-3000:])
        catalog=json.loads(path.read_text())
        hydrostatic=next(element for element in catalog['elements'] if element['type']=='HydrostaticEquilibrium')
        properties={p['name']:p for p in hydrostatic['properties']}
        self.assertEqual(set(properties['coordinateSystem']['choices']),{'elevation','gravityAligned'})
        self.assertEqual(properties['coordinateSystem']['default'],'elevation')
        self.assertFalse(properties['setNames']['required'])
        self.assertIsNone(properties['datumElevation']['units'])

    def test_legacy_vertical_compatibility(self):
        for gravity in ((0.,0.,-9.81), (0.,0.,9.81), (0.,0.,0.)):
            _, old = self.case('legacy-baseline', mode=None, gravity=gravity, baseline=True)
            _, current = self.case('legacy-current', mode=None, gravity=gravity)
            self.compare_fields(old, current)
        _, aligned = self.case('aligned-downward')
        _, legacy = self.case('legacy-downward', mode=None)
        self.compare_fields(legacy, aligned, tolerance=1e-7)
        self.check_analytic(aligned)

    def test_rotated_gravity_and_datum_covariance(self):
        _, reference = self.case('reference-contact-centres', nz=5)
        self.check_analytic(reference)
        for rotation, translation in ((UPWARD,(0.,0.,0.)), (SIDEWAYS,(0.,0.,0.)),
                                      (oblique_rotation(),(0.,0.,0.)), (oblique_rotation(),(12.,-7.,3.))):
            _, transformed = self.case('rotated', rotation=rotation, translation=translation, nz=5)
            self.compare_fields(reference, transformed, tolerance=1e-6)
            self.check_analytic(transformed)

    def test_large_coordinate_contact_precision(self):
        rotation = oblique_rotation()
        translation = transform(rotation, (1e8,0.,0.))
        _, reference = self.case('large-frame-reference', nz=5)
        _, translated = self.case('large-frame-on-contact', nz=5, rotation=rotation, translation=translation)
        self.compare_fields(reference, translated, tolerance=2e-4)
        for key in reference:
            for a,b in zip(reference[key]['phaseVolumeFraction'],translated[key]['phaseVolumeFraction']):
                self.assertAlmostEqual(a,b,delta=1e-12)
        # These true offsets exceed the explicitly bounded contact tie band.
        for offset in (-1e-5,1e-5):
            def shifted(tree):
                equilibrium=tree.find('./FieldSpecifications/HydrostaticEquilibrium')
                value=float(equilibrium.get('phaseContacts').strip('{}'))
                equilibrium.set('phaseContacts',vector((value+offset,)))
            _, values=self.case('large-frame-resolved-offset',nz=5,rotation=rotation,translation=translation,mutate=shifted)
            for key,fields in values.items():
                if key[2] == 5.:
                    self.assertAlmostEqual(fields['phaseVolumeFraction'][0],1. if offset < 0 else 0.,delta=1e-12)
        def thin(tree,gap):
            actual=tuple(float(x)for x in tree.find('Solvers').get('gravityVector').strip('{}').split(','))
            up=tuple(-x/math.hypot(*actual)for x in actual)
            offset=dot(up,translation)
            fluid=tree.find('./Constitutive/InvariantImmiscibleFluid')
            for key,value in {'componentNames':'{C0,C1,C2}','phaseNames':'{gas,oil,water}',
                              'densities':'{500,800,1000}','componentMolarWeight':'{0.044,0.114,0.018}',
                              'viscosities':'{0.001,0.001,0.001}'}.items():fluid.set(key,value)
            relperm=tree.find('./Constitutive/BrooksCoreyRelativePermeability')
            for key,value in {'phaseNames':'{gas,oil,water}','phaseMinVolumeFraction':'{0,0,0}',
                              'phaseRelPermExponent':'{2,2,2}','phaseRelPermMaxValue':'{1,1,1}'}.items():relperm.set(key,value)
            equilibrium=tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            equilibrium.set('componentNames','{C0,C1,C2}')
            equilibrium.set('componentFractionVsElevationTableNames','{gasFraction,oilFraction,waterFraction}')
            equilibrium.set('phaseContacts',vector((5+offset,5+gap+offset)))
            functions=tree.find('Functions')
            for table in list(functions):
                if table.get('name')!='temperature':functions.remove(table)
            coordinates=vector(x+offset for x in (0,4.999,5,5+.9*gap,5+gap,10))
            for name,values in (('gasFraction','{0,0,0,0,1,1}'),('oilFraction','{0,0,1,1,0,0}'),('waterFraction','{1,1,0,0,0,0}')):
                ET.SubElement(functions,'TableFunction',name=name,coordinates=coordinates,values=values)
        gap=1e-6
        # Place a well-shaped cell centre inside the thin contact interval.
        # Sliver cells at 1e8 trigger the existing mesh face-area tolerance,
        # independently of initializer projection or contact selection.
        levels=(0,2,4,6+gap,8,10)
        _, thin_values=self.case('large-frame-resolved-micron-layer',nz=5,rotation=rotation,translation=translation,
                                z_coordinates=levels,mutate=lambda t:thin(t,gap))
        layer=[fields for key,fields in thin_values.items() if abs(key[2]-(5+.5*gap))<1e-6]
        self.assertEqual(len(layer),16)
        for fields in layer:self.assertAlmostEqual(fields['phaseVolumeFraction'][1],1.,delta=1e-12)
        _,log=self.case('large-frame-unresolved-layer',nz=5,rotation=rotation,translation=translation,
                        mutate=lambda t:thin(t,1e-8),success=False)
        self.assertIn('too close to resolve',log)

    def test_zero_gravity_constant_pressure(self):
        _, values = self.case('zero-gravity', gravity=(0.,0.,0.), nz=5)
        self.check_analytic(values, gravity=0.)

    def test_single_phase_and_thermal_coordinates(self):
        def single(tree, thermal=False, compressibility=0.):
            solvers = tree.find('Solvers')
            old = solvers.find('CompositionalMultiphaseFVM')
            replacement = ET.Element('SinglePhaseFVM', name='flow', logLevel='0', discretization='tpfa',
                                     targetRegions='{reservoir}', temperature='350', isThermal='1' if thermal else '0')
            for child in old: replacement.append(copy.deepcopy(child))
            solvers.remove(old); solvers.append(replacement)
            constitutive = tree.find('Constitutive')
            constitutive.remove(constitutive.find('InvariantImmiscibleFluid'))
            constitutive.remove(constitutive.find('BrooksCoreyRelativePermeability'))
            kind = 'ThermalCompressibleSinglePhaseFluid' if thermal else 'CompressibleSinglePhaseFluid'
            fluid = ET.SubElement(constitutive, kind, name='fluid', defaultDensity='1000', defaultViscosity='0.001',
                                  referenceDensity='1000', referencePressure='10000000', compressibility=str(compressibility))
            materials = '{fluid,rock}'
            if thermal:
                fluid.set('referenceTemperature','300'); fluid.set('thermalExpansionCoeff','0'); fluid.set('specificHeatCapacity','4180')
                constitutive.find('CompressibleSolidConstantPermeability').set('solidInternalEnergyModelName','solidEnergy')
                ET.SubElement(constitutive,'SolidInternalEnergy',name='solidEnergy',referenceTemperature='300',
                              referenceInternalEnergy='0',referenceVolumetricHeatCapacity='2000000')
                ET.SubElement(constitutive,'SinglePhaseThermalConductivity',name='conductivity',defaultThermalConductivityComponents='{1,1,1}')
                materials = '{fluid,rock,conductivity}'
            tree.find('./ElementRegions/CellElementRegion').set('materialList',materials)
            equilibrium = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            for key in ('phaseContacts','componentNames','componentFractionVsElevationTableNames'):
                equilibrium.attrib.pop(key)
            if not thermal: equilibrium.attrib.pop('temperatureVsElevationTableName')
        for thermal in (False, True):
            _, base = self.case('single-thermal' if thermal else 'single', mutate=lambda t:single(t,thermal))
            _, rotated = self.case('single-rotated', rotation=oblique_rotation(), mutate=lambda t:single(t,thermal))
            self.assertEqual(base.keys(),rotated.keys())
            for key in base:
                self.assertAlmostEqual(base[key]['pressure'][0],1e7+1000*9.81*(10-key[2]),delta=1e-4)
                self.assertAlmostEqual(base[key]['pressure'][0],rotated[key]['pressure'][0],delta=1e-6)
                if thermal:
                    self.assertAlmostEqual(rotated[key]['temperature'][0],300+5*key[2],delta=1e-9)

        _, compressible = self.case('single-compressible', mutate=lambda t:single(t,False,1e-8))
        _, rotated = self.case('single-compressible-rotated', rotation=oblique_rotation(), mutate=lambda t:single(t,False,1e-8))
        for key in compressible:
            exact = 1e7 - math.log1p(1e-8*1000*9.81*(key[2]-10))/1e-8
            self.assertAlmostEqual(compressible[key]['pressure'][0],exact,delta=1e-3)
            self.assertAlmostEqual(compressible[key]['pressure'][0],rotated[key]['pressure'][0],delta=1e-6)
        def nonconverging(tree):
            single(tree,False,1e-6)
            equilibrium=tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            equilibrium.set('maxNumberOfEquilibrationIterations','1')
            equilibrium.set('equilibrationTolerance','1e-10')
            equilibrium.set('elevationIncrementInHydrostaticPressureTable','10')
        _, log = self.case('single-nonconverging-refused', mutate=nonconverging, success=False)
        self.assertIn('failed to converge',log)

    def test_one_potential_coordinate_target(self):
        def one_level(tree):
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('datumElevation','5')
        _, values = self.case('one-level', nz=1, mutate=one_level)
        for fields in values.values():
            self.assertAlmostEqual(fields['pressure'][0],1e7,delta=1e-6)
            self.assertAlmostEqual(fields['temperature'][0],325.,delta=1e-9)
            self.assertAlmostEqual(fields['phaseVolumeFraction'][0],1.,delta=1e-12)

    def test_capillary_rotation_and_pressure_consistency(self):
        def capillary(tree):
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','0')
            region = tree.find('./ElementRegions/CellElementRegion')
            region.set('materialList','{fluid,rock,relperm,capillary}')
            ET.SubElement(tree.find('Constitutive'),'TableCapillaryPressure', name='capillary', phaseNames='{gas,water}', wettingNonWettingCapPressureTableName='pc')
            ET.SubElement(tree.find('Functions'),'TableFunction', name='pc', coordinates='{0,1}', values='{20000,0}')
        _, reference = self.case('capillary', mutate=capillary)
        _, rotated = self.case('capillary-oblique', rotation=oblique_rotation(), mutate=capillary)
        self.compare_fields(reference, rotated, tolerance=1e-6)
        for (_,_,z), fields in reference.items():
            pc = max(0., min(20000., (1000.-500.)*9.81*(z-5.)))
            expected_water = 1.-pc/20000.
            self.assertAlmostEqual(fields['phaseVolumeFraction'][1], expected_water, delta=2e-5)

    def test_mixed_solver_domains_and_second_mesh(self):
        def mixed(tree):
            mesh = tree.find('Mesh'); mesh.clear()
            ET.SubElement(mesh,'InternalMesh',name='mesh',elementTypes='{C3D8}',xCoords='{0,5,10}',yCoords='{0,10}',zCoords='{0,10}',nx='{2,2}',ny='{4}',nz='{4}',cellBlockNames='{leftCells,rightCells}')
            regions = tree.find('ElementRegions'); regions.clear()
            ET.SubElement(regions,'CellElementRegion',name='leftRegion',cellBlocks='{leftCells}',materialList='{fluid,rock,relperm}')
            ET.SubElement(regions,'CellElementRegion',name='rightRegion',cellBlocks='{rightCells}',materialList='{singleFluid,rock}')
            flow = tree.find('./Solvers/CompositionalMultiphaseFVM')
            flow.set('targetRegions','{leftRegion}')
            single = ET.SubElement(tree.find('Solvers'),'SinglePhaseFVM',name='single',discretization='tpfa',targetRegions='{rightRegion}',temperature='350')
            for child in flow: single.append(copy.deepcopy(child))
            ET.SubElement(tree.find('Constitutive'),'CompressibleSinglePhaseFluid',name='singleFluid',defaultDensity='900',defaultViscosity='0.001',referenceDensity='900',referencePressure='10000000',compressibility='0')
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('objectPath','ElementRegions/leftRegion')
            ET.SubElement(tree.find('FieldSpecifications'),'HydrostaticEquilibrium',name='singleEquilibrium',objectPath='ElementRegions/rightRegion',coordinateSystem='gravityAligned',datumElevation='10',datumPressure='10000000')
        _, serial = self.case('mixed-solver-domains',mutate=mixed)
        for (x,y,z),fields in serial.items():
            if x > 5: self.assertAlmostEqual(fields['pressure'][0],1e7+900*9.81*(10-z),delta=1e-4)
            else: self.check_analytic({(x,y,z):fields})
        if MPIEXEC:
            _, parallel = self.case('mixed-solver-domains-mpi',mutate=mixed,mpi=True)
            for key in serial: self.assertAlmostEqual(serial[key]['pressure'][0],parallel[key]['pressure'][0],delta=1e-6)
        def second_mesh(tree):
            mesh = tree.find('Mesh'); mesh.clear()
            for name,x,block in (('unused','{20,30}','unusedCells'),('active','{0,10}','activeCells')):
                ET.SubElement(mesh,'InternalMesh',name=name,elementTypes='{C3D8}',xCoords=x,yCoords='{0,10}',zCoords='{0,10}',nx='{4}',ny='{4}',nz='{4}',cellBlockNames='{'+block+'}')
            regions = tree.find('ElementRegions'); regions.clear()
            ET.SubElement(regions,'CellElementRegion',name='unusedRegion',meshBody='unused',cellBlocks='{unusedCells}',materialList='{rock}')
            ET.SubElement(regions,'CellElementRegion',name='reservoir',meshBody='active',cellBlocks='{activeCells}',materialList='{fluid,rock,relperm}')
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('targetRegions','{active/reservoir}')
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('objectPath','active/tpfa/ElementRegions/reservoir')
        _, values = self.case('second-mesh-only',mutate=second_mesh,expected_cells=128)
        self.check_analytic({key:value for key,value in values.items() if key[0] < 10})

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
            for mode in (None,'gravityAligned'):
                directory, current = self.case(prefix+'-corrected',mode=mode,nz=10,mutate=lambda t:configure(t,three,True))
                self.assertLess(maximum_phase_flux(current),1e-11)
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
            _, rotated = self.case(prefix+'-rotated',rotation=oblique_rotation(),nz=10,mutate=lambda t:configure(t,three))
            self.compare_fields(current,rotated,tolerance=1e-6)

    def test_capillary_support_and_numerical_refusals(self):
        def capillary(tree, dead_oil=False, initial_phase=False):
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','0')
            phases = '{gas,oil}' if dead_oil else '{gas,water}'
            tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock,relperm,capillary}')
            ET.SubElement(tree.find('Constitutive'),'TableCapillaryPressure',name='capillary',phaseNames=phases,wettingNonWettingCapPressureTableName='pc')
            ET.SubElement(tree.find('Functions'),'TableFunction',name='pc',coordinates='{0,1}',values='{0,20000}' if dead_oil else '{20000,0}')
            eq = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            if dead_oil:
                old = tree.find('./Constitutive/InvariantImmiscibleFluid')
                tree.find('Constitutive').remove(old)
                ET.SubElement(tree.find('Constitutive'),'DeadOilFluid',name='fluid',phaseNames=phases,surfaceDensities='{500,800}',componentMolarWeight='{0.044,0.018}',hydrocarbonFormationVolFactorTableNames='{Bg,Bo}',hydrocarbonViscosityTableNames='{vg,vo}')
                tree.find('./Constitutive/BrooksCoreyRelativePermeability').set('phaseNames',phases)
                eq.set('componentNames','{gas,oil}'); eq.set('equilibrationTolerance','0.001')
                for name,values in (('Bg','{1.27,1,0.7}'),('Bo','{1.09,1,0.9}'),('vg','{0.001,0.001,0.001}'),('vo','{0.001,0.001,0.001}')):
                    ET.SubElement(tree.find('Functions'),'TableFunction',name=name,coordinates='{1000000,10000000,20000000}',values=values)
            if initial_phase:
                eq.attrib.pop('phaseContacts'); eq.set('initialPhaseName','gas')
                tree.find("./Functions/TableFunction[@name='gasFraction']").set('values','{1,1,1,1,1}')
                tree.find("./Functions/TableFunction[@name='waterFraction']").set('values','{0,0,0,0,0}')
        features=json.loads(subprocess.run([GEOS,'--capabilities','--format=json'],capture_output=True,text=True,check=True).stdout)
        if features.get('self-consistent-capillary-initialization',False):
            self.case('pressure-dependent-capillary-supported',mutate=lambda t:capillary(t,dead_oil=True))
        else:
            _, log = self.case('pressure-dependent-capillary-refused',mutate=lambda t:capillary(t,dead_oil=True),success=False)
            self.assertIn('pressure-dependent EOS/capillary coupling',log)
        _, log = self.case('single-phase-capillary-refused',mutate=lambda t:capillary(t,initial_phase=True),success=False)
        self.assertIn('initialPhaseName with active capillarity is not supported',log)
        def unresolvable(tree, single=False, increment='0.00001'):
            eq = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            eq.set('elevationIncrementInHydrostaticPressureTable',increment)
            if single:
                solver = tree.find('./Solvers/CompositionalMultiphaseFVM')
                solver.tag='SinglePhaseFVM'; solver.attrib.pop('useMass')
                old=tree.find('./Constitutive/InvariantImmiscibleFluid'); tree.find('Constitutive').remove(old)
                ET.SubElement(tree.find('Constitutive'),'CompressibleSinglePhaseFluid',name='fluid',defaultDensity='1000',defaultViscosity='0.001',referenceDensity='1000',referencePressure='10000000',compressibility='0')
                tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock}')
                for key in ('phaseContacts','componentNames','componentFractionVsElevationTableNames','temperatureVsElevationTableName'): eq.attrib.pop(key)
        for single in (False, True):
            _, log = self.case('unresolvable-increment-refused',translation=(0.,0.,1e12),mutate=lambda t:unresolvable(t,single),success=False)
            self.assertIn('unresolvable',log)
            _, log = self.case('oversized-table-refused',mutate=lambda t:unresolvable(t,single,'1e-12'),success=False)
            self.assertIn('representable table size',log)

    def test_capillary_entry_pressure_and_inconsistent_mobile_endpoints(self):
        def capillary(tree, impossible=False):
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','0')
            tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock,relperm,capillary}')
            ET.SubElement(tree.find('Constitutive'),'TableCapillaryPressure',name='capillary',phaseNames='{gas,water}',wettingNonWettingCapPressureTableName='pc')
            ET.SubElement(tree.find('Functions'),'TableFunction',name='pc',coordinates='{0.1,1}' if impossible else '{0,1}',values='{1000,0}' if impossible else '{21000,1000}')
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('phaseContacts','{5.4}')
        for mode in (None,'gravityAligned'):
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

    def test_unsupported_obl_refusal(self):
        source = pathlib.Path(__file__).resolve().parents[1]/'inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_obl_3d.xml'
        directory = ROOT/'unsupported-obl'
        directory.mkdir()
        tree = ET.parse(source)
        tree.find('./Solvers/ReactiveCompositionalMultiphaseOBL').set('OBLOperatorsTableFile',str(source.parent/'obl_do_static.txt'))
        events = tree.find('Events')
        events.clear(); events.set('maxTime','1')
        ET.SubElement(events,'PeriodicEvent',name='output',timeFrequency='1',target='/Outputs/vtkOutput')
        ET.SubElement(tree.find('FieldSpecifications'),'HydrostaticEquilibrium',name='unsupportedAligned',objectPath='ElementRegions/Channel',coordinateSystem='gravityAligned',datumElevation='10',datumPressure='10000000')
        tree.write(directory/'input.xml',encoding='utf-8',xml_declaration=True)
        result = subprocess.run([GEOS,'-i',str(directory/'input.xml'),'-o',str(directory)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
        (directory/'run.log').write_text(result.stdout)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('does not implement gravityAligned',result.stdout)
        self.assertFalse(list(directory.glob('*.pvd')))

    def test_three_phase_interfaces(self):
        def three(tree):
            fluid = tree.find('./Constitutive/InvariantImmiscibleFluid')
            for key,value in {'componentNames':'{C0,C1,C2}','phaseNames':'{gas,oil,water}',
                              'densities':'{500,800,1000}','componentMolarWeight':'{0.044,0.114,0.018}',
                              'viscosities':'{0.001,0.001,0.001}'}.items(): fluid.set(key,value)
            relperm = tree.find('./Constitutive/BrooksCoreyRelativePermeability')
            for key,value in {'phaseNames':'{gas,oil,water}','phaseMinVolumeFraction':'{0,0,0}',
                              'phaseRelPermExponent':'{2,2,2}','phaseRelPermMaxValue':'{1,1,1}'}.items(): relperm.set(key,value)
            equilibrium = tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            equilibrium.set('componentNames','{C0,C1,C2}')
            equilibrium.set('componentFractionVsElevationTableNames','{gasFraction,oilFraction,waterFraction}')
            equilibrium.set('phaseContacts','{3,7}')
            functions = tree.find('Functions')
            for table in list(functions):
                if table.get('name') != 'temperature': functions.remove(table)
            for name,values in (('gasFraction','{0,0,0,0,1,1}'),('oilFraction','{0,0,1,1,0,0}'),('waterFraction','{1,1,0,0,0,0}')):
                ET.SubElement(functions,'TableFunction',name=name,coordinates='{0,2.9,3,6.9,7,10}',values=values)
        _, base = self.case('three-phase', nz=5, mutate=three)
        _, rotated = self.case('three-phase-oblique', nz=5, rotation=oblique_rotation(), mutate=three)
        self.compare_fields(base,rotated,tolerance=1e-6)
        for (_,_,z), fields in base.items():
            pressure = 1e7 + 500*9.81*(10-max(z,7)) + 800*9.81*max(0,7-max(z,3)) + 1000*9.81*max(0,3-z)
            self.assertAlmostEqual(fields['pressure'][0],pressure,delta=1e-4)
            phase = 2 if z < 3 else 1 if z < 7 else 0
            self.assertAlmostEqual(fields['phaseVolumeFraction'][phase],1.,delta=1e-12)
            self.assertAlmostEqual(fields['globalCompFraction'][phase],1.,delta=1e-12)

    def test_real_solver_run_and_mpi(self):
        directory, initial = self.case('flow-step', run_solver=True)
        final = records(directory, final=True)
        # The unmodified baseline has the same 1.25 mPa pressure correction.
        # Bound this by 1e-9 of datum pressure; fractions remain at machine scale.
        for key in initial:
            self.assertLessEqual(abs(initial[key]['pressure'][0]-final[key]['pressure'][0]), 1e-2)
            for field in ('temperature','globalCompFraction','phaseVolumeFraction',
                          'fluid_phaseMassDensity','fluid_phaseDensity'):
                for a,b in zip(initial[key][field],final[key][field]):
                    self.assertLessEqual(abs(a-b), 1e-10)
        baseline_directory, _ = self.case('legacy-flow-step', mode=None, baseline=True, run_solver=True)
        self.compare_fields(final, records(baseline_directory, final=True), tolerance=1e-7)
        if MPIEXEC:
            _, parallel = self.case('mpi-oblique', rotation=oblique_rotation(), mpi=True)
            self.compare_fields(initial, parallel, tolerance=1e-6)

    def test_selected_sets_and_refusals(self):
        def select(tree):
            ET.SubElement(tree, 'Geometry')
            ET.SubElement(tree.find('Geometry'),'Box',name='left',xMin='{0,0,0}',xMax='{5,10,10}')
            specs=tree.find('FieldSpecifications')
            specs.find('HydrostaticEquilibrium').set('setNames','{left}')
            for name, field, component, value in (('p','pressure','-1','9000000'),('x0','globalCompFraction','0','1'),('x1','globalCompFraction','1','0')):
                ET.SubElement(specs,'FieldSpecification',name=name,initialCondition='1',setNames='{all}',objectPath='ElementRegions/reservoir',fieldName=field,component=component,scale=value)
        _, selected = self.case('selected-left', mutate=select)
        self.check_analytic(selected, selected=True)
        def overlap(tree):
            select(tree)
            second=copy.deepcopy(tree.find('./FieldSpecifications/HydrostaticEquilibrium'))
            second.set('name','overlap'); second.set('setNames','{all}')
            tree.find('FieldSpecifications').append(second)
        _, log = self.case('overlap-refused', mutate=overlap, success=False)
        self.assertIn('Overlapping hydrostatic initializers',log)
        def missing(tree):
            select(tree)
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('setNames','{left,missing}')
        _, log = self.case('missing-set-refused', mutate=missing, success=False)
        self.assertTrue('missing' in log)
        def mixed_fluids(tree):
            mesh = tree.find('Mesh'); mesh.clear()
            ET.SubElement(mesh,'InternalMesh',name='mesh',elementTypes='{C3D8}',xCoords='{0,5,10}',yCoords='{0,10}',zCoords='{0,10}',
                          nx='{2,2}',ny='{4}',nz='{4}',cellBlockNames='{leftCells,rightCells}')
            regions=tree.find('ElementRegions'); regions.clear()
            ET.SubElement(regions,'CellElementRegion',name='leftRegion',cellBlocks='{leftCells}',materialList='{fluid,rock,relperm}')
            ET.SubElement(regions,'CellElementRegion',name='rightRegion',cellBlocks='{rightCells}',materialList='{fluidB,rock,relperm}')
            tree.find('./Solvers/CompositionalMultiphaseFVM').set('targetRegions','{leftRegion,rightRegion}')
            fluid=copy.deepcopy(tree.find('./Constitutive/InvariantImmiscibleFluid'))
            fluid.set('name','fluidB'); fluid.set('densities','{600,1200}'); tree.find('Constitutive').append(fluid)
            tree.find('./FieldSpecifications/HydrostaticEquilibrium').set('objectPath','ElementRegions/*')
        _, log = self.case('incompatible-fluids-refused', mutate=mixed_fluids, success=False)
        self.assertIn('different fluid models',log)
        def selected_model(tree):
            mixed_fluids(tree)
            ET.SubElement(tree,'Geometry')
            ET.SubElement(tree.find('Geometry'),'Box',name='left',xMin='{0,0,0}',xMax='{5,10,10}')
            specs = tree.find('FieldSpecifications')
            specs.find('HydrostaticEquilibrium').set('setNames','{left}')
            for name,field,component,value in (('p','pressure','-1','9000000'),('x0','globalCompFraction','0','1'),('x1','globalCompFraction','1','0')):
                ET.SubElement(specs,'FieldSpecification',name=name,initialCondition='1',setNames='{all}',objectPath='ElementRegions/*',fieldName=field,component=component,scale=value)
        _, selected_model_values = self.case('selected-model-union',mutate=selected_model)
        self.check_analytic(selected_model_values,selected=True)
        if MPIEXEC:
            _, parallel = self.case('selected-model-union-mpi',mutate=selected_model,mpi=True)
            for key in selected_model_values:
                for field in ('pressure','temperature','globalCompFraction','phaseVolumeFraction'):
                    for a,b in zip(selected_model_values[key][field],parallel[key][field]):
                        self.assertLessEqual(abs(a-b),1e-6)
            _, log = self.case('rank-split-models-refused',mutate=mixed_fluids,mpi=True,success=False)
            self.assertIn('different fluid models',log)
        def mixed_phase_single_declaration(tree):
            eq=tree.find('./FieldSpecifications/HydrostaticEquilibrium')
            eq.attrib.pop('phaseContacts'); eq.set('initialPhaseName','gas')
            tree.find("./Functions/TableFunction[@name='gasFraction']").set('values','{0.5,0.5,0.5,0.5,0.5}')
            tree.find("./Functions/TableFunction[@name='waterFraction']").set('values','{0.5,0.5,0.5,0.5,0.5}')
        _, log = self.case('single-declaration-multiphase-refused',mutate=mixed_phase_single_declaration,success=False)
        self.assertIn('multiple phases',log)
        self.case('legacy-subset-refused', mode=None, mutate=select, success=False)
        self.case('legacy-oblique-refused', mode=None, rotation=oblique_rotation(), success=False)
        self.case('unknown-coordinate-refused', mode='worldNormal', success=False)
        self.case('independent-plane-refused', mutate=lambda t:t.find('./FieldSpecifications/HydrostaticEquilibrium').set('contactNormal','{1,0,1}'), success=False)
        self.case('invalid-temperature-refused', mutate=lambda t:t.find("./Functions/TableFunction[@name='temperature']").set('values','{0,350}'), success=False)
        self.case('invalid-composition-refused', mutate=lambda t:t.find("./Functions/TableFunction[@name='gasFraction']").set('values','{0,0,2,1,1}'), success=False)
        self.case('mixed-phase-declarations-refused', mutate=lambda t:t.find('./FieldSpecifications/HydrostaticEquilibrium').set('initialPhaseName','gas'), success=False)
        self.case('invalid-increment-refused', mutate=lambda t:t.find('./FieldSpecifications/HydrostaticEquilibrium').set('elevationIncrementInHydrostaticPressureTable','0'), success=False)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geos',required=True)
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--mpiexec')
    parser.add_argument('--output',required=True,help='New directory for durable numerical evidence')
    args=parser.parse_args()
    GEOS=str(pathlib.Path(args.geos).absolute())
    BASELINE=str(pathlib.Path(args.baseline).absolute())
    MPIEXEC=args.mpiexec
    ROOT=pathlib.Path(args.output).absolute()
    ROOT.mkdir(parents=True,exist_ok=False)
    unittest.main(argv=[__file__],verbosity=2)
