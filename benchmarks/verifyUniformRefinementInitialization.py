#!/usr/bin/env python3
"""Run full FVM/FEM initialization, solve and restart checks on a sparse-ID mesh.

This correctness harness uses h5py/numpy to inspect checkpoints after GEOS exits;
those packages are test dependencies, not dependencies of the mesh generator.
"""

import argparse
import copy
import json
import os
import shlex
import shutil
import signal
import subprocess
from pathlib import Path
from xml.sax.saxutils import quoteattr
from xml.etree import ElementTree as ET


def make_inputs(directory, cells=(4, 2, 2), domain=(4, 2, 2)):
    nx, ny, nz = cells
    lx, ly, lz = domain
    coordinates = [(lx * x / nx, ly * y / ny, lz * z / nz) for z in range(nz + 1)
                   for y in range(ny + 1) for x in range(nx + 1)]

    def index(x, y, z):
        return (z * (ny + 1) + y) * (nx + 1) + x

    cells = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                cells.append([index(x, y, z), index(x + 1, y, z),
                              index(x + 1, y + 1, z), index(x, y + 1, z),
                              index(x, y, z + 1), index(x + 1, y, z + 1),
                              index(x + 1, y + 1, z + 1), index(x, y + 1, z + 1)])
    base = 9007199254741001
    points_text = " ".join(str(value) for xyz in coordinates for value in xyz)
    connectivity = " ".join(str(value) for cell in cells for value in cell)
    point_ids = " ".join(str(base + 7 * p) for p in range(len(coordinates)))
    cell_ids = " ".join(str(base + 1000 + 13 * c) for c in range(len(cells)))
    offsets = " ".join(str(8 * (c + 1)) for c in range(len(cells)))
    (directory / "coarse.vtu").write_text(f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid><Piece NumberOfPoints="{len(coordinates)}" NumberOfCells="{len(cells)}">
<PointData GlobalIds="pointIds"><DataArray type="Int64" IdType="1" Name="pointIds" format="ascii">{point_ids}</DataArray></PointData>
<CellData GlobalIds="cellIds"><DataArray type="Int64" IdType="1" Name="cellIds" format="ascii">{cell_ids}</DataArray><DataArray type="Int32" Name="attribute" format="ascii">{' '.join('1' for _ in cells)}</DataArray></CellData>
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{points_text}</DataArray></Points>
<Cells><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray><DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray><DataArray type="UInt8" Name="types" format="ascii">{' '.join('12' for _ in cells)}</DataArray></Cells>
</Piece></UnstructuredGrid></VTKFile>
''')
    geometry = f'''<Geometry>
<Box name="left" xMin="{{ -0.001, -0.001, -0.001 }}" xMax="{{ 0.001, {ly + 0.001}, {lz + 0.001} }}"/>
<Box name="right" xMin="{{ {lx - 0.001}, -0.001, -0.001 }}" xMax="{{ {lx + 0.001}, {ly + 0.001}, {lz + 0.001} }}"/>
<Box name="yZero" xMin="{{ -0.001, -0.001, -0.001 }}" xMax="{{ {lx + 0.001}, 0.001, {lz + 0.001} }}"/>
<Box name="zZero" xMin="{{ -0.001, -0.001, -0.001 }}" xMax="{{ {lx + 0.001}, {ly + 0.001}, 0.001 }}"/>
</Geometry>'''
    for kind in ("fvm", "fem"):
        for level in (0, 1, 2):
            solver = "flow" if kind == "fvm" else "solid"
            mesh = f'''<Mesh><VTKMesh name="mesh" file={quoteattr(str(directory / 'coarse.vtu'))} useGlobalIds="1"
scatterMethod="rcb" partitionRefinement="0" uniformRefinement="{level}"/></Mesh>'''
            events = f'''<Events maxTime="2"><PeriodicEvent name="solve" forceDt="1" target="/Solvers/{solver}"/>
<PeriodicEvent name="restart" timeFrequency="1" target="/Outputs/restart"/></Events>
<Outputs><Restart name="restart"/></Outputs>'''
            nonlinear = '<NonlinearSolverParameters newtonTol="1e-8" newtonMaxIter="8"/>'
            linear = '<LinearSolverParameters solverType="gmres" preconditionerType="amg" krylovTol="1e-10" krylovMaxIter="200"/>'
            if kind == "fvm":
                body = f'''<Solvers gravityVector="{{ 0, 0, 0 }}"><SinglePhaseFVM name="flow" logLevel="1" discretization="TPFA" targetRegions="{{ region }}">{nonlinear}{linear}</SinglePhaseFVM></Solvers>
<NumericalMethods><FiniteVolume><TwoPointFluxApproximation name="TPFA"/></FiniteVolume></NumericalMethods>
<ElementRegions><CellElementRegion name="region" cellBlocks="{{ 1_hexahedra }}" materialList="{{ water, rock }}"/></ElementRegions>
<Constitutive><CompressibleSinglePhaseFluid name="water" defaultDensity="1000" defaultViscosity="0.001" referencePressure="0" compressibility="1e-9" viscosibility="0"/><CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="porosity" permeabilityModelName="permeability"/><NullModel name="nullSolid"/><PressurePorosity name="porosity" defaultReferencePorosity="0.1" referencePressure="0" compressibility="1e-9"/><ConstantPermeability name="permeability" permeabilityComponents="{{ 1e-12, 1e-12, 1e-12 }}"/></Constitutive>
<FieldSpecifications><FieldSpecification name="initialPressure" initialCondition="1" objectPath="ElementRegions" fieldName="pressure" scale="1e5" setNames="{{ all }}"/><FieldSpecification name="leftPressure" objectPath="faceManager" fieldName="pressure" scale="1e5" setNames="{{ left }}"/><FieldSpecification name="rightPressure" objectPath="faceManager" fieldName="pressure" scale="0" setNames="{{ right }}"/></FieldSpecifications>'''
            else:
                constraints = "".join(
                    f'<FieldSpecification name="fix{d}" objectPath="nodeManager" fieldName="totalDisplacement" component="{d}" scale="0" setNames="{{ {name} }}"/>'
                    for d, name in enumerate(("left", "yZero", "zZero")))
                body = f'''<Solvers gravityVector="{{ 0, 0, 0 }}"><SolidMechanicsLagrangianFEM name="solid" logLevel="1" discretization="FE1" targetRegions="{{ region }}">{nonlinear}{linear}</SolidMechanicsLagrangianFEM></Solvers>
<NumericalMethods><FiniteElements><FiniteElementSpace name="FE1" order="1"/></FiniteElements></NumericalMethods>
<ElementRegions><CellElementRegion name="region" cellBlocks="{{ 1_hexahedra }}" materialList="{{ rock }}"/></ElementRegions>
<Constitutive><ElasticIsotropic name="rock" defaultDensity="2700" defaultBulkModulus="5.5556e9" defaultShearModulus="4.16667e9"/></Constitutive>
<FieldSpecifications>{constraints}<Traction name="load" objectPath="faceManager" scale="1e6" functionName="loadRamp" direction="{{ 1, 0, 0 }}" setNames="{{ right }}"/></FieldSpecifications>
<Functions><TableFunction name="loadRamp" inputVarNames="{{ time }}" coordinates="{{ 0, 1, 2 }}" values="{{ 0, 1, 2 }}"/></Functions>'''
            (directory / f"{kind}-L{level}.xml").write_text(
                f'<?xml version="1.0"?><Problem>{mesh}{events}{body}{geometry}</Problem>')


def inspect_checkpoint(file, kind, level, rank):
    import h5py
    import numpy as np
    base = "Problem/domain/MeshBodies/mesh/meshLevels/Level0/"
    sub = "ElementRegions/elementRegionsGroup/region/elementSubRegions/1_hexahedra/"
    result = {}
    with h5py.File(file) as data:
        def values(path):
            group = data[base + path]
            raw = group["__values__"][:]
            if "__dimensions__" not in group:
                return raw
            dims = tuple(int(i) for i in group["__dimensions__"][:])
            perm = tuple(int(i) for i in group["__permutation__"][:])
            return raw.reshape(tuple(dims[i] for i in perm)).transpose(
                tuple(perm.index(i) for i in range(len(perm))))

        assert int(data["Problem/Mesh/mesh/uniformRefinement/__values__"][0]) == level
        ghosts = values(sub + "ghostRank")
        if level:
            assert np.all(values(sub + "_geosUniformRootOwner")[ghosts < 0] == rank)
            assert np.all(values(sub + "_geosUniformGeneration") == level)
            assert np.all(values(sub + "_geosUniformRootCellId") > 2**53)
            coarse = np.floor(values(sub + "elementCenter")).astype(np.int64)
            coarse_index = (coarse[:, 2] * 2 + coarse[:, 1]) * 4 + coarse[:, 0]
            expected_roots = 9007199254741001 + 1000 + 13 * coarse_index
            np.testing.assert_array_equal(values(sub + "_geosUniformRootCellId"), expected_roots)
        else:
            assert base + sub + "_geosUniformRootCellId" not in data
        node_position = values("nodeManager/ReferencePosition")
        old_points = np.all(node_position == np.floor(node_position), axis=1)
        old_position = node_position[old_points].astype(np.int64)
        old_index = (old_position[:, 2] * 3 + old_position[:, 1]) * 5 + old_position[:, 0]
        np.testing.assert_array_equal(values("nodeManager/localToGlobalMap")[old_points],
                                      9007199254741001 + 7 * old_index)
        if kind == "fvm":
            position, field, owner = values(sub + "elementCenter"), values(sub + "pressure"), ghosts
        else:
            position = values("nodeManager/ReferencePosition")
            field = values("nodeManager/totalDisplacement")
            owner = values("nodeManager/ghostRank")
        if kind == "fem":
            bulk, shear = 5.5556e9, 4.16667e9
            young = 9 * bulk * shear / (3 * bulk + shear)
            poisson = (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear))
            load_scale = 1
            if "Problem/Functions/loadRamp" in data:
                load_scale = float(data["Problem/Events/time/__values__"][0])
                if int(data["Problem/Events/currentSubEvent/__values__"][0]) > 0:
                    load_scale += float(data["Problem/Events/dt/__values__"][0])
            expected = position * np.array([1, -poisson, -poisson]) * (load_scale * 1e6 / young)
            np.testing.assert_allclose(field, expected, rtol=1e-8, atol=1e-11,
                                       err_msg=f"Analytic uniaxial elastic patch: {file}")
        for i in np.flatnonzero(owner < 0):
            key = tuple(position[i])
            assert key not in result
            result[key] = field[i]
        return result, int(sum(ghosts < 0)), int(sum(ghosts >= 0))


def inspect_restart(reference, restarted, kind):
    """The mesh must be restored exactly; compare the state after continuing."""
    import h5py
    import numpy as np
    base = "Problem/domain/MeshBodies/mesh/meshLevels/Level0/"
    sub = "ElementRegions/elementRegionsGroup/region/elementSubRegions/1_hexahedra/"
    with h5py.File(reference) as a, h5py.File(restarted) as b:
        paths = ["nodeManager/ReferencePosition", "nodeManager/localToGlobalMap",
                 "nodeManager/ghostRank", sub + "localToGlobalMap", sub + "ghostRank"]
        paths += [sub + "_geosUniform" + name for name in (
            "RootCellId", "ParentCellId", "Generation", "ChildOrdinal",
            "RootOwner", "SourceType", "SourceAttribute")]
        for path in paths:
            np.testing.assert_array_equal(a[base + path + "/__values__"][:],
                                          b[base + path + "/__values__"][:], err_msg=path)
        path = sub + "pressure" if kind == "fvm" else "nodeManager/totalDisplacement"
        np.testing.assert_allclose(a[base + path + "/__values__"][:],
                                   b[base + path + "/__values__"][:], rtol=1e-12,
                                   atol=1e-8 if kind == "fvm" else 1e-13,
                                   err_msg=f"Restart continuation state: {restarted}")


def make_well_inputs(directory):
    paths = {}
    for level in (0, 1, 2):
        tree = ET.parse(directory / f"fvm-L{level}.xml")
        root = tree.getroot()
        root.find("Events").set("maxTime", "1")
        root.find("Events/PeriodicEvent").set("target", "/Solvers/coupled")
        flow = root.find("Solvers/SinglePhaseFVM")
        coupled = ET.SubElement(root.find("Solvers"), "SinglePhaseReservoir", name="coupled",
                                flowSolverName="flow", wellSolverName="wells", logLevel="1",
                                initialDt="1", targetRegions="{ region, wellA, wellB }")
        for parameters in list(flow):
            coupled.append(copy.deepcopy(parameters))
        linear = coupled.find("LinearSolverParameters")
        linear.set("solverType", "fgmres")
        linear.set("preconditionerType", "mgr")
        manager = ET.SubElement(root.find("Solvers"), "WellManager", name="wells", logLevel="1",
                                targetRegions="{ wellA, wellB }")
        for name, x, y in (("A", .25, .25), ("B", 3.75, 1.75)):
            region, controls = f"well{name}", f"controls{name}"
            ET.SubElement(root.find("ElementRegions"), "WellElementRegion", name=region, materialList="{ water }")
            well = ET.SubElement(root.find("Mesh/VTKMesh"), "InternalWell", name=f"generator{name}",
                                 wellRegionName=region, wellControlsName=controls, radius="0.01",
                                 numElementsPerSegment="2", polylineSegmentConn="{ { 0, 1 } }",
                                 polylineNodeCoords=f"{{ {{ {x}, {y}, 1.75 }}, {{ {x}, {y}, .25 }} }}")
            ET.SubElement(well, "Perforation", name="perforation", distanceFromHead="1.125")
            control = ET.SubElement(manager, "SinglePhaseWell", name=controls, type="injector", control="BHP")
            ET.SubElement(control, "MaximumBHPConstraint", name="bhp", targetBHP="2e5", referenceElevation="0")
            ET.SubElement(control, "InjectionVolumeRateConstraint", name="rate", volumeRate="1e-8")
        path = directory / f"wells-L{level}.xml"
        tree.write(path, encoding="utf-8", xml_declaration=True)
        paths[level] = path
    for kind, array_name, maximum in (("node", "pointIds", 2**63 - 3),
                                       ("element", "cellIds", 2**63 - 2)):
        mesh = ET.parse(directory / "coarse.vtu")
        array = mesh.find(f".//DataArray[@Name='{array_name}']")
        ids = array.text.split()
        ids[0] = str(maximum)
        array.text = " ".join(ids)
        file = directory / f"well-{kind}-overflow.vtu"
        mesh.write(file, encoding="utf-8", xml_declaration=True)
        tree = ET.parse(paths[0])
        tree.find("Mesh/VTKMesh").set("file", str(file))
        path = directory / f"well-{kind}-overflow.xml"
        tree.write(path, encoding="utf-8", xml_declaration=True)
        paths[kind] = path
    mesh = ET.parse(directory / "coarse.vtu")
    for array_name in ("pointIds", "cellIds"):
        array = mesh.find(f".//DataArray[@Name='{array_name}']")
        array.text = " ".join(str(i) for i in range(len(array.text.split())))
    compact_mesh = directory / "well-compact.vtu"
    mesh.write(compact_mesh, encoding="utf-8", xml_declaration=True)
    tree = ET.parse(paths[0])
    tree.find("Mesh/VTKMesh").set("file", str(compact_mesh))
    compact_input = directory / "well-compact.xml"
    tree.write(compact_input, encoding="utf-8", xml_declaration=True)
    paths["compact"] = compact_input
    return paths


def make_fracture_well_inputs(directory, well_inputs):
    """Split the middle plane into distinct contact sides; keep sparse surface IDs."""
    main = ET.parse(directory / "coarse.vtu")
    piece = main.find(".//Piece")
    coordinates = piece.find("Points/DataArray")
    values = [float(v) for v in coordinates.text.split()]
    points = [tuple(values[i:i + 3]) for i in range(0, len(values), 3)]
    original = len(points)
    point_ids = piece.find("PointData/DataArray")
    ids = [int(v) for v in point_ids.text.split()]
    plane = [p for p, xyz in enumerate(points) if xyz[0] == 2]
    duplicates = {p: original + i for i, p in enumerate(plane)}
    old_ids = ids.copy()
    points.extend(points[p] for p in plane)
    ids.extend(9007199254741001 + 7 * (original + i) for i in range(len(plane)))
    connectivity = piece.find("Cells/DataArray[@Name='connectivity']")
    cells = [int(v) for v in connectivity.text.split()]
    for offset in range(0, len(cells), 8):
        if min(points[p][0] for p in cells[offset:offset + 8]) >= 2:
            cells[offset:offset + 8] = [duplicates.get(p, p) for p in cells[offset:offset + 8]]
    coordinates.text = " ".join(str(v) for xyz in points for v in xyz)
    point_ids.text = " ".join(map(str, ids))
    connectivity.text = " ".join(map(str, cells))
    piece.set("NumberOfPoints", str(len(points)))
    piece.find("CellData/DataArray[@Name='cellIds']").text = " ".join(map(str, range(16)))
    main.write(directory / "fracture-main.vtu", encoding="utf-8", xml_declaration=True)
    index = {(points[p][1], points[p][2]): i for i, p in enumerate(plane)}
    faces = [[index[y, z], index[y + 1, z], index[y + 1, z + 1], index[y, z + 1]]
             for z in range(2) for y in range(2)]
    collocation = [v for p in plane for v in (old_ids[p], ids[duplicates[p]], -1)]
    (directory / "fracture-face.vtu").write_text(f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid>
<Piece NumberOfPoints="9" NumberOfCells="4">
<PointData GlobalIds="pointIds"><DataArray type="Int64" Name="pointIds" format="ascii">{' '.join(map(str, range(1000,1009)))}</DataArray>
<DataArray type="Int64" IdType="1" Name="collocated_nodes" NumberOfComponents="3" format="ascii">{' '.join(map(str, collocation))}</DataArray></PointData>
<CellData GlobalIds="cellIds"><DataArray type="Int64" Name="cellIds" format="ascii">100 101 102 103</DataArray></CellData>
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{' '.join(str(v) for p in plane for v in points[p])}</DataArray></Points>
<Cells><DataArray type="Int64" Name="connectivity" format="ascii">{' '.join(str(p) for face in faces for p in face)}</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">4 8 12 16</DataArray><DataArray type="UInt8" Name="types" format="ascii">9 9 9 9</DataArray></Cells>
</Piece></UnstructuredGrid></VTKFile>''')
    bundle = directory / "fracture-wells.vtm"
    bundle.write_text('''<?xml version="1.0"?><VTKFile type="vtkMultiBlockDataSet" version="1.0">
<vtkMultiBlockDataSet><DataSet name="main" file="fracture-main.vtu"/>
<DataSet name="fracture" file="fracture-face.vtu"/></vtkMultiBlockDataSet></VTKFile>''')
    result = {}
    for level in (0, 1, 2):
        tree = ET.parse(well_inputs[level])
        root = tree.getroot()
        mesh = root.find("Mesh/VTKMesh")
        mesh.set("file", str(bundle))
        mesh.set("faceBlocks", "{ fracture }")
        root.find("Solvers/SinglePhaseFVM").set("targetRegions", "{ region, fault }")
        root.find("Solvers/SinglePhaseReservoir").set("targetRegions", "{ region, fault, wellA, wellB }")
        ET.SubElement(root.find("ElementRegions"), "SurfaceElementRegion", name="fault", faceBlock="fracture",
                      defaultAperture="1e-4", materialList="{ water, fractureRock }")
        constitutive = root.find("Constitutive")
        ET.SubElement(constitutive, "CompressibleSolidParallelPlatesPermeability", name="fractureRock",
                      solidModelName="nullSolid", porosityModelName="fracturePorosity", permeabilityModelName="fracturePerm")
        ET.SubElement(constitutive, "PressurePorosity", name="fracturePorosity", defaultReferencePorosity="1",
                      referencePressure="0", compressibility="0")
        ET.SubElement(constitutive, "ParallelPlatesPermeability", name="fracturePerm")
        path = directory / f"fracture-wells-L{level}.xml"
        tree.write(path, encoding="utf-8", xml_declaration=True)
        result[level] = path
    return result


def inspect_fracture_well_ids(files, level):
    import h5py
    regions_path = "Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup"
    values = {name: set() for name in ("region", "fault", "wellA", "wellB")}
    for file in files:
        with h5py.File(file) as data:
            for name, target in values.items():
                for sub in data[regions_path + "/" + name + "/elementSubRegions"].values():
                    if not isinstance(sub, h5py.Group) or "localToGlobalMap" not in sub:
                        continue
                    ids = sub["localToGlobalMap/__values__"][:]
                    ghosts = sub["ghostRank/__values__"][:]
                    owned = set(int(v) for v in ids[ghosts < 0])
                    assert not owned & target
                    target.update(owned)
    assert len(values["region"]) == 16 * 8**level
    assert len(values["fault"]) == 4 * 4**level
    all_reservoir = values["region"] | values["fault"]
    assert len(all_reservoir) == len(values["region"]) + len(values["fault"])
    well_ids = values["wellA"] | values["wellB"]
    assert len(well_ids) == 4 and not well_ids & all_reservoir
    offset = max(all_reservoir) + 1 if level else 20
    assert well_ids == set(range(offset, offset + 4)), (level, well_ids, offset)
    return {"volume_cells": len(values["region"]), "surface_cells": len(values["fault"]),
            "minimum_well_element_id": min(well_ids), "maximum_reservoir_element_id": max(all_reservoir)}


def inspect_well_ids(files, level, compact=False):
    import h5py
    import numpy as np
    base = "Problem/domain/MeshBodies/mesh/meshLevels/Level0/"
    regions = base + "ElementRegions/elementRegionsGroup/"
    reservoir, wells, nodes = set(), {"wellA": set(), "wellB": set()}, set()
    for file in files:
        with h5py.File(file) as data:
            for name in ("region", "wellA", "wellB"):
                subregions = data[regions + name + "/elementSubRegions"]
                for subregion in subregions.values():
                    if not isinstance(subregion, h5py.Group) or "localToGlobalMap" not in subregion:
                        continue
                    ids = subregion["localToGlobalMap/__values__"][:]
                    ghost = subregion["ghostRank/__values__"][:]
                    target = reservoir if name == "region" else wells[name]
                    owned = set(int(i) for i in ids[ghost < 0])
                    assert not target & owned
                    target.update(owned)
            nodes.update(int(i) for i in data[base + "nodeManager/localToGlobalMap/__values__"][:])
    assert len(reservoir) == 16 * 8**level
    if compact:
        assert reservoir == set(range(16 * 8**level))
    assert all(len(ids) == 2 for ids in wells.values())
    assert not wells["wellA"] & wells["wellB"]
    well_ids = wells["wellA"] | wells["wellB"]
    offset = max(reservoir) + 1 if level else len(reservoir)
    assert min(well_ids) == offset
    assert well_ids == set(range(offset, offset + 4))
    volume_nodes = (4 * 2**level + 1) * (2 * 2**level + 1)**2
    volume_maximum = volume_nodes - 1 if compact else 9007199254741001 + 7 * 44 + volume_nodes - 45
    new_nodes = {i for i in nodes if i > volume_maximum}
    assert new_nodes == set(range(volume_maximum + 1, volume_maximum + 7))
    assert len(nodes) == volume_nodes + 6
    return {"well_elements": len(well_ids), "well_nodes": len(new_nodes),
            "minimum_well_element_id": min(well_ids), "minimum_well_node_id": min(new_nodes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path, help="New directory for inputs, logs and checkpoints")
    parser.add_argument("--launcher", default="mpiexec")
    parser.add_argument("--ranks", nargs="+", type=int, default=[1, 4])
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--wells", action="store_true", help="also solve two wells and check sparse/overflow ID ranges")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without Python -O; this correctness harness uses assertions")
    import numpy as np
    import h5py  # Fail before launching if checkpoint verification is unavailable.
    del h5py
    if 1 not in args.ranks or any(rank < 1 for rank in args.ranks):
        parser.error("--ranks must include 1 and contain only positive values")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    launcher = shlex.split(args.launcher)
    if not launcher:
        parser.error("--launcher must not be empty")
    executable = args.executable.resolve(strict=True)
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    make_inputs(directory)
    well_inputs = make_well_inputs(directory) if args.wells else {}
    fracture_well_inputs = make_fracture_well_inputs(directory, well_inputs) if args.wells else {}
    status = []

    def run(case, ranks, input_file, restart=None, expected_error=None):
        out = directory / case
        out.mkdir()
        command = [str(executable), "-i", str(input_file), "-o", str(out)]
        if restart:
            command += ["-r", str(restart)]
        if ranks > 1 or Path(launcher[0]).name == "srun":
            command = launcher + ["-n", str(ranks)] + command
        log_file = directory / f"{case}.log"
        with log_file.open("w") as log:
            process = subprocess.Popen(command, stdout=log, stderr=log, start_new_session=True)
            try:
                returncode = process.wait(timeout=args.timeout)
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                returncode = "timeout" if isinstance(error, subprocess.TimeoutExpired) else "interrupted"
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
        status.append({"case": case, "status": returncode, "command": command})
        (directory / "status.json").write_text(json.dumps(status, indent=2) + "\n")
        if expected_error:
            assert isinstance(returncode, int) and returncode != 0 and expected_error in log_file.read_text(), case
        else:
            assert returncode == 0, f"{case}: see {log_file} ({returncode})"

    def checkpoint(kind, level, ranks, rank, prefix=None, step=2):
        name = f"{kind}-L{level}"
        return (directory / (prefix or f"{name}-p{ranks}") /
                (name + f"_restart_{step:09}") / f"rank_{rank:07}.hdf5")

    fields = {}
    checks = []
    for ranks in sorted(set(args.ranks)):
        for kind in ("fvm", "fem"):
            for level in (0, 1, 2):
                case = f"{kind}-L{level}-p{ranks}"
                run(case, ranks, directory / f"{kind}-L{level}.xml")
                combined, owned, ghosts = {}, 0, 0
                for rank in range(ranks):
                    field, local, ghost = inspect_checkpoint(checkpoint(kind, level, ranks, rank), kind, level, rank)
                    assert not combined.keys() & field.keys(), case
                    combined.update(field)
                    owned += local
                    ghosts += ghost
                assert owned == 16 * 8**level, case
                expected = 16 * 8**level if kind == "fvm" else (4 * 2**level + 1) * (2 * 2**level + 1)**2
                assert len(combined) == expected, case
                difference = 0
                if ranks == 1:
                    fields[kind, level] = combined
                else:
                    reference = fields[kind, level]
                    assert combined.keys() == reference.keys(), case
                    a = np.asarray([reference[key] for key in sorted(reference)])
                    b = np.asarray([combined[key] for key in sorted(reference)])
                    np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-6 if kind == "fvm" else 1e-10)
                    difference = float(np.max(np.abs(a - b)))
                checks.append({"case": case, "owned_cells": owned, "ghost_cells": ghosts,
                               "unique_fields": len(combined), "maximum_solution_difference": difference})
        for kind in ("fvm", "fem"):
            # The restart event follows the solve. Cycle zero contains the
            # first step; resuming cycle one would only finish the event loop
            # after both solves. Resume zero to exercise an actual second solve.
            restart = directory / f"{kind}-L2-p{ranks}" / f"{kind}-L2_restart_000000000"
            run(f"{kind}-restart-p{ranks}", ranks, directory / f"{kind}-L2.xml", restart)
            assert "Successful nonlinear iterations" in (directory / f"{kind}-restart-p{ranks}.log").read_text()
            for rank in range(ranks):
                reference = checkpoint(kind, 2, ranks, rank)
                resumed = checkpoint(kind, 2, ranks, rank, prefix=f"{kind}-restart-p{ranks}")
                inspect_restart(reference, resumed, kind)
                inspect_checkpoint(resumed, kind, 2, rank)
            run(f"{kind}-changed-level-p{ranks}", ranks, directory / f"{kind}-L1.xml", restart,
                "Restart uniformRefinement differs from the input level count")
        # Simulate an older zero-level checkpoint lacking the new optional
        # wrapper. Copy only our fixture outputs, then remove that field.
        legacy = directory / f"legacy-zero-p{ranks}"
        shutil.copytree(directory / f"fvm-L0-p{ranks}", legacy)
        legacy_restart = legacy / "fvm-L0_restart_000000001"
        import h5py
        for rank in range(ranks):
            with h5py.File(legacy_restart / f"rank_{rank:07}.hdf5", "r+") as data:
                del data["Problem/Mesh/mesh/uniformRefinement"]
        run(f"fvm-legacy-zero-p{ranks}", ranks, directory / "fvm-L0.xml", legacy_restart)
        run(f"fvm-legacy-changed-level-p{ranks}", ranks, directory / "fvm-L1.xml", legacy_restart,
            "Restart uniformRefinement differs from the input level count")
        for rank in range(ranks):
            inspect_checkpoint(checkpoint("fvm", 0, ranks, rank,
                                           prefix=f"fvm-legacy-zero-p{ranks}"), "fvm", 0, rank)
        if args.wells:
            for level in (0, 1, 2):
                case = f"fracture-wells-L{level}-p{ranks}"
                run(case, ranks, fracture_well_inputs[level])
                files = [directory / case / f"fracture-wells-L{level}_restart_000000001" / f"rank_{rank:07}.hdf5"
                         for rank in range(ranks)]
                checks.append({"case": case, **inspect_fracture_well_ids(files, level)})
            for level in (0, 1, 2):
                case = f"wells-L{level}-p{ranks}"
                run(case, ranks, well_inputs[level])
                files = [directory / case / f"wells-L{level}_restart_000000001" / f"rank_{rank:07}.hdf5"
                         for rank in range(ranks)]
                for rank, file in enumerate(files):
                    inspect_checkpoint(file, "fvm", level, rank)
                checks.append({"case": case, **inspect_well_ids(files, level)})
            case = f"wells-compact-L0-p{ranks}"
            run(case, ranks, well_inputs["compact"])
            files = [directory / case / "well-compact_restart_000000001" / f"rank_{rank:07}.hdf5"
                     for rank in range(ranks)]
            checks.append({"case": case, **inspect_well_ids(files, 0, compact=True)})
            for kind in ("node", "element"):
                run(f"well-{kind}-overflow-p{ranks}", ranks, well_inputs[kind],
                    expected_error=f"Global well {kind} ID range exceeds storage" if kind == "node" else None)
                if kind == "element":
                    # Zero refinement retains develop's count-based offset;
                    # near-limit source IDs do not force the wells above them.
                    files = [directory / f"well-{kind}-overflow-p{ranks}" / "well-element-overflow_restart_000000001" /
                             f"rank_{rank:07}.hdf5" for rank in range(ranks)]
                    checks.append({"case": f"well-{kind}-overflow-p{ranks}", **inspect_well_ids(files, 0)})
            legacy = directory / f"legacy-well-width-p{ranks}"
            shutil.copytree(directory / f"wells-L1-p{ranks}", legacy)
            legacy_restart = legacy / "wells-L1_restart_000000000"
            for rank in range(ranks):
                with h5py.File(legacy_restart / f"rank_{rank:07}.hdf5", "r+") as data:
                    regions = data["Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup"]
                    for name in ("wellA", "wellB"):
                        group = regions[name]["elementSubRegions"][name + "UniqueSubRegion"]["nextWellElementIndexGlobal"]
                        old = group["__values__"][:]
                        assert np.all((old >= -(2**31)) & (old < 2**31))
                        del group["__values__"]
                        group.create_dataset("__values__", data=old.astype(np.int32))
            case = f"wells-legacy-width-p{ranks}"
            run(case, ranks, well_inputs[1], legacy_restart)
            files = [directory / case / "wells-L1_restart_000000001" / f"rank_{rank:07}.hdf5"
                     for rank in range(ranks)]
            checks.append({"case": case, **inspect_well_ids(files, 1)})
            for file in files:
                rank = int(file.stem.removeprefix("rank_"))
                reference = directory / f"wells-L1-p{ranks}" / "wells-L1_restart_000000001" / f"rank_{rank:07}.hdf5"
                with h5py.File(file) as data, h5py.File(reference) as original:
                    regions = data["Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup"]
                    for name in ("wellA", "wellB"):
                        values = regions[name]["elementSubRegions"][name + "UniqueSubRegion"]["nextWellElementIndexGlobal/__values__"]
                        assert values.dtype == np.int64
                        np.testing.assert_array_equal(values[:], original[values.name][:])
    (directory / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(f"Passed {len(status)} application runs; checkpoint checks: {directory / 'verification.json'}")


if __name__ == "__main__":
    main()
