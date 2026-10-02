#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only
"""Black-box pre-run contract regressions against an actual compiled GEOS.

Requires Python VTK for authoritative geometry/order comparisons. The normal
simulation reference is produced by the same executable, not a mock mesh.
Example: python3 scripts/testPreRunExports.py --geos /path/to/geosx --mpiexec mpirun
"""
import argparse
import hashlib
import json
import pathlib
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

import vtk

FIXTURE = """<Problem>
  <Mesh><InternalMesh name="mesh" elementTypes="{ C3D8 }"
    xCoords="{ 0, 0.5, 1 }" yCoords="{ 0, 1 }" zCoords="{ 0, 1 }"
    nx="{ 1, 1 }" ny="{ 2 }" nz="{ 2 }" cellBlockNames="{ left, right }"/></Mesh>
  <ElementRegions>
    <CellElementRegion name="leftRegion" cellBlocks="{ left }" materialList="{}"/>
    <CellElementRegion name="rightRegion" cellBlocks="{ right }" materialList="{}"/>
  </ElementRegions>
  <Events maxTime="1"><PeriodicEvent name="output" timeFrequency="1" target="/Outputs/vtk"/></Events>
  <Outputs><VTK name="vtk" plotFileRoot="result" writeGhostCells="1"/></Outputs>
</Problem>"""


def datasets(path):
    result = []

    def visit(node, names):
        if node.tag == "DataSet":
            result.append((tuple(names + [node.get("name")]), path.parent / node.get("file")))
        for child in node:
            visit(child, names + ([node.get("name")] if node.tag == "Block" else []))
    visit(ET.parse(path).getroot(), [])
    return result


def grid_signature(path, identity=False):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return grid_data_signature(reader.GetOutput(), identity)


def grid_data_signature(grid, identity=False):
    assert grid.GetNumberOfCells() > 0
    points = tuple(grid.GetPoint(i) for i in range(grid.GetNumberOfPoints()))
    cells = tuple((grid.GetCellType(i), tuple(grid.GetCell(i).GetPointId(j)
                   for j in range(grid.GetCell(i).GetNumberOfPoints())))
                  for i in range(grid.GetNumberOfCells()))
    if not identity:
        return points, cells
    assert grid.GetPoints().GetDataType() == vtk.VTK_DOUBLE, "authoritative coordinates must survive native reading as float64"
    assert grid.GetFieldData().GetNumberOfArrays() == 0, "pre-run output must have no TIME metadata"
    arrays = []
    for data in (grid.GetPointData(), grid.GetCellData()):
        assert {data.GetArrayName(i) for i in range(data.GetNumberOfArrays())} == {"ghostRank", "localToGlobalMap"}
        arrays.append(tuple((name, tuple(data.GetArray(name).GetTuple1(i)
                         for i in range(data.GetArray(name).GetNumberOfTuples())))
                            for name in ("ghostRank", "localToGlobalMap")))
    return points, cells, tuple(arrays)


def collection_signature(path):
    # Exercise the real composite reader, not only XML paths and individual VTU
    # readers: native clients must reconstruct all named blocks and cells.
    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    assert reader.GetErrorCode() == 0
    result = []
    def visit(block, names):
        assert block.GetNumberOfBlocks() > 0
        for index in range(block.GetNumberOfBlocks()):
            child = block.GetBlock(index)
            assert child is not None, (names, index)
            name = block.GetMetaData(index).Get(vtk.vtkCompositeDataSet.NAME())
            child_names = names + (name,)
            if isinstance(child, vtk.vtkMultiBlockDataSet):
                visit(child, child_names)
            else:
                assert isinstance(child, vtk.vtkUnstructuredGrid), child.GetClassName()
                result.append((child_names, grid_data_signature(child, identity=True)))
    visit(reader.GetOutput(), ())
    assert result
    return result


class PreRunExports(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="geos-prerun-")
        cls.root = pathlib.Path(cls.temp.name)
        cls.deck = cls.root / "input.xml"
        cls.deck.write_text(FIXTURE)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_geos(self, *args, success=True, mpi=False):
        command = [GEOS, *map(str, args)]
        if mpi:
            command = [MPIEXEC, "-np", "2", *command]
        completed = subprocess.run(command, cwd=self.root, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, timeout=180)
        if success:
            self.assertEqual(completed.returncode, 0, (completed.stdout + completed.stderr)[-5000:])
        else:
            self.assertNotEqual(completed.returncode, 0, (completed.stdout + completed.stderr)[-5000:])
        return completed

    def test_capabilities_and_bad_format(self):
        result = self.run_geos("--capabilities", "--format=json")
        capabilities = json.loads(result.stdout)
        self.assertTrue(capabilities["input-catalog"])
        self.assertTrue(capabilities["export-mesh"])
        self.assertEqual(capabilities["inputCatalogScope"], "global")
        result = self.run_geos("--capabilities", "--format=yaml", success=False)
        self.assertEqual(json.loads(result.stdout)["error"]["code"], "invalid-capabilities-arguments")

    def test_global_catalog_repeatability_and_units(self):
        paths = [self.root / f"catalog-{i}.json" for i in range(3)]
        self.run_geos("--input-catalog", paths[0])
        self.run_geos("--input-catalog", paths[1])
        self.run_geos("-i", self.deck, "--input-catalog", paths[2])
        self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
        self.assertEqual(paths[0].read_bytes(), paths[2].read_bytes())
        catalog = json.loads(paths[0].read_bytes())
        self.assertEqual(catalog["scope"], "global")
        self.assertEqual(catalog["mpiSize"], 1)
        elements = {e["path"]: e for e in catalog["elements"]}
        self.assertEqual(len(elements), len(catalog["elements"]))
        # None of these types occur in the input deck: this must remain global metadata.
        for path in ("/Problem/Constitutive/DeadOilFluid", "/Problem/FieldSpecifications/HydrostaticEquilibrium",
                     "/Problem/Solvers/CompositionalMultiphaseFVM/LinearSolverParameters"):
            self.assertIn(path, elements)
        schema_path = self.root / "reference.xsd"
        self.run_geos("-s", schema_path)
        ns = "{http://www.w3.org/2001/XMLSchema}"
        schema = ET.parse(schema_path).getroot()
        types = {t.get("name"): t for t in schema.findall(ns + "complexType")}
        for element in elements.values():
            target = types[element["schemaType"]]
            self.assertEqual({p["name"] for p in element["properties"]},
                             {p.get("name") for p in target.findall(ns + "attribute")})
            choice = target.find(ns + "choice")
            if choice is not None:
                self.assertEqual(element["childChoice"]["minOccurs"], int(choice.get("minOccurs", "1")))
                self.assertEqual(element["childChoice"]["maxOccurs"], choice.get("maxOccurs", "1"))
                expected = [(c.get("name"), c.get("type"), int(c.get("minOccurs", "1")), c.get("maxOccurs", "1"))
                            for c in choice.findall(ns + "element")]
                actual = [(c["name"], c["type"], c["minOccurs"], c["maxOccurs"]) for c in element["children"]]
                self.assertEqual(actual, expected)
        properties = [p for e in elements.values() for p in e["properties"]]
        self.assertTrue(any("gmres" in p["choices"] for p in properties))
        self.assertTrue(any("minimum" in p for p in properties))
        for element in elements.values():
            for prop in element["properties"]:
                self.assertEqual(prop["path"], element["path"] + "/@" + prop["name"])
                self.assertIn(prop["unitsStatus"], ("unknown", "declared"))
                if prop["unitsStatus"] == "unknown":
                    self.assertIsNone(prop["units"])
                else:
                    self.assertTrue(prop["units"])
        # Never overwrite an existing file, even when valid JSON already exists.
        digest = hashlib.sha256(paths[0].read_bytes()).digest()
        result = self.run_geos("--input-catalog", paths[0], success=False)
        diagnostics = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{"error":')]
        self.assertEqual(diagnostics[0]["error"]["code"], "pre-run-export-failed")
        self.assertEqual(digest, hashlib.sha256(paths[0].read_bytes()).digest())

    def test_invalid_deck_and_modes_publish_nothing(self):
        invalid = self.root / "invalid.xml"
        invalid.write_text('<Problem><Mesh><NotAGeosMesh name="bad"/></Mesh></Problem>')
        for flag, suffix in (("--input-catalog", ".json"), ("--export-mesh", ".vtm")):
            output = self.root / ("invalid" + suffix)
            self.run_geos("-i", invalid, flag, output, success=False)
            self.assertFalse(output.exists())
            self.assertFalse(pathlib.Path(str(output) + ".data").exists())
        for flag in ("--input-catalog", "--export-mesh"):
            self.run_geos(flag, success=False)
        for args in (("--export-mesh", self.root / "missing.vtm"),
                     ("-i", self.deck, "--export-mesh", self.root / "bad-extension"),
                     ("-i", self.deck, "-v", "--input-catalog", self.root / "mixed.json")):
            self.run_geos(*args, success=False)
            self.assertFalse(pathlib.Path(args[-1]).exists())

    def compare_mesh(self, mpi):
        suffix = "mpi" if mpi else "serial"
        mesh = self.root / f"mesh-{suffix}.vtm"
        output = self.root / f"normal-{suffix}"
        output.mkdir()
        partitions = ("-x", "2") if mpi else ()
        self.run_geos("-i", self.deck, *partitions, "--export-mesh", mesh, mpi=mpi)
        sidecars = pathlib.Path(str(mesh) + ".data")
        self.assertFalse(list(sidecars.rglob("*.pvd")))
        self.assertFalse(list(sidecars.rglob("*.vtm")))
        self.assertFalse(list(self.root.glob(f".{mesh.name}.partial-*")))
        metadata = json.loads((sidecars / "metadata.json").read_bytes())
        self.assertFalse(metadata["timeLoopEntered"])
        self.assertFalse(metadata["initialConditionsApplied"])
        self.run_geos("-i", self.deck, *partitions, "-o", output, mpi=mpi)
        exported = datasets(mesh)
        self.assertEqual(collection_signature(mesh),
                         [(name, grid_signature(path, identity=True)) for name, path in exported])
        normal = datasets(output / "result/000000.vtm")
        self.assertEqual([name for name, _ in exported], [name for name, _ in normal])
        for (_, a), (_, b) in zip(exported, normal):
            self.assertEqual(grid_signature(a), grid_signature(b))
            grid_signature(a, identity=True)
        before = mesh.read_bytes()
        self.run_geos("-i", self.deck, *partitions, "--export-mesh", mesh, success=False, mpi=mpi)
        self.assertEqual(before, mesh.read_bytes())
        again = self.root / f"mesh-repeat-{suffix}.vtm"
        self.run_geos("-i", self.deck, *partitions, "--export-mesh", again, mpi=mpi)
        self.assertEqual([grid_signature(p, identity=True) for _, p in exported],
                         [grid_signature(p, identity=True) for _, p in datasets(again)])

    def test_mesh_preserves_real64_coordinates(self):
        deck = self.root / "precision.xml"
        deck.write_text(FIXTURE.replace('xCoords="{ 0, 0.5, 1 }"',
                                       'xCoords="{ 100000000, 100000000.5, 100000001 }"'))
        output = self.root / "precision.vtm"
        self.run_geos("-i", deck, "--export-mesh", output)
        self.assertEqual(collection_signature(output),
                         [(name, grid_signature(path, identity=True)) for name, path in datasets(output)])
        xs = set()
        for _, path in datasets(output):
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(str(path))
            reader.Update()
            grid = reader.GetOutput()
            self.assertEqual(grid.GetPoints().GetDataType(), vtk.VTK_DOUBLE)
            xs.update(grid.GetPoint(i)[0] for i in range(grid.GetNumberOfPoints()))
        self.assertEqual(xs, {100000000.0, 100000000.5, 100000001.0})

    def test_mesh_matches_solver_serial(self):
        self.compare_mesh(False)

    def test_mesh_matches_solver_mpi(self):
        if not MPIEXEC:
            self.skipTest("Pass --mpiexec to exercise MPI geometry and ordering")
        self.compare_mesh(True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geos", required=True)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    GEOS = str(pathlib.Path(args.geos).absolute())
    MPIEXEC = args.mpiexec
    unittest.main(argv=[__file__], verbosity=2)
