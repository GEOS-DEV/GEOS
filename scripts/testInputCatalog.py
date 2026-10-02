#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only
"""Independent regression for global input metadata/capability discovery.
Run against the metadata-only PR or a later build. Mesh export is not required.
"""
import argparse
import hashlib
import json
import pathlib
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET


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



class InputCatalog(unittest.TestCase):

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


    def test_rejected_input_and_modes(self):
        invalid = self.root / "invalid.xml"
        invalid.write_text('<Problem><Mesh><NotAGeosMesh name="bad"/></Mesh></Problem>')
        output = self.root / "invalid.json"
        self.run_geos("-i", invalid, "--input-catalog", output, success=False)
        self.assertFalse(output.exists())
        self.run_geos("--input-catalog", success=False)
        self.run_geos("-i", self.deck, "-v", "--input-catalog", output, success=False)
        self.assertFalse(output.exists())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geos", required=True)
    args = parser.parse_args()
    GEOS = str(pathlib.Path(args.geos).absolute())
    MPIEXEC = None
    unittest.main(argv=[__file__], verbosity=2)
