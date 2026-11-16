#!/usr/bin/env python3
"""Standard-library regression checks for thermal input and PFW XML generation."""
import ast
import importlib.util
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
PFW = ROOT / "scripts/preProcessing/materialPointMethodParticleFileWriter"


def load_function(path, name, namespace=None):
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name)
    namespace = {} if namespace is None else namespace
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


class ThermalInputs(unittest.TestCase):
    def test_solver_dictionary_and_xml(self):
        formatter = load_function(PFW / "particleFileWriter.py", "_format_thermal_xml_block")
        value = {"defaultContactResistance": 0.2, "diagnostics": True, "mechanicalCoupling": "Coupled", "couplingTolerance": 1e-10,
                 "heatSources": [{"name": "heater", "particleIDs": [1, 2], "power": 0.5}],
                 "isothermal": [{"name": "bath", "temperature": 300}],
                 "contactResistance": [{"name": "pair", "materials": ["steel", "polymer"], "resistance": 0.3}]}
        root = ET.fromstring(formatter(value))
        self.assertEqual(root.get("enabled"), "1")
        self.assertEqual(root.get("diagnostics"), "1")
        self.assertEqual(root.get("mechanicalCoupling"), "Coupled")
        self.assertEqual(root.get("couplingTolerance"), "1e-10")
        self.assertEqual(root.find("HeatSource").get("particleIDs"), "{ 1, 2 }")
        self.assertEqual(root.find("ContactResistance").get("materials"), "{ steel, polymer }")
        explicit = ET.fromstring(formatter('<Thermal enabled="0" safetyFactor="0.2"/>'))
        self.assertEqual(explicit.get("enabled"), "0")
        self.assertEqual(explicit.get("safetyFactor"), "0.2")
        with self.assertRaises(ValueError):
            formatter("<HeatSource/>")

    def test_material_thermal_properties_survive_formatting(self):
        namespace = {"_required": lambda card, key: card[key], "_format_material_value": str}
        formatter = load_function(PFW / "pfw_materials.py", "_format_material_xml", namespace)
        card = {"name": "solid", "model": "LinearThermalElasticIsotropic", "defaultDensity": 1,
                "thermalSpecificHeat": 3, "thermalConductivity": 4, "thermalHeatCapacitySlope": 0.1,
                "thermalReferenceTemperature": 300, "defaultDrainedLinearTEC": 1e-5, "solverThermalExpansion": 1,
                "plasticHeatFraction": 0.8, "thermalYieldSoftening": 0.01}
        root = ET.fromstring(formatter(card, ["defaultDensity"]))
        self.assertEqual(root.tag, "LinearThermalElasticIsotropic")
        for key in card:
            if key not in ("model",):
                self.assertEqual(root.get(key), str(card[key]))

    def test_all_verification_cases(self):
        spec = importlib.util.spec_from_file_location("thermal_generate", Path(__file__).with_name("generate.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            for case in ("source", "bulk", "groups", "dfg", "materials", "isothermal", "insulated"):
                folder = Path(directory) / case
                xml = module.generate(folder, case, cells=8)
                root = ET.parse(xml).getroot()
                self.assertEqual(root.find(".//Thermal").get("enabled"), "1")
                rows = (folder / "particles.tsv").read_text().splitlines()
                count = len(rows[0].split("\t"))
                self.assertTrue(all(len(row.split("\t")) == count for row in rows[1:]))
                ids = [int(row.split("\t")[0]) for row in rows[1:]]
                self.assertEqual(ids, list(range(1, len(ids) + 1)))
                if case == "dfg":
                    self.assertTrue(all(row.split("\t")[6] == "0" for row in rows[1:]))
                    self.assertTrue(all(float(row.split("\t")[8]) > 0.99 for row in rows[1:]))
                if case == "materials":
                    self.assertEqual(len(root.findall(".//ThermalElasticIsotropic")), 2)


if __name__ == "__main__":
    unittest.main()
