#!/usr/bin/env python3
"""Exercise the actual octant input's mode parser, without importing MPI/PFW.

Only the two environment helper definitions are extracted via AST. A separate
full-writer smoke test is needed to validate XML/particle output. No GEOS is run.
"""
from pathlib import Path
import ast
import os
import unittest
from unittest.mock import patch

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / "scripts/preProcessing/materialPointMethodParticleFileWriter").is_dir())
INPUT = (ROOT / "scripts/preProcessing/materialPointMethodParticleFileWriter/verification/"
         "weakTraceSphericalInclusionOctant/pfw_input_weakTraceSphericalInclusionOctant.py")
TREE = ast.parse(INPUT.read_text(), filename=str(INPUT))
DEFINITIONS = [node for node in TREE.body if isinstance(node, ast.FunctionDef)
               and node.name in ("bool_env", "marker_search_mode_env")]
if len(DEFINITIONS) != 2:
    raise RuntimeError("Expected the actual bool and search-mode environment helpers")
SCOPE = {"os": os}
exec(compile(ast.Module(body=DEFINITIONS, type_ignores=[]), str(INPUT), "exec"), SCOPE)


class OctantSearchOptionTest(unittest.TestCase):
    def mode(self, value=None):
        with patch.dict(os.environ, {}, clear=True):
            if value is not None:
                os.environ["WEAK_TRACE_FAST_MARKER_SEARCH"] = value
            return SCOPE["marker_search_mode_env"]()

    def test_default_is_indexed(self):
        self.assertEqual(self.mode(), 1)

    def test_existing_boolean_spellings(self):
        for token in ("1", "true", "yes", "on", " TRUE "):
            with self.subTest(token=token):
                self.assertEqual(self.mode(token), 1)
        for token in ("0", "false", "no", "off", " OFF "):
            with self.subTest(token=token):
                self.assertEqual(self.mode(token), 0)

    def test_validation_mode(self):
        for token in ("2", " 2 "):
            with self.subTest(token=token):
                self.assertEqual(self.mode(token), 2)

    def test_invalid_modes_rejected(self):
        for token in ("", "3", "-1", "2.0", "auto"):
            with self.subTest(token=token):
                with self.assertRaises(ValueError):
                    self.mode(token)

    def test_solver_keyword_uses_mode_parser(self):
        values = [kw.value for node in ast.walk(TREE) if isinstance(node, ast.Call)
                  for kw in node.keywords if kw.arg == "weakInterfaceTraceFastMarkerSearch"]
        self.assertEqual(len(values), 1)
        self.assertIsInstance(values[0], ast.Call)
        self.assertIsInstance(values[0].func, ast.Name)
        self.assertEqual(values[0].func.id, "marker_search_mode_env")


if __name__ == "__main__":
    unittest.main(verbosity=2)
