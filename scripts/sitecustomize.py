"""Small runtime compatibility hook for the GEOS integrated-test runner.

The Open MPI ATS machine uses one setting for both the maximum MPI rank
capacity and the number of concurrently launched test jobs.  Sanitizer runs
need to keep the former high enough for the full test matrix while bounding
the latter.  When the runner exports GEOS_ITS_ATS_MAX_TESTS, cap only
``maxtests`` after the machine has parsed its normal rank capacity.  The
import hook is lazy so ATS can discover and register its machine modules first.
"""

import builtins
import importlib
import os


def _limit_openmpi_test_jobs() -> None:
    value = os.environ.get("GEOS_ITS_ATS_MAX_TESTS")
    if not value:
        return

    try:
        max_tests = int(value)
    except ValueError:
        return
    if max_tests < 1:
        return

    def patch_module(module) -> None:
        openmpi_machine = getattr(module, "OpenmpiMachine", None)
        if openmpi_machine is None or getattr(openmpi_machine, "_codex_ats_cap", False):
            return

        original_examine_options = openmpi_machine.examineOptions

        def examine_options(self, options):
            original_examine_options(self, options)
            self.maxtests = min(self.maxtests, max_tests)

        openmpi_machine.examineOptions = examine_options
        openmpi_machine._codex_ats_cap = True

    original_import_module = importlib.import_module

    def import_module(name, package=None):
        module = original_import_module(name, package)
        if name in (".machines.openmpi", "machines.openmpi") and package == "ats.atsMachines":
            patch_module(module)
        return module

    importlib.import_module = import_module

    original_import = builtins.__import__

    def import_hook(name, globals=None, locals=None, fromlist=(), level=0):
        module = original_import(name, globals, locals, fromlist, level)
        if name == "ats.atsMachines.machines.openmpi":
            patch_module(module)
        return module

    builtins.__import__ = import_hook


_limit_openmpi_test_jobs()
