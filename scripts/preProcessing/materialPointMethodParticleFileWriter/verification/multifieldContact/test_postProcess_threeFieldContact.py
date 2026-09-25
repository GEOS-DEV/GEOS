#!/usr/bin/env python3
"""Synthetic regression tests for the three-field contact postprocessor."""

from __future__ import annotations

import tempfile
import unittest
import csv
import json
from pathlib import Path

import postProcess_threeFieldContact as post


def synthetic_sample(spec: post.CaseSpec) -> post.Sample:
    if spec.family == "corner":
        projected = ((0.05, 0.05, 0.0), (0.05, 0.0, 0.0), (0.0, 0.05, 0.0))
    elif spec.family == "chain_explicit":
        projected = spec.pre_velocity
    else:
        projected = spec.expected_velocity or spec.pre_velocity
    force = tuple(
        tuple(spec.mass_ratio[body] * (projected[body][component] - spec.pre_velocity[body][component]) for component in range(3))
        for body in range(3)
    )
    if spec.family == "chain_explicit":
        position = ((-0.0625, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0625, 0.0, 0.0))
    else:
        position = ((0.003125, 0.0, 0.0), (0.003125, 0.0, 0.0), (-0.003125, 0.0, 0.0))
    return post.Sample(
        state=0,
        time=0.0,
        mass=spec.mass_ratio,
        active=(1.0, 1.0, 1.0),
        pre=spec.pre_velocity,
        post=projected,
        force=force,
        position=position,
        normal=((1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)),
    )


def write_history(path: Path, spec: post.CaseSpec, sample: post.Sample) -> None:
    row = {"state": sample.state, "time": sample.time}
    for physical_body, raw_field in enumerate(spec.body_to_field):
        row["gridMass_f{0}".format(raw_field)] = sample.mass[physical_body]
        row["gridActive_f{0}".format(raw_field)] = sample.active[physical_body]
        for base, vectors in (
            ("gridUncontactedVelocity", sample.pre),
            ("gridVelocity", sample.post),
            ("gridContactForce", sample.force),
            ("gridSurfacePosition", sample.position),
            ("gridSurfaceNormal", sample.normal),
        ):
            for component_index, component in enumerate(post.COMPONENTS):
                row["{0}_f{1}_{2}".format(base, raw_field, component)] = vectors[physical_body][component_index]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


class EvaluatorTests(unittest.TestCase):
    def test_all_published_oracles_classify_as_expected(self):
        results = []
        for spec in post.CASE_SPECS:
            result = post.evaluate_case(spec, [synthetic_sample(spec)], None, None, post.Tolerances())
            results.append(result)
            expected = "DIAGNOSTIC" if spec.diagnostic else "PASS"
            self.assertEqual(result.status, expected, (spec.key, [(c.name, c.status, c.actual) for c in result.checks]))
        post.compare_permutation(results, post.Tolerances())
        permuted = next(result for result in results if result.spec.key.endswith("Permuted"))
        self.assertEqual(permuted.status, "PASS")
        self.assertEqual(permuted.checks[-1].name, "Permutation equivalence to baseline")
        self.assertEqual(permuted.checks[-1].status, "PASS")

    def test_penetrating_post_velocity_fails_inequality(self):
        spec = post.SPEC_BY_KEY["threeField_sharedInterfaceOverlap"]
        sample = synthetic_sample(spec)
        broken = post.Sample(**{**sample.__dict__, "post": spec.pre_velocity})
        result = post.evaluate_case(spec, [broken], None, None, post.Tolerances())
        self.assertEqual(result.status, "FAIL")
        failed_names = {check.name for check in result.checks if check.status == "FAIL"}
        self.assertIn("Projected normal inequalities", failed_names)

    def test_no_contact_case_rejects_adhesive_impulse(self):
        spec = post.SPEC_BY_KEY["threeField_sharedInterfaceSeparating"]
        sample = synthetic_sample(spec)
        bad_post = ((-0.09, 0.0, 0.0), (-0.05, 0.0, 0.0), (0.095, 0.0, 0.0))
        broken = post.Sample(**{**sample.__dict__, "post": bad_post})
        result = post.evaluate_case(spec, [broken], None, None, post.Tolerances())
        self.assertEqual(result.status, "FAIL")
        self.assertTrue(any(check.name == "No unilateral contact impulse" and check.status == "FAIL" for check in result.checks))

    def test_coupled_solver_diagnostics_are_summarized(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "run.out").write_text(
                "CoupledContactFailureNodeDiagnostics solver=NewtonRaphson nodeGlobalID=8 "
                "activeFields=3 candidatePairs=2 activeConstraints=1 iterations=4 "
                "velocityResidual=2.0e-11 converged=0 regularizedSteps=1 lineSearchReductions=2\n"
                "CoupledContactFailureFieldDiagnostics solver=NewtonRaphson nodeGlobalID=8 "
                "field=0 mass=1.0 surfaceNormal=[1.0,0.0,0.0]\n"
                "CoupledContactFailurePairDiagnostics solver=NewtonRaphson nodeGlobalID=8 pairIndex=0 "
                "fieldA=0 fieldB=2 active=0 activation=inactive bilateral=0 hasGap=1 gap=0.0 "
                "normal=[9.0e-1,4.0e-1,0.0] coulombMargin=0.0\n"
                "CoupledContactFailurePairDiagnostics solver=NewtonRaphson nodeGlobalID=8 pairIndex=1 "
                "fieldA=1 fieldB=2 active=1 activation=gap bilateral=0 hasGap=1 gap=-1.0e-3 "
                "normal=[1.0,0.0,0.0] coulombMargin=-3.0e-12\n"
                "CoupledContactSolverDiagnostics solver=NewtonRaphson solvedNodes=1 "
                "nonconvergedNodes=1 maximumNodeIterations=4 maximumVelocityResidual=2.0e-11 "
                "velocityTolerance=1.0e-12 maximumRegularizedSteps=1 maximumLineSearchReductions=2 "
                "numericalGuardActivations=3 numericallySkippedNodes=0 numericallySkippedPairs=1 "
                "newtonToPGSFallbackNodes=1 solverRollbackNodes=0.\n"
            )
            failures, max_iterations, notes = post.scan_logs(root)
            self.assertEqual(max_iterations, 4)
            self.assertTrue(failures)
            self.assertTrue(
                any("Failure-detail log: nodes=1; fields=1; pairs=2" in note for note in notes)
            )
            self.assertTrue(
                any(
                    "max guard activations=3; skipped nodes=0; skipped pairs=1; "
                    "Newton-to-PGS fallbacks=1; solver rollbacks=0" in note
                    for note in notes
                )
            )

    def test_report_writers_emit_all_formats(self):
        spec = post.SPEC_BY_KEY["threeField_sharedInterfaceOverlap"]
        result = post.evaluate_case(spec, [synthetic_sample(spec)], None, None, post.Tolerances())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            post.write_csv([result], root / "checks.csv")
            post.write_markdown([result], "PASS", root / "report.md")
            post.write_pdf([result], "PASS", root / "report.pdf")
            post.write_tex([result], "PASS", root / "results.tex")
            for filename in ("checks.csv", "report.md", "report.pdf", "results.tex"):
                self.assertGreater((root / filename).stat().st_size, 100 if filename != "results.tex" else 20)

    def test_permuted_csv_is_mapped_back_to_physical_bodies(self):
        spec = post.SPEC_BY_KEY["threeField_sharedInterfacePermuted"]
        # Raw velocity fields are B, C, A because physical A/B/C use 2/0/1.
        raw_body = (1, 2, 0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "history.csv"
            row = {"state": 3, "time": 0.002}
            for raw_field, physical_body in enumerate(raw_body):
                row["gridMass_f{0}".format(raw_field)] = spec.mass_ratio[physical_body]
                row["gridActive_f{0}".format(raw_field)] = 1
                for base, vectors in (
                    ("gridUncontactedVelocity", spec.pre_velocity),
                    ("gridVelocity", spec.expected_velocity),
                    ("gridContactForce", ((0.0, 0.0, 0.0),) * 3),
                    ("gridSurfacePosition", ((1.0, 2.0, 3.0),) * 3),
                    ("gridSurfaceNormal", ((1.0, 0.0, 0.0),) * 3),
                ):
                    for component_index, component in enumerate(post.COMPONENTS):
                        row["{0}_f{1}_{2}".format(base, raw_field, component)] = vectors[physical_body][component_index]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
            sample = post.read_history(path, spec)[0]
            self.assertEqual(sample.mass, spec.mass_ratio)
            self.assertEqual(sample.pre, spec.pre_velocity)
            self.assertEqual(sample.post, spec.expected_velocity)

    def test_suite_cli_compiles_synthetic_histories(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)
            for spec in post.CASE_SPECS:
                write_history(
                    source / "output" / spec.case_id / "three_field_contact_history.csv",
                    spec,
                    synthetic_sample(spec),
                )
            report = source / "compiled_report"
            return_code = post.main(
                [
                    "--source-dir",
                    str(source),
                    "--output-dir",
                    str(report),
                    "--no-visit",
                    "--no-plots",
                ]
            )
            self.assertEqual(return_code, 0)
            payload = json.loads((report / "three_field_contact_summary.json").read_text())
            self.assertEqual(payload["overall_status"], "PASS")
            self.assertEqual(len(payload["cases"]), len(post.CASE_SPECS))


if __name__ == "__main__":
    unittest.main()
