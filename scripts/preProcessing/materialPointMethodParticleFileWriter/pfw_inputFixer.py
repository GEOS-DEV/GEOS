#!/usr/bin/env python3
"""Best-effort updater for obsolete GEOS-MPM particle-file-writer inputs.

The production particleFileWriter.py intentionally assumes that pfw_input files
are already current and explicit.  This helper is a separate migration utility:
run it on an older input file, review the generated *.fixed.py file, and then
run particleFileWriter.py on the fixed file.
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path


CONTACT_GAP = {"0": "Simple", "1": "Implicit", "2": "Softened"}
FTABLE_INTERP = {"0": "Linear", "1": "Cosine", "2": "Smoothstep"}


def _line_replacement(text: str, node: ast.AST, replacement: str) -> tuple[int, int, str]:
    lines = text.splitlines(keepends=True)
    start = node.lineno - 1
    end = getattr(node, "end_lineno", node.lineno)
    indent = re.match(r"\s*", lines[start]).group(0)
    rep_lines = []
    for line in replacement.rstrip("\n").splitlines():
        rep_lines.append(indent + line.lstrip() + "\n")
    return start, end, "".join(rep_lines)


def _expand_solver_attributes(text: str, changes: list[str]) -> str:
    try:
        tree = ast.parse(text)
    except SyntaxError:
        changes.append("could not parse Python AST; skipped solver_attributes expansion")
        return text

    generated_assignments: list[str] = []
    removals: list[tuple[int, int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "solver_attributes" for target in node.targets):
            continue
        if not isinstance(node.value, ast.Dict):
            changes.append("found solver_attributes, but it is not a literal dictionary; left it for manual review")
            continue

        ok = True
        local_assignments = ["# pfw_inputFixer: expanded obsolete solver_attributes into direct PFW keys."]
        for key_node, value_node in zip(node.value.keys, node.value.values):
            if not isinstance(key_node, ast.Constant) or not isinstance(key_node.value, str):
                ok = False
                break
            value_source = ast.get_source_segment(text, value_node)
            if value_source is None:
                ok = False
                break
            local_assignments.append(f'pfw[{key_node.value!r}] = {value_source}')
        if not ok:
            changes.append("found solver_attributes with non-literal keys; left it for manual review")
            continue
        generated_assignments.extend(local_assignments)
        removals.append(_line_replacement(text, node, "# pfw_inputFixer: removed obsolete solver_attributes dictionary."))

    if not generated_assignments:
        return text

    lines = text.splitlines(keepends=True)
    for start, end, replacement in sorted(removals, reverse=True):
        lines[start:end] = [replacement]
    text = "".join(lines)

    # Remove references to the obsolete dictionary from metadata or pfw itself.
    text = re.sub(r"^\s*['\"]solverAttributes['\"]\s*:\s*solver_attributes\s*,?\s*\n", "", text, flags=re.M)
    text = re.sub(r"^\s*['\"]solver_attributes['\"]\s*:\s*solver_attributes\s*,?\s*\n", "", text, flags=re.M)
    text = re.sub(r"^\s*pfw\[['\"]solverAttributes['\"]\]\s*=\s*solver_attributes\s*\n", "", text, flags=re.M)
    text = re.sub(r"^\s*pfw\[['\"]solver_attributes['\"]\]\s*=\s*solver_attributes\s*\n", "", text, flags=re.M)

    lines = text.splitlines(keepends=True)
    insert_after = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*pfw\s*=\s*", line):
            insert_after = i + 1
            break
    assignment_text = "\n".join(generated_assignments) + "\n"
    if insert_after is None:
        lines.insert(0, "pfw = {}\n" + assignment_text)
    else:
        lines.insert(insert_after, assignment_text)
    text = "".join(lines)
    changes.append("expanded solver_attributes into direct pfw assignments")
    return text


def _replace_assignment_enum(text: str, key: str, mapping: dict[str, str], changes: list[str]) -> str:
    pattern = re.compile(
        r"(?P<prefix>pfw\[\s*['\"]" + re.escape(key) + r"['\"]\s*\]\s*=\s*)(?P<value>[012])(?P<suffix>\s*(?:#.*)?$)",
        flags=re.M,
    )

    def repl(match: re.Match[str]) -> str:
        old_value = match.group("value")
        new_value = mapping[old_value]
        changes.append(f"converted {key}={old_value} to {new_value!r}")
        return f'{match.group("prefix")}{new_value!r}{match.group("suffix")}'

    return pattern.sub(repl, text)


def _fix_events(text: str, changes: list[str]) -> str:
    event_names = (
        "Anneal", "BodyForceUpdate", "BoreholePressure", "CohesiveZone",
        "ConfiningPressure", "CrystalHeal", "DeformationUpdate",
        "FrictionCoefficientSwap", "Heal", "InitializeStress",
        "InsertPeriodicContactSurfaces", "MachineSample", "MaterialSwap",
        "PolymerHeal", "ResetDeformationGradient", "TemperatureProfile",
        "TransformParticles",
    )
    event_pattern = re.compile(r"<(?P<tag>" + "|".join(event_names) + r")\b(?P<attrs>[^>]*)>", re.S)

    def repl(match: re.Match[str]) -> str:
        tag = match.group("tag")
        attrs = match.group("attrs")
        original = attrs
        attrs = re.sub(r"\btime\s*=", "startTime=", attrs)
        attrs = re.sub(r"\binterval\s*=", "endTime=", attrs)
        if attrs != original:
            changes.append(f"renamed obsolete event attributes on {tag}")
        return f"<{tag}{attrs}>"

    return event_pattern.sub(repl, text)


def fix_text(text: str) -> tuple[str, list[str]]:
    changes: list[str] = []

    replacements = [
        ("planeSrain", "planeStrain", "renamed planeSrain to planeStrain"),
        ("planeStrian", "planeStrain", "renamed planeStrian to planeStrain"),
        ("ShockEOSVonMisesJ", "MetalShock", "renamed ShockEOSVonMisesJ to MetalShock"),
        ("ShockEOSvonMisesJ", "MetalShock", "renamed ShockEOSvonMisesJ to MetalShock"),
        ("CohesiveZoneReference", "CohesiveZone", "renamed CohesiveZoneReference to CohesiveZone"),
    ]
    for old, new, label in replacements:
        if old in text:
            text = text.replace(old, new)
            changes.append(label)

    before = text
    text = re.sub(r"^\s*pfw\[\s*['\"]useConstantTimeStep['\"]\s*\]\s*=.*\n", "", text, flags=re.M)
    text = re.sub(r"^\s*pfw\[\s*['\"]constantTimeStepValue['\"]\s*\]\s*=.*\n", "", text, flags=re.M)
    if text != before:
        changes.append("removed obsolete constant-time-step PFW controls")

    text = _replace_assignment_enum(text, "contactGapCorrection", CONTACT_GAP, changes)
    text = _replace_assignment_enum(text, "fTableInterpType", FTABLE_INTERP, changes)
    text = _expand_solver_attributes(text, changes)
    text = _fix_events(text, changes)

    header = (
        "# This file was generated by pfw_inputFixer.py.\n"
        "# Review the changes before using it for production calculations.\n"
    )
    if not text.startswith("# This file was generated by pfw_inputFixer.py"):
        text = header + text
    return text, changes


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a best-effort .fixed.py version of an obsolete GEOS-MPM pfw_input file")
    parser.add_argument("input", type=Path, help="input pfw_input_*.py file")
    parser.add_argument("-o", "--output", type=Path, help="output path. Default: <input-stem>.fixed.py")
    parser.add_argument("--dry-run", action="store_true", help="print planned changes without writing the fixed file")
    args = parser.parse_args()

    source = args.input
    text = source.read_text()
    fixed, changes = fix_text(text)
    output = args.output or source.with_name(source.stem + ".fixed" + source.suffix)

    print("pfw_inputFixer summary:")
    if changes:
        for change in changes:
            print("  - " + change)
    else:
        print("  - no automatic changes were required")

    if args.dry_run:
        print(f"dry run: would write {output}")
        return 0

    output.write_text(fixed)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
