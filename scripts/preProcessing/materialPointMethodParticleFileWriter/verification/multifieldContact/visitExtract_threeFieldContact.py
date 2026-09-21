#!/usr/bin/env python3
"""Extract the three-field junction-node history from a GEOS Silo database.

Run this file with VisIt's Python interpreter, not ordinary CPython::

    visit -nowin -cli -s visitExtract_threeFieldContact.py \
        --run-dir <case-run-directory> --output <history.csv>

The resulting CSV is intentionally solver-agnostic.  The companion
``postProcess_threeFieldContact.py`` performs all pass/fail decisions.
"""

from __future__ import print_function

import argparse
import csv
import json
import re
import sys
from pathlib import Path


VECTOR_FIELDS = (
    "gridVelocity",
    "gridUncontactedVelocity",
    "gridContactForce",
    "gridSurfacePosition",
    "gridSurfaceNormal",
)
SCALAR_FIELDS = ("gridMass", "gridActive")
COMPONENTS = ("x", "y", "z")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=".")
    parser.add_argument("--output", default="three_field_contact_history.csv")
    parser.add_argument("--metadata-output", default="")
    parser.add_argument("--database", default="")
    parser.add_argument("--pick", default="0,0,0", help="junction coordinate x,y,z")
    return parser.parse_known_args()[0]


def call_visit(name, *args):
    fn = globals().get(name)
    if fn is None:
        import __main__

        fn = getattr(__main__, name)
    return fn(*args)


def call_visit_keywords(name, **kwargs):
    fn = globals().get(name)
    if fn is None:
        import __main__

        fn = getattr(__main__, name)
    return fn(**kwargs)


def compact(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def metadata_entries(metadata, count_name, get_name):
    result = []
    count_fn = getattr(metadata, count_name, None)
    get_fn = getattr(metadata, get_name, None)
    if count_fn is None or get_fn is None:
        return result
    try:
        count = int(count_fn())
    except Exception:
        return result
    for index in range(count):
        try:
            value = get_fn(index)
            name = str(value.name)
        except Exception:
            continue
        if name:
            result.append(name)
    return result


def discover_database(run_dir, explicit=""):
    if explicit:
        database = str(Path(explicit).expanduser().resolve())
        if not call_visit("OpenDatabase", database):
            raise RuntimeError("could not open database: " + database)
        return database

    silo_dir = run_dir / "siloFiles"
    matches = sorted(path for path in silo_dir.glob("mpm_cpdi_*") if path.is_file()) or sorted(
        path for path in silo_dir.glob("mpm_*") if path.is_file()
    )
    if not matches:
        raise FileNotFoundError("no GEOS Silo files found under " + str(silo_dir))

    if len(matches) == 1:
        database = str(matches[0].resolve())
    else:
        family = run_dir / "three_field_contact.visit"
        family.write_text("\n".join(str(path.resolve()) for path in matches) + "\n")
        database = str(family.resolve())
    if not call_visit("OpenDatabase", database):
        raise RuntimeError("VisIt could not open " + database)
    return database


def resolve_field(names, base, field_index):
    """Resolve a Silo variable while tolerating folder/name variations."""
    wanted_base = compact(base)
    wanted_suffixes = (
        "velocityfield{0}".format(field_index + 1),
        "field{0}".format(field_index + 1),
    )
    candidates = []
    for name in names:
        cname = compact(name)
        if wanted_base not in cname:
            continue
        suffix_score = 0
        for suffix in wanted_suffixes:
            if cname.endswith(suffix):
                suffix_score = 100
                break
        if not suffix_score:
            continue
        score = suffix_score
        if "cellregion1gridfields" in cname:
            score += 30
        if "particle" in cname:
            score -= 100
        if cname.endswith(wanted_base + wanted_suffixes[0]):
            score += 10
        candidates.append((score, len(name), name))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    return candidates[0][2]


def scalar_value(value):
    if isinstance(value, (int, float)):
        return float(value)
    try:
        values = list(value)
    except Exception:
        values = []
    for item in values:
        try:
            return float(item)
        except Exception:
            pass
    try:
        return float(value)
    except Exception:
        return float("nan")


def current_time(state):
    try:
        call_visit("Query", "Time")
        return scalar_value(call_visit("GetQueryOutputValue"))
    except Exception:
        return float(state)


def make_columns(resolved):
    columns = ["state", "time"]
    pick_variables = []
    expression_names = []
    for field_index in range(3):
        for base in SCALAR_FIELDS:
            label = "{0}_f{1}".format(base, field_index)
            columns.append(label)
            if resolved.get(label):
                pick_variables.append(resolved[label])
        for base in VECTOR_FIELDS:
            vector_key = "{0}_f{1}".format(base, field_index)
            vector = resolved.get(vector_key)
            for component_index, component in enumerate(COMPONENTS):
                label = "{0}_{1}".format(vector_key, component)
                columns.append(label)
                if not vector:
                    continue
                expression = "tf_{0}_f{1}_{2}".format(base.lower(), field_index, component)
                call_visit(
                    "DefineScalarExpression",
                    expression,
                    "array_decompose(<{0}>,{1})".format(vector, component_index),
                )
                resolved[label] = expression
                expression_names.append(expression)
                pick_variables.append(expression)
    return columns, pick_variables, expression_names


def main():
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    output = Path(args.output).expanduser()
    if not output.is_absolute():
        output = (run_dir / output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output = Path(args.metadata_output).expanduser() if args.metadata_output else output.with_suffix(".metadata.json")
    if not metadata_output.is_absolute():
        metadata_output = (run_dir / metadata_output).resolve()

    pick = tuple(float(part.strip()) for part in args.pick.split(","))
    if len(pick) != 3:
        raise ValueError("--pick must contain x,y,z")

    database = discover_database(run_dir, args.database)
    metadata = call_visit("GetMetaData", database)
    scalars = metadata_entries(metadata, "GetNumScalars", "GetScalars")
    vectors = metadata_entries(metadata, "GetNumVectors", "GetVectors")

    resolved = {}
    for field_index in range(3):
        for base in SCALAR_FIELDS:
            resolved["{0}_f{1}".format(base, field_index)] = resolve_field(scalars, base, field_index)
        for base in VECTOR_FIELDS:
            resolved["{0}_f{1}".format(base, field_index)] = resolve_field(vectors, base, field_index)

    columns, pick_variables, expression_names = make_columns(resolved)
    mass_plot = resolved.get("gridMass_f0")
    if not mass_plot:
        raise RuntimeError("gridMass_velocityField1 was not found; cannot select the grid mesh")
    call_visit("AddPlot", "Pseudocolor", mass_plot)
    call_visit("DrawPlots")

    try:
        state_count = int(call_visit("TimeSliderGetNStates"))
    except Exception:
        state_count = 1

    rows = []
    for state in range(max(1, state_count)):
        if state_count > 1:
            call_visit("SetTimeSliderState", state)
        try:
            picked = call_visit_keywords("NodePick", coord=pick, vars=tuple(pick_variables))
        except TypeError:
            # Older VisIt wrappers only accept positional arguments here.
            picked = call_visit("NodePick", pick, tuple(pick_variables))
        row = {column: float("nan") for column in columns}
        row["state"] = state
        row["time"] = current_time(state)
        for label, variable in resolved.items():
            if variable and label in row:
                row[label] = scalar_value(picked.get(variable))
        rows.append(row)
        print("Extracted junction node for state {0}/{1}".format(state + 1, state_count))

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    missing = sorted(key for key, value in resolved.items() if "_f" in key and not value)
    metadata_payload = {
        "database": database,
        "junction_coordinate": list(pick),
        "state_count": len(rows),
        "resolved_variables": resolved,
        "missing_variables": missing,
        "scalar_metadata": scalars,
        "vector_metadata": vectors,
    }
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(json.dumps(metadata_payload, indent=2, sort_keys=True) + "\n")
    print("Wrote {0} rows to {1}".format(len(rows), output))
    if missing:
        print("WARNING: missing variables: " + ", ".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main()
