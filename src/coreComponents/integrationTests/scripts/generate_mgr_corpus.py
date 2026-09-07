#!/usr/bin/env python3

"""Generate the offline GEOS MGR strategy corpus for hypredrive.

The output contains one serial (np1) and one four-rank (np4) ASCII layout for
each strategy.  Run this script from any checkout and pass the hypredrive data
directory with --output-dir, for example:

    python3 generate_mgr_corpus.py --output-dir /path/to/hypredrive/data
"""

import argparse
import os
from math import ceil
from pathlib import Path


VARIANTS = (1, 4)


def spec(name, blocks, c_levels, f_relax, g_relax, interp, restrict, coarse,
         coarse_flavor="pressure", amg_f=None, coarse_min=None):
    return {
        "name": name,
        "blocks": blocks,
        "c_levels": c_levels,
        "f_relax": f_relax,
        "g_relax": g_relax,
        "interp": interp,
        "restrict": restrict,
        "coarse": coarse,
        "coarse_flavor": coarse_flavor,
        "amg_f": amg_f or {},
        "coarse_min": coarse_min,
    }


SPECS = [
    spec("singlePhaseReservoirFVM", 3, [[0]], ["ge-inv"], ["none"], [12], [0], [0]),
    spec("thermalSinglePhaseReservoirFVM", 5, [[0, 1]], ["ge-inv"], ["none"], [12], [0], [0], "pressureTemperature"),
    spec("singlePhaseHybridFVM", 2, [[1]], ["l1-jacobi"], ["none"], [2], [0], [0]),
    spec("singlePhaseReservoirHybridFVM", 4, [[0, 1], [1]], ["ge-inv", "jacobi"], ["none", "none"], [12, 2], [0, 0], [0, 0]),
    spec("singlePhasePoromechanics", 4, [[3]], ["amg"], ["none"], [2], [0], [1], amg_f={0: "displacementFiltered"}),
    spec("thermalSinglePhasePoromechanics", 5, [[3, 4]], ["amg"], ["none"], [2], [0], [1], "pressureTemperature", amg_f={0: "displacementFiltered"}),
    spec("hybridSinglePhasePoromechanics", 5, [[3, 4], [4]], ["amg", "l1-jacobi"], ["none", "none"], [2, 2], [0, 0], [1, 0], amg_f={0: "displacementFiltered"}),
    spec("singlePhasePoromechanicsEmbeddedFractures", 7, [[0, 1, 2, 6], [6]], ["ge-piv", "amg"], ["none", "none"], [12, 2], [0, 0], [0, 1], amg_f={1: "displacement"}),
    spec("singlePhasePoromechanicsConformingFractures", 7, [[0, 1, 2, 6], [6]], ["none", "amg"], ["ilu", "none"], [12, 2], [0, 0], [0, 1], amg_f={1: "displacement"}),
    spec("singlePhasePoromechanicsReservoirFVM", 6, [[3, 4, 5], [3]], ["amg", "ge-inv"], ["none", "none"], [2, 12], [0, 0], [1, 0], amg_f={0: "displacementFiltered"}),
    spec("thermalSinglePhasePoromechanicsReservoirFVM", 8, [[3, 4, 5, 6, 7], [3, 4]], ["amg", "ge-inv"], ["none", "none"], [2, 12], [0, 0], [1, 0], "pressureTemperature", amg_f={0: "displacementFiltered"}),
    spec("compositionalMultiphaseFVM", 4, [[0, 1], [0]], ["jacobi", "none"], ["none", "ilu"], [2, 0], [0, 14], [0, 0]),
    spec("compositionalMultiphaseHybridFVM", 6, [[5]], ["ilu"], ["none"], [12], [0], [0]),
    spec("compositionalMultiphaseReservoirFVM", 6, [[0, 1, 2], [0, 1], [0]], ["ge-inv", "jacobi", "none"], ["none", "none", "ilu"], [12, 2, 0], [0, 0, 14], [0, 0, 0]),
    spec("compositionalMultiphaseReservoirHybridFVM", 9, [[0, 1, 2, 3, 4, 5], [3]], ["ge-inv", "ilu"], ["none", "none"], [12, 12], [0, 0], [0, 0]),
    spec("immiscibleMultiphaseFVM", 4, [[0]], ["none"], ["ilu"], [0], [14], [0]),
    spec("reactiveCompositionalMultiphaseOBL", 4, [[0]], ["none"], ["ilu"], [0], [14], [0]),
    spec("thermalCompositionalMultiphaseFVM", 4, [[0, 1, 3], [0, 3]], ["jacobi", "none"], ["blk-gs", "ilu"], [2, 0], [0, 15], [0, 0], "pressureTemperature"),
    spec("thermalCompositionalMultiphaseReservoirFVM", 9, [[0, 1, 2, 3], [0, 1, 3], [0, 3]], ["ge-inv", "jacobi", "none"], ["none", "blk-gs", "ilu"], [12, 2, 0], [0, 0, 15], [0, 0, 0], "pressureTemperature"),
    spec("multiphasePoromechanics", 6, [[3, 4, 5], [3, 4], [3]], ["amg", "jacobi", "none"], ["none", "none", "ilu"], [2, 2, 0], [0, 0, 14], [1, 0, 0], amg_f={0: "displacementFiltered"}),
    spec("multiphasePoromechanicsReservoirFVM", 10, [list(range(3, 10)), [3, 4, 5], [3, 4], [3]], ["amg", "ge-inv", "jacobi", "none"], ["none", "none", "none", "ilu"], [2, 12, 2, 0], [0, 0, 0, 14], [1, 0, 0, 0], amg_f={0: "displacementFiltered"}),
    # Four flow components keep the final pressure/temperature coarse set and
    # leave one F component, so HYPRE's scalar partial column-lumped path is
    # exercised without violating its one-C-point-per-F-block requirement.
    spec("thermalMultiphasePoromechanics", 7, [list(range(3, 7)), [3, 4, 6], [3, 6]], ["amg", "jacobi", "none"], ["none", "none", "ilu"], [2, 2, 0], [0, 0, 15], [1, 0, 0], "pressureTemperature", amg_f={0: "displacementFiltered"}),
    spec("hydrofracture", 4, [[3]], ["amg"], ["none"], [2], [0], [0], amg_f={0: "displacementFiltered"}, coarse_min=1000),
    spec("lagrangianContactMechanics", 6, [[0, 1, 2]], ["l1-jacobi"], ["none"], [12], [0], [0], "displacementFiltered"),
    spec("augmentedLagrangianContactMechanics", 6, [[0, 1, 2]], ["l1-jacobi"], ["none"], [12], [0], [0], "displacementFiltered"),
    spec("lagrangianContactMechanicsBubbleStab", 9, [[0, 1, 2, 6, 7, 8], [0, 1, 2]], ["l1-jacobi", "l1-jacobi"], ["none", "none"], [12, 12], [0, 0], [0, 0], "displacementFiltered"),
    spec("solidMechanicsEmbeddedFractures", 6, [[0, 1, 2]], ["jacobi"], ["none"], [12], [0], [0], "displacementFiltered"),
]


# Keep the on-disk dataset names short, like the existing compflow6k and
# poromech2k datasets.  The numeric suffix is added from the total matrix DOF
# count, not from the number of cells or the number of nonzeros.
SHORT_NAMES = {
    "singlePhaseReservoirFVM": "spres",
    "thermalSinglePhaseReservoirFVM": "tspres",
    "singlePhaseHybridFVM": "sphyb",
    "singlePhaseReservoirHybridFVM": "spreshyb",
    "singlePhasePoromechanics": "spporo",
    "thermalSinglePhasePoromechanics": "tspporo",
    "hybridSinglePhasePoromechanics": "hspporo",
    "singlePhasePoromechanicsEmbeddedFractures": "spporoef",
    "singlePhasePoromechanicsConformingFractures": "spporocf",
    "singlePhasePoromechanicsReservoirFVM": "spporores",
    "thermalSinglePhasePoromechanicsReservoirFVM": "tspporores",
    "compositionalMultiphaseFVM": "cmpf",
    "compositionalMultiphaseHybridFVM": "cmphyb",
    "compositionalMultiphaseReservoirFVM": "cmpres",
    "compositionalMultiphaseReservoirHybridFVM": "cmpreshyb",
    "immiscibleMultiphaseFVM": "immf",
    "reactiveCompositionalMultiphaseOBL": "rcmpobl",
    "thermalCompositionalMultiphaseFVM": "tcmpf",
    "thermalCompositionalMultiphaseReservoirFVM": "tcmpres",
    "multiphasePoromechanics": "mpporo",
    "multiphasePoromechanicsReservoirFVM": "mpporores",
    "thermalMultiphasePoromechanics": "tmpporo",
    "hydrofracture": "hydrofrac",
    "lagrangianContactMechanics": "lcontact",
    "augmentedLagrangianContactMechanics": "alcontact",
    "lagrangianContactMechanicsBubbleStab": "lcontactbs",
    "solidMechanicsEmbeddedFractures": "smef",
}


def append_amg(lines, indent, flavor, coarse_min=None):
    p = " " * indent
    q = " " * (indent + 2)
    r = " " * (indent + 4)
    s = " " * (indent + 6)
    lines.extend([
        f"{p}amg:",
        f"{q}tolerance: 0.0",
        f"{q}max_iter: 1",
        f"{q}print_level: 0",
        f"{q}smoother:",
        f"{r}type: schwarz",
        f"{r}num_levels: 0",
        f"{r}num_sweeps: 1",
        f"{r}ilu:",
        f"{s}reordering: 0",
    ])
    if flavor in ("displacement", "displacementFiltered"):
        lines.extend([
            f"{q}coarsening:",
            f"{r}max_coarse_size: 9",
            f"{r}max_row_sum: 1.0",
            f"{r}strong_th: 0.6",
            f"{r}num_functions: 3",
            f"{r}filter_functions: 0",
            f"{q}relaxation:",
            f"{r}order: 1",
        ])
    elif flavor == "pressureTemperature":
        lines.extend([
            f"{q}aggressive:",
            f"{r}num_levels: 1",
            f"{r}max_nnz_row: 16",
            f"{q}coarsening:",
            f"{r}max_coarse_size: 9",
            f"{r}num_functions: 2",
            f"{q}relaxation:",
            f"{r}order: 1",
        ])
    else:
        lines.extend([
            f"{q}aggressive:",
            f"{r}num_levels: 1",
            f"{r}max_nnz_row: 20",
            f"{r}prolongation_type: 4",
            f"{q}coarsening:",
        ])
        if coarse_min is not None:
            lines.append(f"{r}min_coarse_size: {coarse_min}")
        lines.extend([
            f"{r}max_coarse_size: 9",
            f"{q}relaxation:",
            f"{r}order: 1",
        ])


def append_ilu(lines, indent):
    p = " " * indent
    q = " " * (indent + 2)
    lines.extend([
        f"{p}num_sweeps: 1",
        f"{p}ilu:",
        f"{q}type: 0",
        f"{q}fill_level: 0",
        f"{q}max_iter: 1",
        f"{q}tolerance: 0.0",
        f"{q}reordering: 0",
        f"{q}print_level: 0",
        f"{q}max_row_nnz: 1000",
    ])


def append_relaxation(lines, key, value, indent):
    p = " " * indent
    if value == "amg":
        lines.append(f"{p}{key}:")
        lines.append(f"{p}  num_sweeps: 1")
        append_amg(lines, indent + 2, "__LEVEL__")
        return
    if value == "ilu":
        lines.append(f"{p}{key}:")
        append_ilu(lines, indent + 2)
        return
    if value == "none":
        lines.extend([
            f"{p}{key}:",
            f"{p}  type: \"none\"",
            f"{p}  num_sweeps: 0",
        ])
        return
    lines.append(f"{p}{key}: {value}")


def matrix_text(blocks):
    cells = ceil(1024 / blocks)
    n = cells * blocks
    rows = []
    for cell in range(cells):
        for component in range(blocks):
            row = cell * blocks + component
            entries = []
            for other_cell in range(max(0, cell - 2), min(cells, cell + 3)):
                distance = abs(other_cell - cell)
                for other_component in range(blocks):
                    col = other_cell * blocks + other_component
                    if row == col:
                        value = 30.0 + 0.25 * (component + 1)
                    else:
                        value = -0.08 * (1.0 + 0.02 * (component + other_component + 2)) / (1.0 + distance)
                    entries.append((col, value))
            entries.sort()
            rows.extend(f"{row} {col} {value:.16e}" for col, value in entries)
    return f"0 {n - 1} 0 {n - 1}\n" + "\n".join(rows) + "\n", n, cells, len(rows)


def rhs_text(n):
    return f"0 {n - 1}\n" + "\n".join(f"{i} 1.0" for i in range(n)) + "\n"


def dofmap_text(blocks, n):
    labels = [str(i % blocks) for i in range(n)]
    return f"{n}\n" + "\n".join(labels) + "\n"


def partition_bounds(n, blocks, nranks, rank):
    if n % blocks:
        raise ValueError("matrix size must be an integer number of cell blocks")
    cells = n // blocks
    first_cell = cells * rank // nranks
    last_cell = cells * (rank + 1) // nranks
    return first_cell * blocks, last_cell * blocks


def partition_matrix_text(matrix, n, blocks, nranks, rank):
    start, end = partition_bounds(n, blocks, nranks, rank)
    rows = [
        line for line in matrix.splitlines()[1:]
        if start <= int(line.split()[0]) < end
    ]
    return (
        f"{start} {end - 1} {start} {end - 1}\n"
        + "\n".join(rows)
        + "\n"
    )


def partition_rhs_text(n, blocks, nranks, rank):
    start, end = partition_bounds(n, blocks, nranks, rank)
    return (
        f"{start} {end - 1}\n"
        + "\n".join(f"{i} 1.0" for i in range(start, end))
        + "\n"
    )


def partition_dofmap_text(blocks, n, nranks, rank):
    start, end = partition_bounds(n, blocks, nranks, rank)
    return (
        f"{end - start}\n"
        + "\n".join(str(i % blocks) for i in range(start, end))
        + "\n"
    )


def yaml_text(item):
    lines = [
        "# Offline GEOS MGR strategy case generated from the strategy definition.",
        "linear_system:",
        "  type: ij",
        "  matrix_filename: IJ.out.A",
        "  rhs_filename: IJ.out.b",
        "  dofmap_filename: dofmap.out",
        "  rhs_mode: file",
        "  init_guess_mode: zeros",
        "",
        "solver:",
        "  gmres:",
        "    max_iter: 1",
        "    krylov_dim: 10",
        "    relative_tol: 0.0",
        "    print_level: 0",
        "",
        "preconditioner:",
        "  mgr:",
        "    max_iter: 1",
        "    tolerance: 0.0",
        "    print_level: 0",
        "    non_c_to_f: 1",
        "    nonglk_max_elmts: 1",
        "    pmax: 0",
        "    coarse_th: 1.0e-20",
        f"    num_levels: {len(item['c_levels']) + 1}",
        "    level:",
    ]
    active = list(range(item["blocks"]))
    for level, c_labels in enumerate(item["c_levels"]):
        c_set = set(c_labels)
        f_labels = [label for label in active if label not in c_set]
        lines.append(f"      {level}:")
        lines.append("        f_dofs: [{}]".format(", ".join(str(x) for x in f_labels)))
        f_type = item["f_relax"][level]
        if f_type == "amg":
            lines.extend([
                "        f_relaxation:",
                "          num_sweeps: 1",
            ])
            append_amg(lines, 10, item["amg_f"][level])
        elif f_type == "ilu":
            lines.append("        f_relaxation:")
            append_ilu(lines, 10)
        elif f_type == "none":
            lines.extend([
                "        f_relaxation:",
                "          type: \"none\"",
                "          num_sweeps: 0",
            ])
        else:
            lines.append(f"        f_relaxation: {f_type}")

        g_type = item["g_relax"][level]
        if g_type == "ilu":
            lines.append("        g_relaxation:")
            append_ilu(lines, 10)
        elif g_type == "none":
            lines.extend([
                "        g_relaxation:",
                "          type: \"none\"",
                "          num_sweeps: 0",
            ])
        else:
            lines.append(f"        g_relaxation: {g_type}")
        lines.extend([
            f"        prolongation_type: {item['interp'][level]}",
            f"        restriction_type: {item['restrict'][level]}",
            f"        coarse_level_type: {item['coarse'][level]}",
        ])
        active = c_labels[:]
    lines.extend(["    coarsest_level:"])
    append_amg(lines, 6, item["coarse_flavor"], item["coarse_min"])
    return "\n".join(lines) + "\n"


def dataset_name(item, n):
    try:
        short_name = SHORT_NAMES[item["name"]]
    except KeyError as exc:
        raise RuntimeError(f"missing short name for {item['name']}") from exc
    dof_k = max(1, (n + 500) // 1000)
    return f"{short_name}{dof_k}k"


def case_readme(item, dataset, n, cells, nnz):
    c_text = "; ".join("{" + ", ".join(str(x) for x in labels) + "}" for labels in item["c_levels"])
    return f"""# {dataset}\n\nOffline MGR coverage case corresponding to GEOS strategy `{item['name']}`.\n\n- Total degrees of freedom: {n}\n- Matrix size: {n} x {n}\n- Matrix nonzeros: {nnz}\n- Cell blocks: {cells} with {item['blocks']} dof labels per cell\n- Kept labels by reduction level: {c_text}\n- Variants: `np1/` and `np4/`\n- Run: `cd data/{dataset}/np1 && hypredrive-cli input.yml` or `cd data/{dataset}/np4 && mpiexec -n 4 hypredrive-cli input.yml`\n\nThe `k` suffix in `{dataset}` denotes the total matrix DOF count rounded to\nthousands. The matrix is deliberately small and deterministic: it is a\nblock-banded test operator, not a physics model, and exercises the strategy's\nMGR recipe and dof-map handling offline.\n"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate offline np1/np4 MGR strategy cases."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            os.environ.get("HYPREDRIVE_DATA_DIR", "data")
        ),
        help="hypredrive data directory containing the dataset directories",
    )
    return parser.parse_args()


def main():
    output_dir = parse_args().output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    if len(SPECS) != 27:
        raise RuntimeError(f"expected 27 strategies, found {len(SPECS)}")
    if len(SHORT_NAMES) != len(SPECS):
        raise RuntimeError(f"expected {len(SPECS)} short names, found {len(SHORT_NAMES)}")
    if len(set(SHORT_NAMES.values())) != len(SHORT_NAMES):
        raise RuntimeError("short dataset names must be unique")

    for item in SPECS:
        matrix, n, cells, nnz = matrix_text(item["blocks"])
        dataset = dataset_name(item, n)
        case = output_dir / dataset
        case.mkdir(parents=True, exist_ok=True)

        for nranks in VARIANTS:
            variant = case / f"np{nranks}"
            variant.mkdir(parents=True, exist_ok=True)
            for rank in range(nranks):
                (variant / f"IJ.out.A.{rank:05d}").write_text(
                    partition_matrix_text(matrix, n, item["blocks"], nranks, rank)
                )
                (variant / f"IJ.out.b.{rank:05d}").write_text(
                    partition_rhs_text(n, item["blocks"], nranks, rank)
                )
                (variant / f"dofmap.out.{rank:05d}").write_text(
                    partition_dofmap_text(item["blocks"], n, nranks, rank)
                )
            (variant / "input.yml").write_text(yaml_text(item))

        (case / "README.md").write_text(case_readme(item, dataset, n, cells, nnz))

    overview = [
        "# Offline GEOS MGR strategy datasets",
        "",
        "This directory contains one small, self-contained hypredrive dataset for every",
        "`LinearSolverParameters::MGR::StrategyType` currently defined and dispatched",
        "by GEOS (27 cases). These datasets are intended to be distributed through",
        "the same Zenodo archive as the other hypredrive datasets.",
        "",
        "Each short-named directory contains two ASCII layouts:",
        "",
        "- `np1/`: one-rank files for serial execution",
        "- `np4/`: four rank-local files for MPI execution with four ranks",
        "- `README.md`: exact DOF count and reduction labels for the case",
        "",
        "Each variant contains `IJ.out.A.00000` (and additional rank files in",
        "`np4/`), `IJ.out.b.*`, `dofmap.out.*`, and `input.yml`.",
        "",
        "The matrices are deterministic block-banded operators with approximately",
        "1,024 degrees of freedom each. They are intentionally not physics models:",
        "the offline tests target MGR strategy setup, reductions, smoothers, and",
        "restriction/prolongation paths without requiring a GEOS simulation.",
        "",
        "Run one case from the hypredrive checkout with:",
        "",
        "```bash",
        "cd data/spres1k/np1",
        "hypredrive-cli input.yml",
        "",
        "cd ../np4",
        "mpiexec -n 4 hypredrive-cli input.yml",
        "```",
        "",
        "The output directory is selected with `--output-dir` or the",
        "`HYPREDRIVE_DATA_DIR` environment variable.",
        "",
        "| GEOS strategy | Dataset directory | Total DOFs |",
        "| --- | --- | ---: |",
    ]
    for item in SPECS:
        _, n, _, _ = matrix_text(item["blocks"])
        dataset = dataset_name(item, n)
        overview.append(f"| `{item['name']}` | `{dataset}/` | {n} |")
    (output_dir / "mgr_strategy_cases.md").write_text("\n".join(overview) + "\n")


if __name__ == "__main__":
    main()
