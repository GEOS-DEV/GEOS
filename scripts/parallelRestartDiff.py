import h5py
import numpy as np
import glob
import os
import sys
import itertools
import re
import argparse

# ------------------------------
# Utilities for wildcard expansion
# ------------------------------
def expand_field_paths_per_object(file_name, patterns):
    """
    Expand wildcard patterns and group by base path.
    Raise an error if any specified pattern does not match anything.
    """
    expanded = {}
    with h5py.File(file_name, 'r') as f:
        for pattern in patterns:
            tokens = []
            pos = 0
            has_wildcard = False
            for match in re.finditer(r'\[([^]]+)\]', pattern):
                has_wildcard = True
                tokens.append(pattern[pos:match.start()])
                options = match.group(1).split(',')
                tokens.append(options)
                pos = match.end()
            tokens.append(pattern[pos:])
            option_lists = [t if isinstance(t, list) else [t] for t in tokens]

            matched_any = False
            for combination in itertools.product(*option_lists):
                full_path = ''.join(combination)
                if full_path in f:
                    base_path = extract_base_path(full_path)
                    field_name = full_path.split('/')[-1]
                    expanded.setdefault(base_path, []).append(field_name)
                    matched_any = True

            if not matched_any:
                if has_wildcard:
                    raise ValueError(f"No fields matched pattern: {pattern}")
                else:
                    raise ValueError(f"Specified field path not found in HDF5: {pattern}")

    print(f"[INFO] Expanded {sum(len(v) for v in expanded.values())} fields over {len(expanded)} objects.")
    return expanded


def extract_base_path(field_path):
    return '/'.join(field_path.split('/')[:-1])

# ------------------------------
# Read and Aggregate
# ------------------------------
def infer_coord_dataset_name(base_path):
    """
    Infer the coordinate dataset name based on the container object name.
    """
    container = base_path.split('/')[-1]
    if container.endswith('nodeManager'):
        return 'ReferencePosition'
    elif container.endswith('FaceManager'):
        return 'faceCenter'
    elif container.endswith('EdgeManager'):
        return 'edgeCenter'
    else:
        return 'elementCenter'


def read_and_aggregate_hdf5_series(file_pattern, paths_to_read):
    file_pattern = os.path.expanduser(file_pattern)
    file_list = sorted(glob.glob(file_pattern))
    if not file_list:
        raise RuntimeError(f"No files matched pattern: {file_pattern}")

    aggregated_data = {}

    for base_path, field_list in paths_to_read.items():
        all_rows = []
        for file_name in file_list:
            with h5py.File(file_name, 'r') as f:
                global_index_path = f"{base_path}/localToGlobalMap"
                coord_dataset_name = infer_coord_dataset_name(base_path)
                coord_path = f"{base_path}/{coord_dataset_name}"

                # Global Indices
                group = f[global_index_path]
                dataset = group['__values__']
                global_indices = np.asarray(dataset)
                if '__dimensions__' in dataset.attrs:
                    dims = dataset.attrs['__dimensions__']
                    global_indices = global_indices.reshape(dims)
                elif global_indices.ndim != 1:
                    global_indices = global_indices.reshape((-1,))

                # Coordinates
                group = f[coord_path]
                dataset = group['__values__']
                coords = np.asarray(dataset)
                if '__dimensions__' in group:
                    dims = np.asarray(group['__dimensions__'])
                    coords = coords.reshape(dims)
                elif coords.ndim == 1:
                    coords = coords.reshape((-1, 3))
                if coords.shape[1] != 3:
                    raise ValueError(f"Coordinates should have shape (N,3) but got {coords.shape}")

                # Fields
                field_datas = []
                for path in field_list:
                    field_path = f"{base_path}/{path}"
                    group = f[field_path]
                    dataset = group['__values__']
                    field_data = np.asarray(dataset)

                    if '__dimensions__' in group:
                        dims = np.asarray(group['__dimensions__'])
                        field_data = field_data.reshape(dims)
                        print(f"[INFO] Reshaping field data {field_path} to {dims}")
                    elif field_data.ndim == 1:
                        field_data = field_data.reshape((-1,))
                    else:
                        raise ValueError(f"Unsupported data shape for {field_path}: {field_data.shape}")

                    field_datas.append(field_data)

                # Assemble Rows
                num_entries = global_indices.shape[0]
                for i in range(num_entries):
                    row = [global_indices[i]] + list(coords[i])
                    for field_data in field_datas:
                        if field_data.ndim == 1 or (field_data.shape[1] == 1):
                            # scalar field
                            row.append(field_data[i])
                        else:
                            # vector/tensor field
                            row.extend(field_data[i])
                    all_rows.append(row)


        all_rows = np.array(all_rows)
        sort_order = np.argsort(all_rows[:, 0])
        aggregated_data[base_path] = all_rows[sort_order]

    print(f"[INFO] Aggregated {len(aggregated_data)} object groups.")
    return aggregated_data


# ------------------------------
# Field Name Generation
# ------------------------------
def generate_field_names_per_object(file_name, base_path_to_fields):
    field_names_per_object = {}
    with h5py.File(file_name, 'r') as f:
        for base_path, fields in base_path_to_fields.items():
            field_names = ["globalIndex", "x", "y", "z"]
            for field in fields:
                full_path = f"{base_path}/{field}"
                group = f[full_path]
                dataset = group['__values__']

                if '__dimensions__' in group:
                    dims = np.asarray(group['__dimensions__'])
                else:
                    dims = dataset.shape

                # Correct handling
                if len(dims) == 1:
                    field_names.append(field)
                elif len(dims) == 2 and dims[1] == 1:
                    field_names.append(field)
                else:
                    field_names.extend([f"{field}_{i}" for i in range(dims[1])])

            field_names_per_object[base_path] = field_names
    return field_names_per_object


# ------------------------------
# Printing Utilities
# ------------------------------
def print_aggregated_table(aggregated_data, field_names_per_object, num_rows=None):
    for base_path, array in aggregated_data.items():
        field_names = field_names_per_object.get(base_path)
        if field_names is None:
            raise ValueError(f"No field names available for {base_path}.")

        print(f"\n[INFO] Object: {base_path} ({array.shape[0]} entries)")

        if array.shape[1] != len(field_names):
            raise ValueError(f"Mismatch between array shape and field names for {base_path}.")

        rows_to_print = array.shape[0] if num_rows is None else min(array.shape[0], num_rows)
        header = " | ".join(f"{name:^12}" for name in field_names)
        print(header)
        print("-" * len(header))

        for row in array[:rows_to_print]:
            row_items = []
            for i, value in enumerate(row):
                if i == 0:
                    row_items.append(f"{int(value):12d}")
                else:
                    row_items.append(f"{value:12.5e}")
            print(" | ".join(row_items))





def print_aggregated_diff_table(aggregated1, aggregated2, field_names_per_object, coord_tolerance=1e-8, num_rows=None):
    all_keys = sorted(set(aggregated1.keys()).union(aggregated2.keys()))
    for base_path in all_keys:
        print(f"\n[DIFF] Object: {base_path}")
        array1 = aggregated1.get(base_path)
        array2 = aggregated2.get(base_path)
        field_names = field_names_per_object.get(base_path)

        if field_names is None:
            raise ValueError(f"No field names available for {base_path}.")

        if array1 is None and array2 is None:
            print("[WARN] Missing in BOTH datasets.")
            continue
        if array1 is None:
            print("[WARN] Missing in dataset 1. No comparison possible.")
            continue
        if array2 is None:
            print("[WARN] Missing in dataset 2. No comparison possible.")
            continue

        if array1.shape[1] != len(field_names) or array2.shape[1] != len(field_names):
            raise ValueError(f"Mismatch between array shape and field names for {base_path}.")

        map1 = {int(row[0]): row for row in array1}
        map2 = {int(row[0]): row for row in array2}
        all_indices = sorted(set(map1.keys()).union(map2.keys()))
        rows_to_print = len(all_indices) if num_rows is None else min(len(all_indices), num_rows)

        header = f"{'globalIndex':>12} | {'x1':>12} | {'y1':>12} | {'z1':>12} | {'x2':>12} | {'y2':>12} | {'z2':>12} | {'coord_match':^12}"
        for field in field_names[4:]:
            header += f" | {field + '_1':>12} | {field + '_2':>12} | {field + '_Δ':>12}"
        print(header)
        print("-" * len(header))

        for idx in all_indices[:rows_to_print]:
            row1 = map1.get(idx)
            row2 = map2.get(idx)

            if row1 is None or row2 is None:
                print(f"{idx:12d} | {'MISSING':>12}")
                continue

            coord_diff = np.linalg.norm(row1[1:4] - row2[1:4])
            coord_match = "OK" if coord_diff < coord_tolerance else "MISMATCH"

            line = f"{idx:12d} |"
            for i in range(1, 4):
                line += f" {row1[i]:12.5e} |"
            for i in range(1, 4):
                line += f" {row2[i]:12.5e} |"
            line += f" {coord_match:^12}"

            for i in range(4, len(field_names)):
                val1 = row1[i]
                val2 = row2[i]
                delta = val2 - val1
                line += f" | {val1:12.5e} | {val2:12.5e} | {delta:12.5e}"
            print(line)





# ------------------------------
# Summarizing and Pass/Fail
# ------------------------------
def summarize_aggregated_diff(aggregated1, aggregated2, field_names_per_object, coord_tolerance=1e-8, field_tolerance=1e-8):
    all_keys = sorted(set(aggregated1.keys()).union(aggregated2.keys()))
    overall_summary = {}

    for base_path in all_keys:
        array1 = aggregated1.get(base_path)
        array2 = aggregated2.get(base_path)
        field_names = field_names_per_object.get(base_path)

        if field_names is None:
            raise ValueError(f"No field names available for {base_path}.")

        if array1 is None:
            overall_summary[base_path] = {
                "total_indices": 0,
                "missing_in_1": 1,
                "missing_in_2": 0,
                "coord_mismatches": 0,
                "field_mismatches": 0,
                "max_coord_norm": 0.0,
                "field_max_abs": []
            }
            continue

        if array2 is None:
            overall_summary[base_path] = {
                "total_indices": 0,
                "missing_in_1": 0,
                "missing_in_2": 1,
                "coord_mismatches": 0,
                "field_mismatches": 0,
                "max_coord_norm": 0.0,
                "field_max_abs": []
            }
            continue

        if array1.shape[1] != len(field_names) or array2.shape[1] != len(field_names):
            print( f"[WARN] array1.shape[1] = {array1.shape[1]}.")
            print( f"[WARN] array2.shape[1] = {array2.shape[1]}.")
            print( f"[WARN] len(field_names) = {len(field_names)}.")
            print( f"[WARN] field_names = {field_names}.")
            raise ValueError(f"Mismatch between array shape and field names for {base_path}.")

        map1 = {int(row[0]): row for row in array1}
        map2 = {int(row[0]): row for row in array2}
        all_indices = sorted(set(map1.keys()).union(map2.keys()))

        coords1 = array1[:, 1:4]
        coords2 = array2[:, 1:4]
        max_coord_norm = max(np.max(np.linalg.norm(coords1, axis=1)),
                             np.max(np.linalg.norm(coords2, axis=1)))
        if max_coord_norm == 0.0:
            max_coord_norm = 1.0

        field_max_abs = []
        for i in range(4, len(field_names)):
            max1 = np.max(np.abs(array1[:, i])) if array1.shape[0] > 0 else 0.0
            max2 = np.max(np.abs(array2[:, i])) if array2.shape[0] > 0 else 0.0
            field_max_abs.append(max(max1, max2) or 1.0)

        missing_in_1 = missing_in_2 = coord_mismatches = field_mismatches = 0
        for idx in all_indices:
            row1 = map1.get(idx)
            row2 = map2.get(idx)
            if row1 is None:
                missing_in_1 += 1
                continue
            if row2 is None:
                missing_in_2 += 1
                continue
            coord_diff = np.linalg.norm(row1[1:4] - row2[1:4])
            if coord_diff > coord_tolerance * max_coord_norm:
                coord_mismatches += 1
            for i in range(4, len(field_names)):
                val1, val2 = row1[i], row2[i]
                delta = abs(val2 - val1)

                field_name = field_names[i]
                if 'displacement' in field_name.lower():  # special rule for displacement-like fields
                    scale = max_coord_norm
                else:
                    scale = field_max_abs[i-4]

                if delta > field_tolerance * scale:
                    field_mismatches += 1


        overall_summary[base_path] = {
            "total_indices": len(all_indices),
            "missing_in_1": missing_in_1,
            "missing_in_2": missing_in_2,
            "coord_mismatches": coord_mismatches,
            "field_mismatches": field_mismatches,
            "max_coord_norm": max_coord_norm,
            "field_max_abs": field_max_abs
        }

    return overall_summary



def check_pass_fail(overall_summary):
    overall_pass = True
    for base_path, summary in overall_summary.items():
        print(f"\n[CHECK] Object: {base_path}")
        fail = False
        if summary['missing_in_1']:
            print(f"[FAIL] missing_in_1: {summary['missing_in_1']}")
            fail = True
        if summary['missing_in_2']:
            print(f"[FAIL] missing_in_2: {summary['missing_in_2']}")
            fail = True
        if summary['coord_mismatches']:
            print(f"[FAIL] coord_mismatches: {summary['coord_mismatches']}")
            fail = True
        if summary['field_mismatches']:
            print(f"[FAIL] field_mismatches: {summary['field_mismatches']}")
            fail = True
        if not fail:
            print("[PASS] No mismatches found.")
        else:
            overall_pass = False

    return overall_pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two sets of parallel HDF5 restart files.")

    parser.add_argument("--baseline", required=True, help="path to baseline (e.g., '~/Downloads/old/rank_*.hdf5')")
    parser.add_argument("--restart", required=True, help="path to restart (e.g., '~/Downloads/new/rank_*.hdf5')")
    parser.add_argument("--fields", required=True, action="append", help="Field paths with optional wildcards (can specify multiple).")

    args = parser.parse_args()

    baseline_file_pattern = args.baseline
    restart_file_pattern = args.restart
    raw_field_patterns = args.fields
    
    one_example_file = glob.glob(os.path.expanduser(baseline_file_pattern))[0]

    base_path_to_fields = expand_field_paths_per_object(one_example_file, raw_field_patterns)

    print("\nExpanded base_path_to_fields:")
    for base_path, fields in base_path_to_fields.items():
        print(f"  {base_path} -> {fields}")

    field_names_per_object = generate_field_names_per_object(one_example_file, base_path_to_fields)

    aggregated1 = read_and_aggregate_hdf5_series(baseline_file_pattern, base_path_to_fields)
    aggregated2 = read_and_aggregate_hdf5_series(restart_file_pattern, base_path_to_fields)


    overall_summary = summarize_aggregated_diff(aggregated1, aggregated2, field_names_per_object)

    print("\nSummary of Differences:")
    for base_path, summary in overall_summary.items():
        print(f"\n[Summary for Object: {base_path}]")
        for key, value in summary.items():
            if isinstance(value, list):
                for field, val in zip(field_names_per_object[base_path][4:], value):
                    print(f"  Max abs {field:>15}: {val:.5e}")
            else:
                print(f"{key:>20}: {value}")

    passed = check_pass_fail(overall_summary)


    if not passed:
        print_aggregated_diff_table(aggregated1, aggregated2, field_names_per_object)
        sys.exit(1)
