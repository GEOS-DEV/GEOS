import h5py
import numpy as np
import glob
import os

def read_and_aggregate_hdf5_series(file_pattern, paths_to_read, global_index_path, coord_path):
    """
    Reads a series of HDF5 files (one per rank), aggregates into a unified 2D array:
    [globalIndex, x, y, z, field1, field2, ...]
    """

    # Expand ~ and glob files
    file_pattern = os.path.expanduser(file_pattern)
    file_list = sorted(glob.glob(file_pattern))
    if not file_list:
        raise RuntimeError(f"No files matched pattern: {file_pattern}")

    print(f"[INFO] Found {len(file_list)} files matching {file_pattern}")
    for file_name in file_list:
        print(f"[INFO] Loading file: {file_name}")

    # Initialize storage
    all_global_indices = []
    all_coords = []
    all_fields = {path: [] for path in paths_to_read}

    # Loop over all files
    for file_name in file_list:
        with h5py.File(file_name, 'r') as f:

            # --- Global Indices ---
            if global_index_path not in f:
                raise ValueError(f"Global index group {global_index_path} not found in {file_name}")

            global_group = f[global_index_path]
            if '__values__' not in global_group:
                raise ValueError(f"__values__ dataset not found inside global index group {global_index_path}")

            global_dataset = global_group['__values__']
            global_indices = np.asarray(global_dataset)
            all_global_indices.append(global_indices)


            # Reshape if needed (usually 1D, but safer to be general)
            if '__dimensions__' in global_dataset.attrs:
                dims = global_dataset.attrs['__dimensions__']
                global_indices = global_indices.reshape(dims)
            elif global_indices.ndim != 1:
                global_indices = global_indices.reshape((-1,))


            # --- Coordinates ---
            if coord_path not in f:
                raise ValueError(f"Coordinate path {coord_path} not found in {file_name}")
            coord_group = f[coord_path]
            if '__values__' not in coord_group:
                raise ValueError(f"__values__ dataset not found under {coord_path} in {file_name}")
            coords_dataset = coord_group['__values__']
            coords = np.asarray(coords_dataset)

            # Try reshaping properly
            if '__dimensions__' in coords_dataset.attrs:
                dims = coords_dataset.attrs['__dimensions__']
                coords = coords.reshape(dims)
            elif '__dimensions__' in coord_group.attrs:
                dims = coord_group.attrs['__dimensions__']
                coords = coords.reshape(dims)
            elif coords.ndim == 1:
                coords = coords.reshape((-1, 3))
            if coords.shape[1] != 3:
                raise ValueError(f"Coordinates should have shape (N,3) but got {coords.shape}")

            all_coords.append(coords)

            # --- Fields ---
            for path in paths_to_read:
                if path not in f:
                    raise ValueError(f"Field group {path} not found in {file_name}")

                group = f[path]
                if '__values__' not in group:
                    raise ValueError(f"__values__ dataset not found inside group {path} in {file_name}")
                
                dataset = group['__values__']
                field_data = np.asarray(dataset)

                # Try reshaping properly
                if '__dimensions__' in dataset.attrs:
                    dims = dataset.attrs['__dimensions__']
                    field_data = field_data.reshape(dims)
                elif field_data.ndim == 1:
                    field_data = field_data.reshape((-1,))  # probably 1D fields like pressure, saturation

                all_fields[path].append(field_data)


    # --- Concatenate all ranks ---
    all_global_indices = np.concatenate(all_global_indices, axis=0)
    all_coords = np.concatenate(all_coords, axis=0)
    for path in all_fields:
        all_fields[path] = np.concatenate(all_fields[path], axis=0)

    # --- Sort by global index ---
    sort_order = np.argsort(all_global_indices)
    all_global_indices = all_global_indices[sort_order]
    all_coords = all_coords[sort_order]
    for path in all_fields:
        all_fields[path] = all_fields[path][sort_order]

    # --- Assemble final array ---
    columns = [all_global_indices[:, np.newaxis],  # (N,1)
               all_coords]                         # (N,3)
    for path in paths_to_read:
        columns.append(all_fields[path][:, np.newaxis])  # (N,1)

    final_array = np.hstack(columns)

    print(f"[INFO] Aggregated final array shape: {final_array.shape}")

    return final_array

def print_aggregated_table(array, field_names, num_rows=10):
    """
    Pretty-prints an aggregated array with headers.

    Parameters:
    - array (np.ndarray): The 2D numpy array to print
    - field_names (list of str): List of column names
    - num_rows (int): Number of rows to print (default: 10)
    """

    if array.shape[1] != len(field_names):
        raise ValueError(f"Mismatch: array has {array.shape[1]} columns but {len(field_names)} field names provided.")

    # Print header
    header_line = " | ".join(f"{name:^12}" for name in field_names)
    print(header_line)
    print("-" * len(header_line))

    # Print each row
    for row in array[:num_rows]:
        row_items = []
        for i, value in enumerate(row):
            if i == 0:
                row_items.append(f"{int(value):12d}")  # globalIndex as integer
            else:
                row_items.append(f"{value:12.5e}")      # floating-point for others
        print(" | ".join(row_items))


def print_aggregated_diff_table(array1, array2, field_names, tolerance=1e-8, num_rows=None):
    """
    Pretty-prints a side-by-side comparison of two aggregated arrays.

    Parameters:
    - array1, array2 (np.ndarray): Aggregated arrays to compare
    - field_names (list of str): List of column names
    - tolerance (float): Coordinate comparison tolerance
    - num_rows (int or None): Number of rows to print (default: all)
    """

    if array1.shape[1] != len(field_names) or array2.shape[1] != len(field_names):
        raise ValueError("Mismatch between array shape and number of field names.")

    # Build lookup by globalIndex
    map1 = {int(row[0]): row for row in array1}
    map2 = {int(row[0]): row for row in array2}

    # Find union of all globalIndices
    all_indices = sorted(set(map1.keys()).union(map2.keys()))

    if num_rows is None:
        num_rows = len(all_indices)

    # Print headers
    header = f"{'globalIndex':>12} | {'x1':>12} | {'y1':>12} | {'z1':>12} | {'x2':>12} | {'y2':>12} | {'z2':>12} | {'coord_match':^12}"
    for field in field_names[4:]:  # Skip globalIndex, x, y, z
        header += f" | {field + '_1':>12} | {field + '_2':>12} | {field + '_Δ':>12}"
    print(header)
    print("-" * len(header))

    # Print rows
    for idx in all_indices[:num_rows]:
        row1 = map1.get(idx)
        row2 = map2.get(idx)

        if row1 is None:
            print(f"{idx:12d} | {'(missing)':>12} | {'(missing)':>12} | {'(missing)':>12} | {row2[1]:12.5e} | {row2[2]:12.5e} | {row2[3]:12.5e} | {'MISSING 1':^12}")
            continue
        if row2 is None:
            print(f"{idx:12d} | {row1[1]:12.5e} | {row1[2]:12.5e} | {row1[3]:12.5e} | {'(missing)':>12} | {'(missing)':>12} | {'(missing)':>12} | {'MISSING 2':^12}")
            continue

        # Compare coordinates
        coord_diff = np.linalg.norm(row1[1:4] - row2[1:4])
        coord_match = "OK" if coord_diff < tolerance else "MISMATCH"

        line = f"{idx:12d} |"
        for i in range(1, 4):
            line += f" {row1[i]:12.5e} |"
        for i in range(1, 4):
            line += f" {row2[i]:12.5e} |"
        line += f" {coord_match:^12}"

        # Compare fields
        for i in range(4, len(field_names)):
            val1 = row1[i]
            val2 = row2[i]
            delta = val2 - val1
            line += f" | {val1:12.5e} | {val2:12.5e} | {delta:12.5e}"

        print(line)

def summarize_aggregated_diff(array1, array2, field_names, coord_tolerance=1e-8, field_tolerance=1e-8):
    """
    Summarizes differences between two aggregated arrays with scaled tolerances.

    Parameters:
    - array1, array2 (np.ndarray): Aggregated arrays
    - field_names (list of str): List of column names
    - coord_tolerance (float): Relative tolerance for coordinates (default 1e-8)
    - field_tolerance (float): Relative tolerance for fields (default 1e-8)

    Returns:
    - dict: summary counts
    """

    if array1.shape[1] != len(field_names) or array2.shape[1] != len(field_names):
        raise ValueError("Mismatch between array shape and field names length.")

    # Build maps by globalIndex
    map1 = {int(row[0]): row for row in array1}
    map2 = {int(row[0]): row for row in array2}

    all_indices = sorted(set(map1.keys()).union(map2.keys()))

    # Compute maximum norm of coordinates
    coords1 = array1[:, 1:4]
    coords2 = array2[:, 1:4]
    max_coord_norm = max(np.max(np.linalg.norm(coords1, axis=1)),
                         np.max(np.linalg.norm(coords2, axis=1)))
    if max_coord_norm == 0.0:
        max_coord_norm = 1.0  # prevent division by zero

    # Compute maximum absolute value for each field
    field_max_abs = []
    for i in range(4, len(field_names)):
        max1 = np.max(np.abs(array1[:, i])) if array1.shape[0] > 0 else 0.0
        max2 = np.max(np.abs(array2[:, i])) if array2.shape[0] > 0 else 0.0
        field_max_abs.append(max(max1, max2) or 1.0)  # fallback to 1.0 to avoid division by zero

    # Initialize counters
    missing_in_1 = 0
    missing_in_2 = 0
    coord_mismatches = 0
    field_mismatches = 0

    # Now check all rows
    for idx in all_indices:
        row1 = map1.get(idx)
        row2 = map2.get(idx)

        if row1 is None:
            missing_in_1 += 1
            continue
        if row2 is None:
            missing_in_2 += 1
            continue

        # Check coordinates
        coord_diff = np.linalg.norm(row1[1:4] - row2[1:4])
        if coord_diff > coord_tolerance * max_coord_norm:
            coord_mismatches += 1

        # Check fields
        for i in range(4, len(field_names)):
            val1 = row1[i]
            val2 = row2[i]
            delta = abs(val2 - val1)
            max_val = field_max_abs[i-4]
            if delta > field_tolerance * max_val:
                field_mismatches += 1

    return {
        "total_indices": len(all_indices),
        "missing_in_1": missing_in_1,
        "missing_in_2": missing_in_2,
        "coord_mismatches": coord_mismatches,
        "field_mismatches": field_mismatches,
        "max_coord_norm": max_coord_norm,
        "field_max_abs": field_max_abs
    }



# Example usage:
if __name__ == "__main__":
    file_pattern1 = "~/Downloads/baseline_integratedTests-pr3624-11053-ae011c7/singlePhaseFlow/incompressible_pebi3d_02/0to1_restart_000000001/rank_000000*.hdf5"
    file_pattern2 = "~/Downloads/baseline_integratedTests-pr3626-11174-d100cd2/singlePhaseFlow/incompressible_pebi3d_02/0to1_restart_000000001/rank_000000*.hdf5"
    paths_to_compare = ["/Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup/Domain/elementSubRegions/pentagonalPrisms/pressure"]
    global_index_path = "/Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup/Domain/elementSubRegions/pentagonalPrisms/localToGlobalMap"  # or None
    coord_path = "/Problem/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup/Domain/elementSubRegions/pentagonalPrisms/elementCenter"  # fallback

    aggregated1 = read_and_aggregate_hdf5_series(file_pattern1, paths_to_compare, global_index_path, coord_path)
    aggregated2 = read_and_aggregate_hdf5_series(file_pattern2, paths_to_compare, global_index_path, coord_path)

    # Suppose you have:
    field_names = ["globalIndex", "x", "y", "z", "pressure"]

#    print_aggregated_diff_table(aggregated1, aggregated2, field_names)

    summary = summarize_aggregated_diff(aggregated1, aggregated2, field_names)

    print("\nSummary of Differences:")
    for key, value in summary.items():
        print(f"{key:>20}: {value}")

