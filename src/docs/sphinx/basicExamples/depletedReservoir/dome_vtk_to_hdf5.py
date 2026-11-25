r"""
Convert the GEOS VTK output of the dome example into a single HDF5 file.

GEOS writes one ``.vtm`` per output time (plus a ``.pvd`` index) and, inside
each ``.vtm``, one ``.vtu`` per MPI rank and per region.  This script reads that
tree and stores everything the documentation needs in one compact, rank
independent file:

    /times/seconds        (num_times,)      output times [s]
    /times/years          (num_times,)      output times [yr]  (1 yr = 31557600 s)
    /mesh/center          (num_cells, 3)    cell centres  x, y, z   [m]
    /mesh/region          (num_cells,)      index in the alphabetically sorted region names
                                            (0 = CAPROCK, 1 = CLAY, 2 = SAND; attribute
                                            ``region_names`` lists the names)
    /fields/<name>        (num_times, num_cells)   one dataset per scalar field, NaN in the
                                            cells where the field does not exist
                                            (e.g. the inactive CAPROCK)
                                            With --field-times a dataset only stores some
                                            times; its attribute ``time_index`` gives the
                                            matching rows of /times.
                                            With --round the attribute ``precision`` gives
                                            the (lossy) precision applied.

Field names
-----------
Scalar arrays keep their VTK name (``pressure``).  Multi-component arrays are
split into one dataset per component, named ``<array>_<component>``:

    phaseVolumeFraction_gas, phaseVolumeFraction_water
    globalCompFraction_CH4, globalCompFraction_CO2, globalCompFraction_H2S, ...

    The component labels are taken from the VTK file when it carries them.  When
    it does not, they are rebuilt from ``--phases`` according
    to the number of components of the array.

Usage
-----
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd                 # writes dome_soreide_whitson_bench.h5
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd -o dome_results.h5
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --list-fields
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --times-years 0 5 8 13 15.5 18 20
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --fields pressure 'phaseVolumeFraction*'

Generating dome_soreide_whitson_bench.h5 (the file read by plot_dome.py)
-------------------------------------------------------------------------
1. Run the dome example with GEOS; it writes ``<case>/VTK.OUTPUT.pvd``.
2. Run (replace ``VTK.OUTPUT.pvd`` with its path if needed)::

    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd -o dome_soreide_whitson_bench.h5 \
        --fields pressure phaseVolumeFraction globalCompFraction phaseCompFraction \
        --keep pressure phaseVolumeFraction_gas globalCompFraction_CH4 'phaseCompFraction_gas/CO2' \
        --times-years 0 5 8 13 13.5 15.5 18 20 \
        --field-times pressure=0,5,8,13,15.5,18,20 phaseVolumeFraction_gas=0,5,8,13,15.5,18,20 \
                      globalCompFraction_CH4=0,5,8,13,15.5,18,20 'phaseCompFraction_gas/CO2=13.5,15.5,18,20' \
        --round pressure=bits12 phaseVolumeFraction_gas=1e-4 globalCompFraction_CH4=1e-4 \
                'phaseCompFraction_gas/CO2=1e-4'

   * ``--fields`` / ``--keep``: only the four datasets used by the figures;
   * ``--times-years``: the times shown in plot_dome.py (``TIMES_ALL`` and ``TIMES_CO2``);
   * ``--field-times``: each field only at the times it is plotted (lossless);
   * ``--round``: LOSSY precision reduction (pressure: relative error < 2**-12,
     fractions: absolute error <= 5e-5), invisible in the figures.  Drop this
     option to keep the full float32 precision (about 3.5 MB instead of 1 MB).

   If a figure is added or its times change in plot_dome.py, update ``--keep``,
   ``--times-years`` and ``--field-times`` accordingly.
3. Copy the file next to plot_dome.py, or point ``DOME_HDF5`` (or ``--hdf5``) to it.

Requirements: numpy, h5py, vtk (``pip install numpy h5py vtk``).
"""
import argparse
import datetime
import fnmatch
import os
import sys
import xml.etree.ElementTree as xml_tree
import numpy
import vtk
from vtk.util import numpy_support

SECONDS_PER_YEAR = 31557600.0
DEFAULT_FIELDS = ("pressure", "phaseVolumeFraction", "globalCompFraction")
DEFAULT_PHASES = ("gas", "water")

# --------------------------------------------------------------------------
# Index files
# --------------------------------------------------------------------------
def read_pvd(pvd_path):
    """
    Parse a ParaView Data (PVD) file to extract simulation time steps.

    Parameters:
        pvd_path (str): The file path to the .pvd file.

    Returns:
        list: A time-sorted list of tuples, where each tuple contains:
              - (float) time_in_seconds
              - (str) absolute_vtm_path
    """
    root_element = xml_tree.parse(pvd_path).getroot()
    base_directory = os.path.dirname(os.path.abspath(pvd_path))
    time_steps = []
    
    for dataset in root_element.iter("DataSet"):
        time_value = float(dataset.attrib["timestep"])
        file_path = os.path.normpath(os.path.join(base_directory, dataset.attrib["file"]))
        time_steps.append((time_value, file_path))
        
    time_steps.sort(key=lambda step: step[0])
    return time_steps


def region_files(vtm_path):
    """
    Map region names to their corresponding list of .vtu files.

    Parameters:
        vtm_path (str): The file path to the .vtm file.

    Returns:
        dict: A dictionary mapping region names to a list of absolute file paths 
              for the .vtu files belonging to that region.
    """
    root_element = xml_tree.parse(vtm_path).getroot()
    base_directory = os.path.dirname(os.path.abspath(vtm_path))
    region_to_files_map = {}

    def walk_blocks(block, is_in_cell_regions, current_path):
        for child in block:
            if child.tag == "Block":
                block_name = child.attrib.get("name", "")
                walk_blocks(child, is_in_cell_regions or block_name == "CellElementRegion",
                            current_path + [block_name])
            elif child.tag == "DataSet" and is_in_cell_regions:
                # Expected structure: .../CellElementRegion/<REGION>/rank_i
                region_name = current_path[-1]
                file_path = os.path.normpath(os.path.join(base_directory, child.attrib["file"]))
                region_to_files_map.setdefault(region_name, []).append(file_path)

    top_block = root_element.find("vtkMultiBlockDataSet")
    if top_block is None:
        raise RuntimeError(f"{vtm_path} is not a valid vtkMultiBlockDataSet file")
    walk_blocks(top_block, False, [])
    
    return region_to_files_map

def read_vtu(file_path):
    """
    Read an Unstructured Grid (.vtu) file using VTK.

    Parameters:
        file_path (str): The file path to the .vtu file.

    Returns:
        vtkUnstructuredGrid: The loaded VTK unstructured grid data.
    """
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(file_path)
    vtu_reader.Update()
    return vtu_reader.GetOutput()

def cell_centers(unstructured_grid):
    """
    Compute and return the cell centers for a given grid.

    Parameters:
        unstructured_grid (vtkUnstructuredGrid): The mesh grid to process.

    Returns:
        numpy.ndarray: A (num_cells, 3) array of cell center coordinates, compressed to float32.
    """
    center_filter = vtk.vtkCellCenters()
    center_filter.SetInputData(unstructured_grid)
    center_filter.Update()
    centers_array = numpy_support.vtk_to_numpy(center_filter.GetOutput().GetPoints().GetData())
    return centers_array.astype(numpy.float32)

def component_labels(array_base_name, num_components, vtk_data_array, phases):
    """
    Generate readable names for the components of a multi-component array.

    Parameters:
        array_base_name (str): The base name of the array (e.g., 'pressure').
        num_components (int): The number of components in the array.
        vtk_data_array (vtkDataArray): The VTK array holding component metadata.
        phases (list): List of expected phase names.

    Returns:
        list: A list of strings representing the label for each component.
    """
    labels = []
    for component_index in range(num_components):
        component_name = vtk_data_array.GetComponentName(component_index)
        labels.append(component_name if component_name else None)
        
    if all(labels):
        return [label if label.startswith(array_base_name + "_") else f"{array_base_name}_{label}" for label in labels]
        
    if num_components == len(phases):
        inferred_names = list(phases)
    else:
        inferred_names = [str(component_index) for component_index in range(num_components)]
        
    return [f"{array_base_name}_{name}" for name in inferred_names]


def read_region(region_files_list, phases, fields_to_keep):
    """
    Read and concatenate all the rank files of one specific region.

    Parameters:
        region_files_list (list): A list of absolute file paths to .vtu files.
        phases (list): List of phase names used for component naming.
        fields_to_keep (list): Array patterns specifying which fields to export.

    Returns:
        tuple: A tuple containing:
               - (numpy.ndarray) centers: The sorted cell centers.
               - (dict) fields: Dictionary mapping field names to 1D numpy arrays of sorted field data.
    """
    data_chunks = []
    
    for file_name in sorted(region_files_list):
        unstructured_grid = read_vtu(file_name)
        num_cells = unstructured_grid.GetNumberOfCells()
        if num_cells == 0:
            continue
            
        valid_cells_mask = numpy.ones(num_cells, dtype=bool)                    # Used to drop ghost cells
        cell_data = unstructured_grid.GetCellData()
        
        for ghost_name, ownership_check in (("vtkGhostType", lambda attr: attr == 0),
                                            ("ghostRank", lambda attr: attr < 0)):
            ghost_array = cell_data.GetArray(ghost_name)
            if ghost_array is not None:
                valid_cells_mask &= ownership_check(numpy_support.vtk_to_numpy(ghost_array).reshape(num_cells, -1)[:, 0])
                
        region_data = {
            "centers": cell_centers(unstructured_grid)[valid_cells_mask], 
            "fields": {}
        }
        
        for data_index in range(cell_data.GetNumberOfArrays()):
            vtk_array = cell_data.GetAbstractArray(data_index)
            if not isinstance(vtk_array, vtk.vtkDataArray):
                continue
                
            array_base_name = vtk_array.GetName()
            selected_pattern = None
            for pattern in fields_to_keep:
                if pattern in array_base_name:
                    selected_pattern = pattern
                    break
            if selected_pattern is None:
                continue
                
            numpy_data = numpy_support.vtk_to_numpy(vtk_array).astype(numpy.float64).reshape(num_cells, -1)[valid_cells_mask]
            
            if numpy_data.shape[1] == 1:
                region_data["fields"][selected_pattern] = numpy_data[:, 0]
            else:
                component_labels_list = component_labels(selected_pattern, numpy_data.shape[1], vtk_array, phases)
                for comp_index, label in enumerate(component_labels_list):
                    region_data["fields"][label] = numpy_data[:, comp_index]
                    
        data_chunks.append(region_data)

    if not data_chunks:
        return numpy.empty((0, 3), dtype=numpy.float32), {}

    combined_centers = numpy.concatenate([chunk["centers"] for chunk in data_chunks])
    all_field_names = set().union(*[chunk["fields"].keys() for chunk in data_chunks])
    
    combined_fields = {
        key: numpy.concatenate([chunk["fields"].get(key, numpy.full(len(chunk["centers"]), numpy.nan)) for chunk in data_chunks]) 
        for key in all_field_names
    }

    # Canonical order: sort by rounded (x, y, z)
    rounded_centers = numpy.round(combined_centers, 3)
    sort_order = numpy.lexsort((rounded_centers[:, 2], rounded_centers[:, 1], rounded_centers[:, 0]))
    
    return combined_centers[sort_order], {k: v[sort_order] for k, v in combined_fields.items()}


def read_step(vtm_file_path, regions, phases, fields_to_keep):
    """
    Read all requested data for a single output time step.

    Parameters:
        vtm_file_path (str): The file path to the .vtm file for the time step.
        regions (list): List of regions to process.
        phases (list): List of expected phase names.
        fields_to_keep (list): Array patterns specifying which fields to export.

    Returns:
        tuple: A tuple containing concatenated arrays across all regions:
               - (numpy.ndarray) centers
               - (numpy.ndarray) regions
               - (dict) fields mapping to concatenated 1D data arrays
    """
    region_to_files_map = region_files(vtm_file_path)

    centers_list, region_ids_list, fields_list = [], [], []
    for region_index, region_name in enumerate(regions):
        region_centers, region_fields = read_region(
            region_to_files_map.get(region_name, []), phases, fields_to_keep
        )
        centers_list.append(region_centers)
        region_ids_list.append(numpy.full(len(region_centers), region_index, dtype=numpy.int8))
        fields_list.append(region_fields)
        
    cell_counts = [len(centers) for centers in centers_list]
    unique_field_names = sorted(set().union(*[fields.keys() for fields in fields_list]))
    
    combined_fields = {}
    for field_name in unique_field_names:
        combined_fields[field_name] = numpy.concatenate(
            [fields.get(field_name, numpy.full(count, numpy.nan)) for fields, count in zip(fields_list, cell_counts)]
        )
        
    return (numpy.concatenate(centers_list), numpy.concatenate(region_ids_list), combined_fields)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def select_times(available_steps, requested_years):
    """
    Filter the available simulation steps to match requested output times in years.

    Parameters:
        available_steps (list): List of all available (time_in_seconds, file_path) tuples.
        requested_years (list): List of specific years to export (can be None or empty).

    Returns:
        list: A subset of available_steps that are closest to the requested_years.
    """
    if not requested_years:
        return available_steps
        
    times_in_years = numpy.array([step[0] for step in available_steps]) / SECONDS_PER_YEAR
    picked_steps = []
    
    for target_year in requested_years:
        closest_index = int(numpy.argmin(numpy.abs(times_in_years - target_year)))
        if abs(times_in_years[closest_index] - target_year) > 0.26:   # outputs are typically every 0.5 yr
            print(f"warning: no output within 0.25 yr of {target_year} yr "
                  f"(closest: {times_in_years[closest_index]:.2f} yr)", file=sys.stderr)
        if available_steps[closest_index] not in picked_steps:
            picked_steps.append(available_steps[closest_index])
            
    return sorted(picked_steps, key=lambda step: step[0])


def parse_field_options(entries, convert):
    """
    Parse ``PATTERN=VALUE`` entries into a list of (glob pattern, converted value).
    """
    options = []
    for entry in entries or []:
        pattern, separator, value = entry.rpartition("=")
        if not separator or not pattern:
            sys.exit(f"invalid option '{entry}': expected PATTERN=VALUE")
        try:
            options.append((pattern, convert(value)))
        except ValueError as error:
            sys.exit(f"invalid option '{entry}': {error}")
    return options


def match_field_option(field_name, options):
    """
    Return the value of the first option whose pattern matches the field name, or None.
    """
    for pattern, value in options:
        if fnmatch.fnmatchcase(field_name, pattern):
            return value
    return None


def parse_precision(spec):
    """
    Validate a precision spec: ``bitsN`` (keep N float32 mantissa bits) or a rounding step.
    """
    if spec.startswith("bits"):
        bits = int(spec[4:])
        if not 1 <= bits <= 23:
            raise ValueError("number of mantissa bits must be in [1, 23]")
    elif float(spec) <= 0.0:
        raise ValueError("rounding step must be positive")
    return spec


def reduce_precision(values, spec):
    """
    Lossy precision reduction that makes the data compress much better.

    ``bitsN`` keeps N mantissa bits (relative error < 2**-N); a number rounds the
    values to a multiple of that step (absolute error <= step / 2).  NaN is kept.
    """
    values = numpy.asarray(values, dtype=numpy.float32)
    if spec is None:
        return values
    if spec.startswith("bits"):
        mask = numpy.uint32((0xFFFFFFFF << (23 - int(spec[4:]))) & 0xFFFFFFFF)
        reduced = (values.view(numpy.uint32) & mask).view(numpy.float32)
    else:
        step = float(spec)
        reduced = (numpy.round(values / step) * step).astype(numpy.float32)
    return numpy.where(numpy.isnan(values), numpy.float32(numpy.nan), reduced)


def main(arguments=None):
    """
    Main entry point for the script.

    Parameters:
        arguments (list, optional): Command line arguments. Defaults to sys.argv.

    Returns:
        int: Exit status code (0 for success).
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pvd", help="GEOS VTK index file (e.g. VTK.OUTPUT.pvd)")
    parser.add_argument("-o", "--output", default="dome_soreide_whitson_bench.h5",
                        help="Output HDF5 file path (default: %(default)s)")
    parser.add_argument("--fields", nargs="+", default=list(DEFAULT_FIELDS),
                        help="array-name patterns to export; use '*' for everything "
                             "(default: %(default)s)")
    parser.add_argument("--phases", nargs="+", default=list(DEFAULT_PHASES),
                        help="Phase names for multi-component labeling (default: %(default)s)")
    parser.add_argument("--times-years", nargs="+", type=float, default=None,
                        help="only export the outputs closest to these times [yr]")
    parser.add_argument("--keep", nargs="+", default=None,
                        help="only store the datasets whose final name matches one of these "
                             "glob patterns (e.g. pressure 'phaseCompFraction_gas/CO2')")
    parser.add_argument("--field-times", nargs="+", default=None, metavar="PATTERN=T1,T2,...",
                        help="only store the outputs closest to these times [yr] for the matching "
                             "datasets (e.g. 'phaseCompFraction_gas/CO2=13.5,15.5,18,20')")
    parser.add_argument("--round", nargs="+", default=None, metavar="PATTERN=SPEC",
                        help="LOSSY: reduce the precision of the matching datasets; SPEC is bitsN "
                             "(keep N mantissa bits) or a rounding step (e.g. pressure=bits12 "
                             "'phaseVolumeFraction_gas=1e-4')")
    parser.add_argument("--list-fields", action="store_true",
                        help="print the arrays found in the first output and exit")
    parsed_args = parser.parse_args(arguments)
    field_time_options = parse_field_options(
        parsed_args.field_times, lambda value: [float(year) for year in value.split(",")])
    round_options = parse_field_options(parsed_args.round, parse_precision)

    simulation_steps = select_times(read_pvd(parsed_args.pvd), parsed_args.times_years)
    if not simulation_steps:
        sys.exit("no time step found in " + parsed_args.pvd)

    fields_to_keep = ["*"] if parsed_args.list_fields else parsed_args.fields
    print(f"reading {len(simulation_steps)} time steps from {parsed_args.pvd}")
    
    # Using meaningful variable names instead of the *0 variables
    regions = sorted(region_files(simulation_steps[0][1]).keys())
    initial_centers, initial_regions, initial_fields = read_step(
        simulation_steps[0][1], regions, parsed_args.phases, fields_to_keep
    )
    
    if parsed_args.list_fields:
        print(f"{len(initial_centers)} cells; fields in the first output:")
        for field_name in sorted(initial_fields):
            defined_count = int(numpy.sum(~numpy.isnan(initial_fields[field_name])))
            print(f"  {field_name:45s} defined in {defined_count} cells")
        return 0

    import h5py
    extracted_field_names = sorted(initial_fields)
    if parsed_args.keep:
        extracted_field_names = [name for name in extracted_field_names
                                 if any(fnmatch.fnmatchcase(name, pattern) for pattern in parsed_args.keep)]
    if not extracted_field_names:
        sys.exit("none of the requested fields was found; try --list-fields")
        
    num_times = len(simulation_steps)
    num_cells = len(initial_centers)
    times_in_years = numpy.array([step[0] for step in simulation_steps]) / SECONDS_PER_YEAR

    # Map: field -> {index in /times: row in the field dataset}
    field_rows = {}
    for field_name in extracted_field_names:
        requested_years = match_field_option(field_name, field_time_options)
        if requested_years is None:
            field_rows[field_name] = {index: index for index in range(num_times)}
            continue
        time_indices = sorted({int(numpy.argmin(numpy.abs(times_in_years - year))) for year in requested_years})
        for year in requested_years:
            if numpy.min(numpy.abs(times_in_years - year)) > 0.26:
                print(f"warning: {field_name}: no exported output within 0.25 yr of {year} yr", file=sys.stderr)
        field_rows[field_name] = {index: row for row, index in enumerate(time_indices)}
    precisions = {name: match_field_option(name, round_options) for name in extracted_field_names}

    with h5py.File(parsed_args.output, "w") as hdf5_file:
        hdf5_file.attrs["title"] = "GEOS dome example - CO2 injection into a depleted reservoir"
        hdf5_file.attrs["created"] = datetime.datetime.now().isoformat(timespec="seconds")
        hdf5_file.attrs["phases"] = list(parsed_args.phases)
        hdf5_file.attrs["seconds_per_year"] = SECONDS_PER_YEAR
        hdf5_file.attrs["num_cells"] = num_cells
        
        times_array = numpy.array([step[0] for step in simulation_steps])
        time_group = hdf5_file.create_group("times")
        time_group.create_dataset("seconds", data=times_array)
        time_group.create_dataset("years", data=times_array / SECONDS_PER_YEAR)
        
        # We apply high compression to mesh datasets. VTI-like meshes will compress down significantly.
        mesh_group = hdf5_file.create_group("mesh")
        mesh_group.create_dataset("center", data=initial_centers, compression="gzip", compression_opts=9)
        mesh_group.create_dataset("region", data=initial_regions, compression="gzip", compression_opts=9)
        mesh_group.attrs["region_names"] = list(regions)
        
        fields_group = hdf5_file.create_group("fields")
        field_datasets = {}
        for field_name in extracted_field_names:
            dataset = fields_group.create_dataset(field_name, shape=(len(field_rows[field_name]), num_cells),
                                                  dtype="f4", chunks=(1, num_cells), compression="gzip",
                                                  compression_opts=4, shuffle=True, fillvalue=numpy.nan)
            if len(field_rows[field_name]) != num_times:
                dataset.attrs["time_index"] = numpy.array(sorted(field_rows[field_name]), dtype=numpy.int32)
            if precisions[field_name] is not None:
                dataset.attrs["precision"] = precisions[field_name]
            field_datasets[field_name] = dataset

        def write_fields(time_index, step_fields):
            for field_name in extracted_field_names:
                row = field_rows[field_name].get(time_index)
                if row is not None and field_name in step_fields:
                    field_datasets[field_name][row, :] = reduce_precision(step_fields[field_name],
                                                                          precisions[field_name])

        write_fields(0, initial_fields)

        for time_index in range(1, num_times):
            print(f"  [{time_index + 1:3d}/{num_times}] t = {simulation_steps[time_index][0] / SECONDS_PER_YEAR:7.2f} yr",
                  flush=True)
            step_centers, _, step_fields = read_step(
                simulation_steps[time_index][1], regions, parsed_args.phases, fields_to_keep
            )
            
            if len(step_centers) != num_cells or not numpy.allclose(step_centers, initial_centers, atol=1e-6):
                sys.exit(f"cell layout of {simulation_steps[time_index][1]} differs from the first "
                         "output; the mesh must be fixed in time")

            write_fields(time_index, step_fields)

    print(f"wrote {parsed_args.output}: {num_times} times, {num_cells} cells, {len(extracted_field_names)} fields")
    return 0


if __name__ == "__main__":
    sys.exit(main())
