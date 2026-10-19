"""
Figures of the "CO2 injection into a depleted reservoir" (dome) example.

Every ``plot_*`` function below creates one figure from the HDF5 file written
by ``dome_vtk_to_hdf5.py``.  They are called by the documentation through the
matplotlib ``plot`` directive of Sphinx, e.g.

    .. plot:: ../inputFiles/.../dome/plot_dome.py plot_pressure

so that no image is stored in the repository: the pictures (including the
rock-type map of the dome) are regenerated from the simulation output every
time the documentation is built.

Where things are read from
--------------------------
* results   : ``dome_results.hdf5`` next to this script, or the file given by
              the environment variable ``DOME_HDF5``;
* wells     : ``dome_mesh_bench.xml`` (``InternalWell`` polylines);
* schedule  : ``dome_soreide_whitson_base.xml`` (the three ``*.RATE`` tables)
              and the gas-water contact (``phaseContacts``).

Reading the wells, the schedule and the contact from the input files keeps the
figures consistent with the input if the example is modified.

Stand-alone use
---------------
    python plot_dome.py                       # write every figure as a PNG
    python plot_dome.py plot_pressure --hdf5 run/dome_results.hdf5 --outdir figs
"""
import argparse
import csv
import os
import re
import sys
import warnings
import xml.etree.ElementTree as ET

import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Set global default font size
plt.rcParams.update({'font.size': 16})

# --------------------------------------------------------------------------
# Directories Configuration
# --------------------------------------------------------------------------
geos_directory = "."
if "sphinx" in sys.modules:
    geos_directory = "../../../../.."
data_directory = os.path.normpath(os.path.abspath(os.path.join(geos_directory, "src/docs/sphinx/basicExamples/depletedReservoir")))
case_directory = os.path.normpath(os.path.abspath(os.path.join(geos_directory, "inputFiles/compositionalMultiphaseWell")))

# --------------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------------
HDF5_NAME = "dome_soreide_whitson_bench.h5"
MESH_XML = "dome_mesh_bench.xml"
PROPS_XML = "dome_properties.xml"
BASE_XML = "dome_soreide_whitson_base.xml"

# (time in years, caption).  The captions follow the schedule of the example:
# 0-8 yr production (+ water injection from 5 yr), 8-13 yr shut-in,
# 13-18 yr gas injection, 18-20 yr shut-in.
TIMES_ALL = [
    (0.0, "initial state"),
    (5.0, "end of depletion"),
    (8.0, "end of water flooding"),
    (13.0, "end of shut-in, start of gas injection"),
    (15.5, "2.5 yr into gas injection"),
    (18.0, "end of gas injection"),
    (20.0, "end of simulation"),
]
TIMES_CO2 = [
    (13.5, "0.5 yr into gas injection"),
    (15.5, "2.5 yr into gas injection"),
    (18.0, "end of gas injection"),
    (20.0, "end of simulation"),
]

ROCK_COLORS = ["#00aaff", "#ff0000", "#00ff00"]   # SAND, CLAY, CAPROCK
ROCK_NAMES = ["SAND", "CLAY", "CAPROCK"]
CAPROCK_GRAY = "#d9d9d9"

# Names of the datasets written by dome_vtk_to_hdf5.py (first match wins)
ALIASES = {
    "pressure": ["pressure"],
    "gas_saturation": ["phaseVolumeFraction_gas", "phaseVolumeFraction_0"],
    "CH4": ["globalCompFraction_CH4"],
    "CO2": ["phaseCompFraction_gas,CO2"],
}


# --------------------------------------------------------------------------
# Input-file helpers
# --------------------------------------------------------------------------
_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def extract_numbers(text):
    """
    Extract all floating point numbers from a string.

    Args:
        text (str): The string containing numerical values.

    Returns:
        numpy.ndarray: An array of the extracted float values.
    """
    return numpy.array([float(value) for value in _NUM.findall(text)])


def read_wells(path=None):
    """
    Read the perforation mid-point (x, y, z) of every well from the mesh file.

    Args:
        path (str, optional): The path to the mesh XML file. Defaults to MESH_XML in case_directory.

    Returns:
        dict: A dictionary mapping well names to their (x, y, z) coordinate tuples.
    """
    file_path = path if path else os.path.join(case_directory, MESH_XML)
    wells_dictionary = {}
    for well_element in ET.parse(file_path).getroot().iter("InternalWell"):
        points = extract_numbers(well_element.attrib["polylineNodeCoords"]).reshape(-1, 3)
        perforation = well_element.find("Perforation")
        distance = float(perforation.attrib["distanceFromHead"]) if perforation is not None else 0.0
        vector = points[1] - points[0]
        wells_dictionary[well_element.attrib["wellRegionName"]] = tuple(points[0] + vector / numpy.linalg.norm(vector) * distance)
    return wells_dictionary


def read_gwc(path=None):
    """
    Read the initial gas-water contact depth from the base XML file.

    Args:
        path (str, optional): The path to the base XML file. Defaults to BASE_XML in case_directory.

    Returns:
        float: The gas-water contact depth.
    """
    file_path = path if path else os.path.join(case_directory, BASE_XML)
    for element in ET.parse(file_path).getroot().iter("HydrostaticEquilibrium"):
        return float(extract_numbers(element.attrib["phaseContacts"])[0])
    raise ValueError("Gas-water contact not found in XML")


def read_table(name, path=None):
    """
    Read the coordinates and values of a TableFunction from the base input file.

    Args:
        name (str): The name of the TableFunction to extract.
        path (str, optional): The path to the base XML file. Defaults to BASE_XML in case_directory.

    Returns:
        tuple: A tuple containing (coordinates_array, values_array).
    """
    file_path = path if path else os.path.join(case_directory, BASE_XML)
    for table_element in ET.parse(file_path).getroot().iter("TableFunction"):
        if table_element.attrib.get("name") == name:
            return extract_numbers(table_element.attrib["coordinates"]), extract_numbers(table_element.attrib["values"])
    raise KeyError(name)


def get_rate_at_time(coordinates, values, time_point):
    """
    Get the interpolated value from a table with piecewise constant interpolation.

    Args:
        coordinates (numpy.ndarray): The x-coordinates (times) of the table.
        values (numpy.ndarray): The y-values (rates) of the table.
        time_point (float): The specific time to interpolate the value for.

    Returns:
        float: The interpolated rate.
    """
    index = numpy.searchsorted(coordinates, time_point, side="right") - 1
    return 0.0 if index < 0 else values[index]


def read_schedule(path=None):
    """
    Read the rate schedules for the producer, water injector, and gas injector.

    Args:
        path (str, optional): The path to the base XML file. Defaults to BASE_XML in case_directory.

    Returns:
        dict: A dictionary mapping 'P', 'W', 'G' to their respective (coordinates, values) tables.
              Returns None if the files or tables cannot be found.
    """
    try:
        return {"P": read_table("PRODUCER.RATE", path),
                "W": read_table("WATER.INJECTOR.RATE", path),
                "G": read_table("GAS.INJECTOR.RATE", path)}
    except (KeyError, OSError):
        return None


PHASE_LABELS = {frozenset(): "shut-in", frozenset("P"): "depletion",
                frozenset("PW"): "water flooding", frozenset("G"): "gas injection"}
PHASE_COLORS = {"shut-in": "#eeeeee", "depletion": "#fddbc7",
                "water flooding": "#d1e5f0", "gas injection": "#d9f0d3"}


def operating_phases(schedule, time_end):
    """
    Determine the operating phases over time based on active rate tables.

    Args:
        schedule (dict): The dictionary of well rate tables.
        time_end (float): The final time point to consider in seconds.

    Returns:
        list: A list of tuples (time_start, time_end, phase_label).
    """
    time_cuts = {0.0, time_end}
    for coordinates, _ in schedule.values():
        time_cuts.update(coord for coord in coordinates if 0.0 < coord < time_end)
    sorted_cuts = sorted(time_cuts)
    phases = []
    for time_start, time_end_segment in zip(sorted_cuts[:-1], sorted_cuts[1:]):
        time_midpoint = 0.5 * (time_start + time_end_segment)
        active_wells = frozenset(key for key, (coords, vals) in schedule.items() if get_rate_at_time(coords, vals, time_midpoint) > 0)
        phases.append((time_start, time_end_segment, PHASE_LABELS.get(active_wells, "+".join(sorted(active_wells)))))
    return phases


# --------------------------------------------------------------------------
# Results access
# --------------------------------------------------------------------------
def get_hdf5_path():
    """
    Determine the location of the HDF5 results file.

    Returns:
        str or None: The path to the HDF5 file if found, otherwise None.
    """
    candidates = [os.environ.get("DOME_HDF5"), os.path.join(data_directory, HDF5_NAME), os.path.join(os.getcwd(), HDF5_NAME)]
    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate
    return None


class Results:
    """Thin reader that exposes the quasi-2D (x, z) section of the dome."""

    def __init__(self, path, target_y_coordinate=None):
        """
        Initialize the Results reader from an HDF5 file.

        Args:
            path (str): Path to the HDF5 file.
            target_y_coordinate (float, optional): The target y-coordinate to slice. Defaults to the median y.
        """
        import h5py
        try:
            self.h5 = h5py.File(path, "r")
        except Exception as e:
            raise RuntimeError(f"Failed to open HDF5 file at {path} {e}") from e
            
        # Handle updated time keys and conversion
        self.time_years = self.h5["times"]["years"][()]

        mesh_group = self.h5["mesh"]
        self.center = mesh_group["center"][:]
        
        # Handle potential missing 'size' and new 'region' key
        self.size = mesh_group["size"][:] if "size" in mesh_group else None
        self.region = mesh_group["region_id"][:] if "region_id" in mesh_group else mesh_group["region"][:]
        
        y_slices = numpy.unique(self.center[:, 1].round(3))
        if target_y_coordinate is None:
            target_y_coordinate = y_slices[len(y_slices) // 2]
        self.selected_y = y_slices[numpy.argmin(numpy.abs(y_slices - target_y_coordinate))]
        self.mask = numpy.abs(self.center[:, 1] - self.selected_y) < 1e-3
        center_x = self.center[self.mask, 0].round(3)
        center_z = self.center[self.mask, 2].round(3)
        self.unique_x, self.unique_z = numpy.unique(center_x), numpy.unique(center_z)
        self._index_x = numpy.searchsorted(self.unique_x, center_x)
        self._index_z = numpy.searchsorted(self.unique_z, center_z)
        self.edges_x, self.edges_z = calculate_edges(self.unique_x), calculate_edges(self.unique_z)

    def get_time_index(self, target_time_years):
        """
        Get the closest time index in the results for a target time.

        Args:
            target_time_years (float): The target time in years.

        Returns:
            int: The index of the closest output time step.
        """
        index = int(numpy.argmin(numpy.abs(self.time_years - target_time_years)))
        if abs(self.time_years[index] - target_time_years) > 0.3:
            warnings.warn(f"no output near {target_time_years} yr; using {self.time_years[index]:.2f} yr")
        return index

    def get_dataset(self, field_key):
        """
        Retrieve a dataset by key, accounting for aliases.

        Args:
            field_key (str): The logical name of the field.

        Returns:
            h5py.Dataset: The HDF5 dataset for the field.
        """
        for name in ALIASES.get(field_key, [field_key]):
            path = name.split(",")
            data_set = self.h5["fields"][path[0]]
            for p in path[1:]:
                data_set = data_set[p]
            if data_set is not None:
                return data_set
        raise KeyError(f"none of {ALIASES.get(field_key, [field_key])} found in the HDF5 file; "
                       "run dome_vtk_to_hdf5.py --list-fields on the VTK output")

    def map_to_grid(self, values_1d):
        """
        Map a 1D array of values to the 2D grid slice.

        Args:
            values_1d (numpy.ndarray): 1D array of values corresponding to the active cells.

        Returns:
            numpy.ndarray: A 2D array of the values on the grid.
        """
        grid = numpy.full((len(self.unique_z), len(self.unique_x)), numpy.nan)
        grid[self._index_z, self._index_x] = values_1d[self.mask]
        return grid

    def get_field_2d(self, field_key, time_index):
        """
        Get a 2D slice of a specific field at a specific time index.

        Args:
            field_key (str): The field to retrieve.
            time_index (int): The time step index.

        Returns:
            numpy.ndarray: 2D array of the requested field.
        """
        return self.map_to_grid(self.get_dataset(field_key)[time_index, :].astype(float))

    def get_region_2d(self):
        """
        Get a 2D slice of the region IDs.

        Returns:
            numpy.ndarray: 2D array of the region IDs.
        """
        return self.map_to_grid(self.region.astype(float))


def calculate_edges(centers):
    """
    Calculate grid edges from cell centers.
    """
    if len(centers) == 1:
        return numpy.array([centers[0] - 0.5, centers[0] + 0.5])
    edges = numpy.empty(len(centers) + 1)
    edges[1:-1] = 0.5 * (centers[1:] + centers[:-1])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def load_results():
    """
    Load the Results object from the configured HDF5 file path.

    Returns:
        Results or None: The Results object, or None if the file could not be found.
    """
    file_path = get_hdf5_path()
    if file_path is None:
        warnings.warn(f"{HDF5_NAME} not found: run dome_vtk_to_hdf5.py on the GEOS "
                      "output or set DOME_HDF5")
        return None
    return Results(file_path, target_y_coordinate=list(read_wells().values())[0][1])

def create_placeholder_figure(text="Results file not found"):
    """
    Create a placeholder matplotlib figure when data is missing.
    """
    figure, axis = plt.subplots(figsize=(8, 1.5))
    axis.text(0.5, 0.5, f"{text}\n(run dome_vtk_to_hdf5.py on the GEOS output "
            f"to create {HDF5_NAME})", ha="center", va="center")
    axis.axis("off")
    return figure

# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def decorate_axes(axis, results_data, wells, gas_water_contact, add_labels=False):
    """
    Draw interfaces, gas-water contact, and wells on one matplotlib axis.
    """
    label_color = "white"

    axis.set_facecolor(CAPROCK_GRAY)                    # inactive cells (CAPROCK)
    axis.contour(results_data.unique_x, results_data.unique_z, results_data.get_region_2d(), levels=[0.5, 1.5], colors="k",
               linewidths=0.4, alpha=0.5)
    axis.axhline(gas_water_contact, color=label_color, ls="--", lw=0.7)
    gas_injector = wells.get("GAS.INJECTOR", wells.get("WATER.INJECTOR"))
    producer_well = wells.get("PRODUCER")

    for well_node, marker_style, label_name in ((gas_injector, "v", "injectors"), (producer_well, "^", "producer")):
        if well_node is not None:
            axis.plot(well_node[0], well_node[2], marker=marker_style, ms=12, mfc="k", mec="w", mew=0.8,
                    zorder=5, clip_on=False)
    if add_labels:
        axis.annotate("gas-water contact", (results_data.edges_x[-1], gas_water_contact), xytext=(-4, 2),
                    textcoords="offset points", ha="right", va="bottom", color=label_color)
        if gas_injector is not None:
            axis.annotate("injectors", (gas_injector[0], gas_injector[2]), xytext=(6, 4),
                        textcoords="offset points", color=label_color)
        if producer_well is not None:
            axis.annotate("producer", (producer_well[0], producer_well[2]), xytext=(-6, 4),
                        textcoords="offset points", ha="right", color=label_color)
    axis.set_xlim(results_data.edges_x[0], results_data.edges_x[-1])
    axis.set_ylim(results_data.edges_z[0], results_data.edges_z[-1])
    axis.set_yticks([-4680, -4400, -4000])

def get_fraction_limits(arrays):
    """
    Calculate appropriate visual limits for fraction data.
    """
    max_value = numpy.nanmax([numpy.nanmax(array_data) for array_data in arrays])
    max_value = 1.0 if not numpy.isfinite(max_value) or max_value <= 0 else min(1.0, numpy.ceil(max_value * 20) / 20)
    return 0.0, max_value

def plot_stacked_fields(field_key, target_times, colormap_name, colorbar_label, use_diverging_norm=False, show_as_change=False,
           is_fraction=False):
    """
    Create one figure with key times stacked vertically sharing one color scale.
    """
    results_data = load_results()
    if results_data is None:
        return create_placeholder_figure()
    wells, gas_water_contact = read_wells(), read_gwc()
    time_indices = [results_data.get_time_index(time_val) for time_val, _ in target_times]
    arrays_to_plot = [results_data.get_field_2d(field_key, index) for index in time_indices]
    
    if field_key == "pressure":
        arrays_to_plot = [array_data / 1e6 for array_data in arrays_to_plot]                  # Pa -> MPa
    if show_as_change:
        reference_array = results_data.get_field_2d(field_key, results_data.get_time_index(0.0)) / 1e6
        arrays_to_plot = [array_data - reference_array for array_data in arrays_to_plot]

    if use_diverging_norm:
        max_abs_value = numpy.nanmax([numpy.nanmax(numpy.abs(array_data)) for array_data in arrays_to_plot]) or 1.0
        plot_arguments = dict(norm=TwoSlopeNorm(0.0, -max_abs_value, max_abs_value))
    elif is_fraction:
        vmin, vmax = get_fraction_limits(arrays_to_plot)
        plot_arguments = dict(vmin=vmin, vmax=vmax)
    else:
        plot_arguments = dict(vmin=numpy.nanmin([numpy.nanmin(array_data) for array_data in arrays_to_plot]),
                  vmax=numpy.nanmax([numpy.nanmax(array_data) for array_data in arrays_to_plot]))

    num_subplots = len(target_times)
    figure, axes_list = plt.subplots(num_subplots, 1, figsize=(10, 1.8 * num_subplots + 0.8), sharex=True,
                             layout="constrained")
    axes_list = numpy.atleast_1d(axes_list)
    for subplot_index, (axis, (time_val, caption), array_data, time_index) in enumerate(zip(axes_list, target_times, arrays_to_plot, time_indices)):
        color_mesh = axis.pcolormesh(results_data.edges_x, results_data.edges_z, array_data, cmap=colormap_name, shading="flat",
                            rasterized=True, **plot_arguments)
        decorate_axes(axis, results_data, wells, gas_water_contact, add_labels=(subplot_index == 0))
        axis.set_title(f"t = {results_data.time_years[time_index]:g} yr: {caption}", loc="left", pad=2)
    axes_list[-1].set_xlabel("x (m)")
    figure.supylabel("elevation z (m)")
    figure.colorbar(color_mesh, ax=list(axes_list), label=colorbar_label, pad=0.01, shrink=0.8, aspect=30)
    plt.show()
    return figure

# --------------------------------------------------------------------------
# Public figures (one per ``.. plot::`` directive)
# --------------------------------------------------------------------------
def plot_rock_types():
    """
    Rock-type map of the dome, rebuilt from the cell regions of the output.
    """
    results_data = load_results()
    if results_data is None:
        return create_placeholder_figure()
    wells, gas_water_contact = read_wells(), read_gwc()
    figure, axis = plt.subplots(figsize=(10, 4.2), layout="constrained")
    axis.pcolormesh(results_data.edges_x, results_data.edges_z, results_data.get_region_2d(), cmap=ListedColormap(ROCK_COLORS),
                  vmin=-0.5, vmax=2.5, shading="flat", rasterized=True)
    decorate_axes(axis, results_data, wells, gas_water_contact, add_labels=True)
    axis.set_xlabel("x (m)")
    axis.set_ylabel("elevation z (m)")
    axis.legend(handles=[Patch(color=color_code, label=name) for color_code, name in zip(ROCK_COLORS, ROCK_NAMES)]
              + [Line2D([], [], color="k", ls="--", lw=0.7, label="gas-water contact")],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, frameon=False)
    plt.show()
    return figure


def plot_pressure():
    """
    Absolute pressure (MPa) at the key times.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    return plot_stacked_fields("pressure", TIMES_ALL, "jet", "pressure (MPa)")


def plot_pressure_change():
    """
    Pressure change with respect to the initial state (MPa).

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    return plot_stacked_fields("pressure", TIMES_ALL[1:], "jet",
                  "pressure change from t = 0 (MPa)", use_diverging_norm=True, show_as_change=True)


def plot_gas_saturation():
    """
    Gas saturation at the key times.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    return plot_stacked_fields("gas_saturation", TIMES_ALL, "jet", "gas saturation (-)",
                  is_fraction=True)


def plot_co2_fraction():
    """
    CO2 global mole fraction during and after the gas injection.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    return plot_stacked_fields("CO2", TIMES_CO2, "jet", "CO$_2$ gas mole fraction (-)",
                  is_fraction=True)


def plot_ch4_fraction():
    """
    CH4 global mole fraction at the key times.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    return plot_stacked_fields("CH4", TIMES_ALL, "jet", "CH$_4$ global mole fraction (-)",
                  is_fraction=True)

def _read_well_csv(filename):
    """
    Helper to read well rates and BHP from a GEOS control CSV.
    Converts time to years and pressure to MPa.
    """
    path = os.path.join(data_directory, filename)
    if not os.path.isfile(path):
        return None
    
    times, bhps, gas_rates, water_rates = [], [], [], []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        if not headers:
            return None
        
        # Make lookup case-insensitive and robust against trailing spaces
        lower_headers = {h: h.strip().lower() for h in headers}
        
        # Identify columns
        t_key = next((h for h in headers if 'time' in lower_headers[h]), None)
        bhp_key = next((h for h in headers if 'bottom hole pressure' in lower_headers[h] or 'bhp' in lower_headers[h]), None)
        gas_key = next((h for h in headers if 'phase0 surface volumetric' in lower_headers[h] and 'rate' in lower_headers[h]), None)
        water_key = next((h for h in headers if 'phase1 surface volumetric' in lower_headers[h] and 'rate' in lower_headers[h]), None)
        
        for row in reader:
            try:
                times.append(float(row[t_key]) if t_key and row[t_key] else 0.0)
                bhps.append(float(row[bhp_key]) if bhp_key and row[bhp_key] else 0.0)
                gas_rates.append(numpy.abs(float(row[gas_key])) if gas_key and row[gas_key] else 0.0)
                water_rates.append(numpy.abs(float(row[water_key])) if water_key and row[water_key] else 0.0)
            except ValueError:
                pass # Skip rows with problematic data
                
    return {
        'time': numpy.array(times) / 31557600.0, # seconds to years
        'bhp': numpy.array(bhps) / 1e6,          # Pa to MPa
        'gas_rate': numpy.array(gas_rates),
        'water_rate': numpy.array(water_rates)
    }

def plot_well(well_name):
    """
    Well surface volumetric rates read from CSV files.
    Plots producer and gas injector (gas rates) on left axis.
    Plots water injector (water rate) on right axis.
    """
    fig, ax_p = plt.subplots(figsize=(10, 5), layout="constrained")
   
    data_well = _read_well_csv(f"{well_name}.csv")  
    if not data_well:
        return create_placeholder_figure("Well control CSV files not found")
        
    is_producer = "PRODUCER" in well_name.upper()
    is_gas = "GAS" in well_name.upper()
    ax_w, ax_g = None, None
    if is_producer or not is_gas:
        ax_w = ax_p.twinx()
    if is_producer or is_gas:
        ax_g = ax_p.twinx()

    # Plotting gas rates on the left axis
    ax_p.plot(data_well['time'], data_well['bhp'], color='k', lw=2, label="BHP")
    if ax_g:
        ax_g.plot(data_well['time'], data_well['gas_rate'], color='r', lw=2, label="Gas Rate")
    if ax_w:
        ax_w.plot(data_well['time'], data_well['water_rate'], color='b', lw=2, label="Water Rate")
        
    ax_p.set_xlabel("time (years)")
    ax_p.set_ylabel("BHP (MPa)")
    
    if ax_w:
        ax_w.set_ylabel("water rate (m$^3$/s)", color='b')
        ax_w.tick_params(axis='y', labelcolor='b')
    if ax_g:
        ax_g.set_ylabel("gas rate (m$^3$/s)", color='r')
        ax_g.tick_params(axis='y', labelcolor='r')
    if ax_w and ax_g:
        ax_w.spines['right'].set_position(('outward', 60))

    # Combine legends from both axes
    lines_1, labels_1 = ax_p.get_legend_handles_labels()
    lines_2, labels_2 = ax_w.get_legend_handles_labels() if ax_w else ([], [])
    lines_3, labels_3 = ax_g.get_legend_handles_labels() if ax_g else ([], [])
    ax_p.legend(lines_1 + lines_2 + lines_3, labels_1 + labels_2 + labels_3, loc=None)
    
    ax_p.grid(alpha=0.3, lw=0.4)
    plt.show()
    return fig

def plot_well_producer():
    return plot_well("PRODUCER.CONTROL")

def plot_well_water_injector():
    return plot_well("WATER.INJECTOR.CONTROL")

def plot_well_gas_injector():
    return plot_well("GAS.INJECTOR.CONTROL")

def get_field_scale(name, path=None):
    """
    Get the ``scale`` of a FieldSpecification of the properties file.

    Args:
        name (str): Name of the field specification.
        path (str, optional): The path to the properties XML file. Defaults to PROPS_XML in case_directory.

    Returns:
        float: The scale value.
    """
    file_path = path if path else os.path.join(case_directory, PROPS_XML)
    for element in ET.parse(file_path).getroot().iter("FieldSpecification"):
        if element.attrib.get("name") == name:
            return float(element.attrib["scale"])
    raise KeyError(name)


def plot_saturation_functions():
    """
    Relative permeabilities and capillary pressures of SAND and CLAY.

    Built from the tables of the input file, no simulation output needed.  The
    clay capillary curve is a J-function; it is converted to pressure with
    J * sigma * sqrt(porosity / permeability), using the surface tension of the
    model, the clay porosity and its permeability in the X direction.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    surface_tension = None
    for element in ET.parse(os.path.join(case_directory, BASE_XML)).getroot().iter("JFunctionCapillaryPressure"):
        surface_tension = float(element.attrib["wettingNonWettingSurfaceTension"])
    
    porosity = get_field_scale("POROSITY.CLAY")
    permeability_x = get_field_scale("PERMX.CLAY")
    j_function_to_pascals = surface_tension * numpy.sqrt(porosity / permeability_x)

    figure, (axis_1, axis_2) = plt.subplots(1, 2, figsize=(14, 3.6), layout="constrained")
    for rock_type, line_style in (("SAND", "-"), ("CLAY", "--")):
        sat_water, rel_perm_water = read_table(f"KRW.{rock_type}")
        sat_gas, rel_perm_gas = read_table(f"KRG.{rock_type}")
        axis_1.plot(sat_water, rel_perm_water, line_style, color="tab:blue", label=f"water, {rock_type}")
        axis_1.plot(1.0 - sat_gas, rel_perm_gas, line_style, color="tab:red", label=f"gas, {rock_type}")
        sat_water_pc, cap_pressure = read_table(f"PCGW.{rock_type}")
        axis_2.plot(sat_water_pc, cap_pressure * (j_function_to_pascals if rock_type == "CLAY" else 1.0) / 1e6, line_style,
                color="k", label=rock_type)
    axis_1.set_xlabel("water saturation (-)")
    axis_1.set_ylabel("relative permeability (-)")
    axis_1.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    axis_2.set_xlabel("water saturation (-)")
    axis_2.set_ylabel("gas-water\ncapillary pressure (MPa)")
    axis_2.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    for target_axis in (axis_1, axis_2):
        target_axis.set_xlim(0, 1)
        target_axis.grid(alpha=0.3, lw=0.4)
    return figure

ALL_FIGURES = ["plot_rock_types", "plot_saturation_functions",
               "plot_well_producer", "plot_well_water_injector", "plot_well_gas_injector",
               "plot_pressure", "plot_pressure_change", "plot_gas_saturation",
               "plot_co2_fraction", "plot_ch4_fraction"]

def run_command_line_interface(argv=None):
    """
    Parse arguments and run the figure generation from the command line.
    """
    argument_parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    argument_parser.add_argument("figures", nargs="*", help=f"any of {ALL_FIGURES} (default: all)")
    argument_parser.add_argument("--hdf5", help="results file (default: %s)" % HDF5_NAME)
    argument_parser.add_argument("--outdir", default=".")
    argument_parser.add_argument("--dpi", type=int, default=150)
    arguments = argument_parser.parse_args(argv)
    
    if arguments.hdf5:
        os.environ["DOME_HDF5"] = arguments.hdf5
    os.makedirs(arguments.outdir, exist_ok=True)
    
    for figure_name in arguments.figures or ALL_FIGURES:
        figure = globals()[figure_name]()
        output_path = os.path.join(arguments.outdir, f"{figure_name}.png")
        figure.savefig(output_path, dpi=arguments.dpi, bbox_inches="tight")
        plt.close(figure)
        print("wrote", output_path)

# The Sphinx plot directive runs this file with __name__ == "__main__", so the
# command-line entry point is disabled when Sphinx is loaded.
if __name__ == "__main__" and "sphinx" not in sys.modules:
    import matplotlib
    matplotlib.use("Agg")
    run_command_line_interface()
