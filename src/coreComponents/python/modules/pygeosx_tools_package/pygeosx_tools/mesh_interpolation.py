
import numpy as np
from scipy import stats, interpolate


def apply_to_bins(fn, position, value, bins, collapse_edges=True):
    """
    Apply a function to values that are located within a series of bins
    Note: if a bin is empty, this function will fill a nan value

    Args:
        fn (function): Function that takes a single scalar or array input
        position (np.ndarray): A 1D list/array describing the location of each sample
        value (np.ndarray): A 1D list/array of values at each location
        bins (np.ndarray): The bin edges for the position data
        collapse_edges (bool): Controls the behavior of edge-data (default=True)

    Returns:
        np.ndarray: an array of function results for each bin
    """
    # Sort values into bins
    Nr = len(bins) + 1
    Ibin = np.digitize(position, bins)

    if collapse_edges:
        Nr -= 2
        Ibin -= 1
        Ibin[Ibin == -1] = 0
        Ibin[Ibin == Nr] = Nr-1

    # Apply functions to bins
    binned_values = np.zeros(Nr)
    for ii in range(Nr):
        tmp = (Ibin == ii)
        if np.sum(tmp):
            binned_values[ii] = fn(value[tmp])
        else:
            # Empty bin
            binned_values[ii] = np.NaN

    return binned_values


def extrapolate_nan_values(x, y, slope_scale=0.0):
    """
    Fill in any nan values in two 1D arrays by extrapolating

    Args:
        x (np.ndarray): 1D list/array of positions
        y (np.ndarray): 1D list/array of values
        slope_scale (float): value to scale the extrapolation slope (default=0.0)
        
    Returns:
        np.ndarray: The input array with nan values replaced by extrapolated data
    """
    Inan = np.isnan(y)
    reg = stats.linregress(x[~Inan], y[~Inan])
    y[Inan] = reg[0] * x[Inan] * slope_scale + reg[1]
    return y


def get_random_realization(x, bins, value, rand_fill=0, rand_scale=0, slope_scale=0):
    """
    Get a random realization for a noisy signal with a set of bins

    Args:
        x (np.ndarray): 1D list/array of positions
        bins (np.ndarray): 1D list/array of bin edges
        value (np.ndarray): 1D list/array of values
        rand_fill (float): The standard deviation to use where data is not defined (default=0)
        rand_scale (float): Value to scale the standard deviation for the realization (default=0)
        slope_scale (float): Value to scale the extrapolation slope (default=0.0)

    Returns:
        np.ndarray: An array containing the random realization
    """
    y_mean = apply_to_bins(np.mean, x, value, bins)
    y_std = apply_to_bins(np.std, x, value, bins)

    # Extrapolate to fill the upper/lower bounds
    x_mid = bins[:-1] + 0.5 * (bins[1] - bins[0])
    y_mean = extrapolate_nan_values(x_mid, y_mean, slope_scale)
    y_std[np.isnan(y_std)] = rand_fill

    # Add a random perturbation to the target value to match missing high/lows
    y_final = y_mean + (rand_scale * y_std * np.random.randn(len(y_mean)))
    return y_final


def get_realizations(x, bins, targets):
    """
    Get random realizations for noisy signals on target bins

    Args:
        x (np.ndarray): 1D list/array of positions
        bins (np.ndarray): 1D list/array of bin edges
        targets (dict): Dict of geosx target keys, inputs to get_random_realization

    Returns:
        dict: Dictionary of random realizations

    """
    results = {}
    for ka in targets.keys():
        results[ka] = get_random_realization(x, bins, **targets[ka])
    return results


def get_surface_interpolator(fname, z_origin=2700.0):
    """
    Parse a surface definition file

    Args:
        fname (str): Path of the target file
        z_origin (float): Vertical offset of a target mesh file (default=2700.0)

    Returns:
        tuple: Surface interpolators (linear-depth, nearest-depth, linear-x, nearest-x)
    """
    nodes = []
    with open(fname, 'r') as f:
        for line in f:
            if ('#' not in line):
                tmp = line.split(',')
                if ('.' in line):
                    nodes.append([float(x) for x in tmp])

    # Convert values
    nodes = np.array(nodes)
    x = np.ascontiguousarray(np.squeeze(nodes[:, 2])) / 3.281
    y = np.ascontiguousarray(np.squeeze(nodes[:, 0])) / 3.281
    z = (np.ascontiguousarray(np.squeeze(nodes[:, 1])) - z_origin) / -3.281

    # Generate four interpolators to reconstruct the surface
    interp_a = interpolate.LinearNDInterpolator((x, y), z)
    interp_b = interpolate.NearestNDInterpolator((x, y), z)
    interp_c = interpolate.NearestNDInterpolator((x, y), x)
    interp_d = interpolate.NearestNDInterpolator((x, y), y)

    return interp_a, interp_b, interp_c, interp_d


def build_conforming_tables(property_interps, surface_interps, axes, extrapolate_surfaces=True):
    """
    Build a set of conforming tables property, surface interpolators

    Args:
        property_interps a dict of 1D property interpolators with depth
        surface_interps a dict of interpolators for surfaces
        axes a list of x,y,z axis definitions
        extrapolate_surfaces a flag to indicate whether surface elevations should be extrapolated (default=True)

    Returns:
        tuple: Conforming tables, surface elevations at grid, grid values
    """
    # Sort the surfaces by depth at the origin
    surfaces = list(surface_interps.keys())
    z_pilot = [surface_interps[ka][1](0.0, 0.0) for ka in surfaces]
    Ia = np.argsort(z_pilot)
    z_pilot = np.array(z_pilot)[Ia]
    surfaces = np.array(surfaces)[Ia]

    # Interpolate values onto a new grid
    x_grid, y_grid, z_grid = axes
    zo = np.concatenate([[-1e6], z_pilot, [0.0]], axis=0)
    X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')
    conforming_tables = {k: np.zeros(np.shape(X)) for k in property_interps.keys()}
    surface_elevations = np.zeros((len(x_grid), len(y_grid), len(surfaces)))

    for ii in range(len(x_grid)):
        for jj in range(len(y_grid)):
            for kk, s in zip(range(len(surfaces)), surfaces):
                # Find the current layer positions
                tmp_z = surface_interps[s][0](x_grid[ii], y_grid[jj])
                if np.isnan(tmp_z):
                    # The point is out of the surface xy range
                    tmp_z = surface_interps[s][1](x_grid[ii], y_grid[jj])
                    if extrapolate_surfaces:
                        # Find nearest point, add in the average delta_z from other layers
                        nearest_x = surface_interps[s][2](x_grid[ii], y_grid[jj])
                        nearest_y = surface_interps[s][3](x_grid[ii], y_grid[jj])
                        other_delta_z = []
                        for kk in range(0, len(surfaces)):
                            tmp_dz = surface_interps[s][0](x_grid[ii], y_grid[jj]) - surface_interps[s][0](nearest_x, nearest_y)
                            if (not np.isnan(tmp_dz)):
                                other_delta_z.append(tmp_dz)
                        if len(other_delta_z):
                            tmp_z += np.mean(other_delta_z)
                # Save the final surface elevation
                surface_elevations[ii, jj, kk] = tmp_z

            # Map the origin/current surface elevations
            za = np.concatenate(([-1e6], np.squeeze(surface_elevations[ii, jj, :]), [0.0]), axis=0)
            zb = interpolate.interp1d(za, zo, kind='linear')(z_grid)

            # Interplate properties
            for ka in property_interps:
                conforming_tables[ka][ii, jj, :] = property_interps[ka](zb)

    return conforming_tables, surface_elevations, (X, Y, Z)


