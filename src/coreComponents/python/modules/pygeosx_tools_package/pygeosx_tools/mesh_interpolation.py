import numpy as np
from scipy import stats


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
        Ibin[Ibin == Nr] = Nr - 1

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


def get_random_realization(x,
                           bins,
                           value,
                           rand_fill=0,
                           rand_scale=0,
                           slope_scale=0):
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
    for k, t in targets.items():
        results[k] = get_random_realization(x, bins, **t)
    return results
