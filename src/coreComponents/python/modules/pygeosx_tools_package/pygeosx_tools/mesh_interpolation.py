
import numpy as np
from scipy import stats, interpolate


def apply_to_bins(fn, position, value, bins, collapse_edges=True):
  """
  @brief apply a function to values that are located within a series of bins
  @param fn a function that takes a single scalar or array input
  @param position an array describing the location of each sample
  @param value an array of values at each location
  @param bins the bin edges for the position data
  @param collapse_edges a flag that controls the behavior of edge-data (default=True)
  @return binned_values an array of function results for each bin
          Note: if a bin is empty, this function will fill a nan value
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
  @brief fill in any nan values in two 1D arrays by extrapolating
  @param x the position array
  @param y the value array
  @param slope_scale the value to scale the extrapolation slope (default=0.0)
  @return y the value array with nan values replaced by extrapolated data
  """
  Inan = np.isnan(y)
  reg = stats.linregress(x[~Inan], y[~Inan])
  y[Inan] = reg[0] * x[Inan] * slope_scale + reg[1]
  return y


def get_random_realization(x, bins, value=0.0, rand_fill=0, rand_scale=0, slope_scale=0):
  """
  @brief get a random realization for a noisy signal with a set of bins
  @param x the location of data samples
  @param bins the bin edges for the position data
  @param value the value of the data samples
  @param rand_fill the expected standard deviation for missing bins (default=0)
  @param rand_scale the value to scale the standard deviation for the realization (default=0)
  @param slope_scale the value to scale the extrapolation slope (default=0.0)
  @return y_final an ndarary containing the random realization
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


def get_realizations(z, bins, targets):
  """
  @brief get random realizations for noisy signals on target bins
  @param z position information
  @param targets dict of target values to get realizations
  """
  results = {}
  for ka in targets.keys():
    results[ka] = get_random_realization(z, bins, **targets[ka])
  return results


def get_surface_interpolator(fname, z_origin=2700.0):
  """
  @brief parse a surface definition file
  @param fname the path of the target file
  @param z_origin the vertical offset of the file (default=2700.0)
  @return interp_a a linear interpolator for depth
  @return interp_b a nearest-neighbor interpolator for depth
  @return interp_c a nearest-neighbor interpolator for x
  @return interp_d a nearest-neighbor interpolator for y
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
  @brief build a set of conforming tables property, surface interpolators
  @param property_interps a dict of 1D property interpolators with depth
  @param surface_interps a dict of interpolators for surfaces
  @param axes a list of x,y,z axis definitions
  @param extrapolate_surfaces a flag to indicate whether surface elevations should be extrapolated (default=True)
  @return conforming_tables a dict containing tables that conform to the surfaces
  @return surface_elevations a dict containing surface elevations in the xy plane
  @return grid the table grid
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


