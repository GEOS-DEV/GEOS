
import numpy as np
from scipy import stats, interpolate


def apply_to_bins(fn, position, value, bins, collapse_edges=True):
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


def extrapolate_nan_values(x, y, slope_scale):
  Inan = np.isnan(y)
  reg = stats.linregress(x[~Inan], y[~Inan])
  y[Inan] = reg[0] * x[Inan] * slope_scale + reg[1]
  return y


def bin_extrapolate_values(x, x_mid, bins, value=0.0, rand_fill=0, rand_scale=0, slope_scale=0):
  y_mean = apply_to_bins(np.mean, x, value, bins)
  y_std = apply_to_bins(np.std, x, value, bins)

  # Extrapolate to fill the upper/lower bounds
  y_mean = extrapolate_nan_values(x_mid, y_mean, slope_scale)
  y_std[np.isnan(y_std)] = rand_fill

  # Add a random perturbation to the target value to match missing high/lows
  y_final = y_mean + (rand_scale * y_std * np.random.randn(len(y_mean)))
  return y_final


def mesh_bin_values(z,
                    bins,
                    targets):
  """
  @brief Bin values onto the target mesh
  @param z location of the grid nodes
  @param targets dict of target values to conform
  """
  z_mid = bins[:-1] + 0.5 * (bins[1] - bins[0])
  results = {}
  for ka in targets.keys():
    results[ka] = bin_extrapolate_values(z, z_mid, bins,
                                         **targets[ka])
  return results


def get_surface_interpolator(fname,
                             z_origin=2700.0):
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


def build_conforming_tables(property_interps, surface_interps, x_grid, y_grid, z_grid, extrapolate_surfaces=True):
  # Sort the surfaces by depth at the origin
  surfaces = list(surface_interps.keys())
  z_pilot = [surface_interps[ka][1](0.0, 0.0) for ka in surfaces]
  Ia = np.argsort(z_pilot)
  z_pilot = np.array(z_pilot)[Ia]
  surfaces = np.array(surfaces)[Ia]

  # Interpolate values onto a new grid
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


