
import numpy as np


def parse_las(fname, variable_start='~CURVE INFORMATION', body_start='~A'):
  variable_names = []
  variable_units = []
  variable_values = []

  with open(fname) as f:
    file_location = 0
    for line in f:
      # Preamble
      if (file_location == 0):
        if variable_start in line:
          file_location += 1

      # Variable definitions
      elif (file_location == 1):
        if ('#' not in line):
          if body_start in line:
            file_location += 1
          else:
            tmp = line.split()
            variable_names.append(tmp[0])
            variable_units.append(tmp[1])

      # Body
      else:
        variable_values.append([float(x) for x in line.split()])

  # Format results
  results = {}
  variable_values = np.array(variable_values)

  for ii in range(0, len(variable_names)):
    results[variable_names[ii]] = {'units': variable_units[ii],
                                   'values': np.squeeze(variable_values[:, ii])}

  return results


def convert_E_nu_to_K_G(E, nu):
  K = E / (3.0 * (1 - 2.0 * nu))
  G = E / (2.0 * (1 + nu))
  return K, G


def estimate_shmin(z, rho, nu):
  k = nu / (1.0 - nu)
  sigma_h = k * rho * 9.81 * z
  return sigma_h

