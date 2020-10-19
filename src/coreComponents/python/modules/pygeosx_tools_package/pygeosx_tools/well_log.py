
import numpy as np


def parse_las(fname, variable_start='~CURVE INFORMATION', body_start='~A'):
  """
  @brief parse an las format log file
  @param fname the path to the log file
  @param variable_start a string that indicates the start of variable header information (default = '~CURVE INFORMATION')
  @param variable_start a string that indicates the start of the log body (default = '~A')
  @return a dict containing the values and unit definitions for each variable in the log
  """
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
  """
  @brief convert young's modulus and poisson's ratio to bulk and shear modulus
         the input/output parameters can be either scalars or arrays
  @param E the young's modulus
  @param nu the poisson's ratio
  @return K the bulk modulus
  @return G the shear modulus
  """
  K = E / (3.0 * (1 - 2.0 * nu))
  G = E / (2.0 * (1 + nu))
  return K, G


def estimate_shmin(z, rho, nu):
  """
  @brief estimate the minimum horizontal stress using the poisson's ratio
         the input/output parameters can be either scalars or arrays
  @param z the depth
  @param rho the density
  @return nu the poisson's ratio
  @return sigma_h the minimum horizontal stress
  """
  k = nu / (1.0 - nu)
  sigma_h = k * rho * 9.81 * z
  return sigma_h

