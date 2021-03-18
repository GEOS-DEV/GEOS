
import numpy as np


def parse_las(fname, variable_start='~CURVE INFORMATION', body_start='~A'):
    """
    Parse an las format log file

    Args:
        fname (str): Path to the log file
        variable_start (str): A string that indicates the start of variable header information (default = '~CURVE INFORMATION')
        body_start (str): a string that indicates the start of the log body (default = '~A')

    Returns:
        np.ndarray: a dict containing the values and unit definitions for each variable in the log
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
    Convert young's modulus and poisson's ratio to bulk and shear modulus

    Args:
        E (float, np.ndarray): Young's modulus
        nu (float, np.ndarray): Poisson's ratio

    Returns:
        tuple: bulk modulus, shear modulus with same size as inputs
    """
    K = E / (3.0 * (1 - 2.0 * nu))
    G = E / (2.0 * (1 + nu))
    return K, G


def estimate_shmin(z, rho, nu):
    """
    Estimate the minimum horizontal stress using the poisson's ratio

    Args:
        z (float, np.ndarray): Depth
        rho (float, np.ndarray): Density
        nu (float, np.ndarray): Poisson's ratio

    Returns:
        float: minimum horizontal stress
    """
    k = nu / (1.0 - nu)
    sigma_h = k * rho * 9.81 * z
    return sigma_h

