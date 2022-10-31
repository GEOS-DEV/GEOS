import numpy as np
import re


def parse_las(fname, variable_start='~C', body_start='~A'):
    """
    Parse an las format log file

    Args:
        fname (str): Path to the log file
        variable_start (str): A string that indicates the start of variable header information (default = '~CURVE INFORMATION')
        body_start (str): a string that indicates the start of the log body (default = '~A')

    Returns:
        np.ndarray: a dict containing the values and unit definitions for each variable in the log
    """
    results = {}
    variable_order = []

    # The expected format of the varible definition block is:
    # name.units code:description
    variable_regex = re.compile('\s*([^\.^\s]*)\s*(\.[^ ]*) ([^:]*):(.*)')

    with open(fname) as f:
        file_location = 0
        for line in f:
            line = line.split('#')[0]
            if line:
                # Preamble
                if (file_location == 0):
                    if variable_start in line:
                        file_location += 1

                # Variable definitions
                elif (file_location == 1):
                    # This is not a comment line
                    if body_start in line:
                        file_location += 1
                    else:
                        match = variable_regex.match(line)
                        if match:
                            variable_order.append(match[1])
                            results[match[1]] = {
                                'units': match[2][0:],
                                'code': match[3],
                                'description': match[4],
                                'values': []
                            }
                        else:
                            # As a fall-back use the full line
                            variable_order.append(line[:-1])
                            results[line[:-1]] = {'units': '', 'code': '', 'description': '', 'values': []}

                # Body
                else:
                    for k, v in zip(variable_order, line.split()):
                        results[k]['values'].append(float(v))

    # Convert values to numpy arrays
    for k in results:
        results[k]['values'] = np.array(results[k]['values'])

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
