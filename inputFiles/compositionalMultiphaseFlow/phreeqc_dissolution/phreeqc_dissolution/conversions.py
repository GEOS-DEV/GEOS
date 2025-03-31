def bar2atm(input_pressure):
    return input_pressure / 1.01325


def bar2pa(input_pressure):
    return input_pressure * 100000


def atm2bar(input_pressure):
    return input_pressure * 1.01325


def ml_min2ft3_d(input_rate):
    return input_rate / 19.664


def convert_rate(input_rate):
    """
    Input in ml/min;
    Output in m3/day.
    """
    return input_rate * 60 * 24 / 100 / 100 / 100


def convert_composition(component_stream, E):
    import numpy as np
    element_stream = np.zeros(E.shape[0])
    for i in range(E.shape[0]):
        element_stream[i] = np.divide(np.sum(np.multiply(E[i], component_stream)),
                                      np.sum(np.multiply(E, component_stream)))
    return element_stream


def correct_composition(composition, comp_min):
    import numpy as np
    mask = np.zeros(len(composition))
    for i in range(len(composition)):
        if composition[i] == 0:
            mask[i] = 1
    factor = np.count_nonzero(mask)
    composition = np.multiply(composition, 1 - factor * comp_min)
    composition += mask * comp_min
    return composition[:-1]


def calculate_injection_stream(q_water, q_co2, temperature, pressure_bar):
    import CoolProp.CoolProp as CP

    # Set up constants
    molar_mass_water = 0.018016     # kg/mol
    molar_mass_co2 = 0.04401        # kg/mol

    # Evaluate ratio
    ratio_co2 = 1
    ratio_water = q_water / q_co2

    # Convert state values
    pressure = bar2pa(pressure_bar)       # Pa

    # Get and densities
    rho_water = CP.PropsSI('D', 'T', temperature, 'P', pressure, 'Water')
    rho_co2 = CP.PropsSI('D', 'T', temperature, 'P', pressure, 'CarbonDioxide')

    # Calculated masses, assume 1 fraction to be 1 m3
    mass_water = ratio_water * rho_water   # kg
    mass_co2 = ratio_co2 * rho_co2         # kg

    # Calculate moles
    mole_water = mass_water / molar_mass_water      # mole
    mole_co2 = mass_co2 / molar_mass_co2            # mole
    return mole_water, mole_co2


def get_mole_fractions(mole_water, mole_co2):
    mole_total = mole_water + mole_co2
    mole_fraction_water = mole_water / mole_total
    mole_fraction_co2 = mole_co2 / mole_total
    return mole_fraction_water, mole_fraction_co2
