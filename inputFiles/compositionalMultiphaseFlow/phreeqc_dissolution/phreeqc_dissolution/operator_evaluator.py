from darts.physics.base.operators_base import OperatorsBase
from phreeqc_dissolution.conversions import bar2pa
from darts.engines import *
import CoolProp.CoolProp as CP
import os.path as osp
import numpy as np

physics_name = osp.splitext(osp.basename(__file__))[0]

# Reservoir operators working with the following state:
# state: (pressure, overall mineral molar fractions, overall fluid molar fractions)
class my_own_acc_flux_etor(OperatorsBase):
    def __init__(self, input_data, properties):
        super().__init__(properties, thermal=properties.thermal)
        # Store your input parameters in self here, and initialize other parameters here in self
        self.input_data = input_data
        self.min_z = input_data.min_z
        self.temperature = input_data.temperature
        self.exp_w = input_data.exp_w
        self.exp_g = input_data.exp_g
        self.kin_fact = input_data.kin_fact
        self.property = properties
        self.counter = 0

        # Operator order
        self.ACC_OP = 0  # accumulation operator - ne
        self.FLUX_OP = self.ACC_OP + self.ne  # flux operator - ne * nph
        self.UPSAT_OP = self.FLUX_OP + self.ne * self.nph  # saturation operator (diffusion/conduction term) - nph
        self.GRAD_OP = self.UPSAT_OP + self.nph  # gradient operator (diffusion/conduction term) - ne * nph
        self.KIN_OP = self.GRAD_OP + self.ne * self.nph  # kinetic operator - ne

        # extra operators
        self.GRAV_OP = self.KIN_OP + self.ne  # gravity operator - nph
        self.PC_OP = self.GRAV_OP + self.nph  # capillary operator - nph
        self.PORO_OP = self.PC_OP + self.nph  # porosity operator - 1
        self.ENTH_OP = self.PORO_OP + 1  # enthalpy operator - nph
        self.TEMP_OP = self.ENTH_OP + self.nph  # temperature operator - 1
        self.PRES_OP = self.TEMP_OP + 1
        self.n_ops = self.PRES_OP + 1

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def get_overall_composition(self, state):
        if self.thermal:
            z = state[1:-1]
        else:
            z = state[1:]
        z = np.append(z, 1 - np.sum(z[self.property.flash_ev.fc_mask[:-1]]))
        return z

    def evaluate(self, state, values):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # state and values numpy vectors:
        state_np = state.to_numpy()
        values_np = values.to_numpy()

        # pore pressure
        pressure = state_np[0]
        # get overall molar composition
        z = self.get_overall_composition(state_np)

        # call flash:
        nu_v, x, y, rho_phases, kin_state, _, _ = self.property.flash_ev.evaluate(state_np)
        nu_s = state_np[1]
        nu_v = nu_v * (1 - nu_s) # convert to overall molar fraction
        nu_a = 1 - nu_v - nu_s

        # molar densities in kmol/m3
        rho_a, rho_v = rho_phases['aq'], rho_phases['gas']
        rho_s = self.property.density_ev['solid'].evaluate(pressure) / self.property.Mw['Solid']

        # viscosities
        mu_a = CP.PropsSI('V', 'T', self.temperature, 'P|liquid', bar2pa(pressure), 'Water') * 1000
        try:
            mu_v = CP.PropsSI('V', 'T', self.temperature, 'P|gas', bar2pa(pressure), 'CarbonDioxide') * 1000
        except ValueError:
            mu_v = 0.05
        
        # Get saturations
        if nu_v > 0:
            sv = nu_v / rho_v / (nu_v / rho_v + nu_a / rho_a + nu_s / rho_s)
            sa = nu_a / rho_a / (nu_v / rho_v + nu_a / rho_a + nu_s / rho_s)
            ss = nu_s / rho_s / (nu_v / rho_v + nu_a / rho_a + nu_s / rho_s)
        else:
            sv = 0
            sa = nu_a / rho_a / (nu_a / rho_a + nu_s / rho_s)
            ss = nu_s / rho_s / (nu_a / rho_a + nu_s / rho_s)

        # Need to normalize to get correct Brook-Corey relative permeability
        sa_norm = sa / (sv + sa)
        sv_norm = sv / (sv + sa)

        kr_a = self.property.rel_perm_ev['liq'].evaluate(sa_norm)
        kr_v = self.property.rel_perm_ev['gas'].evaluate(sv_norm)

        # all properties are in array, and can be separate
        self.x = np.array([y, x])
        self.rho_m = np.array([rho_v, rho_a])
        self.kr = np.array([kr_v, kr_a])
        self.mu = np.array([mu_v, mu_a])
        self.compr = self.property.rock_compr_ev.evaluate(pressure)
        self.sat = np.array([sv_norm, sa_norm])

        # Densities
        rho_t = rho_a * sa + rho_s * ss + rho_v * sv
        rho_f = rho_a * sa_norm + rho_v * sv_norm

        # Kinetic reaction rate
        kin_rate = self.property.kinetic_rate_ev.evaluate(kin_state, ss, rho_s, self.min_z, self.kin_fact)

        nc = self.property.nc
        nph = 2
        ne = nc

        """ CONSTRUCT OPERATORS HERE """
        values_np[:] = 0.
        
        """ Alpha operator represents accumulation term: """
        values_np[self.ACC_OP] = z[0] * rho_t
        values_np[self.ACC_OP + 1:self.ACC_OP + nc] = (1 - ss) * z[1:] * rho_f

        """ Beta operator represents flux term: """
        for j in range(nph):
            values_np[self.FLUX_OP + j * self.ne:self.FLUX_OP + j * self.ne + self.nc] = self.x[j] * self.rho_m[j] * self.kr[j] / self.mu[j]

        """ Gamma operator for diffusion (same for thermal and isothermal) """
        shift = ne + ne * nph
        for j in range(nph):
            values_np[self.UPSAT_OP + j] = self.compr * self.sat[j]

        """ Chi operator for diffusion """
        dif_coef = np.array([0, 1, 1, 1, 1]) * 5.2e-10 * 86400
        for i in range(nc):
            for j in range(nph):
                values_np[self.GRAD_OP + i * nph + j] = dif_coef[i] * self.rho_m[j] * self.x[j][i]
                # values[shift + ne * j + i] = 0

        """ Delta operator for reaction """
        for i in range(ne):
            values_np[self.KIN_OP + i] = self.input_data.stoich_matrix[i] * kin_rate

        """ Gravity and Capillarity operators """
        # E3-> gravity
        for i in range(nph):
            values_np[self.GRAV_OP + i] = 0

        # E4-> capillarity
        for i in range(nph):
            values_np[self.PC_OP + i] = 0

        # E5_> porosity
        values_np[self.PORO_OP] = 1 - ss

        # values[shift + 3 + 2 * nph + 1] = kin_state['SR']
        # values[shift + 3 + 2 * nph + 2] = kin_state['Act(H+)']

        return 0

# Operators required for initialization, to convert given volume fraction to molar one
# state: (pressure, overall mineral volume fractions, fluid molar fractions)
class my_own_comp_etor(my_own_acc_flux_etor):
    def __init__(self, input_data, properties):
        super().__init__(input_data, properties)  # Initialize base-class
        self.fluid_mole = 1
        self.counter = 0
        self.props_name = ['z_solid']

    def evaluate(self, state, values):
        state_np = state.to_numpy()
        values_np = values.to_numpy()
        pressure = state_np[0]
        ss = state_np[1] # volume fraction in initialization

        # initial flash
        _, _, _, _, _, fluid_volume, _ = self.property.flash_ev.evaluate(state_np)

        # evaluate molar fraction
        solid_volume = fluid_volume * ss / (1 - ss)         # m3
        solid_mole = solid_volume * self.property.density_ev['solid'].evaluate(pressure) / self.property.Mw['Solid']
        nu_s = solid_mole / (solid_mole + self.fluid_mole)
        values_np[0] = nu_s

        return 0

class my_own_rate_evaluator(operator_set_evaluator_iface):
    # Simplest class existing to mankind:
    def __init__(self, properties, temperature):
        # Initialize base-class
        super().__init__()
        self.property = properties
        self.temperature = temperature

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        # Composition vector and pressure from state:
        state_np = state.to_numpy()
        values_np = values.to_numpy()
        pressure = state_np[0]

        # zc = np.append(state_np[2:], 1 - np.sum(state_np[1:]))
        # Perform Flash procedure here:
        vap, x, y, rho_phases, _, _, _ = self.property.flash_ev.evaluate(state_np)

        # Note: officially three phases are present now
        rho_w = rho_phases['aq']
        mu_w = CP.PropsSI('V', 'T', self.temperature, 'P|liquid', bar2pa(pressure), 'Water') * 1000        # Pa * s

        rho_g = rho_phases['gas']

        try:
            mu_g = CP.PropsSI('V', 'T', self.temperature, 'P|gas', bar2pa(pressure), 'CarbonDioxide') * 1000       # Pa * s
        except ValueError:
            mu_g = 16.14e-6 * 1000     # Pa * s, for 50 C

        # Easiest example, constant volumetric phase rate:
        values[0] = 0   # vapor phase
        values[1] = 1 / mu_w    # liquid phase

        return 0

class my_own_property_evaluator(operator_set_evaluator_iface):
    def __init__(self, input_data, properties):
        # Initialize base-class
        super().__init__()
        self.input_data = input_data
        self.property = properties
        self.props_name = ['z' + prop for prop in properties.flash_ev.phreeqc_species]

    def evaluate(self, state, values):
        state_np = state.to_numpy()
        values_np = values.to_numpy()
        _, _, _, _, _, _, molar_fractions = self.property.flash_ev.evaluate(state_np)
        values_np[:molar_fractions.size] = molar_fractions

        return 0
