from darts.models.cicd_model import CICDModel
from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer
from darts.physics.properties.flash import ConstantK
from darts.physics.properties.basic import ConstFunc, PhaseRelPerm
from darts.physics.properties.density import DensityBasic

class Model(CICDModel):
    def __init__(self, n_points=50):
        # Call base class constructor
        super().__init__()
        self.n_obl_points = n_points

        self.set_physics()

        self.initial_values = {self.physics.vars[0]: 50,
                               self.physics.vars[1]: 0.1,
                               self.physics.vars[2]: 0.2
                               }

    def set_physics(self):
        """Physical properties"""
        self.zero = 1e-8
        # Create property containers:
        components = ['CO2', 'C1', 'H2O']
        phases = ['gas', 'oil']
        thermal = 0
        Mw = [44.01, 16.04, 18.015]

        property_container = PropertyContainer(phases_name=phases, components_name=components,
                                               Mw=Mw, min_z=self.zero / 10, temperature=1.)

        """ properties correlations """
        property_container.flash_ev = ConstantK(len(components), [4, 2, 1e-1], self.zero)
        property_container.density_ev = dict([('gas', DensityBasic(compr=1e-3, dens0=200)),
                                              ('oil', DensityBasic(compr=1e-5, dens0=600))])
        property_container.viscosity_ev = dict([('gas', ConstFunc(0.05)),
                                                ('oil', ConstFunc(0.5))])
        property_container.rel_perm_ev = dict([('gas', PhaseRelPerm("gas")),
                                               ('oil', PhaseRelPerm("oil"))])

        """ Activate physics """
        self.physics = Compositional(components, phases, self.timer,
                                     n_points=self.n_obl_points, min_p=1, max_p=300, min_z=self.zero/10, max_z=1-self.zero/10)
        self.physics.add_property_region(property_container)
        self.engine = self.physics.init_physics(platform='cpu')
        return
