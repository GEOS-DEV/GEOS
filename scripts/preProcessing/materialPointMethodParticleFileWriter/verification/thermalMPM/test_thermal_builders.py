#!/usr/bin/env python3
"""Pure-Python structural tests for the PFW case construction, with GEOS stubs."""
from __future__ import annotations
import sys
import types
import unittest
import xml.etree.ElementTree as ET

from thermal_cases import CASES, RC_GLOBAL, RC_GROUP, RC_MATERIAL, RC_COMBINED

class FakeBox:
    def __init__(self, name, x0, x1, *, mat, group, **kwargs):
        self.name, self.x0, self.x1 = name, x0, x1
        self.mat, self.group = mat, group


def fake_set_tracers(pfw, points, variables, **kwargs):
    pfw['_tracers'] = (points, variables, kwargs)


class TestCaseInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        module=types.ModuleType('pfw_geometryObjects')
        module.box=FakeBox
        module2=types.ModuleType('pfw_tracerPoints')
        module2.set_tracers=fake_set_tracers
        cls.previous={k:sys.modules.get(k) for k in ('pfw_geometryObjects','pfw_tracerPoints','thermal_pfw')}
        sys.modules['pfw_geometryObjects']=module
        sys.modules['pfw_tracerPoints']=module2
        sys.modules.pop('thermal_pfw',None)
        from thermal_pfw import build
        cls.build=staticmethod(build)

    @classmethod
    def tearDownClass(cls):
        for key,original in cls.previous.items():
            if original is None:sys.modules.pop(key,None)
            else:sys.modules[key]=original

    def test_all_case_inputs(self):
        for name in CASES:
            with self.subTest(case=name):
                pfw=self.build(name)
                self.assertEqual(pfw['thermal']['enabled'],1)
                self.assertTrue(pfw['thermal']['diagnostics'])
                self.assertEqual(len(pfw['materials']),max(o.mat for o in pfw['objects'])+1)
                root=ET.fromstring('<Constitutive>'+pfw['materialPropertyString']+'</Constitutive>')
                self.assertEqual(len(root),len(pfw['materials']))
                self.assertEqual([a.attrib['name'] for a in root],pfw['materials'])
                self.assertTrue(pfw['_tracers'][0])
                self.assertIn('temperature',pfw['_tracers'][1])
                for key in ('contactResistance','heatSources','isothermal'):
                    self.assertIsInstance(pfw['thermal'][key],list)

    def test_contact_precedence(self):
        for case,expected in (
             ('contact_global',(RC_GLOBAL,[])),
             ('contact_material',(RC_GLOBAL,[RC_MATERIAL])),
             ('contact_group',(RC_GLOBAL,[RC_MATERIAL,RC_GROUP])),
             ('contact_combined',(RC_GLOBAL,[RC_MATERIAL,RC_GROUP,RC_COMBINED]))):
            with self.subTest(case=case):
                vals=self.build(case)['thermal']
                self.assertEqual(vals['defaultContactResistance'],expected[0])
                self.assertEqual([r['resistance'] for r in vals['contactResistance']],expected[1])

    def test_particle_selection(self):
        heat=self.build('source_local')['thermal']['heatSources']
        self.assertEqual(len(heat),1)
        self.assertEqual(heat[0]['targetRegions'],['ParticleRegion2'])
        bath=self.build('dirichlet_transient')['thermal']['isothermal']
        self.assertEqual([r['targetRegions'] for r in bath],
                         [['ParticleRegion1'],['ParticleRegion3']])

    def test_nonlinear_thermal_models(self):
        for name,attr in [('capacity_linear','thermalHeatCapacitySlope'),
                          ('conductivity_linear','thermalConductivitySlope')]:
            root=ET.fromstring(self.build(name)['materialPropertyString']) if name=='capacity_linear' else None
            if root is None:
                root=ET.fromstring('<M>'+self.build(name)['materialPropertyString']+'</M>')[0]
            self.assertEqual(root.tag,'LinearThermalElasticIsotropic')
            self.assertIn(attr,root.attrib)

    def test_dfg_inputs(self):
        pfw=self.build('contact_dfg')
        self.assertEqual(pfw['damageFieldPartitioning'],1)
        self.assertEqual(pfw['needsNeighborList'],1)
        self.assertNotEqual(pfw['objects'][0].damage,pfw['objects'][1].damage)


if __name__=='__main__':unittest.main()
