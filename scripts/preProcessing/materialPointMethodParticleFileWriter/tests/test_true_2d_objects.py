"""Native shape, legacy extrusion and physical Weibull-volume contracts.

The single-rank MPI adapter is for constructor imports only; no MPI claim.
"""
import contextlib
import io
from pathlib import Path
import sys
import types
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
if 'mpi4py' not in sys.modules:
    try:
        import mpi4py
    except ImportError:
        comm=types.SimpleNamespace(Get_rank=lambda:0, Get_size=lambda:1, bcast=lambda v,root=0:v,
                                   allreduce=lambda v,op=None:v, Barrier=lambda:None)
        sys.modules['mpi4py']=types.SimpleNamespace(MPI=types.SimpleNamespace(COMM_WORLD=comm))
import pfw_geometryObjects as geom

class NativeObjectTests(unittest.TestCase):
    def same_queries(self, shape, xy, skin=.08):
        reference=None
        for z in (-1e6, -.001, 0, .001, 1e6):
            point=[*xy,z]
            result=np.r_[shape.isInterior(point,skin),shape.getSurfaceNormal(point),shape.getSurfacePosition(point)]
            self.assertTrue(np.all(np.isfinite(result)))
            self.assertEqual(result[3],0);self.assertEqual(result[6],0)
            if reference is None: reference=result
            else: np.testing.assert_array_equal(result,reference)
        return reference

    def test_disk_annulus_and_flags(self):
        disk=geom.disk('disk',[1.,2.],1.)
        self.assertEqual(disk.isInterior([1.,2.,0],0),0)
        result=self.same_queries(disk,[1.95,2.])
        self.assertGreater(result[0],0)
        np.testing.assert_allclose(result[1:4],[1,0,0])
        np.testing.assert_allclose(result[4:],[.05,0,0])
        annulus=geom.disk('annulus',[0,0],2,ri=1)
        result=self.same_queries(annulus,[1.02,0])
        np.testing.assert_allclose(result[1:4],[-1,0,0])
        np.testing.assert_allclose(result[4:],[-.02,0,0])
        outer_only=geom.disk('outer',[0,0],2,ri=1,flaggedSurfaces=[True,False])
        self.assertEqual(outer_only.isInterior([1.02,0,0],.1),0)
        np.testing.assert_allclose(outer_only.getSurfaceNormal([1.02,0,0]),[1,0,0])
        no_flags=geom.disk('none',[0,0],1,flaggedSurfaces=[False,False])
        np.testing.assert_array_equal(no_flags.getSurfaceNormal([.99,0,0]),np.zeros(3))
        for bad in (dict(r=-1),dict(r=1,ri=2),dict(r=float('nan'))):
            with self.assertRaises(ValueError):geom.disk('bad',[0,0],**bad)

    def test_rectangle_polygon_foam_and_bed_are_z_independent(self):
        shapes=[geom.rectangle('r',[0,0],[1,1]),
                geom.polygon('p',[[0,0],[1,0],[1,1],[0,1]])]
        with contextlib.redirect_stdout(io.StringIO()):
            shapes.append(geom.poissonDiskFoam('foam',[0,0],[1,1],poreDiameter=.1,porosity=.05,
                                              periodic=[False,False],dim=2,seed=17))
            bed=geom.packedSphericalBed('bed',[0,0],[1,1],materials=[dict(mat=0,volumeFraction=.08,meanDiameter=.15,stdDiameter=0)],
                                        dim=2,periodic=[False,False],method='rsa',seed=17)
        shapes.append(bed)
        for shape in shapes:
            for xy in ([.02,.5],[.3,.2],[.5,.5]):
                self.same_queries(shape,xy)

    def test_legacy_extrusions_ignore_end_caps_and_preserve_metadata(self):
        results=[]
        for depth in (1e-8,.01,100):
            cylinder=geom.cylinder('old',[0,0,-depth],[0,0,depth],1.,mat=7,group=3)
            box=geom.box('oldbox',[0,0,-depth],[1,1,depth],flaggedSurfaces=[True,False,True,False,True,True])
            wrapped=geom.union('union',cylinder,box)
            geom.prepare_plane_strain_objects([wrapped,cylinder])
            self.assertEqual((cylinder.mat,cylinder.group),(7,3))
            self.assertEqual(box.dim,2)
            np.testing.assert_array_equal(box.flaggedSurfaces,[True,False,False,True])
            results.append(self.same_queries(cylinder,[.95,0]))
            self.same_queries(box,[.02,.5])
        for result in results[1:]:np.testing.assert_array_equal(result,results[0])
        tilted=geom.cylinder('tilted',[0,0,0],[1,0,1],1.)
        geom.prepare_plane_strain_objects([tilted])
        self.assertNotIn('_planeStrainSection',vars(tilted))

class WeibullThicknessTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(174)
        self.sites=rng.uniform(.05,.95,(12,2))
        self.shape=geom.rectangle('solid',[0.,0.],[1.,1.])

    def check_binding(self,obj):
        areas,draws,reference,modulus,floor,target=obj._nativeWeibull
        for thickness in (.001,1.,3.,.001):
            geom.bind_plane_strain_thickness([obj,obj],thickness)
            expected=lambda a,r:((reference/np.maximum(np.asarray(a)*thickness,floor))*np.log(r)/np.log(.5))**(1/modulus)
            actual=getattr(obj,target)
            if target=='cellSliceStrengthScale':
                for a,r,value in zip(areas,draws,actual):np.testing.assert_allclose(value,expected(a,r),rtol=1e-14)
            else:np.testing.assert_allclose(actual,expected(areas,draws),rtol=1e-14)
        self.assertEqual(obj.domainThickness,.001)

    def test_all_builtin_2d_weibull_wrappers(self):
        with contextlib.redirect_stdout(io.StringIO()):
            wrappers=[geom.voronoiWeibullBoxWrapper('v',self.shape,[0,0],[1,1],.2,2.,5.,8,.01,vpts=self.sites,dim=2),
                geom.layeredStrengthScaleBoxWrapper('l',self.shape,[0,0,0],[1,1,0],matDir=[1.,0.,0.],
                    weibullReferenceVolume=2.,weibullModulus=5.,weibullSeed=8,weibullLayerThickness=.2,
                    weibullMinVolume=.01,nIntegrationPoints=4,dim=2),
                geom.voronoiWrapperWithLayeredStrengthScale('vl',self.shape,[0,0,0],[1,1,0],.2,2.,.01,2.,5.,8,.2,
                    nIntegrationPoints=4,vpts=self.sites,dim=2),
                geom.layeredVoronoiWeibullWrapper('lv',self.shape,[0,0,0],[1,1,0],.2,2.,5.,8,.01,.2,
                    nIntegrationPoints=4,matDirs=np.eye(3),vpts=self.sites,dim=2)]
        for wrapper in wrappers:self.check_binding(wrapper)

    def test_nested_custom_hook_and_conflicting_thickness(self):
        calls=[]
        class Custom:
            def setPlaneStrainThickness(self,t):calls.append(t)
        shared=Custom();geom.bind_plane_strain_thickness([shared,shared],.3)
        self.assertEqual(calls,[.3])
        with contextlib.redirect_stdout(io.StringIO()):
            wrapper=geom.layeredStrengthScaleBoxWrapper('l',self.shape,[0,0,0],[1,1,0],matDir=[1.,0.,0.],
                    weibullLayerThickness=.2,weibullMinVolume=.01,dim=2,domainThickness=.3)
        with self.assertRaisesRegex(ValueError,'disagrees'):geom.bind_plane_strain_thickness([wrapper],.4)
        for bad in (0,-1,float('nan')):
            with self.assertRaises(ValueError):geom.bind_plane_strain_thickness([],bad)

if __name__=='__main__':unittest.main()
