"""Run the actual renderer selection/camera code against a VisIt API adapter."""
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pfw_visit_render as render

class Visit2DTests(unittest.TestCase):
    def test_native_metadata_does_not_use_tensor_components(self):
        for dimensions,expected in (([2,2],True),([3,3],False),([2,3],False),([],False)):
            md=types.SimpleNamespace(GetNumMeshes=lambda:len(dimensions),
                                     GetMeshes=lambda i:types.SimpleNamespace(spatialDimension=dimensions[i]))
            with patch.object(render,'call_visit',return_value=md):
                self.assertEqual(render.native_xy_database('case'),expected)

    def test_camera_fit_uses_view2d_and_keeps_aspect(self):
        view=types.SimpleNamespace();calls=[]
        def api(name,*args):
            calls.append(name)
            if name=='GetView2D':return view
        with patch.object(render,'call_visit',side_effect=api):
            render.configure_fitted_view('xy',(-1,1,-.5,.5),1200,600,1.1,False,native_2d=True)
        x0,x1,y0,y1=view.windowCoords
        self.assertAlmostEqual((x1-x0)/(y1-y0),2)
        self.assertLess(x0,-1);self.assertGreater(x1,1)
        self.assertEqual(calls,['GetView2D','SetView2D'])
        with self.assertRaises(ValueError):
            render.configure_fitted_view('xz',(-1,1,-1,1,0,0),600,600,1.1,False,native_2d=True)

    def test_invalid_extents_reset_view(self):
        with patch.object(render,'call_visit') as api:
            render.configure_fitted_view2d((0,float('nan'),0,1),600,600,1.1)
            api.assert_called_once_with('ResetView')

if __name__=='__main__':unittest.main()
