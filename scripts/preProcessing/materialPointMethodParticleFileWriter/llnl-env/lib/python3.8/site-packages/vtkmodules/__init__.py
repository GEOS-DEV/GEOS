r"""
Currently, this package is experimental and may change in the future.
"""
from __future__ import absolute_import
import sys
import importlib.util

from pathlib import Path

from importlib.abc import MetaPathFinder
from importlib.machinery import ExtensionFileLoader

LOADING_STACK = []

def find_lib_path(paths, vtk_module_name):
    for ext in [".pyd", ".so"]:
        for base_path in paths:
            for f in Path(base_path).glob(f"{vtk_module_name}.*{ext}"):
                return str(f.resolve())


class VTKMetaHook(MetaPathFinder):
    """Attach a custom  loaded for vtk native library loading to defer loading of pure python dependencies"""
    def find_spec(self, fullname, path, target=None):
        if fullname.startswith("vtkmodules.vtk"):
            vtk_module_name = fullname.split(".")[1]
            module_path = find_lib_path(path, vtk_module_name)
            if module_path is None:
                return None

            LOADING_STACK.append(fullname)
            return importlib.util.spec_from_file_location(fullname, module_path, loader=VTKLoader(fullname, module_path))

        return None


class VTKLoader(ExtensionFileLoader):
    """Flush any pending dependency load once initialize() phase is done"""
    def exec_module(self, module):
        super().exec_module(module)

        # Process pending dependencies only if the module match the first load request
        if len(LOADING_STACK) and LOADING_STACK[0] == module.__name__:
            LOADING_STACK.clear()
            on_vtk_module_init_completed()



# Register our hook for vtk library loader
sys.meta_path.insert(0, VTKMetaHook())


def _windows_dll_path():
    import os
    _vtk_python_path = './vtkmodules'
    _vtk_dll_path = 'bin'
    # Compute the DLL path based on the location of the file and traversing up
    # the installation prefix to append the DLL path.
    _vtk_dll_directory = os.path.dirname(os.path.abspath(__file__))
    # Loop while we have components to remove.
    while _vtk_python_path not in ('', '.', '/'):
        # Strip a directory away.
        _vtk_python_path = os.path.dirname(_vtk_python_path)
        _vtk_dll_directory = os.path.dirname(_vtk_dll_directory)
    _vtk_dll_directory = os.path.join(_vtk_dll_directory, _vtk_dll_path)
    if os.path.exists(_vtk_dll_directory):
        # We never remove this path; it is required for VTK to work and there's
        # no scope where we can easily remove the directory again.
        _ = os.add_dll_directory(_vtk_dll_directory)

    # Build tree support.
    try:
        from . import _build_paths

        # Add any paths needed for the build tree.
        for path in _build_paths.paths:
            if os.path.exists(path):
                _ = os.add_dll_directory(path)
    except ImportError:
        # Relocatable install tree (or non-Windows).
        pass


# CPython 3.8 added behaviors which modified the DLL search path on Windows to
# only search "blessed" paths. When importing SMTK, ensure that SMTK's DLLs are
# in this set of "blessed" paths.
if sys.version_info >= (3, 8) and sys.platform == 'win32':
    _windows_dll_path()


#------------------------------------------------------------------------------
# this little trick is for static builds of VTK. In such builds, if
# the user imports this Python package in a non-statically linked Python
# interpreter i.e. not of the of the VTK-python executables, then we import the
# static components importer module.
def _load_vtkmodules_static():
    if 'vtkmodules_vtkCommonCore' not in sys.builtin_module_names:
        import _vtkmodules_static

#_load_vtkmodules_static()


#------------------------------------------------------------------------------
# list the contents
__all__ = [
        'vtkCommonCore',
    'vtkWebCore',
    'vtkCommonMath',
    'vtkCommonTransforms',
    'vtkCommonDataModel',
    'vtkCommonExecutionModel',
    'vtkIOCore',
    'vtkImagingCore',
    'vtkIOImage',
    'vtkIOXMLParser',
    'vtkIOXML',
    'vtkCommonMisc',
    'vtkFiltersCore',
    'vtkRenderingCore',
    'vtkRenderingContext2D',
    'vtkRenderingFreeType',
    'vtkRenderingSceneGraph',
    'vtkRenderingVtkJS',
    'vtkIOExport',
    'vtkWebGLExporter',
    'vtkInteractionStyle',
    'vtkFiltersGeneral',
    'vtkFiltersSources',
    'vtkInteractionWidgets',
    'vtkViewsCore',
    'vtkViewsInfovis',
    'vtkCommonComputationalGeometry',
    'vtkCommonSystem',
    'vtkFiltersCellGrid',
    'vtkIOCellGrid',
    'vtkIOLegacy',
    'vtkDomainsChemistry',
    'vtkRenderingHyperTreeGrid',
    'vtkRenderingUI',
    'vtkRenderingOpenGL2',
    'vtkRenderingContextOpenGL2',
    'vtkRenderingVolume',
    'vtkImagingMath',
    'vtkRenderingVolumeOpenGL2',
    'vtkViewsContext2D',
    'vtkImagingColor',
    'vtkTestingRendering',
    'vtkSerializationManager',
    'vtkRenderingVolumeAMR',
    'vtkPythonContext2D',
    'vtkParallelCore',
    'vtkRenderingParallel',
    'vtkRenderingVRModels',
    'vtkRenderingVR',
    'vtkRenderingMatplotlib',
    'vtkRenderingLabel',
    'vtkRenderingLOD',
    'vtkRenderingLICOpenGL2',
    'vtkRenderingImage',
    'vtkRenderingExternal',
    'vtkRenderingCellGrid',
    'vtkIOXdmf2',
    'vtkIOVeraOut',
    'vtkIOVPIC',
    'vtkIOTecplotTable',
    'vtkIOTRUCHAS',
    'vtkIOSegY',
    'vtkIOParallelXML',
    'vtkIOLSDyna',
    'vtkIOParallelLSDyna',
    'vtkIOExodus',
    'vtkIOParallelExodus',
    'vtkIOPLY',
    'vtkIOPIO',
    'vtkIOMovie',
    'vtkIOOggTheora',
    'vtkIOOMF',
    'vtkIONetCDF',
    'vtkIOMotionFX',
    'vtkIOGeometry',
    'vtkIOParallel',
    'vtkIOMINC',
    'vtkIOImport',
    'vtkIOIOSS',
    'vtkIOHDF',
    'vtkIOH5part',
    'vtkIOH5Rage',
    'vtkIOGeoJSON',
    'vtkIOFLUENTCFF',
    'vtkIOVideo',
    'vtkIOFDS',
    'vtkIOInfovis',
    'vtkIOExportPDF',
    'vtkRenderingGL2PSOpenGL2',
    'vtkIOExportGL2PS',
    'vtkIOEngys',
    'vtkIOEnSight',
    'vtkIOERF',
    'vtkIOCityGML',
    'vtkIOChemistry',
    'vtkIOCesium3DTiles',
    'vtkIOCONVERGECFD',
    'vtkIOCGNSReader',
    'vtkIOAsynchronous',
    'vtkIOAMR',
    'vtkInteractionImage',
    'vtkInfovisLayout',
    'vtkImagingStencil',
    'vtkImagingStatistics',
    'vtkImagingGeneral',
    'vtkImagingOpenGL2',
    'vtkImagingMorphological',
    'vtkImagingFourier',
    'vtkIOSQL',
    'vtkRenderingAnnotation',
    'vtkImagingHybrid',
    'vtkGeovisCore',
    'vtkFiltersTopology',
    'vtkFiltersTensor',
    'vtkFiltersSelection',
    'vtkFiltersSMP',
    'vtkFiltersPython',
    'vtkFiltersProgrammable',
    'vtkFiltersModeling',
    'vtkFiltersPoints',
    'vtkFiltersStatistics',
    'vtkFiltersParallelStatistics',
    'vtkFiltersImaging',
    'vtkFiltersExtraction',
    'vtkFiltersGeometry',
    'vtkFiltersHybrid',
    'vtkFiltersHyperTree',
    'vtkFiltersTexture',
    'vtkFiltersParallel',
    'vtkFiltersParallelImaging',
    'vtkFiltersParallelDIY2',
    'vtkFiltersTemporal',
    'vtkFiltersGeometryPreview',
    'vtkFiltersGeneric',
    'vtkFiltersFlowPaths',
    'vtkFiltersAMR',
    'vtkDomainsChemistryOpenGL2',
    'vtkCommonPython',
    'vtkChartsCore',
    'vtkCommonColor',
    'vtkImagingSources',
    'vtkInfovisCore',
    'vtkAcceleratorsVTKmCore',
    'vtkAcceleratorsVTKmDataModel',
    'vtkAcceleratorsVTKmFilters',
    'vtkFiltersVerdict',
    'vtkFiltersReduction',
    'all',
    'gtk',
    'numpy_interface',
    'qt',
    'test',
    'tk',
    'util',
    'wx',

]
#------------------------------------------------------------------------------
# get the version
__version__ = "9.4.1"

#------------------------------------------------------------------------------
# describe import dependencies to properly define Python @override
MODULE_MAPPER = {
    "vtkCommonDataModel": [
        "vtkmodules.util.data_model",
    ],
    "vtkCommonExecutionModel": [
        "vtkmodules.util.execution_model",
    ],
}
LOADED_MODULES = set()
PENDING_LOADED_MODULES = set()

def register_vtk_module_dependencies(vtk_module_name, *import_names):
    """Method to call for registering external override on vtkmodule load"""
    MODULE_MAPPER.setdefault(vtk_module_name, []).extend(import_names)

    # If already loaded let's make sure we import it now
    if vtk_module_name in LOADED_MODULES:
        for import_name in import_names:
            importlib.import_module(import_name)


def on_vtk_module_init(module_name):
    """Automatically called by vtkmodule when they are loaded"""
    if module_name in LOADED_MODULES:
        return

    PENDING_LOADED_MODULES.add(module_name)


def on_vtk_module_init_completed():
    pending = list(PENDING_LOADED_MODULES)
    PENDING_LOADED_MODULES.clear()

    for module_name in pending:
        LOADED_MODULES.add(module_name)
        for import_name in MODULE_MAPPER.get(module_name, []):
            importlib.import_module(import_name)
