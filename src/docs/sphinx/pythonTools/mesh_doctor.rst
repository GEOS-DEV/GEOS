``mesh_doctor``
---------------

``mesh_doctor`` is a ``python`` executable that can be used through the command line to perform various checks, validations, and tiny fixes to the ``vtk`` mesh that are meant to be used in ``geos``.
``mesh_doctor`` is organized as a collection of modules with their dedicated sets of options.
The current page will introduce those modules, but the details and all the arguments can be retrieved by using the ``--help`` option for each module.

Modules
^^^^^^^

To list all the modules available through ``mesh_doctor``, you can simply use the ``--help`` option, which will list all available modules as well as a quick summary.

.. command-output:: python mesh_doctor.py --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

Then, if you are interested in a specific module, you can ask for its documentation using the ``mesh_doctor module_name --help`` pattern.
For example

.. command-output:: python mesh_doctor.py collocated_nodes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``mesh_doctor`` loads its module dynamically.
If a module can't be loaded, ``mesh_doctor`` will proceed and try to load other modules.
If you see a message like

.. code-block:: bash

    [1970-04-14 03:07:15,625][WARNING] Could not load module "collocated_nodes": No module named 'vtkmodules'

then most likely ``mesh_doctor`` could not load the ``collocated_nodes`` module, because the ``vtk`` python package was not found.
Thereafter, the documentation for module ``collocated_nodes`` will not be displayed.
You can solve this issue by installing the dependencies of ``mesh_doctor`` defined in its ``requirements.txt`` file (``python -m pip install -r requirements.txt``).

Here is a list and brief description of all the modules available.

``collocated_nodes``
""""""""""""""""""""

Displays the neighboring nodes that are closer to each other than a prescribed threshold.
It is not uncommon to define multiple nodes for the exact same position, which will typically be an issue for ``geos`` and should be fixed.

.. command-output:: python mesh_doctor.py collocated_nodes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``element_volumes``
"""""""""""""""""""

Computes the volumes of all the cells and displays the ones that are below a prescribed threshold.
Cells with negative volumes will typically be an issue for ``geos`` and should be fixed.

.. command-output:: python mesh_doctor.py element_volumes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``fix_elements_orderings``
""""""""""""""""""""""""""

It sometimes happens that an exported mesh does not abide by the ``vtk`` orderings.
The ``fix_elements_orderings`` module can rearrange the nodes of given types of elements.
This can be convenient if you cannot regenerate the mesh.

.. command-output:: python mesh_doctor.py fix_elements_orderings --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``generate_cube``
"""""""""""""""""

This module conveniently generates cubic meshes in ``vtk``.
It can also generate fields with simple values.
This tool can also be useful to generate a trial mesh that will later be refined or customized.

.. command-output:: python mesh_doctor.py generate_cube --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``generate_fractures``
""""""""""""""""""""""

For a conformal fracture to be defined in a mesh, ``geos`` requires the mesh to be split at the faces where the fracture gets across the mesh.
The ``generate_fractures`` module will split the mesh and generate the multi-block ``vtk`` files.

.. command-output:: python mesh_doctor.py generate_fractures --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``generate_global_ids``
"""""""""""""""""""""""

When running ``geos`` in parallel, `global ids` can be used to refer to data across multiple ranks.
The ``generate_global_ids`` can generate `global ids` for the imported ``vtk`` mesh.

.. command-output:: python mesh_doctor.py generate_global_ids --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``non_conformal``
"""""""""""""""""

This module will detect elements which are close enough (there's a user defined threshold) but which are not in front of each other (another threshold can be defined).
`Close enough` can be defined in terms or proximity of the nodes and faces of the elements.
The angle between two faces can also be precribed.
This module can be a bit time consuming.

.. command-output:: python mesh_doctor.py non_conformal --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``self_intersecting_elements``
""""""""""""""""""""""""""""""

Some meshes can have cells that auto-intersect.
This module will display the elements that have faces intersecting.

.. command-output:: python mesh_doctor.py self_intersecting_elements --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``supported_elements``
""""""""""""""""""""""

``geos`` supports a specific set of elements.
Let's cite the standard elements like `tetrahedra`, `wedges`, `pyramids` or `hexahedra`.
But also prismes up to 11 faces.
``geos`` also supports the generic ``VTK_POLYHEDRON``/``42`` elements, which are converted on the fly into one of the elements just described.

The ``supported_elements`` check will validate that no unsupported element is included in the input mesh.
It will also verify that the ``VTK_POLYHEDRON`` cells can effectively get converted into a supported type of element.

.. command-output:: python mesh_doctor.py supported_elements --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``Using mesh_doctor in paraview``
""""""""""""""""""""""""""""""""""

Using mesh_doctor as a programmable filter
____________________________________________

To use ``mesh_doctor`` in Paraview as a python programmable filter, a python package install is required first in Paraview python resolved
path. Paraview is storing its python ressources under its *lib/pythonX.X* depending on the paraview version, *e.g* Paraview 5.11 is working
with python 3.9. As a results the following command will install ``mesh_doctor`` package into Paraview resolved path.

.. command-output:: python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps --upgrade --target /path/to/Paraview/lib/python3.9/ mesh_doctor

.. note::
    ``pip`` is installing the ``mesh_doctor`` package from the test.pypi repo, which is intended to test package deployment.
    Once stabilized and ``mesh_doctor`` uploaded onto the main package repo, this should be dropped out.

Once the installation done, the */path/to/Paraview/lib/pythonX.X* should holds ``mesh_doctor`` package content, *i.e.* ``checks`` and ``parsing``.
Then launching ``Paraview`` and loading our *mesh.vtu*, as an example, we will design a *Programmable python filter* relying on *element_volumes* from
``mesh_doctor``. Add such a filter pipelined after the mesh reader, in the script section paste the following,

.. code-block:: python
    :linenos:

    mesh = inputs[0].VTKObject
    tol = 1.2e-6

    from checks import element_volumes
    import vtk

    res = element_volumes.__check(mesh, element_volumes.Options(tol))
    #print(res)
    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)
    for cell_index, volume in res.element_volumes:
        ids.InsertNextValue(cell_index)

    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
    selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
    selectionNode.SetSelectionList(ids)
    selection = vtk.vtkSelection()
    selection.AddNode(selectionNode)
    extracted = vtk.vtkExtractSelection()
    extracted.SetInputDataObject(0, mesh)
    extracted.SetInputData(1, selection)
    extracted.Update()
    print("There are {} cells under {} m3 vol".format(extracted.GetOutput().GetNumberOfCells(), tol))
    output.ShallowCopy(extracted.GetOutput())

Here we rely on ``pyvtk`` interface more than on Paraview adaptation, for legacy and reusability reasons. This is the reason
for the full ``import vtk`` instead of ``from paraview import vtk``, the `vtkSelectionNode` being fully wrapped in paraview
and not accessible otherwise.

On line 7, we leverage ``mesh_doctor`` package to provide us with pairs of `(index,volumes)` of cells with volumes lower
than tolerance `tol`. As input of *Programmable Python Filter* is wrapped in a `dataset_adapter.UnstructuredGrid`, we rely on
the copy of the inital VTKObject `inputs[0].VTKObject` to ensure consistency with our ``pyvtk`` workflow.

What follows is ``pyvtk`` steps in oder to convert into input struct and extract from the original mesh this list of cells.
Eventually, the `extracted` selection is shallow-copied to the output and then accessible in ``Paraview``. An helper print
is left and should be reported in *Output Message* of ``Paraview`` (and in launching terminal if exist).

Using mesh_doctor as a paraview plugins
____________________________________________

Another way of leveraging ``mesh_doctor`` in ``Paraview`` is to wrap it in a python plugin that would be loadable through the
``Paraview`` interface under **Tools | Manage Plugins/Extensions** and **Load New** looking for ``mesh_doctor-pvplugin.py``.
(see `Paraview How To <https://www.paraview.org/Wiki/ParaView/Plugin_HowTo#Using_Plugins>`_  for more details).

The file ``mesh_doctor-pvplugin.py`` is located under the ``geosx_mesh_doctor`` module in GEOS. Once the plugin loaded and a mesh opened,
it should appear in filter list as *Mesh Doctor(GEOS)*. It displays a parameter value box allowing the user to enter the volume he wants as
threshold to select cells based on ``element_volumes`` capability. Once applied, it extracts selected set of cells as a new unstructured grid.