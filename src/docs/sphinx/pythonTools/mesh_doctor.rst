``mesh_doctor``
---------------

``mesh_doctor`` is a ``python`` executable that can be used through the command line to perform various checks, validations and tiny fixes to the ``vtk`` mesh that are meant to be used in ``geos``.
``mesh_doctor`` is organized as a collection of modules with their dedicated set of options.
The current page will introduce those modules, but the details and all the arguments can be retrieved by using the ``--help`` option for each module.

Modules
^^^^^^^

To list all the modules available through ``mesh_doctor`` you can simply use the ``--help`` option which will provide all the modules as well as a quick summuray.

.. command-output:: python mesh_doctor.py --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

.. program-output:: pwd

Then, it you are interested in a specific module, you can ask for its documentation using the ``mesh_doctor module_name --help`` pattern.
For example

.. command-output:: python mesh_doctor.py collocated_nodes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``mesh_doctor`` loads its module dynamically.
If a module can't be loaded, ``mesh_doctor`` will proceed and try to load the other modules.
If you see a message like

.. code-block:: bash

    [1969-07-21 02:56:15,625][WARNING] Could not load module "collocated_nodes": No module named 'vtkmodules'

then most likely ``mesh_doctor`` could not load the ``collocated_nodes`` modules because the ``vtk`` modules was not found.
Consistently, the documentation for module ``collocated_nodes`` will not be displayed.
You can solve this issue by installing the dependencies of ``mesh_doctor`` defined in its ``requirements.txt`` file (``python -m pip install -r requirements.txt``).

Here is now a summury of all the modules available.

``collocated_nodes``
""""""""""""""""""""

Displays the nodes that are closer than a prescribed threshold.
It is not uncommon to find multiple nodes in the exact same position.
This will typically be an issue for ``geos`` and should be fixed.

.. program-output:: python mesh_doctor.py collocated_nodes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``element_volumes``
"""""""""""""""""""

Computes the volumes of all the cells and displays the ones that are above a prescribed threshold.
Cells with negative volumes will typically be an issue for ``geos`` and should be fixed.

.. command-output:: python mesh_doctor.py element_volumes --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``fix_elements_orderings``
""""""""""""""""""""""""""

It can happen that an exported mesh does not abide by the ``vtk`` orderings.
The ``fix_elements_orderings`` module can rearrange the nodes of given types of elements.
This can be convenient if you cannot regenerate the mesh.

.. program-output:: python mesh_doctor.py fix_elements_orderings --help
   :cwd: ../../../oreComponents/python/modules/geosx_mesh_doctor

``generate_cube``
"""""""""""""""""

This module conveniently generates cubic meshes in ``vtk``.
It can also generate fields with simple values.
This tool can also be useful to generate a first mesh that will be refined or customized.

.. program-output:: python mesh_doctor.py generate_cube --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``generate_fractures``
""""""""""""""""""""""

For a conformal to be defined in a mesh, ``geos`` requires the mesh to be split at the faces where the fracture gets across the mesh.
The ``generate_fractures`` module will split the mesh and generate the multi-block ``vtk`` files.

.. program-output:: python mesh_doctor.py generate_fractures --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``generate_global_ids``
"""""""""""""""""""""""

When running ``geos`` in parallel, using global ids can be used to refer to data across the ranks.
The ``generate_global_ids`` can generate global ids for the input ``vkt`` mesh.

.. program-output:: python mesh_doctor.py generate_global_ids --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``non_conformal``
"""""""""""""""""

This module will detect elements close enough (there's a user defined threshold) but are not in front of each other (another threshold can be defined).
This module can be a little time consuming.

.. program-output:: python mesh_doctor.py non_conformal --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``self_intersecting_elements``
""""""""""""""""""""""""""""""

Some meshes can have elements that auto-intersect.
This module will display the elements that have faces intersecting.

.. program-output:: python mesh_doctor.py self_intersecting_elements --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor

``supported_elements``
""""""""""""""""""""""

``geos`` supports a specific set of elements.
Let's cite the standard elements like `tetrahedra`, `wedges`, `pyramids` or `hexahedra`.
But also prismes up to 11 faces.
The ``supported_elements`` check will validate that no unsupported element is defined in the input mesh.
Also, ``geos`` supports the generic ``VTK_POLYHEDRON``/``42`` elements, which are converted on the fly into one of the elements described above.
The ``supported_elements`` check will also verify that those ``VTK_POLYHEDRON`` cells can effectively get converted.

.. program-output:: python mesh_doctor.py supported_elements --help
   :cwd: ../../../coreComponents/python/modules/geosx_mesh_doctor
