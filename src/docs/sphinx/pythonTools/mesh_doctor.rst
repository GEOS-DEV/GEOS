``mesh_doctor``
---------------

``mesh_doctor`` is a ``python`` executable that can be used through the command line to perform various checks, validations and tiny fixes to the ``vtk`` mesh that are meant to be used in ``geos``.
``mesh_doctor`` is organized as a collection of modules with their dedicated set of options.
The current page will introduce those modules, but the details and all the arguments can be retrieved by using the ``--help`` option for each module.

Modules
^^^^^^^

To list all the modules available through ``mesh_doctor`` you can simply use the ``--help`` option which will provide all the modules as well as a quick summuray.

.. program-output:: mesh_doctor --help

Then, it you are interested in a specific module, you can ask for its documentation using the ``mesh_doctor module_name --help`` pattern.
For example

.. program-output:: mesh_doctor collocated_nodes --help

``mesh_doctor`` loads its module dynamically.
If a module can't be loaded, ``mesh_doctor`` will proceed and try to load the other modules.
If you see a message like

.. code-block:: bash

    [2023-09-15 14:21:25,625][WARNING] Could not load module "collocated_nodes": No module named 'vtkmodules'

then most likely ``mesh_doctor`` could not load the ``collocated_nodes`` modules because the ``vtk`` modules was not found.
Consistently, the documentation for module ``collocated_nodes`` will not be displayed.
You can solve this issue by installing the dependencies of ``mesh_doctor`` defined in its ``requirements.txt`` file (``python -m pip install -r requirements.txt``).

Here is now a summury of all the modules available.

``collocated_nodes``
""""""""""""""""""""
.. program-output:: mesh_doctor collocated_nodes --help


``element_volumes``
"""""""""""""""""""
.. program-output:: mesh_doctor element_volumes --help


``fix_elements_orderings``
"""""""""""""""""""""""""
.. program-output:: mesh_doctor fix_elements_orderings --help


``generate_cube``
"""""""""""""""""
.. program-output:: mesh_doctor generate_cube --help


``generate_fractures``
""""""""""""""""""""""
.. program-output:: mesh_doctor generate_fractures --help


``generate_global_ids``
"""""""""""""""""""""""
.. program-output:: mesh_doctor generate_global_ids --help


``non_conformal``
"""""""""""""""""
.. program-output:: mesh_doctor non_conformal --help


``self_intersecting_elements``
"""""""""""""""""""""""""""""
.. program-output:: mesh_doctor self_intersecting_elements --help


``supported_elements``
""""""""""""""""""""""

``geos`` supports a specific set of elements.
Let's cite the standard elements like `tetrahedra`, `wedges`, `pyramids` or `hexahedra`.
But also prismes up to 11 faces.
The ``supported_elements`` check will validate that no unsupported element is defined in the input mesh.
Also, ``geos`` supports the generic ``VTK_POLYHEDRON``/``42`` elements, which are converted on the fly into one of the elements described above.
The ``supported_elements`` check will also verify that those ``VTK_POLYHEDRON`` cells can effectively get converted.

.. program-output:: mesh_doctor supported_elements --help
