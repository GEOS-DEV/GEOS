.. _pfw_geometryObjects:

############################################
pfw_geometryObjects
############################################


Overview
========================
The ``pfw_geometryObjects.py`` script defines the geometry object framework used by the Particle File Writer (PFW) to populate GEOS-MPM material particle files. 
Rather than requiring users to explicitly list every particle location, PFW works by sampling the simulation domain and asking a collection of geometry 
objects whether each candidate particle lies inside a material region, near a surface, or within a special feature such as a cohesive interface. 
This makes it possible to build complex simulation setups from reusable geometric building blocks.

A key part of this workflow is the ``make_objects`` command defined in the input file. Instead of constructing all geometry directly inside ``particleFileWriter.py``, 
users provide either an ``objects`` list or a ``make_objects`` function that returns the relevant geometry objects for the simulation. This separation is useful 
because geometry construction can become expensive for complex cases, especially when Voronoi tessellations, many particles, or many distinct regions are involved. 
By isolating geometry construction in ``make_objects``, the workflow is more modular, easier to debug, and better suited for rank-aware object generation in parallel runs.

This file also defines common defaults, validation utilities, transform operations, and set operations that geometry objects rely on. 
Together, these components provide the interface that allows PFW to assign material IDs, particle types, surface information, and other 
spatially varying particle properties in a consistent way.

- :ref:`Geometry Object Framework <geometry_object_framework>`

  - :ref:`What Geometry Objects Must Provide <what_geometry_objects_must_provide>`

- :ref:`Default Property Definitions <default_geometry_object_properties>`

- :ref:`MPI and Parallel Context <mpi_and_parallel_context>`

- :ref:`Type Checking Utilities <type_checking_utilities>`

  - :ref:`Scalar Type Checking <scalar_type_checking>`
  - :ref:`Array Type Checking <array_type_checking>`

- :ref:`Utility Functions <utility_functions>`

  - :ref:`Logging and File Utilities <logging_and_file_utilities>`
  - :ref:`Spatial Utilities <spatial_utilities>`
  - :ref:`Basis Construction Utilities <basis_construction_utility_function>`

- :ref:`Transforms <transforms>`

  - :ref:`Homogeneous Transform Matrices <homogeneous_transform_matrices>`
  - :ref:`Transform Wrapper <transform_wrapper>`

- :ref:`Set Operations <set_operations>`

  - :ref:`SetOperation Base Class <setoperation_base_class>`
  - :ref:`Union <union_operation>`
  - :ref:`Intersection <intersection_operation>`
  - :ref:`Difference <difference_operation>`

- :ref:`A Simple Example <simple_example>`

  - :ref:`Example 1 <example_1>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _geometry_object_framework:

Geometry Object Framework
-------------------------

This section explains the role of geometry objects in the Particle File Writer workflow and describes the interface that user-defined geometry objects are expected to support.

.. _why_make_objects_exists:

Why ``make_objects`` Exists
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Particle File Writer does not directly define simulation geometry itself. Instead, it reads an input file and expects that file to provide either:

- an explicit ``objects`` list, or
- a ``make_objects`` function that returns the geometry objects for the run

This design is important because geometry construction can be one of the most expensive parts of the setup process. For simple cases, users may define the objects directly. For more complicated cases, especially when geometry generation involves many regions, procedural construction, or expensive operations such as Voronoi tessellations, it is more efficient and more maintainable to generate those objects in a dedicated function.

Using ``make_objects`` also makes it possible to construct geometry in a rank-aware way. In large simulations, each MPI rank may only need the subset of objects relevant to its assigned domain slice. This reduces setup cost and avoids unnecessary duplication of complex geometry on every rank.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``make_objects`` separates geometry construction from particle writing.
- This improves modularity and makes complex geometries easier to manage.
- Rank-aware object generation can significantly reduce setup overhead for large simulations.
- Users can still provide a direct ``objects`` list for simpler cases.


:ref:`Back to top <pfw_geometryObjects>`

=====

.. _what_geometry_objects_must_provide:

What Geometry Objects Must Provide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At minimum, a geometry object must provide the functionality needed for the Particle File Writer to evaluate each candidate particle location. 

For every point in the domain, PFW queries the geometry object to determine:

- whether the point lies inside the object, outside the object, or near a surface
- which material the particle should belong to
- which contact group the particle should be assigned to
- what particle type and additional properties (e.g., damage, porosity, temperature, surface information) should be assigned to the particle

The geometry object answers these questions through its methods and attributes. 
Based on these answers, PFW decides whether to create a particle at that location and how to assign its properties.
In practice, this means that geometry objects must implement methods such as:

- ``isInterior(pt, skinDepth)``
- ``getSubregions()``

and may optionally provide methods or attributes such as:

- ``getVelocity(pt)``
- ``getDamage(pt)``
- ``getPorosity(pt)``
- ``getSurfaceNormal(pt)``
- ``getSurfacePosition(pt)``

These methods allow PFW to assign particle properties consistently during particle generation. Addition, many objects also provide optional methods or attributes used to assign particle properties.

:ref:`Back to top <pfw_geometryObjects>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

.. _default_geometry_object_properties:

Default Geometry Object Properties
=================================

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Parameter
     - Type
     - Description

   * - ``_defaultMat``
     - lookup
     - Index into ``pfw["materials"]``. ``0`` corresponds to the first material.

   * - ``_defaultGroup``
     - label
     - Contact group identifier. Used to distinguish interacting bodies (e.g., substrate vs. indentor).

   * - ``_defaultVelocity``
     - vector
     - Initial particle velocity in ``[x, y, z]`` (mm/us). Default is zero (at rest).

   * - ``_defaultDamage``
     - scalar
     - Damage state of the particle. ``0`` = intact, ``1`` = fully damaged.

   * - ``_defaultPorosity``
     - scalar
     - Initial porosity. ``0`` = fully dense, values > 0 indicate void space.

   * - ``_defaultStrengthScale``
     - scalar
     - Multiplier applied to material strength. ``1.0`` means no scaling.

   * - ``_defaultTemperature``
     - scalar
     - Initial temperature in Kelvin.

   * - ``_defaultSurfaceFlag``
     - flag
     - Surface classification: ``0`` = internal, ``1`` = fully damaged, ``2`` = surface, ``3`` = cohesive surface.

   * - ``_defaultSurfaceNormal``
     - vector
     - Surface normal direction. A zero vector disables explicit normals (solver computes them implicitly).

   * - ``_defaultSurfacePosition``
     - vector
     - Position of the surface relative to the particle, used in contact calculations.

   * - ``_defaultMatDir``
     - tensor
     - Material orientation tensor. (identity = no preferred direction / isotropic default)

   * - ``_defaultSurfaceTraction``
     - vector
     - Applied surface traction. Default is zero which disables explicit normals (solver computes implicitly)

   * - ``_defaultParticleType``
     - enum
     - Particle interpolation type: ``0`` = SinglePoint, ``1`` = SinglePoint B-spline, ``2`` = CPDI, ``3`` = CPTI, ``4`` = CPDI2. ***See note below.

   * - ``_defaultCZTag``
     - label
     - Cohesive zone (CZ) identifier. ``0`` indicates no cohesive behavior.

.. code-block:: python

   _defaultMat = 0  
   _defaultGroup = 0  
   _defaultVelocity = np.array([0.0, 0.0, 0.0])  
   _defaultDamage = 0.0  
   _defaultPorosity = 0.0  
   _defaultStrengthScale = 1.0  
   _defaultTemperature = 300.0  
   _defaultSurfaceFlag = 2  
   _defaultSurfaceNormal = np.array([0.0, 0.0, 0.0])   
   _defaultSurfacePosition = np.array([0.0, 0.0, 0.0]) 
   _defaultMatDir = np.array([[1.0, 0.0, 0.0],
                             [0.0, 1.0, 0.0],
                             [0.0, 0.0, 1.0]])    # tensor: identity here
   _defaultSurfaceTraction = np.array([0.0, 0.0, 0.0])    
   _defaultParticleType = 2     #see note below
   _defaultCZTag = 0  

The Particle File Writer usually checks for object methods such as ``getVelocity`` or ``getDamage`` first, and falls back to fixed attributes when available. This allows both simple constant-property objects and more advanced spatially varying objects to use the same interface.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``isInterior`` and ``getSubregions`` are the core requirements.
- Many other property getters are optional, depending on which particle fields are requested.
- Geometry objects can support either constant properties or spatially varying callable behavior.
- Consistency of this interface is what allows PFW to treat many different object types uniformly.

.. raw:: html

   <h4><u>*** Some Notes on Particle type</u></h4>

Particle interpolation methods define how particles interact with the background grid in MPM. 

.. list-table:: Particle Interpolation Types
   :widths: 20 30 50
   :header-rows: 1

   * - ID
     - Name
     - Description

   * - ``0``
     - SinglePoint
     - Simplest approach. Particle treated as a point; fastest but least accurate, prone to noise and grid-crossing errors. computationally inexpensive but can suffer from numerical artifacts during large deformation.

   * - ``1``
     - SinglePoint B-spline
     - Improves upon SinglePoint by using smoother interpolation functions, reducing noise and improving stability. Uses smoother B-spline interpolation; reduces numerical noise compared to SinglePoint.

   * - ``2``
     - CPDI
     - Convected Particle Domain Interpolation; represents each particle as a finite domain that deforms with the material, significantly improving accuracy and robustness for large-deformation problems.

   * - ``3``
     - CPTI
     - Convected Particle Tetrahedral Interpolation; extends CPDI to use a tetrahedral representation for improved geometric fidelity.

   * - ``4``
     - CPDI2
     - Higher-order extension of CPDI that further improves accuracy and stability by better capturing particle domain deformation, particularly under extreme distortion.



:ref:`Back to top <pfw_geometryObjects>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _mpi_and_parallel_context:

MPI and Parallel Context
------------------------

This section initializes the MPI communication context used by ``pfw_geometryObjects.py``.

The script imports ``MPI`` from ``mpi4py`` and defines:

- ``comm`` as the MPI communicator
- ``g_rank`` as the current process rank
- ``num_ranks`` as the total number of MPI processes

These values are used by utility functions and geometry-generation workflows that need rank-aware behavior, such as logging or distributed object creation.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``g_rank`` is used instead of a local rank variable so it can be accessed throughout the module.
- Geometry construction may use rank-aware logic even before the particle file itself is written.
- MPI-aware logging is supported through rank-specific output files.

:ref:`Back to top <pfw_geometryObjects>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _type_checking_utilities:

Type Checking Utilities
-----------------------

This section defines helper utilities for validating scalar and array inputs used by geometry objects.

.. _scalar_type_checking:

Scalar Type Checking
^^^^^^^^^^^^^^^^^^^^

The ``type_check_scalar`` function verifies that a scalar variable has an allowed type. Optionally, it can also allow callable inputs.

.. code-block:: python

   def type_check_scalar(s, varName, scalarType, canBeFunc=False):

If the provided value does not match one of the allowed scalar types, the function raises a ``TypeError`` with a descriptive message. When ``canBeFunc=True``, callable objects are also accepted.

.. raw:: html

   <h4><u>Notes</u></h4>

- This is useful for geometry parameters that may be either fixed values or functions of position.
- Error messages include the variable name and the expected type list.
- This helps catch malformed object definitions early.

:ref:`Back to top <pfw_geometryObjects>`

=====

.. _array_type_checking:

Array Type Checking
^^^^^^^^^^^^^^^^^^^

The ``type_check_array`` function validates array-like inputs such as vectors or small tensors.

.. code-block:: python

   def type_check_array(arr, varName, size, scalarType, canBeFunc=False):

This function checks:

- whether the input is a list or NumPy array
- whether its size matches the expected dimension or one of several allowed dimensions
- whether each entry has an allowed scalar type

If ``canBeFunc=True``, callable inputs are also accepted.

.. raw:: html

   <h4><u>Notes</u></h4>

- This is useful for validating coordinates, dimensions, directions, and other structured geometry inputs.
- The ``size`` argument may be a single integer or a list/tuple of allowed sizes.
- As with scalar checking, this helps users catch bad object definitions before particle generation begins.

:ref:`Back to top <pfw_geometryObjects>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">




.. _utility_functions:

Utility Functions
-----------------

This section collects helper functions used across the geometry-object framework for logging, file handling, spatial checks, periodic mapping, interpolation, and local basis construction.


=====

.. _logging_and_file_utilities:

Logging and File Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^

The module includes several helper functions for logging and file inspection:

- ``log2file(msg)``
  Writes a rank-tagged message to a rank-specific MPI log file.

- ``print2file(file_name, text)``
  Appends a line of text to a standard file.

- ``countFileLines(fname)``
  Counts the number of lines in a file efficiently using buffered binary reads.

- ``fOffsetLineNum(fname, line_number)``
  Returns the byte offset associated with a given line number using buffered reads.

- ``fileOffsetFromLineNum(fname, line_number)``
  Computes a character offset for a given line number using standard text iteration.

These functions support debugging, logging, and large-file operations used elsewhere in the PFW workflow.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``log2file`` uses MPI-safe per-rank output.
- ``countFileLines`` is used in other parts of the PFW workflow to inspect generated files.
- Offset utilities are useful when working with large particle files or rank-local outputs.

:ref:`Back to top <pfw_geometryObjects>`

=====

.. _spatial_utilities:

Spatial Utilities
^^^^^^^^^^^^^^^^^

The spatial utility functions (functions that operate on positions in space) support common geometric checks and periodic-domain handling:

- ``mapToRange(x, xmin, xmax)``
  Maps coordinates into a periodic range.

- ``inside_box(x, x0, dx, periodic)``
  Checks whether a point lies inside a box.

- ``outside(x, x0, dx, periodic)``
  Checks whether a point lies outside a box, with periodic handling where requested.

- ``smoothStep3(r, rf, rf0)``
  Provides a smooth cubic transition between two radii, useful for interpolation or smoothed geometric transitions.

These helpers are widely applicable to many geometry definitions.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``mapToRange`` is especially useful for periodic domains.
- ``inside_box`` and ``outside`` support quick bounding checks before more expensive geometry tests.
- ``smoothStep3`` provides a smooth, bounded transition useful for soft region definitions.

:ref:`Back to top <pfw_geometryObjects>`


=====

.. _basis_construction_utility_function:

Basis Construction Utility Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section provides a utility function for constructing a local coordinate system from a single direction vector. The function ``generate_orthonormal_axes(v)`` constructs an orthonormal basis from a given nonzero vector.

.. code-block:: python

   def generate_orthonormal_axes(v):

The input vector defines the first axis, and the remaining two axes are constructed to be orthogonal and normalized.

This allows for defining local coordinate systems, anisotropic material directions, or rotated geometry frames.

.. raw:: html

   <h4><u>Notes</u></h4>

- The input vector must be nonzero.
- The output basis is right-handed.
- This utility is primarily used when geometry or material behavior depends on a preferred direction.



:ref:`Back to top <pfw_geometryObjects>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _transforms:

Transforms
----------

This section defines geometric transforms and transform wrappers that can be applied to existing objects. 
A transform is a mathematical operation that changes the position, orientation, or size of a geometry object without altering its intrinsic shape. 
Common transforms include translation (moving), rotation, reflection, and scaling.
**Homogeneous transform matrices** define these operations mathematically, while **transform wrappers** use them 
to evaluate geometry queries in transformed coordinate systems without modifying the underlying object.

=====

.. _homogeneous_transform_matrices:

Homogeneous Transform Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The functions ``translate``, ``scale``, ``reflect``, and ``rotate`` generate 4x4 homogeneous transformation matrices for common geometric operations.

- ``translate(dx)``
  Creates a translation matrix.

- ``scale(ds)``
  Creates a scaling matrix. A scalar input is expanded isotropically.

- ``reflect(a0, x0)``
  Creates a reflection transform about a plane with normal ``a0`` passing through ``x0``.

- ``rotate(alpha, a0, x0)``
  Creates a rotation transform about axis ``a0`` through point ``x0``.

These matrices provide a consistent way to compose geometric transformations.

.. raw:: html

   <h4><u>Notes</u></h4>

- The transforms are expressed in homogeneous coordinates.
- Reflection and rotation are defined about user-specified geometric references.
- These helpers are typically used through wrapper objects rather than directly by the particle writer.

=====

.. _transform_wrapper:

Transform Wrapper
^^^^^^^^^^^^^^^^^

- The ``transform`` class wraps an existing geometry object and applies a geometric transform to it.

- Transform wrappers are useful for creating rotated, translated, reflected, or scaled copies of existing shapes without redefining the object.

- The wrapper computes and stores the **inverse transform**, which is used to map query points into the object's original coordinate system.

  - An inverse transform reverses a transformation, mapping a point from the transformed coordinate system back into the original coordinate system.
  - Geometry objects are defined in their original coordinate system, so they cannot directly evaluate transformed points.
  - Instead of modifying the object, the wrapper applies the inverse transform to each query point before passing it to the base object.

- As a result, geometric queries such as ``isInterior`` are evaluated in the object's local coordinate system:

  - A point in the transformed space is mapped back using the inverse transform
  - The base object evaluates the query using its original definition

- The wrapper delegates the following operations to the underlying object after transforming the query point:

  - ``isInterior``
  - ``getSurfaceNormal``
  - ``getSurfacePosition``
  - ``getGroup``
  - ``getMatDir``
  - ``getDamage``
  - ``getPorosity``
  - ``getTemperature``
  - ``getSurfaceTraction``

- ``xMin`` and ``xMax`` currently default to infinite bounds, so transformed objects may not provide tight x-range filtering by default.

.. code-block:: python

   class transform(BaseWrapper):

This design allows users to reuse a single base geometry definition while generating multiple transformed configurations.



:ref:`Back to top <pfw_geometryObjects>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _set_operations:

Set Operations
--------------

This section defines composite geometry objects built from two sub-objects using Boolean-style set operations.

.. _setoperation_base_class:

SetOperation Base Class
^^^^^^^^^^^^^^^^^^^^^^^

The ``SetOperation`` base class provides the common structure for combining two geometry objects, ``subObjA`` and ``subObjB``.
The ``SetOperation`` class defines a base interface for combining geometry objects. Child classes such as ``union``, ``intersection``, and ``difference`` 
inherit from the SetOperations class and implement the corresponding set-theoretic behavior.

.. code-block:: python

   class SetOperation(Geometry):

It stores the two sub-objects, defines whether inherited defaults should come from object A or B, and provides shared implementations for many property-access methods such as:

- ``getSurfaceNormal``
- ``getSurfacePosition``
- ``getGroup``
- ``getDamage``
- ``getStrengthScale``
- ``getMatDir``
- ``getSubregions``

Derived classes override ``isInterior`` and may refine other behaviors depending on the set operation.

.. raw:: html

   <h4><u>Notes</u></h4>

- This base class reduces duplication among union/intersection/difference operations.
- Property inheritance can default to either sub-object A or B.
- Composite objects can still participate in the normal PFW geometry workflow.

=====

.. _union_operation:

Union (A ∪ B)
^^^^^^^^^^^^^^^

The ``union`` class creates a geometry object representing the union of two sub-objects.

.. code-block:: python

   class union(SetOperation):

A point is considered inside the union if it lies inside either object. Surface normals and positions are chosen based on which contributing object defines the relevant surface at that point.

.. raw:: html

   <h4><u>Notes</u></h4>

- Useful for building larger composite regions from simpler primitives.
- If both sub-objects contribute a surface at the same point, the surface normal may be averaged.
- Interior/surface classification depends on the flags returned by the sub-objects.

=====

.. _intersection_operation:

Intersection (A ∩ B)
^^^^^^^^^^^^^^^^^^^^

The ``intersection`` class creates a geometry object representing the overlap of two sub-objects.

.. code-block:: python

   class intersection(SetOperation):

A point is considered inside the intersection only if it lies inside both objects. Surface information is selected by comparing the surface positions returned by the two sub-objects.

.. raw:: html

   <h4><u>Notes</u></h4>

- Useful for trimming one object by another or defining overlap regions.
- Surface selection is based on which surface is geometrically closer.
- ``getCZTag`` is currently hardcoded to return the value from sub-object A.

=====

.. _difference_operation:

Difference (A - B)
^^^^^^^^^^^^^^^^^^

The ``difference`` class represents one object with another object removed.

.. code-block:: python

   class difference(SetOperation):

This is intended to define regions such as holes, voids, or cut-outs, although the current implementation appears less complete than ``union`` and ``intersection``.

.. raw:: html

   <h4><u>Notes</u></h4>

- This operation is conceptually useful for defining subtracted geometry.
- The current implementation should be reviewed carefully before relying on it in production cases.
- Unlike the other set operations, its use of member names suggests it may need cleanup or validation.



.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


:ref:`Back to top <pfw_geometryObjects>`


.. _simple_example:

A Simple Example
----------------

Below is an example of how to define and use a simple box geometry object within the PFW workflow.

.. _example_1:

Example 1: Defining and Using the ``box`` Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Class definition**

.. code-block:: python

   #############################################
   class box(Geometry):
       """
       Geometry object for creating a grid-aligned box defined by two corners.
       """

       def __init__(self, 
                    name, 
                    x0, 
                    x1, 
                    vel=_defaultVelocity, 
                    mat=_defaultMat, 
                    group=_defaultGroup,
                    particleType=_defaultParticleType,
                    dim=3, 
                    flaggedSurfaces=[True, True, True, True, True, True]):

           # __init__ is your constructor:
           # it runs when the object is created and initializes its properties
           #---------------------------------------------------------------------------
           super().__init__(name,
                            vel=vel,
                            mat=mat,
                            group=group,
                            particleType=particleType)

           # super().__init__ calls the base ``Geometry`` class so shared properties (material, group, velocity, etc.) are initialized correctly.
           #---------------------------------------------------------------------------

           # Validate your inputs
           type_check_scalar(dim, "Dimension", int)
           type_check_array(x0, "x0", dim, floatType)
           type_check_array(x1, "x1", dim, floatType)
           #---------------------------------------------------------------------------


           # flaggedSurfaces controls which faces of the box are treated as "active surfaces"
           # For a 3D box: [xmin, ymin, zmin, xmax, ymax, zmax]
           # True = surface is considered, False = ignored
           type_check_array(flaggedSurfaces, "Flagged surfaces", 2*dim, bool)
           #---------------------------------------------------------------------------


           # Store geometry data on the object
           # assigns things like ``self.xmin``, ``self.xmax``, etc. defines the actual geometry and stores it on the object.

           self.dim = dim                          # number of spatial dimensions (usually 3)
           self.x0 = np.asarray(x0[:self.dim])     # first corner
           self.x1 = np.asarray(x1[:self.dim])     # opposite corner

           # Compute min/max bounds of the box
           self.xmin = np.minimum(self.x0, self.x1)
           self.xmax = np.maximum(self.x0, self.x1)

           # Store which faces are active
           self.flaggedSurfaces = np.asarray(flaggedSurfaces[:2*self.dim])
           #---------------------------------------------------------------------------


       def isInterior(self, pt, skinDepth):
           # REQUIRED method:
           # tells PFW whether a point is inside (0), on surface (2), or outside (-1)

           # note**** : getSubregions() is also required by the PFW workflow, but it is not shown here because it is 
           # typically implemented in the base class (Geometry) or inherited through wrappers.
           # see the 'What Geometry Objects Must Provide' section above

           x = np.asarray(pt[:self.dim])

           if np.all(np.logical_and(x >= self.xmin, x < self.xmax)):

               # Check if near a flagged surface
               s = np.hstack((x <= self.xmin + skinDepth,
                              x >= self.xmax - skinDepth))

               if np.any(np.logical_and(s, self.flaggedSurfaces)):
                   return _defaultSurfaceFlag   # surface particle
               else:
                   return 0                    # interior particle

           return -1                            # outside object
           #---------------------------------------------------------------------------


       def getSurfaceNormal(self, pt):
           # OPTIONAL (but commonly used):
           # returns outward normal direction of nearest surface

           x = np.asarray(pt[:self.dim])
           m = np.zeros((2*self.dim))
           m[np.logical_not(self.flaggedSurfaces)] = np.inf

           dx = np.hstack((self.xmin - x, self.xmax - x)) + m
           minI = np.argmin(np.absolute(dx))

           d = minI % self.dim
           s = -1 if int(math.floor(minI / self.dim) == 0) else 1

           surfaceNormal = np.array([0.0, 0.0, 0.0])
           surfaceNormal[d] = s

           return surfaceNormal
           #---------------------------------------------------------------------------


       def getSurfacePosition(self, pt):
           # OPTIONAL:
           # returns distance to nearest surface (used in contact)

           x = np.asarray(pt[:self.dim])
           m = np.zeros((2*self.dim))
           m[np.logical_not(self.flaggedSurfaces)] = np.inf

           dx = np.hstack((self.xmin - x, self.xmax - x)) + m
           minI = np.argmin(np.absolute(dx))

           d = minI % self.dim

           surfacePosition = np.array([0.0, 0.0, 0.0])
           surfacePosition[d] = dx[minI]

           return surfacePosition
           #---------------------------------------------------------------------------


       def xMin(self):
           # OPTIONAL optimization:
           # gives lower x-bound for fast filtering
           return self.xmin[0]


       def xMax(self):
           # OPTIONAL optimization:
           # gives upper x-bound for fast filtering
           return self.xmax[0]



**Usage**

The following is defined in your :doc:`pfw_input file <pfw_input>`:

.. code-block:: python

   bar = geom.box(
       'bar',
       [0.0, -domainWidth/2, -domainWidth/2],
       [domainLength, domainWidth/2, domainWidth/2],
       vel=[0.0, 0.0, 0.0],
       mat=0,
       group=0
   )

   pfw["objects"] = [bar]


:ref:`Back to top <pfw_geometryObjects>`