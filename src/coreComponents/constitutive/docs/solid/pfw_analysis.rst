.. _pfw_analysis:

############################################
PFW Analysis
############################################

Overview
========================

The ``pfw_analysis`` script provides a structured workflow for analyzing GEOS Material Point Method (MPM) simulation outputs generated using the Particle File Writer (PFW) framework. 
It integrates utilities for locating valid simulation directories, importing reaction and box-average history data, and organizing these data into consistent field objects. 
The script includes post-processing tools for filtering, trimming, and subsampling datasets, as well as functions for computing derived quantities such as strain, stress, and pressure. 
Visualization helpers and configurable analysis options support flexible plotting and data export, while the ``MPMJob`` and ``DataObj`` classes provide a unified interface for managing simulation metadata and field data. 
Together, these components enable efficient, reproducible analysis of MPM simulation results.


Overview
========================

- **Data Import/Export Functions**
  
  - :ref:`has_mpm_file <has_mpm_file>`
  - :ref:`format_file_paths <format_file_paths>`
  - :ref:`read_from_reaction_file <read_from_reaction_file>`
  - :ref:`read_from_box_average_file <read_from_box_average_file>`
  - :ref:`write_data_to_csv_file <write_data_to_csv>`
  - :ref:`write_data_to_console_file <write_data_to_console>`


- **Data Post Processing Classes**
  
  - :ref:`Trim <Trim>`
  - :ref:`RemoveNonMonotonicEntries <RemoveNonMonotonicEntries>`
  - :ref:`MedianFilter <MedianFilter>`
  - :ref:`SubSample <SubSample>`
  - :ref:`compute_domain_strain <compute_domain_strain>`
  - :ref:`compute_domain_stress <compute_domain_stress>`
  - :ref:`compute_pressure <compute_pressure>`

- **Data Visualization and Plotting Function**
  
  - :ref:`lighten_color <lighten_color>`

- **Supporting Classes and Script Execution**
  
  - :ref:`AnalysisOptions <AnalysisOptions>`
  - :ref:`MPMJob <MPMJob>`
  - :ref:`DataObj <DataObj>`
  - :ref:`main_execution_block <main_execution_block>`







Data Import/Export
==================

The ``pfw_analysis`` script includes a set of helper functions for locating valid MPM run directories, reading simulation output files, and exporting processed data. 
These routines support the post-processing workflow by identifying directories that contain ``pfw_input`` files, loading reaction-force and box-average history data, 
and writing formatted results either to a CSV file or directly to the console.


=====

.. raw:: html

   <h3><u>has_mpm_file function</u></h3>

.. _has_mpm_file:

The ``has_mpm_file`` function checks whether a directory contains a Python input file matching the naming convention used for Particle File Writer runs.
This function is used to determine whether a directory should be treated as a valid MPM run directory during recursive file discovery.


.. code-block:: python

   def has_mpm_file(path):

This function loops over the contents of the directory and searches for any file whose name matches the regular expression:

.. code-block:: text

   .*pfw_input_.*\.py

If such a file is found, the function returns ``True``. Otherwise, it returns ``False``.


Parameters
^^^^^^^^^^

``path``
  Path to the directory being checked.

Returns
^^^^^^^

``True``
  If the directory contains a ``pfw_input`` Python file.

``False``
  If no matching file is found.


:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>format_file_paths function</u></h3>

.. _format_file_paths:


The ``format_file_paths`` function searches through one or more run locations and collects directories that appear to contain valid MPM input files.
- It is intended to recursively identify simulation directories that contain Particle File Writer input files.

..- In the current code, the recursive call uses ``formatFilePaths(...)`` rather than ``format_file_paths(...)``. If ``formatFilePaths`` the author has not seen this defined elsewhere in the code, and mig should be corrected.


.. code-block:: python

   def format_file_paths(runLocations, paths):

For each combination of ``runLocation`` and ``path``, the function constructs a candidate directory path. If that path is a directory, it checks whether the directory contains a valid ``pfw_input`` file by calling ``has_mpm_file``.

- If a valid MPM input file is found, the directory is added to the returned file list.
- If not, the function continues searching recursively through the contents of that directory.


Parameters
^^^^^^^^^^

``runLocations``
  A list of base directory paths to search.

``paths``
  A list of relative paths or subpaths to evaluate beneath each run location.


Returns
^^^^^^^

A list of directory paths corresponding to valid MPM run directories.

:ref:`Back to top <pfw_analysis>`

=====


.. raw:: html

   <h3><u>read_from_reaction_file function</u></h3>

.. _read_from_reaction_file:


The ``read_from_reaction_file`` function reads reaction history data from a comma-separated output file and stores each column as a ``DataObj``.

.. code-block:: python

   def read_from_reaction_file(filename):


This function uses ``numpy.genfromtxt`` to load the contents of the file using a comma delimiter. The first row is assumed to contain column headers, so the returned data begins from row index 1.

The file is expected to contain columns in the following order:

#. Time    : simulation time micro seconds
#. ``F00`` : the deformation gradient in x direction
#. ``F11`` : the deformation gradient in y direction
#. ``F22`` : the deformation gradient in z direction
#. ``Lx``  : current length of domain in x direction
#. ``Ly``  : current length of domain in y direction
#. ``Lz``  : current length of domain in z direction
#. ``Rxm`` : normal component of the reaction at faces in x-
#. ``Rxp`` : normal component of the reaction at faces in x+
#. ``Rym`` : normal component of the reaction at faces in y-
#. ``Ryp`` : normal component of the reaction at faces in y-
#. ``Rzm`` : normal component of the reaction at faces in z-
#. ``Rzp`` : normal component of the reaction at faces in z-
#. ``L00`` : grid velocity gradient in x direction
#. ``L11`` : grid velocity gradient in y direction
#. ``L22`` : grid velocity gradient in z direction

The function creates a ``DataObj`` for each column, although only a subset of these objects is returned.



Parameters
^^^^^^^^^^

``filename``
  Path to the reaction history CSV file.


Returns
^^^^^^^

The function returns:
``time``, ``F00``, ``F11``, ``F22``, ``Lx``, ``Ly``, ``Lz``, ``Rxm``, ``Rxp``, ``Rym``, ``Ryp``, ``Rzm``, ``Rzp``

Notes
^^^^^

- Although ``L00``, ``L11``, and ``L22`` are created inside the function, they are not currently included in the return statement.
- This function assumes the reaction file follows a fixed column ordering.

:ref:`Back to top <pfw_analysis>`

=====


.. raw:: html

   <h3><u>read_from_box_average_file function</u></h3>


.. _read_from_box_average_file:




The ``read_from_box_average_file`` function reads box-averaged simulation history data from a comma-separated output file and stores each column as a ``DataObj``.

.. code-block:: python

   def read_from_box_average_file(filename):



This function uses ``numpy.genfromtxt`` to load the file using a comma delimiter. The first row is assumed to contain headers, so the data begins from row index 1.

The file contains columns corresponding to:

#. Time    : simulation time micro seconds
#. ``Sxx`` : normal stress in the x direction
#. ``Syy`` : normal stress in the y direction
#. ``Szz`` : normal stress in the z direction
#. ``Sxy`` : shear stress in the x direction
#. ``Syz`` : shear stress in the y direction
#. ``Sxz`` : shear stress in the z direction
#. Density : sample density in mm mg\ :sup:`-3`
#. Damage  : a binary indicating undamaged (0) or damaged (1)
#. Internal energy : units of mJ mm\ :sup:`-3`
#. Kinetic energy : units of mJ mm\ :sup:`-3`
#. ``epxx`` : normal plastic strain in the x direction
#. ``epyy`` : normal plastic strain in the y direction 
#. ``epzz`` : normal plastic strain in the z direction
#. ``epyz`` : shear plastic strain in the yz direction
#. ``epxz`` : shear plastic strain in the xz direction
#. ``epxy`` : shear plastic strain in the xy direction
#. Material volume mm:sup:`3`

Each column is converted into a ``DataObj`` and returned for downstream analysis or plotting.



Parameters
^^^^^^^^^^

``filename``
  Path to the box-average history CSV file.


Returns
^^^^^^^

``time``, ``sxx``, ``syy``, ``szz``, ``sxy``, ``syz``, ``sxz``, ``density``, ``damage``, ``internalEnergy``, ``kineticEnergy``, ``epxx``, ``epyy``, ``epzz``, ``epyz``, ``epxz``, ``epxy``, ``matVol``

Notes
^^^^^

This function assumes the box-average file follows a fixed column order and that all expected columns are present.

:ref:`Back to top <pfw_analysis>`

=====


.. raw:: html

   <h3><u>write_data_to_csv function</u></h3>


.. _write_data_to_csv:


The ``write_data_to_csv`` function writes a list of ``DataObj`` fields to a comma-separated text file.

.. code-block:: python

   def write_data_to_csv(filename, data_array):


This function first checks that all fields in ``data_array`` contain the same number of entries. It then writes:

- a header row containing the field names
- one data row for each entry in the dataset

Each value is formatted according to the ``format`` attribute stored in the corresponding ``DataObj``.

Parameters
^^^^^^^^^^

``filename``
  Name of the output CSV file.

``data_array``
  A list of ``DataObj`` instances to be written.


Returns
^^^^^^^

This function does not return a value. It writes formatted data directly to the specified file.

Notes
^^^^^

- All fields must have the same length.
- If a field length does not match the expected number of entries, the function raises an assertion error.
- Output formatting is controlled by the ``format`` attribute of each ``DataObj``.

:ref:`Back to top <pfw_analysis>`

=====


.. raw:: html

   <h3><u>write_data_to_console function</u></h3>

.. _write_data_to_console:


The ``write_data_to_console`` function prints a list of ``DataObj`` fields to standard output using the same formatting logic as ``write_data_to_csv``.

.. code-block:: python

   def write_data_to_console(filename, data_array):


This function verifies that all fields have the same number of entries, prints a header row containing the field names, and then prints each row of formatted data to the console.

Parameters
^^^^^^^^^^

``filename``
  This argument is present in the function signature but is not used inside the function.

``data_array``
  A list of ``DataObj`` instances to be printed.



Returns
^^^^^^^

This function does not return a value. It prints the formatted data directly to standard output.

Notes
^^^^^

- The ``filename`` argument is currently unused and may be retained only for interface consistency.
- All fields must contain the same number of entries.
- Formatting is controlled by the ``format`` attribute of each ``DataObj``.

:ref:`Back to top <pfw_analysis>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


Data Post Processing
====================

The ``pfw_analysis`` script includes a collection of data post-processing utilities for modifying, filtering, and deriving simulation output quantities.
These routines support downstream analysis by trimming datasets, removing invalid or repeated entries, filtering noisy signals, subsampling large arrays, and computing derived mechanical quantities such as strain, stress, and pressure.


=====

.. raw:: html

   <h3><u>Trim class</u></h3>

.. _Trim:

The ``Trim`` class truncates an input array to a specified length.
This can be useful when only the first portion of a dataset is needed for analysis or plotting.

.. code-block:: python

   class Trim:
     def __init__(self, length):
       self.length = length

     def postprocess(self, x):
       return x[:self.length]

The ``postprocess`` method returns the first ``length`` entries of the input array.

Parameters
^^^^^^^^^^

``length``
  Number of entries to retain from the beginning of the dataset.

Returns
^^^^^^^

The ``postprocess`` method returns a truncated version of the input array containing only the first ``length`` entries.

Notes
^^^^^

- This operation preserves the original ordering of the data.
- If ``length`` exceeds the size of the input array, the full array is returned.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>RemoveNonMonotonicEntries class</u></h3>

.. _RemoveNonMonotonicEntries:

The ``RemoveNonMonotonicEntries`` class removes entries from a dataset for which a reference array does not increase monotonically.
This is useful when simulation output contains repeated or decreasing independent-variable values, such as duplicated or non-increasing time entries.

.. code-block:: python

   class RemoveNonMonotonicEntries:
     def __init__(self, x_in):
       self.x_in = x_in
       self.maxX = 0.0
       self.mask = np.ones(len(self.x_in), dtype=bool)
       for ii,t in enumerate(self.x_in):
         if (t<=self.maxX):
           self.mask[ii] = False
         else:
           self.maxX = t

     def postprocess(self, x):
       return x[self.mask,...]

During initialization, the class constructs a Boolean mask based on the entries of ``x_in``.
Entries that do not strictly increase relative to the maximum previously encountered value are marked for removal.
The ``postprocess`` method applies this mask to another array.

Parameters
^^^^^^^^^^

``x_in``
  Reference array used to determine which entries are monotonic and should be retained.

Returns
^^^^^^^

The ``postprocess`` method returns the input array with all entries corresponding to non-monotonic values in ``x_in`` removed.

Notes
^^^^^

- This class is commonly used to clean data before interpolation, filtering, or plotting.
- The same mask can be applied to multiple dependent-variable arrays, provided they share the same indexing as ``x_in``.
- The monotonicity check is strictly increasing, so repeated values are removed.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>MedianFilter class</u></h3>

.. _MedianFilter:

The ``MedianFilter`` class applies a median filter to a one-dimensional dataset in order to reduce local noise while preserving large-scale trends.
This can be useful for smoothing noisy simulation signals without introducing as much blurring as a moving average filter.

.. code-block:: python

   class MedianFilter:
     def __init__(self, window_size):
       self.window_size = window_size

     def postprocess(self, x):
       x = np.array(x)
       N = len(x)
       x_out = np.copy(x)
       for n in range(self.window_size+1, N-self.window_size):
         x_out[n] = np.median(x[(n-self.window_size):(n+self.window_size+1)])
       return x_out

The ``postprocess`` method replaces each interior value with the median over a local window centered on that entry.

Parameters
^^^^^^^^^^

``window_size``
  Number of neighboring entries on each side of a point to include in the median window.

Returns
^^^^^^^

The ``postprocess`` method returns a filtered copy of the input array.

Notes
^^^^^

- Boundary entries near the beginning and end of the array are not modified.
- The filtering window has a total width of ``2 * window_size + 1``.
- This implementation operates on one-dimensional data.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>SubSample class</u></h3>

.. _SubSample:

The ``SubSample`` class reduces the number of samples in a dataset.
This is useful for downsampling large arrays prior to plotting, storage, or comparison with lower-resolution datasets.

.. code-block:: python

   class SubSample:
     def __init__(self, method="average", stride=None, numSamples=None):
       self.method = method
       self.stride = stride
       self.numSamples = numSamples

     def postprocess(self, x):
       if len(x) < self.numSamples:
         return x

       numX = len(x)
       if self.stride == None:
         self.stride = int(np.floor( numX / self.numSamples ))

       if self.method == "average":
         return np.average(np.pad(x, (0, self.stride - numX % self.stride), mode='constant', constant_values=x[-1]).reshape(-1, int(self.stride)), axis=1)

       if self.method == "nearest":
         return x[np.round(np.linspace(0, numX-1, num=self.numSamples)).astype(int)]

The ``postprocess`` method supports two subsampling methods:

- ``"average"``: averages values over blocks of size ``stride``
- ``"nearest"``: selects approximately evenly spaced existing entries

If ``stride`` is not provided, it is computed automatically from the total number of points and the requested number of samples.

Parameters
^^^^^^^^^^

#. ``method``:  Method used for subsampling. Supported options are ``"average"`` and ``"nearest"``.
#. ``stride``:  Number of entries per subsampling block. If not provided, it is computed automatically.
#. ``numSamples``:  Target number of output samples.

Returns
^^^^^^^

The ``postprocess`` method returns a subsampled array.

Notes
^^^^^

- If the input dataset already contains fewer entries than ``numSamples``, it is returned unchanged.
- For the ``"average"`` method, the array is padded at the end as needed using the final value.
- For the ``"nearest"`` method, output values are selected directly from the original array.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>compute_domain_strain function</u></h3>

.. _compute_domain_strain:

The ``compute_domain_strain`` function computes normal domain strains from the diagonal components of the deformation gradient tensor.
It can return either engineering strain or logarithmic strain, along with the Jacobian determinant.

.. code-block:: python

   def compute_domain_strain(F00, F11, F22, engineeringStrain=False):

If ``engineeringStrain`` is set to ``True``, the function computes:

.. math::

   e_{xx} = F_{00} - 1,\quad
   e_{yy} = F_{11} - 1,\quad
   e_{zz} = F_{22} - 1

Otherwise, it computes logarithmic strain:

.. math::

   e_{xx} = \log(F_{00}),\quad
   e_{yy} = \log(F_{11}),\quad
   e_{zz} = \log(F_{22})

The Jacobian is computed as:

.. math::

   J = F_{00} F_{11} F_{22}

Parameters
^^^^^^^^^^

``F00``
  Deformation gradient component in the x direction.

``F11``
  Deformation gradient component in the y direction.

``F22``
  Deformation gradient component in the z direction.

``engineeringStrain``
  Boolean flag indicating whether engineering strain should be used. If ``False``, logarithmic strain is computed.

Returns
^^^^^^^

The function returns:
``exx``, ``eyy``, ``ezz``, ``J``

Notes
^^^^^

- The returned quantities are stored as ``DataObj`` instances.
- ``J`` represents the ratio of deformed volume to reference volume.
- Logarithmic strain is often preferred for large-deformation analysis.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>compute_domain_stress function</u></h3>

.. _compute_domain_stress:

The ``compute_domain_stress`` function computes normal reaction-based stresses from the initial cross-sectional areas, deformation gradient components, and boundary reaction forces.
It can return either engineering stress or true stress depending on the value of ``engineeringStress``.

.. code-block:: python

   def compute_domain_stress(Axx0, Ayy0, Azz0, F00, F11, F22, Rxm, Rxp, Rym, Ryp, Rzm, Rzp, engineeringStress=False):

The function begins with the reference cross-sectional areas:

- ``Ax = Axx0``
- ``Ay = Ayy0``
- ``Az = Azz0``

If ``engineeringStress`` is ``False``, the current areas are computed using the deformation gradient components:

.. math::

   A_x = A_{xx0} F_{11} F_{22}

.. math::

   A_y = A_{yy0} F_{00} F_{22}

.. math::

   A_z = A_{zz0} F_{00} F_{11}

The reaction-based stresses are then computed from the boundary reactions and corresponding areas.

Parameters
^^^^^^^^^^

#. ``Axx0``:  Reference cross-sectional area normal to the x direction.
#. ``Ayy0``:  Reference cross-sectional area normal to the y direction.
#. ``Azz0``:  Reference cross-sectional area normal to the z direction.
#. ``F00``:  Deformation gradient component in the x direction.
#. ``F11``:  Deformation gradient component in the y direction.
#. ``F22``:  Deformation gradient component in the z direction.
#. ``Rxm``:  Reaction force at the x- boundary.
#. ``Rxp``:  Reaction force at the x+ boundary.
#. ``Rym``:  Reaction force at the y- boundary.
#. ``Ryp``:  Reaction force at the y+ boundary.
#. ``Rzm``:  Reaction force at the z- boundary.
#. ``Rzp``:  Reaction force at the z+ boundary.
#. ``engineeringStress``: Boolean flag indicating whether reference-area stress should be used. If ``False``, the stress is based on the current deformed area.

Returns
^^^^^^^

The function returns:
``Rsxm``, ``Rsxp``, ``Rsxx``, ``Rsym``, ``Rsyp``, ``Rsyy``, ``Rszm``, ``Rszp``, ``Rszz``

Notes
^^^^^

- The returned quantities are stored as ``DataObj`` instances.
- The averaged normal stresses are computed as half the difference between opposite-face reactions divided by the corresponding area.
- When ``engineeringStress`` is ``False``, the function computes true stress using updated cross-sectional areas.
- The current implementation should be reviewed carefully for axis consistency in the stress calculations and returned object names.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>compute_pressure function</u></h3>

.. _compute_pressure:

The ``compute_pressure`` function computes pressure from the three normal stress components.

.. code-block:: python

   def compute_pressure(bsxx, bsyy, bszz):

The pressure is computed as:

.. math::

   p = -\frac{1}{3}(b\sigma_{xx} + b\sigma_{yy} + b\sigma_{zz})

where compressive stress corresponds to positive pressure.

Parameters
^^^^^^^^^^

#. ``bsxx``:  Normal stress component in the x direction.
#. ``bsyy``:  Normal stress component in the y direction.
#. ``bszz``:  Normal stress component in the z direction.

Returns
^^^^^^^

The function returns:
``Pressure``

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


Data Visualization and Plotting
===============================

The ``pfw_analysis`` script includes helper utilities for improving the visualization of simulation results.
These routines support plotting workflows by modifying color properties to enhance clarity and distinguishability between multiple datasets.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>lighten_color function</u></h3>

.. _lighten_color:


The ``lighten_color`` function modifies a given color to produce a lighter version.
This is useful when plotting multiple curves or datasets where variations of a base color are needed for visual distinction.


.. code-block:: python

   def lighten_color(color, amount=0.5):


This function converts the input color to HLS (Hue, Lightness, Saturation) space, adjusts the lightness component, and converts it back to RGB format.

The lightening operation is performed by scaling the distance from maximum lightness:

.. math::

   L_{\text{new}} = 1 - \text{amount} \cdot (1 - L)

where ``L`` is the original lightness value.


Parameters
^^^^^^^^^^

``color``
  Input color. This can be provided as:

  - a matplotlib color string (e.g., ``'g'``)
  - a hexadecimal color string (e.g., ``'#F034A3'``)
  - an RGB tuple (e.g., ``(0.3, 0.55, 0.1)``)

``amount``
  Float controlling the degree of lightening.

  - ``amount = 0`` → no change  
  - ``amount = 1`` → maximum lightening  
  - typical values range between ``0.3`` and ``0.7``  


Returns
^^^^^^^

The function returns an RGB tuple corresponding to the lightened color.


Notes
^^^^^

- The function uses ``matplotlib.colors`` and ``colorsys`` for color conversion.
- If the input color is a named matplotlib color, it is first converted to its corresponding RGB value.
- This function is commonly used to generate visually distinct variations of a base color when plotting multiple datasets.
- The hue and saturation are preserved while only the lightness component is modified.


Examples
^^^^^^^^

.. code-block:: python

   lighten_color('g', 0.3)
   lighten_color('#F034A3', 0.6)
   lighten_color((0.3, 0.55, 0.1), 0.5)


:ref:`Back to top <pfw_analysis>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

 
Supporting Classes and Script Execution
==========================

The ``pfw_analysis`` script also contains several supporting classes and execution utilities that organize analysis settings, store simulation metadata, register data fields, and manage command-line execution.
These components provide the structural framework that allows the import, post-processing, and plotting routines to operate on a simulation job in a consistent way.


=====

.. raw:: html

   <h3><u>AnalysisOptions class</u></h3>

.. _AnalysisOptions:

The ``AnalysisOptions`` class stores user-configurable settings that control how analysis results are read, filtered, written, and plotted.

.. code-block:: python

   class AnalysisOptions:
       def __init__(self):

This class defines a collection of Boolean flags and scalar parameters used to configure post-processing behavior, including whether reaction forces are plotted, whether box-average data are read, whether filtered stresses are used, and whether plots are displayed interactively or saved to file.

Parameters
^^^^^^^^^^

This class does not take any input arguments beyond ``self``.
All options are initialized internally with default values.

Attributes
^^^^^^^^^^

``plotReactions``
  If ``True``, reaction results are plotted. If ``False``, results are only written to the console.

``readBoxSums``
  If ``True``, the ``boxAverageHistory.csv`` file is read.

``livePlot``
  If ``True``, the plotting window remains open for interactive viewing.

``writePostprocessFile``
  If ``True``, post-processed data are written to file.

``useEngineeringStrain``
  If ``True``, engineering strain is used instead of logarithmic strain.

``flipCompression``
  If ``True``, stress is flipped in sign for compression simulations.

``filterStresses``
  If ``True``, stress spikes above a threshold are filtered from the reaction data.

``stressMaxThreshold``
  Maximum allowed stress threshold.

``stressMinThreshold``
  Minimum allowed stress threshold.

``enforceAxesLimits``
  If ``True``, manually specified axis limits are used in plots.

``yStressMax``
  Upper limit for the stress axis in plots.

``yStressMin``
  Lower limit for the stress axis in plots.

``filterByMedian``
  If ``True``, median filtering is applied to the data.

``windowSize``
  Window size used for the median filter.

Returns
^^^^^^^

This class does not return a value.
It creates an object containing analysis configuration parameters.

Notes
^^^^^

- This class serves as a centralized container for analysis settings.
- The default values are selected for typical reaction-history and box-average post-processing workflows.
- Threshold and filtering options may need to be adjusted for different material models or loading conditions.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>MPMJob class</u></h3>

.. _MPMJob:

The ``MPMJob`` class represents a single MPM simulation job directory and provides methods for gathering metadata, reading output files, registering fields, and computing derived quantities.

.. code-block:: python

   class MPMJob:
     def __init__(self, job_dir_path):

When initialized, this class stores the job directory path, constructs the expected input filename, initializes a field dictionary, and gathers metadata from the job input file and output files.

Parameters
^^^^^^^^^^

``job_dir_path``
  Path to the simulation job directory.

Attributes
^^^^^^^^^^

``job_dir_path``
  Path to the simulation directory.

``job_name``
  Name of the simulation job, inferred from the directory name.

``job_input_file``
  Expected Particle File Writer input filename of the form ``pfw_input_<job_name>.py``.

``fields``
  Dictionary of registered ``DataObj`` fields associated with the simulation.

Methods
^^^^^^^

``gather_job_metadata``
  Reads job metadata from the input file and available output files.

``read_reaction_file``
  Reads the ``reactionHistory.csv`` file and registers its fields.

``read_from_box_average_file``
  Reads the ``boxAverageHistory.csv`` file and registers its fields.

``registerField``
  Adds a ``DataObj`` instance to the internal field dictionary.

``applyPostProcess``
  Applies a post-processing filter to one field or to all registered fields.

``compute_domain_strain``
  Computes derived strain fields and the Jacobian from the deformation gradient.

``compute_domain_stress``
  Computes reaction-based stresses from boundary reaction forces.

Returns
^^^^^^^

This class does not return a value.
It creates an object representing a simulation job and its associated metadata and fields.

Notes
^^^^^

- The class dynamically imports the job input file in order to access variables defined in the corresponding ``pfw_input`` script.
- Initial domain dimensions, sample dimensions, grid resolution, particles-per-cell, and cross-sectional areas are inferred during metadata gathering.
- The class checks whether reaction and box-average history files are present before attempting to read them.
- Metadata are printed to the console when the object is initialized.
- Some stress-component calculations in ``compute_domain_stress`` should be reviewed carefully for axis consistency.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>gather_job_metadata method</u></h3>

The ``gather_job_metadata`` method reads simulation metadata from the Particle File Writer input file and available job outputs.

.. code-block:: python

   def gather_job_metadata(self):

This method verifies that the expected input file exists, imports it as a Python module, reads the ``pfw`` dictionary, and extracts geometric, discretization, and material information from the job definition.

The method registers metadata such as:

- density
- presence of reaction and box-average files
- domain dimensions
- sample dimensions
- periodicity
- grid resolution
- particles per cell
- cell spacing
- initial cross-sectional areas
- initial sample volume
- particle count
- plane strain setting

Notes
^^^^^

- This method is called automatically during initialization of the ``MPMJob`` object.
- Default values are used when optional job attributes are not explicitly defined in the imported input file.
- The number of particles is obtained by counting the lines in the particle file.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>read_reaction_file method</u></h3>

The ``read_reaction_file`` method reads the reaction-history CSV file associated with the job and registers each column as a ``DataObj``.

.. code-block:: python

   def read_reaction_file(self):

This method reads ``reactionHistory.csv`` from the job directory and stores the resulting quantities in the ``fields`` dictionary.

Returns
^^^^^^^

This method does not return a value.
It registers the reaction-history fields in the ``MPMJob`` object.

Notes
^^^^^

- The reaction file must exist, or the method raises an assertion error.
- The registered field names include ``Time``, ``F00``, ``F11``, ``F22``, ``Lx``, ``Ly``, ``Lz``, ``Rxm``, ``Rxp``, ``Rym``, ``Ryp``, ``Rzm``, ``Rzp``, ``L00``, ``L11``, and ``L22``.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>read_from_box_average_file method</u></h3>

The ``read_from_box_average_file`` method reads the box-average history CSV file associated with the job and registers each column as a ``DataObj``.

.. code-block:: python

   def read_from_box_average_file(self):

This method reads ``boxAverageHistory.csv`` from the job directory and stores the resulting quantities in the ``fields`` dictionary.

Returns
^^^^^^^

This method does not return a value.
It registers the box-average fields in the ``MPMJob`` object.

Notes
^^^^^

- The box-average file must exist, or the method raises an assertion error.
- Registered field names are prefixed with ``B`` to distinguish them from reaction-history fields, for example ``BSxx``, ``BDensity``, and ``BIE``.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>registerField method</u></h3>

The ``registerField`` method adds a field to the internal field dictionary.

.. code-block:: python

   def registerField(self, field):

Parameters
^^^^^^^^^^

``field``
  A ``DataObj`` instance to be stored in the job's field dictionary.

Returns
^^^^^^^

This method does not return a value.

Notes
^^^^^

- The field is stored using ``field.name`` as the dictionary key.
- Registering a field with an existing name will overwrite the previous entry.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>applyPostProcess method</u></h3>

The ``applyPostProcess`` method applies a post-processing filter to one registered field or to all fields.

.. code-block:: python

   def applyPostProcess(self, fieldname, filter):

Parameters
^^^^^^^^^^

``fieldname``
  Name of the field to be post-processed, or ``"all"`` to apply the filter to every registered field.

``filter``
  An object with a ``postprocess`` method, such as ``Trim``, ``MedianFilter``, ``SubSample``, or ``RemoveNonMonotonicEntries``.

Returns
^^^^^^^

This method does not return a value.

Notes
^^^^^

- If ``fieldname`` is ``"all"``, the filter is applied to all entries in the ``fields`` dictionary.
- Otherwise, the filter is applied only to the selected field.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>compute_domain_strain method</u></h3>

The ``compute_domain_strain`` method computes domain strain fields and the Jacobian from the stored deformation-gradient components.

.. code-block:: python

   def compute_domain_strain(self, engineeringStrain=False):

This method uses the registered ``F00``, ``F11``, and ``F22`` fields to compute ``exx``, ``eyy``, ``ezz``, and ``J``, then stores them as new ``DataObj`` fields.

Parameters
^^^^^^^^^^

``engineeringStrain``
  If ``True``, engineering strain is computed. Otherwise, logarithmic strain is computed.

Returns
^^^^^^^

This method does not return a value.
It registers the derived strain quantities in the job's field dictionary.

Notes
^^^^^

- The computation follows the same logic as the standalone ``compute_domain_strain`` function.
- Derived fields are stored under the names ``exx``, ``eyy``, ``ezz``, and ``J``.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>compute_domain_stress method</u></h3>

The ``compute_domain_stress`` method computes reaction-based stresses from stored reaction-force and deformation-gradient fields.

.. code-block:: python

   def compute_domain_stress(self, engineeringStress=False):

This method uses the registered reaction fields and the initial cross-sectional areas stored in the ``MPMJob`` object to compute reaction stresses and store them as new ``DataObj`` fields.

Parameters
^^^^^^^^^^

``engineeringStress``
  If ``True``, stresses are computed using the reference areas. Otherwise, updated deformed areas are used.

Returns
^^^^^^^

This method does not return a value.
It registers the derived stress quantities in the job's field dictionary.

Notes
^^^^^

- The computation follows the same logic as the standalone ``compute_domain_stress`` function.
- Derived fields are stored under names such as ``Rsxx``, ``Rsyy``, and ``Rszz``.
- The current implementation should be reviewed carefully for axis consistency in some stress-component calculations.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>DataObj class</u></h3>

.. _DataObj:

The ``DataObj`` class stores a named dataset together with formatting information and any processed version of that dataset.

.. code-block:: python

   class DataObj:
       def __init__(self, name, data, format=".6f"):

This class is used throughout ``pfw_analysis`` to represent imported or derived scalar fields in a consistent way.

Parameters
^^^^^^^^^^

``name``
  Name of the field.

``data``
  Underlying dataset, typically a NumPy array.

``format``
  String format specifier used when writing the data to file or console output.

Attributes
^^^^^^^^^^

``name``
  Name of the field.

``data``
  Original unprocessed data.

``format``
  String format specifier used for output formatting.

``processedData``
  Current processed version of the data. Initially this is identical to ``data``.

Methods
^^^^^^^

``applyPostProcess``
  Applies a post-processing filter and updates ``processedData``.

``getData``
  Returns the current processed data.

Returns
^^^^^^^

This class does not return a value.
It creates a field object used throughout the analysis workflow.

Notes
^^^^^

- The original raw data are preserved in ``data``.
- Post-processing operations modify only ``processedData``.
- This class provides a simple interface for storing both imported and derived quantities.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>applyPostProcess method for DataObj</u></h3>

The ``applyPostProcess`` method updates the processed version of the dataset using a filter object.

.. code-block:: python

   def applyPostProcess(self, filter):

Parameters
^^^^^^^^^^

``filter``
  An object with a ``postprocess`` method that accepts the current processed data and returns the modified result.

Returns
^^^^^^^

This method does not return a value.

Notes
^^^^^

- The processed data are updated in place.
- This design allows multiple post-processing operations to be chained sequentially.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>getData method</u></h3>

The ``getData`` method returns the current processed version of the dataset.

.. code-block:: python

   def getData(self):

Returns
^^^^^^^

The current processed data array.

Notes
^^^^^

- If no post-processing has been applied, this method returns the original data.

:ref:`Back to top <pfw_analysis>`

=====

.. raw:: html

   <h3><u>Main execution block</u></h3>

.. _main_execution_block:

The main execution block allows ``pfw_analysis.py`` to be run directly from the command line for job analysis and plotting.

.. code-block:: python

   if __name__ == "__main__":

This block uses ``argparse`` to define command-line options, validates the provided job directory, creates an ``MPMJob`` object, and performs selected analysis actions such as reading reaction data and plotting results.

Command-line arguments
^^^^^^^^^^^^^^^^^^^^^^

``jobdir``
  Location of the job to analyze.

``-i``, ``--interactive``
  Display plots interactively.

``-c``, ``--console``
  Write output to the console.

``-x``, ``--xyz``
  List of spatial directions to plot.

``-e``, ``--export``
  Write output to CSV.

``-s``, ``--save``
  Save plot output to a PNG file.

``-p``, ``--plot``
  Select which job data to plot.

``-f``, ``--fields``
  List of field names to output.

``-o``, ``--output``
  Base name of output files.

Returns
^^^^^^^

This block does not return a value.
It executes the requested analysis workflow when the script is run directly.

Notes
^^^^^

- The job directory must exist or the script aborts with an assertion error.
- The current implementation includes a simple plotting branch for reaction data.
- This block serves as the command-line entry point for the script.