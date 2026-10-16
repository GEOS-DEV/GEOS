.. _ContributingInputFilesDoc:

************************
Contributing Input Files
************************
As part of the development cycle, functional input files should be introduced 
into the repository in order to provide 1) testing, 2) examples of proper use of
the code.
The input files should be modularized such that the majority of the input
resides in a base file, and the variations in input are contained in specific
files for integrated/smoke testing, and benchmarking. 
For example if we had a single input file named ``myProblem.xml``, we would break 
it up into ``myProblem_base.xml``, ``myProblem_smoke.xml``, and 
``myProblem_benchmark.xml``.
Each of the ``smoke/benchmark`` files should include the ``base`` file using an
include block as follows:

.. code-block:: xml

     <Included>
       <File name="./myProblem_base.xml"/>
     </Included>
     
The files should be placed in the appropriate application specific subdirectory
under the ``GEOS/inputFiles`` directory. 
For example, the ``beamBending`` problem input files reside in the 
``inputFiles/solidMechanics`` directory. 
To turn the new files into integrated tests, register the ``smoke`` deck in the ``.ats``
file that sits in the same directory (for example ``inputFiles/solidMechanics/beamBending.ats``),
or create a new ``.ats`` file there if the directory has none.
The ``.ats`` files are found automatically by ``inputFiles/main.ats``.
See :ref:`IntegratedTests` for the syntax of these files and for how to generate the corresponding baselines.
