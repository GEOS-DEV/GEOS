.. _kgdToughnessDominated:


#######################################################
Toughness dominated KGD problem
#######################################################

**Input file**

The xml input files for this test case are located at:

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Base.xml

and

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Example.xml

the corresponding integrated test with coarser mesh and smaller injection time is defined by the xml:

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Smoke.xml

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------



.. math::
   a = b

where
- :math:`E` is the Young's modulus

In this example, we focus our attention on the ``Solvers``, ``Mesh`` tags.

-----------------------------------------------------------
Solver
-----------------------------------------------------------

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- SPHINX_SOLVER -->
  :end-before: <!-- SPHINX_END -->



The parameters used in the simulation are summarized in the following table.

  +----------------+-----------------------+------------------+-------------------+
  | Symbol         | Parameter             | Units            | Value             |
  +================+=======================+==================+===================+
  | :math:`E`      | Young's modulus       | [Pa]             | 10\ :sup:`4`      |
  +----------------+-----------------------+------------------+-------------------+

---------------------------------
Inspecting results
---------------------------------

This plot .

.. figure:: mesh.png
   :align: center
   :width: 1000
   :figclass: align-center


.. figure:: propagation.gif
   :align: center
   :width: 1000
   :figclass: align-center



.. plot::

  

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the toughness dominated KGD example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
