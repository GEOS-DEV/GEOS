.. _AdvancedExampleFreeFreeBeamModes:


#####################################################################
Free-Free Beam: Modal Analysis versus Euler-Bernoulli Beam Theory
#####################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example computes the vibration modes of a free-free steel beam with the ``Modal`` time integration option of the :ref:`SolidMechanicsLagrangianFEM` solver, and compares them with the Euler-Bernoulli beam theory.
The beam is 10 m long and has a square cross-section of 0.2 m by 0.2 m.
No displacement boundary condition is applied, so the beam is a free body.
Its stiffness matrix :math:`\mathbf{K}` is singular, and the first six modes are the rigid-body modes at zero frequency.
They are followed by the elastic modes of the beam: bending, torsion and axial modes.

The beam theory gives the frequencies of the bending modes as

.. math::
   f_n = \frac{(\beta_n L)^2}{2 \pi L^2} \sqrt{\frac{E I}{\rho A}},

where :math:`L` is the length, :math:`E` the Young modulus, :math:`\rho` the density, :math:`A` the area of the cross-section and :math:`I` its second moment of area.
The numbers :math:`\beta_n L = 4.7300, 7.8532, 10.9956, \dots` are the roots of :math:`\cos(x)\cosh(x) = 1`.
The cross-section is square, so every bending frequency is double: the beam bends in two perpendicular planes with the same frequency.
The mode shape is

.. math::
   \phi_n(x) = \cosh(\beta_n x) + \cos(\beta_n x) - \sigma_n \left( \sinh(\beta_n x) + \sin(\beta_n x) \right),
   \qquad \sigma_n = \frac{\cosh(\beta_n L) - \cos(\beta_n L)}{\sinh(\beta_n L) - \sin(\beta_n L)}.

GEOS normalizes its modes with :math:`\boldsymbol{\phi}^T \mathbf{M} \boldsymbol{\phi} = 1`.
For the beam, this means :math:`\rho A \int_0^L \phi^2 \, dx = 1`.
The first axial frequency is :math:`f = \frac{1}{2L} \sqrt{E/\rho}`.

**Input files**

This example uses three GEOS xml files:

.. code-block:: console

  inputFiles/solidMechanics/modalFreeFreeBeam_base.xml
  inputFiles/solidMechanics/modalFreeFreeBeam_arnoldi_smoke.xml
  inputFiles/solidMechanics/modalFreeFreeBeam_block_smoke.xml

The base file has the mesh, the material and the outputs.
The two other files differ by the eigensolver settings.
A Python script that gives the beam theory solutions and plots the comparison is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/solidMechanics/FreeFreeBeamModes/FreeFreeBeamModes_vs_EulerBernoulli.py

------------------------------------------------------------------
Mesh
------------------------------------------------------------------

The internal mesh generator builds the beam with 100 by 4 by 4 trilinear hexahedra, that is 7 575 unknowns:

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/modalFreeFreeBeam_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MODAL_BEAM_MESH -->
    :end-before: <!-- SPHINX_MODAL_BEAM_MESH_END -->

------------------------------------------------------------------
Material
------------------------------------------------------------------

The material is an isotropic elastic steel, in the International System of Units:

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/modalFreeFreeBeam_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MODAL_BEAM_MATERIAL -->
    :end-before: <!-- SPHINX_MODAL_BEAM_MATERIAL_END -->

------------------------------------------------------------------
Solver
------------------------------------------------------------------

The solver block selects the ``Modal`` option:

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/modalFreeFreeBeam_arnoldi_smoke.xml
    :language: xml
    :start-after: <!-- SPHINX_MODAL_BEAM_SOLVER -->
    :end-before: <!-- SPHINX_MODAL_BEAM_SOLVER_END -->

The solver computes ``modalNumModes`` modes, the six rigid-body modes included.
The shift ``modalShiftFrequency`` is a negative frequency.
The shifted matrix :math:`\mathbf{K} + \alpha \mathbf{M}`, with :math:`\alpha = (2 \pi f)^2`, is then positive definite although the beam is a free body.
The ``arnoldi`` eigensolver solves one linear system with this matrix for each Krylov vector.
The linear solver tolerance ``krylovTol`` is two orders of magnitude smaller than ``modalTolerance``.
The default ``modalCompletenessCheck`` verifies that no copy of a repeated eigenvalue is missed.
The second input file uses ``modalBlockSize="6"``, which is the multiplicity of the rigid-body modes, instead of this check.
The attribute ``modalDeflateRigidBodyModes="1"`` removes the rigid-body modes from the eigensolve, which lowers the cost of a free body.

The solver runs once.
It prints a table with the frequency, the residual and the participation factors of each mode in the log.
It saves the mode shapes in the nodal fields ``modeShape1``, ``modeShape2``, and so on.
A ``VTK`` output writes them.

------------------------------------------------------------------
Collecting the mode shapes
------------------------------------------------------------------

The integrated test collects the first and the second bending mode along the axis of the beam with ``PackCollection`` tasks.
The modes 7 and 8 are the two polarizations of the first bending mode.
The modes 9 and 10 are the two polarizations of the second one.

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/modalFreeFreeBeam_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MODAL_BEAM_TASKS -->
    :end-before: <!-- SPHINX_MODAL_BEAM_TASKS_END -->

------------------------------------------------------------------
A comparison between GEOS results and beam theory
------------------------------------------------------------------

The figure uses a run with ``modalNumModes="20"``.
The results are in the files ``FreeFreeBeamFrequencies.txt`` and ``FreeFreeBeamModeShapes.txt``.
The left plot compares the bending frequencies.
Each GEOS value is the mean of a pair of modes.
The right plot compares the mass-normalized mode shapes.
The polarization of a bending mode in the plane of the cross-section is arbitrary.
The script projects each computed mode on the beam theory shape.

.. plot:: docs/sphinx/advancedExamples/validationStudies/solidMechanics/FreeFreeBeamModes/FreeFreeBeamModes_vs_EulerBernoulli.py

The shapes agree with the beam theory.
The GEOS frequencies are 5 % higher for the first bending mode.
This is the effect of the shear locking of the standard trilinear hexahedron, which has only four elements across the thickness.
The difference falls with the mesh: it is 13 %, 7.5 % and 5 % with 60, 80 and 100 elements along the beam.
The ratio of the frequencies of the second and the first bending mode is 2.749 in GEOS and 2.756 in the beam theory.
The axial mode at 252.4 Hz agrees with the closed form to :math:`10^{-4}`.
The first six frequencies are zero up to rounding errors (below 0.001 Hz).

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Eigensolvers**

The ``lobpcg`` eigensolver applies only a preconditioner and finds the lowest modes.
It is much faster for a free body, but it needs ``modalDeflateRigidBodyModes="1"`` and a good preconditioner.
A single multigrid cycle is not enough for this slender beam.
See :ref:`SolidMechanicsLagrangianFEM` for the description of both eigensolvers.

**Integrated test**

The two decks run as integrated tests on one, two and four MPI ranks.
They compare the mode shapes with the beam theory with the curve checker.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
