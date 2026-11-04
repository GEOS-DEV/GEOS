.. _AddingNewSolver:

Adding a new Physics Solver
###########################

In this tutorial, you will learn how to construct a new GEOS Physics Solver class.
We will use *LaplaceFEM* solver, computing the solution of the Laplace problem in
a specified material, as a starting point.

.. math::
   :nowrap:

   \begin{eqnarray*}
   D^* \Delta X = f \quad \mbox{in\;} \Omega \\
   X = X^g \quad \mbox{on\;} \Gamma_g
   \end{eqnarray*}

It is advised to read :ref:`XML_and_classes` preliminary to this tutorial.
The goal of this document is to explain how to develop a new solver that solves
Laplace's equation with a constant diffusion coefficient that is specified by users in the XML input.

For readability, member functions in the text will be referenced by their names but their
arguments will be omitted.

*LaplaceFEM* overview
=====================
The *LaplaceFEM* solver can be found in ``./src/coreComponents/physicsSolvers/simplePDE/``.
Let us inspect declarations in ``LaplaceFEM.hpp`` and implementations in ``LaplaceFEM.cpp``
before diving into specifying a new solver class that meets our needs.

Declaration file (reference)
----------------------------
The included header is ``physicsSolvers/simplePDE/LaplaceBaseH1.hpp`` which declares the base class ``LaplaceBaseH1``, shared by all Laplace solvers. Moreover, ``physicsSolvers/simplePDE/LaplaceBaseH1.hpp`` includes the following headers:

 - ``common/format/EnumStrings.hpp`` which includes facilities for enum-string conversion (useful for reading enum values from input);
 - ``physicsSolvers/PhysicsSolverBase.hpp`` which declares the abstraction class shared by all physics solvers.
 - ``fieldSpecification/FieldSpecificationManager.hpp`` which declares a manager used to access and to set field on the discretized domain.

Let us jump forward to the class enum and variable as they contain the data used
specifically in the implementation of *LaplaceFEM*.

class enums and variables (reference)
`````````````````````````````````````
The class exhibits two (protected) member variables:

 - ``m_fieldName`` which stores the name of the diffused variable (*e.g.* the temperature) as a `string`;
 - ``m_timeIntegrationOption`` an `enum` value allowing to dispatch with respect to the transient treatment.

``TimeIntegrationOption`` is an `enum` specifying the transient treatment which can be chosen
respectively between *SteadyState* and *ImplicitTransient* depending on whether we are interested in
the transient state.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_TIMEINTOPT
   :end-before: //END_SPHINX_INCLUDE_TIMEINTOPT

In order to register an enumeration type with the Data Repository and have its value read from input,
we must define stream insertion/extraction operators. This is a common task, so GEOS provides
a facility for automating it. Upon including ``common/format/EnumStrings.hpp``, we can call the following macro
at the namespace scope (in this case, right after the ``LaplaceBaseH1`` class definition is complete):

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGENUM
   :end-before: //END_SPHINX_INCLUDE_REGENUM

Once explained the main variables and enum, let us start reading through the different member functions:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_BEGINCLASS
   :end-before: //END_SPHINX_INCLUDE_BEGINCLASS

Start looking at the class *LaplaceFEM* constructor and destructor declarations
shows the usual `string` ``name`` and `Group*` pointer to ``parent`` that are required
to build the global file-system like structure of GEOS (see :ref:`GroupPar` for details).
It can also be noted that the nullary constructor is deleted on purpose to avoid compiler
automatic generation and user misuse.

The next method ``catalogName()`` is static and returns the key to be added to the *Catalog* for this type of solver
(see :ref:`ObjectCatalogPar` for details). It has to be paired with the following macro in the implementation file.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTER
   :end-before: //END_SPHINX_INCLUDE_REGISTER

Finally, the member function ``registerDataOnMesh()`` is declared in the ``LaplaceBaseH1`` class as

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTERDATAONMESH
   :end-before: //END_SPHINX_INCLUDE_REGISTERDATAONMESH

It is used to assign fields onto the discretized mesh object and
will be further discussed in the :ref:`Implementation` section.

The next block consists in solver interface functions. These member functions set up
and specialize every time step from the system matrix assembly to the solver stage.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_SOLVERINTERFACE
   :end-before: //END_SPHINX_INCLUDE_SOLVERINTERFACE

Furthermore, the following functions are inherited from the base class.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_SOLVERINTERFACE
   :end-before: //END_SPHINX_INCLUDE_SOLVERINTERFACE

Eventually, ``applyDirichletBCImplicit()`` is the working specialized member functions called
when ``applyBoundaryConditions()`` is called in this particular class override.

Browsing the base class ``PhysicsSolverBase``, it can be noted that most of the solver interface functions are called during
either ``PhysicsSolverBase::linearImplicitStep()`` or ``PhysicsSolverBase::nonlinearImplicitStep()`` depending on the solver strategy chosen.

Another central member function is ``postInputInitialization()``, inherited from ``Group`` and overridden in ``PhysicsSolverBase``.
It is called once the input has been read from the XML file, and is the place where a solver can check or post-process its input values.
``LaplaceFEM`` does not need to override it: as shown in the :ref:`Implementation` section, its input values are registered
directly on its member variables, which are therefore filled automatically when the XML file is parsed.

Let us focus on a ``struct`` that plays an important role: the *viewKeyStruct* structure.

*viewKeyStruct* structure (reference)
`````````````````````````````````````

This embedded instantiated structure is a common pattern shared by all solvers.
It stores ``dataRepository::ViewKey`` type objects that are used as binding data
between the input XML file and the source code.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_VIEWKEY
   :end-before: //END_SPHINX_INCLUDE_VIEWKEY

We can check that in the *LaplaceFEM* companion integratedTest

.. literalinclude:: ../../../../../inputFiles/simplePDE/Laplace_base.xml
   :language: xml
   :start-after: <Solvers>
   :end-before: </Solvers>

In the following section, we will see where this binding takes place.

.. _Implementation:

Implementation File (reference)
-------------------------------
Switching to implementation, we will focus on few implementations, leaving details
to other tutorials. The ``LaplaceFEM`` constructor is implemented as follows.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_CONSTRUCTOR
   :end-before: //END_SPHINX_INCLUDE_CONSTRUCTOR

As we see, it calls the ``LaplaceBaseH1`` constructor, that is implemented as follows.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_CONSTRUCTOR
   :end-before: //END_SPHINX_INCLUDE_CONSTRUCTOR

Checking out the constructor, we can see that the use of a ``registerWrapper<T>(...)``
allows us to register the key value from the `enum` ``viewKeyStruct`` defining them as:

 - ``InputFlags::OPTIONAL`` if they are optional and can be provided;
 - ``InputFlags::REQUIRED`` if they are required and will throw error if not;

and their associated descriptions for auto-generated docs.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTERDATAONMESH
   :end-before: //END_SPHINX_INCLUDE_REGISTERDATAONMESH

``registerDataOnMesh()`` is browsing all subgroups in the mesh ``Group`` object and
for all nodes in the sub group:

 - register the observed field under the chosen ``m_fieldName`` key;
 - apply a default value;
 - set the output verbosity level (here ``PlotLevel::LEVEL_0``);
 - set the field associated description for auto generated docs.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_ASSEMBLY
   :end-before: //END_SPHINX_INCLUDE_ASSEMBLY

``assembleSystem()`` will be our core focus as we want to change the diffusion coefficient from its
implicit unit value to a user-defined value read from the XML file. This method is in charge of constructing
the FEM system matrix in parallel. For each mesh level targeted by the solver, it retrieves the degree-of-freedom numbers
of the nodes and builds a ``LaplaceFEMKernelFactory``. The generic function ``finiteElement::regionBasedKernelApplication``
then launches the kernel on all the target cell sub-regions, for the finite element space of the solver discretization.
The element-level computation is implemented in the kernel class ``LaplaceFEMKernel`` (file ``LaplaceFEMKernels.hpp``), whose
``quadraturePointKernel()`` adds the contribution of each quadrature point to the element matrix:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEMKernels.hpp
   :language: c++
   :start-at: void quadraturePointKernel(
   :end-before: /**

As we spotted the place to change in a code to get a user-defined diffusion coefficient into the game, let us jump
to writing our new *LaplaceDiffFEM* solver.

Start doing your own Physics solver
===================================
As we will extend *LaplaceFEM* capabilities, we will derive publicly from it, and write a kernel
deriving from ``LaplaceFEMKernel``.

.. note::

  The code snippets below illustrate the procedure; they are not compiled as part of GEOS.
  Always refer to the current sources of ``LaplaceFEM`` for the exact signatures.

Declaration File
----------------

As there is only few places where we have to change, the whole declaration file is reported below and
commented afterwards.

.. code-block:: c++

  #include "physicsSolvers/simplePDE/LaplaceFEM.hpp"

  namespace geos
  {

  class LaplaceDiffFEM : public LaplaceFEM
  {
  public:

    LaplaceDiffFEM() = delete;

    LaplaceDiffFEM( string const & name,
                    Group * const parent );

    static string catalogName() { return "LaplaceDiffFEM"; }

    string getCatalogName() const override { return catalogName(); }

    virtual void
    assembleSystem( real64 const time,
                    real64 const dt,
                    DomainPartition & domain,
                    DofManager const & dofManager,
                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                    arrayView1d< real64 > const & localRhs ) override;

    struct viewKeyStruct : public LaplaceFEM::viewKeyStruct
    {
      static constexpr char const * diffusionCoefficient() { return "diffusionCoefficient"; }
    };

  protected:

    real64 m_diffusion;

  };

  } /* namespace geos */

We intend to have a user-defined diffusion coefficient, we then need a ``real64`` class variable ``m_diffusion``
to store it.

Consistently with *LaplaceFEM*, we will also delete the nullary constructor and declare a constructor with the same arguments for
forwarding to the ``Group`` base class. Another mandatory step is to define the static ``catalogName()`` method (and the virtual
``getCatalogName()``), which gives the XML tag of the new solver and is used to register it in the catalog.

Then as mentioned in :ref:`Implementation`, the diffusion coefficient is used when assembling the matrix coefficient. Hence
we will have to override the ``assembleSystem()`` function as detailed below.

Moreover, if we want to introduce a new binding between the input XML and the code we will have to work on the
``struct viewKeyStruct`` and on the constructor.

Our new solver ``viewKeyStruct`` will have its own structure inheriting from the *LaplaceFEM* one to have the ``timeIntegrationOption``
and ``fieldName`` keys. It will also define a ``diffusionCoefficient`` key to be bound to the user defined homogeneous coefficient on one hand
and to our ``m_diffusion`` class variable on the other.


Implementation File
-------------------
As we have seen in :ref:`Implementation`, the first place where to implement a new register from XML input is
in the constructor. Registering the wrapper with a pointer to ``m_diffusion`` binds the XML attribute to the member variable:
its value is read and stored automatically when the input file is parsed, so that there is no need to override ``postInputInitialization()``.
The ``diffusionCoefficient`` entry will be asked as a required input. If not provided, the error thrown will ask for it.

.. code-block:: c++

  LaplaceDiffFEM::LaplaceDiffFEM( string const & name,
                                  Group * const parent ):
    LaplaceFEM( name, parent ),
    m_diffusion( 0.0 )
  {
    registerWrapper( viewKeyStruct::diffusionCoefficient(), &m_diffusion ).
      setInputFlag( dataRepository::InputFlags::REQUIRED ).
      setDescription( "Uniform diffusion coefficient for the Laplace equation" );
  }

  REGISTER_CATALOG_ENTRY( PhysicsSolverBase, LaplaceDiffFEM, string const &, Group * const )

Now that we have required, read and bound the user-defined diffusion value to a variable, we can use it in the construction of our
matrix. The element-level computation takes place in the kernel, so we first derive a new kernel from ``LaplaceFEMKernel``,
which receives the diffusion coefficient as an additional constructor argument and uses it in ``quadraturePointKernel()``.
Since the kernels are launched through templates on their exact type, defining ``quadraturePointKernel()`` in the derived kernel
is enough for it to be used instead of the base class one.

.. code-block:: c++
  :emphasize-lines: 24,44,49

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            typename FE_TYPE >
  class LaplaceDiffFEMKernel : public LaplaceFEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >
  {
  public:
    using Base = LaplaceFEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;
    using typename Base::StackVariables;
    using Base::maxNumTestSupportPointsPerElem;

    LaplaceDiffFEMKernel( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          localIndex const targetRegionIndex,
                          SUBREGION_TYPE const & elementSubRegion,
                          FE_TYPE const & finiteElementSpace,
                          CONSTITUTIVE_TYPE & inputConstitutiveType,
                          arrayView1d< globalIndex const > const inputDofNumber,
                          globalIndex const rankOffset,
                          CRSMatrixView< real64, globalIndex const > const inputMatrix,
                          arrayView1d< real64 > const inputRhs,
                          real64 const inputDt,
                          string const fieldName,
                          real64 const diffusion ):
      Base( nodeManager, edgeManager, faceManager, targetRegionIndex, elementSubRegion,
            finiteElementSpace, inputConstitutiveType, inputDofNumber, rankOffset,
            inputMatrix, inputRhs, inputDt, fieldName ),
      m_diffusion( diffusion )
    {}

    GEOS_HOST_DEVICE
    inline
    void quadraturePointKernel( localIndex const k,
                                localIndex const q,
                                StackVariables & stack ) const
    {
      GEOS_UNUSED_VAR( k );
      real64 dNdX[ maxNumTestSupportPointsPerElem ][ 3 ];
      real64 const detJ = FE_TYPE::calcGradN( q, stack.xLocal, stack.feStack, dNdX );
      for( localIndex a = 0; a < stack.numRows; ++a )
      {
        for( localIndex b = 0; b < stack.numCols; ++b )
        {
          stack.localJacobian[ a ][ b ] += m_diffusion * LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ;
        }
      }
    }

  protected:
    real64 const m_diffusion;
  };

  using LaplaceDiffFEMKernelFactory = finiteElement::KernelFactory< LaplaceDiffFEMKernel,
                                                                    arrayView1d< globalIndex const > const,
                                                                    globalIndex const,
                                                                    CRSMatrixView< real64, globalIndex const > const,
                                                                    arrayView1d< real64 > const,
                                                                    real64 const,
                                                                    string const,
                                                                    real64 const >;

The overridden ``assembleSystem()`` is then identical to the one of *LaplaceFEM*, except that it builds the new kernel factory
and passes ``m_diffusion`` to it. Note that ``m_fieldName`` is a protected member of ``LaplaceBaseH1``, so it is directly accessible
from the derived solver.

.. code-block:: c++
  :emphasize-lines: 17-18

  void LaplaceDiffFEM::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & regionNames )
    {
      NodeManager & nodeManager = mesh.getNodeManager();
      string const dofKey = dofManager.getKey( m_fieldName );
      arrayView1d< globalIndex const > const &
      dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofKey );

      LaplaceDiffFEMKernelFactory kernelFactory( dofIndex, dofManager.rankOffset(), localMatrix, localRhs, dt,
                                                 m_fieldName, m_diffusion );

      string const dummyString = "dummy";
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< >,
                                      constitutive::NullModel,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              dummyString,
                                                              kernelFactory );
    } );
  }

This completes the implementation of our new solver *LaplaceDiffFEM*.

Note: For consistency do not forget to change LaplaceFEM to LaplaceDiffFEM in the guards comments

Last steps
==========

After assembling both declarations and implementations for our new solver, the final steps go as:

 - add the new header files to the ``simplePDESolvers_headers`` list of ``physicsSolvers/simplePDE/CMakeLists.txt``;
 - add the new source file to the ``simplePDESolvers_sources`` list of the same file;
 - check that Doxygen comments are properly set in our solver class;
 - uncrustify it to match the code style by going to the build folder and running the command: ``make uncrustify_style``;
 - write unit tests for each new features in the solver class;
 - write an input file and register it as an integrated test for the solver class (see :ref:`ContributingInputFilesDoc` and :ref:`IntegratedTests`).
