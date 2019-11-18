.. _AddingNewSolver:

Adding a new Physics Solver
##################################################

In this tutorial, you will learn how to construct a new GEOSX Physics Solver class.
We will use *LaplaceFEM* solver, computing the solution of the Laplace problem in
a specified material, as a starting point.

.. math::
   :nowrap:

   \begin{eqnarray*}
   D^* \Delta X = f \quad \mbox{in\;} \Omega \\
   X = X^g \quad \mbox{on\;} \Gamma_g
   \end{eqnarray*}

It is advised to read :ref:`XML_and_classes` preliminary to this tutorial.
The point here is to develop a new solver which will read constant diffusion coeff from the XML input.

For readability, member functions in the text will be referenced by their names but their
arguments will be ommited.


*LaplaceFEM* overview
=============================
The *LaplaceFEM* solver can be found in ``./src/coreComponents/physicsSolvers/simplePDE/``.
Let us inspect declarations in ``LaplaceFEM.hpp`` and implementations in ``LaplaceFEM.cpp``
before diving into specifying a new solver class that meets our needs.

Declaration file
-----------------
The four included headers are:
 - ``physicsSolver/SolverBase.hpp`` which declares the abstraction class shared by all physics solvers
 - ``managers/FieldSpecification/FieldSpecificationManager.hpp`` which declares a manager used to access and to set field on the discretized domain.
 - ``linearAlgebra/DofManager.hpp`` which declares a manager used to handle degree of freedom indexation and manipulation
 - ``linearAlgebra/interfaces/InterfaceTypes.hpp`` which declares an interface to linear solvers and linear algebra librairies.

Right after that, a struct is designed to store the maximal stable timestep.
   .. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
      :language: c++
      :start-after: //START_SPHINX_INCLUDE_00
      :end-before: //END_SPHINX_INCLUDE_00

Some forward declarations are following :
 - ``FieldSpecificationBase`` is forward declared as we will use it to set the field on the boundary,
 - ``FiniteElementBase`` is forward declared as we will use finite element dicretization,
 - ``DomainPartition`` is forward declared as we will use MPI parallelism and its associated domain decomposition

Let us jump forward to the class enum and variable as they contain the data used
specifically in the implementation of *LaplaceFEM*.

class enums and variables
``````````````````````````````
The class exhibits three member variables:
 - ``m_fieldName`` which stores the name of the diffused variable (*e.g.* the temperature) as a `string`
 - ``m_stabletdt`` type-specialized struct which stores the maximum timestep as a `real64`
 - ``m_timeIntegrationOption`` an `enum` value allowing to dispatch with respect to the transient treatment

``timeIntegrationOption`` is an `enum` specifying the transient treatment which can be chosen respectively
between *SteadyState* , *ImplicitTransient* and *ExplicitTransient* depending on wheter we are interested
in the transient state and wheter we want it to be discretized as a backward or a forward Euler scheme.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

Once explained the main variables and enum, let us start reading through the different member functions:

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: //END_SPHINX_INCLUDE_02

Start looking at the class *LaplaceFEM* constructor and desctructor declarations
shows the usual `string` ``name`` and `Group*` pointer to ``parent`` that are required
to build the global file-system like structure of GEOSX. (see :ref:`GroupPar` for details).
It can also be noted that the nullary constructor is deleted on purpose to avoid compiler
automatic generation and user misuse.

The next method ``CatalogName()`` is static and return the key to be added to the *Catalog* for this type of solver
(see :ref:`ObjectCatalogPar` for details). It has to be paired with the following macro in the implementation file

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

The following member function ``RegisterDataOnMesh()`` is used to assign fields onto the discretized mesh object and
will be further discussed in the :ref:`Implementation` section.

The next block consists in solver interface functions. These member functions are meant to set up
and specialized every step towards the system matrix assembly and the solver stage.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_03
   :end-before: //END_SPHINX_INCLUDE_03


Eventually, ``ApplyDirichletBC_implicit()`` is the working specialized member functions called
as ``ApplyBoundaryConditions()`` is called in this particular class override.

Browsing the Mother class ``SolverBase``, it can be noted that most of the solver interface functions are called during
either ``SolverBase::LinearImplicitStep()`` or ``SolverBase::NonLinearImplicitStep()`` depending on the solver strategy chosen.

The next couple of methods are inlined getters whose name are self-explanatory.

Switching to protected members, ``PostProcessInput()`` is a central member function and
will be called by ``Group`` object after input is read from XML entry file.
It will set and dispatch solver variables from the Mother class ``BaseSolver`` to the most dervied class.
For *LaplaceFEM*, it will allow us to set the right time integration scheme based ont the XML value
as will be further explored in the next :ref:`Implementation` section.

Let us jump back to focus on a ``struct`` which play also an important role,

*viewKeyStruct* structure
`````````````````````````````

This embedded instantiated structure is a common pattern shared by all solvers.
It stores ``dataRepository::ViewKey`` type objects that are used as binding data
between the input XML file and the source code.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_04
   :end-before: //END_SPHINX_INCLUDE_04

We can check that in the *LaplaceFEM* companion integratedTest

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/integratedTests/10x10x10_LaplaceFEM.xml
   :language: xml
   :start-after: <Solvers>
   :end-before: </Solvers>

In the following section, we sill see where this binding takes place.

.. _Implementation:

Implementation File
--------------------
Switching to implementation, we will focus on few implementations, leaving details
to other tutorials.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

Checking out the constructor, we can see that the use of a ``registerWrapper<T>(...)``
allow us to register the key value from the enum ``viewKeyStruct`` defining them as :

 - ``InputFlags::OPTIONAL`` if they are optional and can be provided
 - ``InputFlags::REQUIRED`` if they are required and will throw error if not

and their associated descriptions for auto-generated docs.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: //END_SPHINX_INCLUDE_02

``RegisterDataOnMesh()`` is browsing all subgroups in the mesh ``Group`` object and
for all nodes in the sub group

 - register the observed field under the chosen ``m_fieldName`` key
 - apply a default value
 - set the output verbosity level (here ``PlotLevel::LEVEL_0``)
 - set the field associated description for auto generated docs

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_03
   :end-before: //END_SPHINX_INCLUDE_03

``PostProcessInput()`` will ensure all dispatches and assignements of all read values from the base class defined
to the deepest derived class. In the Mother class *BaseSolver*, it will set the gravity vector value and the linear solver parameters.
In *LaplaceFEM* implementation, it will assign the ``m_timeIntegrationOption``
to the read value and throw a runtime error if it not among the `enum` values.

.. literalinclude:: ../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_04
   :end-before: //END_SPHINX_INCLUDE_04

``AssembleSystem()`` will be our core focus as we want to change the diffusion coefficient from its
hard coded value to a XML read user-defined value. One can see that this method is in charge of constructing
in a parallel fashion the FEM system matrix. Bringing ``nodeManager`` and ``ElementRegionManager`` from domain local
``MeshLevel`` object together with ``FiniteElementDiscretizationManager`` from the ``NumericalMethodManager``, it uses
nodes embedded loops on degrees of freedom in a local index embedded loops to fill a matrix and a rhs containers.

As we spotted the place to change in a code to get a user-defined diffusion coefficient into the game, let us jump
to writing our new *LaplaceDiffFEM* solver.

Note: we might want to remove final keyword from ``PostProcessInput()`` as it will prevent you from overriding it

Start doing your own Physic solver
===================================
As we will extend *LaplaceFEM* capabilities, we will derive publicly from it

Declaration File
-----------------

As there is only few places where we have to change, the whole declaration file is reported below and
commented afterwards.

.. code-block:: c++

                #include "physicsSolvers/simplePDE/LaplaceFEM.hpp"

                namespace geosx
                {

                class LaplaceDiffFEM : public LaplaceFEM
                {
                public:

                LaplaceDiffFEM() = delete;

                LaplaceDiffFEM( const std::string& name,
                                Group * const parent );

                virtual ~LaplaceDiffFEM() override;

                static string CatalogName() { return "LaplaceDiffFEM"; }


                virtual void
                AssembleSystem( real64 const time,
                                real64 const dt,
                                DomainPartition * const domain,
                                DofManager const & dofManager,
                                ParallelMatrix & matrix,
                                ParallelVector & rhs ) override;


                struct viewKeyStruct : public LaplaceFEM::viewKeyStruct
                {
                   dataRepository::ViewKey diffusionCoeff = {"diffusionCoeff"};
                }  laplaceDiffFEMViewKeys;

                protected:
                virtual void PostProcessInput() override final;


                private:
                   real64 m_diffusion;

                };


We intend to have a user-defined diffusion coefficient, we then need a `real64` class variable ``m_diffusion``
to store it.

Consistently with *LaplaceFEM*, we will also elete the nullary constructor and declare a constructor with the same arguments for
forwarding to `Group` master class. Another mandatory step is to override the static ``CatalogName()`` method to properly
register any data from the new solver class.

Then as mentionned in :ref:`Implementation`, the diffusion coefficient is used when assembling the matrix coefficient. Hence
we will have to override the ``AssembleSystem()`` function as detailed below.

Moreover, if we want to introduce a new binding between the input XML and the code we will have to work on the three
``struct viewKeyStruct`` , ``PostProcessInput()`` and the constructor.

Our new solver ``viewKeyStruct`` will have its own structure inheriting from the *LaplaceFEM* one to have the ``timeIntegrationOption``
and ``fieldName`` field. It will also create a ``diffusionCoeff`` field to be bound to the user defined homogeneous coefficient on one hand
and to our ``m_diffusion`` class variable on the other.


Implementation File
---------------------
As we have seen in :ref:`Implementation`, the first place where to implement a new register from XML input is
in the constructor. The ``diffusionCoeff`` entry we have defined in the ``laplaceDiffFEMViewKeys``
will then be asked as a required input. If not provided, the error thrown will ask for it described asked
an "input uniform diffusion coeff for the laplace equation".

.. code-block:: c++

                LaplaceDiffFEM::LaplaceDiffFEM( const std::string& name,
                        Group * const parent ):
                LaplaceFEM( name, parent ), m_diffusion(0.0)
                {

                 registerWrapper<string>(laplaceDiffFEMViewKeys.diffusionCoeff.Key())->
                 setInputFlag(InputFlags::REQUIRED)->
                 setDescription("input uniform diffusion coeff for the laplace equation");
                }

Another important spot for binding the value of the XML read parameter to our ``m_diffusion`` is in ``PostProcessInput()``

.. code-block:: c++

     void LaplaceDiffFEM::PostProcessInput()
     {
       LaplaceFEM::PostProcessInput();

       string sDiffCoeff = this->getReference<string>(laplaceDiffFEMViewKeys.diffusionCoeff);
       this->m_diffusion = std::stof(sDiffCoeff);

    }

Now that we have required, read and bind the user-defined diffusion value to a variable, we can use it in the construction of our
matrix into the overrriden ``AssembleSystem()``.

.. code-block:: c++
  :emphasize-lines: 17-19

      // begin element loop, skipping ghost elements
      for( localIndex k=0 ; k<elementSubRegion->size() ; ++k )
      {
        if(elemGhostRank[k] < 0)
        {
          element_rhs = 0.0;
          element_matrix = 0.0;
          for( localIndex q=0 ; q<n_q_points ; ++q)
          {
            for( localIndex a=0 ; a<numNodesPerElement ; ++a)
            {
              elemDofIndex[a] = dofIndex[ elemNodes( k, a ) ];

              for( localIndex b=0 ; b<numNodesPerElement ; ++b)
              {

                element_matrix(a,b) += detJ[k][q] *
                                       m_diffusion *
                                     + Dot( dNdX[k][q][a], dNdX[k][q][b] );
              }

            }
          }
          matrix.add( elemDofIndex, elemDofIndex, element_matrix );
          rhs.add( elemDofIndex, element_rhs );
        }
      }


This complete the implementation of our new solver *LaplaceDiffFEM*.

Nonetheless, the compiler should complain that ``m_fieldName`` is privately as inherited from *LaplaceFEM*. One should then either promote ``m_fieldName`` to protected
or add a getter in *LaplaceFEM* class to correct the error. The getter option has been chosen and the fix in our solver is then:

.. code-block:: c++

 array1d<globalIndex> const & dofIndex =
    nodeManager->getReference< array1d<globalIndex> >( dofManager.getKey( getFieldName() ) );


Note: For consistency do not forget to change LaplaceFEM to LaplaceDiffFEM in the guards comments

Last steps
===========

After assembling both declarations and implementations for our new solver, the final steps go as:

 - Add declarations to parent CMakeLists.txt (here add to ``physicsSolvers_headers`` )
 - Add implementations to parent CMakeLists.txt (here add to ``physicsSolvers_sources``)
 - Check that Doxygen comments are properly set in our solver class
 - Uncrustify it to match the code style
 - Write unit tests for each new features in the solver class
 - Write an integratedTests for the solver class
