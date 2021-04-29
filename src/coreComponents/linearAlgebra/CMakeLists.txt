#
# Specify all headers
#
set( linearAlgebra_headers
     common/common.hpp
     interfaces/dense/BlasLapackFunctions.h
     interfaces/dense/BlasLapackLA.hpp
     interfaces/InterfaceTypes.hpp
     common/LinearOperator.hpp
     interfaces/MatrixBase.hpp
     interfaces/VectorBase.hpp
     solvers/BicgstabSolver.hpp
     solvers/BlockPreconditioner.hpp
     solvers/CgSolver.hpp
     solvers/GmresSolver.hpp
     solvers/KrylovSolver.hpp
     solvers/KrylovUtils.hpp
     common/LinearSolverBase.hpp
     common/PreconditionerBase.hpp
     solvers/PreconditionerIdentity.hpp
     solvers/PreconditionerJacobi.hpp
     solvers/PreconditionerBlockJacobi.hpp
     solvers/SeparateComponentPreconditioner.hpp
     utilities/Arnoldi.hpp
     utilities/BlockOperatorView.hpp
     utilities/BlockOperatorWrapper.hpp
     utilities/BlockOperator.hpp
     utilities/BlockVectorView.hpp
     utilities/BlockVectorWrapper.hpp
     utilities/BlockVector.hpp
     utilities/LAIHelperFunctions.hpp
     utilities/LinearSolverParameters.hpp
     utilities/LinearSolverResult.hpp
     utilities/InverseNormalOperator.hpp
     utilities/NormalOperator.hpp
     utilities/TransposeOperator.hpp
     common/traits.hpp
     DofManager.hpp
     DofManagerHelpers.hpp )

#
# Specify all sources
#
set( linearAlgebra_sources
     interfaces/dense/BlasLapackLA.cpp
     solvers/BicgstabSolver.cpp
     solvers/BlockPreconditioner.cpp
     solvers/CgSolver.cpp
     solvers/GmresSolver.cpp
     solvers/KrylovSolver.cpp
     solvers/SeparateComponentPreconditioner.cpp
     DofManager.cpp )

set( dependencyList mesh blas lapack RAJA)


list( APPEND linearAlgebra_headers
      interfaces/direct/SuiteSparse.hpp
      interfaces/direct/SuperLUDist.hpp )

list( APPEND linearAlgebra_sources
      interfaces/direct/SuiteSparse.cpp
      interfaces/direct/SuperLUDist.cpp )

list( APPEND dependencyList suitesparse superlu_dist )

if( ENABLE_TRILINOS )

    list( APPEND linearAlgebra_headers
          interfaces/trilinos/EpetraExport.hpp
          interfaces/trilinos/EpetraMatrix.hpp
          interfaces/trilinos/EpetraVector.hpp
          interfaces/trilinos/EpetraUtils.hpp
          interfaces/trilinos/TrilinosPreconditioner.hpp
          interfaces/trilinos/TrilinosSolver.hpp
          interfaces/trilinos/TrilinosInterface.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/trilinos/EpetraExport.cpp
          interfaces/trilinos/EpetraMatrix.cpp
          interfaces/trilinos/EpetraVector.cpp
          interfaces/trilinos/TrilinosPreconditioner.cpp
          interfaces/trilinos/TrilinosSolver.cpp
          interfaces/trilinos/TrilinosInterface.cpp )

    list( APPEND dependencyList trilinos )

endif()

if( ENABLE_HYPRE )

    list( APPEND linearAlgebra_headers
          interfaces/hypre/HypreExport.hpp
          interfaces/hypre/HypreMGRStrategies.hpp
          interfaces/hypre/HypreMatrix.hpp
          interfaces/hypre/HypreVector.hpp
          interfaces/hypre/HyprePreconditioner.hpp
          interfaces/hypre/HypreSolver.hpp
          interfaces/hypre/HypreInterface.hpp
          interfaces/hypre/HypreUtils.hpp 
          interfaces/hypre/mgrStrategies/HypreMGRSinglePhasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/HypreMGRHybridSinglePhasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/HypreMGRMultiphasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/HypreMGRCompositionalMultiphaseFVM.hpp
          interfaces/hypre/mgrStrategies/HypreMGRCompositionalMultiphaseReservoir.hpp
          interfaces/hypre/mgrStrategies/HypreMGRCompositionalMultiphaseHybridFVM.hpp
          interfaces/hypre/mgrStrategies/HypreMGRLagrangianContactMechanics.hpp)

    list( APPEND linearAlgebra_sources
          interfaces/hypre/HypreExport.cpp
          interfaces/hypre/HypreMatrix.cpp
          interfaces/hypre/HypreVector.cpp
          interfaces/hypre/HyprePreconditioner.cpp
          interfaces/hypre/HypreSolver.cpp
          interfaces/hypre/HypreInterface.cpp )

    list( APPEND dependencyList hypre )

endif()

if( ENABLE_PETSC )

    list( APPEND linearAlgebra_headers
          interfaces/petsc/PetscExport.hpp
          interfaces/petsc/PetscVector.hpp
          interfaces/petsc/PetscMatrix.hpp
          interfaces/petsc/PetscPreconditioner.hpp
          interfaces/petsc/PetscSolver.hpp
          interfaces/petsc/PetscInterface.hpp
          interfaces/petsc/PetscUtils.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/petsc/PetscExport.cpp
          interfaces/petsc/PetscVector.cpp
          interfaces/petsc/PetscMatrix.cpp
          interfaces/petsc/PetscPreconditioner.cpp
          interfaces/petsc/PetscSolver.cpp
          interfaces/petsc/PetscInterface.cpp )

    list( APPEND dependencyList petsc )

endif()

if( ENABLE_OPENMP )
    list( APPEND dependencyList openmp )
endif()

if ( ENABLE_CUDA )
    list( APPEND dependencyList cuda )
endif()

if ( ENABLE_MKL )
    list( APPEND dependencyList mkl )
endif()

if( ENABLE_MPI )
    list( APPEND dependencyList mpi )
endif()



blt_add_library( NAME                  linearAlgebra
                 SOURCES               ${linearAlgebra_sources}
                 HEADERS               ${linearAlgebra_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( linearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

add_subdirectory(unitTests)

geosx_add_code_checks( PREFIX linearAlgebra )
