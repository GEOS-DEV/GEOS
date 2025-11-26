###############################################################################
Code Rules
###############################################################################

.. contents:: Table of Contents
   :depth: 3
   :local:
   :backlinks: none

============
Introduction
============

This document outlines the coding rules and best practices that must be followed when contributing to GEOS. 
These rules ensure code consistency, performance, maintainability, and reliability across the entire codebase, 
especially considering our HPC context with CPU/GPU execution environments.

All code must comply with these rules before being merged into the develop branch. Violations should be caught 
during code review, and CI tests enforce many of these rules automatically.

=======
Logging
=======

Output Stream Rules
-------------------

**Rule: No Direct Output Streams**
  Never use ``std::cout``, ``std::cerr``, ``printf``, or any direct output functions in production code.
  Always use the GEOS logging macros defined in ``Logger.hpp``.

  **Rationale:** Direct output streams bypass MPI rank control and log level filtering, leading to 
  uncontrolled output in parallel runs.

  .. code-block:: cpp

    // WRONG
    std::cout << "Value is " << value << std::endl;
    printf("Error: %s\n", message);
    std::cerr << "Debug: " << variable << "\n";

    // CORRECT
    GEOS_LOG( "Value is " << value );
    GEOS_LOG_RANK_0( "Value is " << value );  // Only on rank 0
    GEOS_LOG_LEVEL( logInfo::Convergence, "Value is " << value );

**Rule: Avoid Log Spam**
  Use appropriate log levels and conditional logging to avoid flooding output, especially in loops
  and frequently called functions.

  **Rationale:** Excessive logging degrades performance and makes important messages hard to find.

  .. code-block:: cpp

    // WRONG - logs every iteration
    for( localIndex iter = 0; iter < 1000; ++iter )
    {
      GEOS_LOG( "Iteration " << iter );
      // ...
    }

    // CORRECT - conditional logging
    for( localIndex iter = 0; iter < 1000; ++iter )
    {
      GEOS_LOG_RANK_0_IF( iter % 100 == 0, 
                          GEOS_FMT( "Iteration {}: residual = {:.2e}", iter, residual ) );
      // ...
    }

GPU Kernel Logging Restrictions
--------------------------------

**Rule: No Logging in Device Code**
  Never use logging macros in RAJA kernels or any code that may execute on GPU in production builds.

  **Rationale:** Logging in device code causes cache reservations, produces suboptimal code, and 
  will fail compilation on most GPU architectures.

  .. code-block:: cpp

    // WRONG - will fail on GPU
    forAll< policy >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      GEOS_LOG( "Processing element " << i );  // COMPILATION ERROR on GPU!
      array[i] *= 2.0;
    });

    // CORRECT - log before/after kernel
    GEOS_LOG_RANK_0( "Starting kernel execution for " << size << " elements" );
    forAll< policy >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      array[i] *= 2.0;
    });
    GEOS_LOG_RANK_0( "Kernel execution complete" );

    // OK for temporary debugging only (remove before merge)
    #ifdef DEBUG_KERNEL
    forAll< policy >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      if( i == debugIndex )
        printf( "Debug: element %d value %f\n", i, array[i] );
      array[i] *= 2.0;
    });
    #endif

Structured Data Output
----------------------

**Rule: Use Proper Formatting Tools**
  For structured data output, use specialized formatting tools instead of manual formatting.

  **Rationale:** Ensures consistent, maintainable, and correctly aligned output.

  .. code-block:: cpp

    // WRONG - manual formatting
    GEOS_LOG( "Iter   Residual   Time" );
    GEOS_LOG( iter << "    " << residual << "    " << time );

    // CORRECT - use TableFormatter
    TableFormatter table;
    table.addColumn( "Iteration" ).addColumn( "Residual" ).addColumn( "Time [s]" );
    table.addRow( iter, residual, elapsedTime );
    GEOS_LOG( table.toString() );

================
Error Management
================

Error Handling Hierarchy
------------------------

**Rule: Use Appropriate Error Macros**
  Choose the correct error macro based on severity and recoverability.

  **Rationale:** Proper error handling enables graceful degradation and helps debugging.

  .. admonition:: Error Macro Selection Guide

    - ``GEOS_ERROR``: Unrecoverable errors requiring immediate termination
    - ``GEOS_THROW``: Errors that should propagate with context
    - ``GEOS_WARNING``: Non-critical issues that don't stop execution
    - ``GEOS_ASSERT``: Debug-only checks (compiled out in release)

Fatal Errors
------------

**Rule: Use GEOS_ERROR for Unrecoverable States**
  Use ``GEOS_ERROR`` when the program cannot continue safely.

  **Rationale:** Immediate termination prevents data corruption and cascading failures.

  .. code-block:: cpp

    // CORRECT - unrecoverable dimension mismatch
    GEOS_ERROR_IF( matrix.numRows() != rhs.size(),
                   GEOS_FMT( "Matrix-vector dimension mismatch: {} != {}", 
                            matrix.numRows(), rhs.size() ) );
    
    // CORRECT - using comparison macros
    GEOS_ERROR_IF_NE_MSG( status, SUCCESS, 
                          "Solver failed with status code " << status );

    // CORRECT - null pointer check
    GEOS_ERROR_IF( ptr == nullptr,
                   "Critical pointer is null in " << __FUNCTION__ );

Recoverable Errors
------------------

**Rule: Use GEOS_THROW for Propagating Errors**
  Use ``GEOS_THROW`` when parent code might handle the error.

  **Rationale:** Allows error context to be added at different levels.

  .. code-block:: cpp

    // CORRECT - input validation
    GEOS_THROW_IF( dt <= 0.0, 
                   GEOS_FMT( "Invalid timestep: {} <= 0", dt ),
                   InputError );

    // CORRECT - with data context
    GEOS_THROW_IF( !materialModel,
                   "Material model not found",
                   std::runtime_error,
                   getWrapperDataContext( viewKeyStruct::constitutiveNameString() ) );

  .. warning::
    **Never catch and recover from exceptions in normal flow.** Exceptions should only be used
    for error propagation due to performance implications and GPU incompatibility.

Non-Critical Issues
-------------------

**Rule: Use GEOS_WARNING for Non-Blocking Issues**
  Use ``GEOS_WARNING`` for issues that should be reported but don't require stopping execution.

  **Rationale:** Allows users to identify potential problems without disrupting computation.

  .. code-block:: cpp

    // CORRECT - performance warning
    GEOS_WARNING_IF( numIterations > 0.8 * maxIterations,
                     GEOS_FMT( "Solver approaching iteration limit: {}/{}", 
                              numIterations, maxIterations ) );

    // CORRECT - deprecation notice
    GEOS_WARNING( "Method computeOld() is deprecated, use computeNew() instead" );

=====================
Parallel/MPI Rules
=====================

Collective Operation Consistency
---------------------------------

**Rule: All Ranks Must Call Collective Operations**
  Every MPI collective operation must be called by all ranks in the communicator with compatible parameters.

  **Rationale:** Mismatched collective calls cause deadlocks or undefined behavior.

  .. code-block:: cpp

    // WRONG - conditional collective causes DEADLOCK
    void processData( real64 const localValue )
    {
      if( MpiWrapper::commRank() == 0 )  
      {
        real64 globalSum = MpiWrapper::sum( localValue );  // DEADLOCK! Only rank 0 calls
      }
    }

    // WRONG - different operations on different ranks
    void inconsistentOperation( real64 const value )
    {
      if( MpiWrapper::commRank() == 0 )
        real64 result = MpiWrapper::sum( value );      // sum on rank 0
      else
        real64 result = MpiWrapper::max( value );      // max on others - DEADLOCK!
    }

    // CORRECT - all ranks participate
    void processDataCorrectly( real64 const localValue )
    {
      // All ranks call sum
      real64 globalSum = MpiWrapper::sum( localValue );
      
      // Only rank 0 uses the result
      if( MpiWrapper::commRank() == 0 )
      {
        GEOS_LOG( "Global sum: " << globalSum );
      }
    }

    // CORRECT - consistent operations with different data
    void consistentOperation( array1d< real64 > const & data )
    {
      real64 localMax = data.empty() ? 0.0 : LvArray::math::max( data.begin(), data.end() );
      real64 globalMax = MpiWrapper::max( localMax );  // All ranks call with their local value
    }

Communication Batching
----------------------

**Rule: Batch Small Communications**
  Combine multiple small messages into larger ones to reduce communication overhead.

  **Rationale:** MPI has significant per-message overhead; batching improves performance.

  .. code-block:: cpp

    // WRONG - many small communications
    void sendDataInefficiently( arrayView1d< real64 const > const & data,
                                int const destination )
    {
      for( localIndex i = 0; i < data.size(); ++i )
      {
        MpiWrapper::send( &data[i], 1, destination, tag + i );  // N messages!
      }
    }
    
    // CORRECT - single batched communication
    void sendDataEfficiently( arrayView1d< real64 const > const & data,
                              int const destination )
    {
      MpiWrapper::send( data.data(), data.size(), destination, tag );  // 1 message
    }

    // CORRECT - pack multiple arrays
    void sendMultipleArrays( arrayView1d< real64 const > const & array1,
                            arrayView1d< real64 const > const & array2 )
    {
      buffer_type buffer;
      buffer.pack( array1.size() ).pack( array1 ).pack( array2.size() ).pack( array2 );
      MpiWrapper::send( buffer.data(), buffer.size(), destination, tag );
    }

Non-Blocking Communication
--------------------------

**Rule: Overlap Communication with Computation**
  Use non-blocking communication to hide latency when possible.

  **Rationale:** Enables computation-communication overlap, improving parallel efficiency.

  .. code-block:: cpp

    // WRONG - sequential communication and computation
    void sequentialGhostUpdate( DomainPartition & domain )
    {
      sendReceiveGhosts( domain );       // Blocking communication
      computeInteriorCells( domain );     // Computation waits
      computeBoundaryCells( domain );
    }

    // CORRECT - overlapped communication and computation
    void overlappedGhostUpdate( DomainPartition & domain )
    {
      // Start non-blocking ghost exchange
      MPI_Request ghostRequests[2];
      startGhostExchange( domain, ghostRequests );
      
      // Compute interior while ghosts are in flight
      computeInteriorCells( domain );
      
      // Wait for ghosts before boundary computation
      MpiWrapper::waitAll( 2, ghostRequests );
      
      // Now safe to compute boundary cells
      computeBoundaryCells( domain );
    }

======================================
Precision and Numerical Considerations
======================================

Floating-Point Comparisons
--------------------------

**Rule: Never Use Direct Equality for Floating-Point**
  Always use tolerance-based comparisons for floating-point numbers.

  **Rationale:** Floating-point arithmetic has inherent rounding errors.

  .. code-block:: cpp

    // WRONG - exact comparison
    if( value == 0.0 ) { }
    if( a == b ) { }
    
    // CORRECT - absolute tolerance for zero comparison
    constexpr real64 machinePrecision = std::numeric_limits< real64 >::epsilon();
    constexpr real64 tolerance = 100.0 * machinePrecision;
    
    if( LvArray::math::abs( value ) < tolerance ) 
    { 
      // value is effectively zero
    }
    
    // CORRECT - relative tolerance for general comparison
    bool areEqual( real64 const a, real64 const b )
    {
      real64 const diff = LvArray::math::abs( a - b );
      real64 const magnitude = LvArray::math::max( LvArray::math::abs( a ), 
                                                    LvArray::math::abs( b ) );
      return diff <= tolerance * magnitude;
    }

    // CORRECT - combined absolute and relative tolerance
    bool robustEqual( real64 const a, real64 const b,
                      real64 const relTol = 1e-10,
                      real64 const absTol = 1e-14 )
    {
      return LvArray::math::abs( a - b ) <= relTol * LvArray::math::max( 
               LvArray::math::abs( a ), LvArray::math::abs( b ) ) + absTol;
    }

Division Safety
---------------

**Rule: Check Denominators Before Division**
  Always verify denominators are not zero or near-zero before division.

  **Rationale:** Prevents NaN/Inf propagation and numerical instabilities.

  .. code-block:: cpp

    // WRONG - unprotected division
    real64 const normalizedResidual = residual / initialResidual;
    real64 const strainRate = velocityGradient / thickness;

    // CORRECT - protected division
    real64 computeNormalizedResidual( real64 const residual, 
                                      real64 const initialResidual )
    {
      if( initialResidual > machinePrecision )
        return residual / initialResidual;
      else
        return residual;  // or return a flag indicating special case
    }

    // CORRECT - with error reporting
    real64 safeDivide( real64 const numerator, 
                       real64 const denominator,
                       string const & context )
    {
      GEOS_ERROR_IF( LvArray::math::abs( denominator ) < machinePrecision,
                     GEOS_FMT( "Division by zero in {}: denominator = {:.2e}", 
                              context, denominator ) );
      return numerator / denominator;
    }

Overflow Prevention
-------------------

**Rule: Check for Potential Overflow**
  Validate that operations won't exceed type limits, especially for index calculations.

  **Rationale:** Overflow leads to undefined behavior and memory corruption.

  .. code-block:: cpp

    // WRONG - unchecked multiplication
    globalIndex totalDofs = numElements * numNodesPerElement * numDofsPerNode;

    // CORRECT - overflow checking
    globalIndex computeTotalDofs( globalIndex const numElements,
                                  localIndex const numNodesPerElement,
                                  localIndex const numDofsPerNode )
    {
      globalIndex const maxValue = std::numeric_limits< globalIndex >::max();
      
      // Check first multiplication
      GEOS_ERROR_IF( numElements > maxValue / numNodesPerElement,
                     "DOF calculation would overflow: elements * nodes" );
      
      globalIndex const temp = numElements * numNodesPerElement;
      
      // Check second multiplication
      GEOS_ERROR_IF( temp > maxValue / numDofsPerNode,
                     "DOF calculation would overflow: total DOFs" );
      
      return temp * numDofsPerNode;
    }

    // CORRECT - capped increment for iterative methods
    real64 applyIncrement( real64 const current, real64 const increment )
    {
      constexpr real64 maxRelativeChange = 10.0;
      
      if( LvArray::math::abs( increment ) > maxRelativeChange * LvArray::math::abs( current ) )
      {
        GEOS_WARNING( GEOS_FMT( "Large increment capped: {:.2e} -> {:.2e}", 
                               increment, maxRelativeChange * current ) );
        return current + LvArray::math::signum( increment ) * maxRelativeChange * 
               LvArray::math::abs( current );
      }
      return current + increment;
    }

Conditioning Monitoring
-----------------------

**Rule: Monitor Numerical Conditioning**
  Check and report conditioning issues in linear systems and algorithms.

  **Rationale:** Poor conditioning leads to inaccurate results and convergence failures.

  .. code-block:: cpp

    // CORRECT - condition number monitoring
    class LinearSystemSolver
    {
      void checkConditioning( Matrix const & matrix )
      {
        real64 const conditionNumber = estimateConditionNumber( matrix );
        
        GEOS_WARNING_IF( conditionNumber > 1e10,
                        GEOS_FMT( "Poor conditioning detected: κ = {:.2e}", conditionNumber ) );
        
        GEOS_ERROR_IF( conditionNumber > 1e16,
                      "Matrix is numerically singular" );
        
        // Adjust solver based on conditioning
        if( conditionNumber > 1e8 )
        {
          m_solver.useRobustPreconditioner();
          m_solver.increaseTolerance( 10.0 );
        }
      }
    };

======================
Performance Guidelines  
======================

Loop Optimization
-----------------

**Rule: Hoist Loop Invariants**
  Move computations that don't change during iterations outside the loop.

  **Rationale:** Reduces redundant computations, especially critical in GPU kernels with limited registers.

  .. code-block:: cpp

    // WRONG - recomputing invariants
    forAll< policy >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      real64 const factor = viscosity * deltaTime / density;  // Computed size times!
      velocity[i] += factor * acceleration[i];
    });

    // CORRECT - precompute invariants
    real64 const factor = viscosity * deltaTime / density;  // Computed once
    forAll< policy >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      velocity[i] += factor * acceleration[i];
    });

    // CORRECT - complex invariant hoisting
    void computeStress( arrayView3d< real64 > const & stress,
                        arrayView3d< real64 const > const & strain )
    {
      // Hoist material property calculations
      real64 const lambda = bulkModulus - 2.0 * shearModulus / 3.0;
      real64 const twoMu = 2.0 * shearModulus;
      
      forAll< policy >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        // Only element-specific calculations in loop
        real64 const volumetricStrain = strain[k][0][0] + strain[k][1][1] + strain[k][2][2];
        stress[k][0][0] = lambda * volumetricStrain + twoMu * strain[k][0][0];
        stress[k][1][1] = lambda * volumetricStrain + twoMu * strain[k][1][1];
        stress[k][2][2] = lambda * volumetricStrain + twoMu * strain[k][2][2];
      });
    }

Memory Access Patterns
----------------------

**Rule: Optimize Memory Access for Cache and Coalescing**
  Access memory sequentially and ensure coalesced access on GPUs.

  **Rationale:** Memory bandwidth is often the bottleneck; proper access patterns can improve performance 10x or more.

  .. code-block:: cpp

    // WRONG - strided access (poor cache usage)
    void transposeNaive( arrayView2d< real64 > const & output,
                         arrayView2d< real64 const > const & input )
    {
      forAll< policy >( input.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        for( localIndex j = 0; j < input.size( 1 ); ++j )
        {
          output[j][i] = input[i][j];  // Strided write!
        }
      });
    }

    // CORRECT - tiled transpose for better cache usage
    void transposeTiled( arrayView2d< real64 > const & output,
                        arrayView2d< real64 const > const & input )
    {
      constexpr localIndex tileSize = 32;
      
      forAll< policy >( input.size( 0 ) / tileSize, 
                       [=] GEOS_HOST_DEVICE ( localIndex const tileI )
      {
        for( localIndex tileJ = 0; tileJ < input.size( 1 ) / tileSize; ++tileJ )
        {
          // Process tile with better locality
          for( localIndex i = 0; i < tileSize; ++i )
          {
            for( localIndex j = 0; j < tileSize; ++j )
            {
              localIndex const globalI = tileI * tileSize + i;
              localIndex const globalJ = tileJ * tileSize + j;
              output[globalJ][globalI] = input[globalI][globalJ];
            }
          }
        }
      });
    }

Kernel Fusion
-------------

**Rule: Combine Related Kernels**
  Fuse multiple simple kernels to reduce memory traffic and launch overhead.

  **Rationale:** Memory bandwidth is precious; each kernel launch has overhead.

  .. code-block:: cpp

    // WRONG - multiple passes over data
    void updateParticlesSeparate( arrayView1d< R1Tensor > const & velocity,
                                  arrayView1d< R1Tensor > const & position,
                                  arrayView1d< R1Tensor const > const & acceleration,
                                  real64 const dt )
    {
      // First kernel - update velocity
      forAll< policy >( velocity.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        velocity[i] += acceleration[i] * dt;
      });
      
      // Second kernel - update position (requires updated velocity)
      forAll< policy >( position.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        position[i] += velocity[i] * dt;
      });
    }

    // CORRECT - single fused kernel
    void updateParticlesFused( arrayView1d< R1Tensor > const & velocity,
                               arrayView1d< R1Tensor > const & position,
                               arrayView1d< R1Tensor const > const & acceleration,
                               real64 const dt )
    {
      forAll< policy >( velocity.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        velocity[i] += acceleration[i] * dt;
        position[i] += velocity[i] * dt;  // Single pass through memory
      });
    }

==========
Data Types
==========

Type Consistency
----------------

**Rule: Use GEOS Type Aliases**
  Always use GEOS-defined type aliases instead of native C++ types.

  **Rationale:** Ensures portability, consistency, and allows global type changes.

  .. list-table:: Type Correspondence
     :header-rows: 1
     :widths: 25 25 50

     * - GEOS Type
       - Native Type
       - Usage
     * - ``real64``
       - ``double``
       - Default floating-point
     * - ``real32``  
       - ``float``
       - Single precision when memory-limited
     * - ``integer``
       - ``int``
       - General integer values
     * - ``localIndex``
       - ``GEOS_LOCALINDEX_TYPE``
       - Indexing within MPI rank
     * - ``globalIndex``
       - ``GEOS_GLOBALINDEX_TYPE``
       - Global indexing across ranks
     * - ``size_t``
       - ``std::size_t``
       - Size and count types
     * - ``string``
       - ``std::string``
       - String data
     * - ``string_view``
       - ``std::string_view``
       - Non-owning string references

  .. code-block:: cpp

    // WRONG - using native types
    double pressure = 101325.0;
    int numIterations = 0;
    std::vector< double > values;

    // CORRECT - using GEOS types
    real64 pressure = 101325.0;
    integer numIterations = 0;
    array1d< real64 > values;

Container Usage
---------------

**Rule: Use GEOS Array Classes**
  Never use STL containers for numerical data. Use GEOS array classes for GPU portability.

  **Rationale:** STL containers are not GPU-portable and lack necessary memory management features.

  .. code-block:: cpp

    // WRONG - STL containers
    std::vector< double > forces;
    double matrix[100][100];
    std::array< double, 3 > vector;

    // CORRECT - GEOS arrays
    array1d< real64 > forces;
    array2d< real64 > matrix( 100, 100 );
    R1Tensor vector = { 0.0, 0.0, 0.0 };
    
    // CORRECT - views for kernel arguments
    arrayView1d< real64 > const forcesView = forces.toView();
    arrayView2d< real64 > const matrixView = matrix.toView();
    
    // CORRECT - stack arrays for small, fixed-size data
    stackArray1d< real64, 10 > localValues;

======================
Physics-Specific Rules
======================

Unit Consistency
----------------

**Rule: Document and Verify Physical Units**
  Always document units in function signatures and verify dimensional consistency.

  **Rationale:** Unit errors are a common source of bugs in physics simulations.

  .. code-block:: cpp

    /**
     * @brief Compute pressure from equation of state
     * @param[in] density fluid density [kg/m³]
     * @param[in] bulkModulus bulk modulus [Pa]
     * @param[in] referenceDensity reference density [kg/m³]
     * @return pressure [Pa]
     */
    GEOS_HOST_DEVICE
    real64 computePressure( real64 const density,        // [kg/m³]
                           real64 const bulkModulus,     // [Pa]
                           real64 const referenceDensity ) // [kg/m³]
    {
      // Dimensional analysis: [Pa] = [Pa] * [dimensionless]
      return bulkModulus * ( density / referenceDensity - 1.0 );
    }

    // CORRECT - unit conversion utilities
    namespace units
    {
      constexpr real64 convertPressure( real64 const value, 
                                        PressureUnit const from,
                                        PressureUnit const to )
      {
        real64 const pascals = ( from == PressureUnit::PSI ) ? value * 6894.76 : value;
        return ( to == PressureUnit::PSI ) ? pascals / 6894.76 : pascals;
      }
    }

Physical Bounds Checking
------------------------

**Rule: Enforce Physical Constraints**
  Validate that state variables remain within physically meaningful bounds.

  **Rationale:** Unphysical values indicate errors and can cause solver failures.

  .. code-block:: cpp

    class ConstitutiveModel
    {
      void updateState( real64 const pressure,
                       real64 const temperature,
                       real64 const porosity )
      {
        // Check physical bounds
        GEOS_ERROR_IF( pressure < -1e5, 
                       GEOS_FMT( "Pressure {:.2e} Pa below cavitation limit", pressure ) );
        
        GEOS_ERROR_IF( temperature < 0.0,
                       GEOS_FMT( "Absolute temperature {:.2f} K is negative", temperature ) );
        
        GEOS_ERROR_IF( porosity < 0.0 || porosity > 1.0,
                       GEOS_FMT( "Porosity {:.3f} outside [0,1]", porosity ) );
        
        // Warn about unusual but valid values
        GEOS_WARNING_IF( temperature > 1000.0,
                        GEOS_FMT( "High temperature {:.0f} K may be outside model validity", 
                                  temperature ) );
      }
    };

Conservation Verification
-------------------------

**Rule: Verify Conservation Laws**
  Check that mass, momentum, and energy are conserved within tolerance.

  **Rationale:** Conservation violations indicate bugs or numerical issues.

  .. code-block:: cpp

    class FlowSolver
    {
      void verifyMassConservation( DomainPartition const & domain )
      {
        real64 const totalMassInitial = computeTotalMass( domain, 0 );
        real64 const totalMassCurrent = computeTotalMass( domain, currentTime );
        real64 const netFlux = computeCumulativeFlux( domain );
        
        real64 const imbalance = LvArray::math::abs( 
          ( totalMassCurrent - totalMassInitial ) - netFlux );
        
        real64 const relativeTolerance = 1e-10;
        real64 const scale = LvArray::math::max( 
          LvArray::math::abs( totalMassInitial ), 
          LvArray::math::abs( netFlux ) );
        
        GEOS_ERROR_IF( imbalance > relativeTolerance * scale,
                       GEOS_FMT( "Mass conservation violated: imbalance = {:.2e} kg ({:.2e}%)",
                                imbalance, 100.0 * imbalance / scale ) );
      }
    };

=====================
Solver-Specific Rules
=====================

Convergence Criteria
--------------------

**Rule: Implement Robust Convergence Checks**
  Use both absolute and relative tolerances with appropriate scaling.

  **Rationale:** Single tolerance types fail for problems with varying scales.

  .. code-block:: cpp

    class ConvergenceCriteria
    {
      bool checkConvergence( real64 const residualNorm,
                            real64 const initialResidualNorm,
                            real64 const solutionNorm,
                            real64 const updateNorm ) const
      {
        // Absolute tolerance for near-zero problems
        if( residualNorm < m_absoluteTolerance )
        {
          GEOS_LOG_LEVEL( logInfo::Convergence, 
                         "Converged: absolute tolerance" );
          return true;
        }
        
        // Relative reduction from initial
        if( initialResidualNorm > 0 && 
            residualNorm < m_relativeTolerance * initialResidualNorm )
        {
          GEOS_LOG_LEVEL( logInfo::Convergence,
                         "Converged: relative reduction" );
          return true;
        }
        
        // Solution-scaled tolerance for nonlinear problems
        if( solutionNorm > 0 && 
            residualNorm < m_relativeTolerance * solutionNorm )
        {
          GEOS_LOG_LEVEL( logInfo::Convergence,
                         "Converged: solution-scaled" );
          return true;
        }
        
        // Update-based convergence for Newton methods
        if( updateNorm < m_updateTolerance * solutionNorm )
        {
          GEOS_LOG_LEVEL( logInfo::Convergence,
                         "Converged: small update" );
          return true;
        }
        
        return false;
      }
    };

Linear Solver Configuration
---------------------------

**Rule: Adapt Solver to Problem Characteristics**
  Choose solver and preconditioner based on matrix properties.

  **Rationale:** Optimal solver choice can improve performance by orders of magnitude.

  .. code-block:: cpp

    void configureLinearSolver( LinearSolver & solver,
                                localIndex const matrixSize,
                                real64 const estimatedConditionNumber,
                                bool const isSymmetric )
    {
      // Small problems: direct solver
      if( matrixSize < 1000 )
      {
        solver.parameters().setDirect();
        GEOS_LOG_LEVEL( logInfo::LinearSolver, "Using direct solver for small system" );
      }
      // Symmetric positive definite: CG
      else if( isSymmetric && estimatedConditionNumber < 1e6 )
      {
        solver.parameters().setSolver( LinearSolverParameters::SolverType::CG );
        solver.parameters().setPreconditioner( LinearSolverParameters::PreconditionerType::AMG );
        GEOS_LOG_LEVEL( logInfo::LinearSolver, "Using CG+AMG for SPD system" );
      }
      // Ill-conditioned or nonsymmetric: GMRES
      else
      {
        solver.parameters().setSolver( LinearSolverParameters::SolverType::GMRES );
        solver.parameters().setKrylovSubspaceSize( 50 );
        
        if( estimatedConditionNumber > 1e10 )
        {
          solver.parameters().setPreconditioner( LinearSolverParameters::PreconditionerType::ILU );
          solver.parameters().setILUFill( 2 );
          GEOS_LOG_LEVEL( logInfo::LinearSolver, "Using GMRES+ILU(2) for ill-conditioned system" );
        }
        else
        {
          solver.parameters().setPreconditioner( LinearSolverParameters::PreconditionerType::AMG );
          GEOS_LOG_LEVEL( logInfo::LinearSolver, "Using GMRES+AMG for general system" );
        }
      }
    }

========================
Memory Management Rules
========================

View Lifetime Management
------------------------

**Rule: Never Outlive Parent Arrays**
  Ensure array views don't outlive their parent arrays.

  **Rationale:** Dangling views cause segmentation faults and undefined behavior.

  .. code-block:: cpp

    // WRONG - returning view of local array
    arrayView1d< real64 > getDanglingView()
    {
      array1d< real64 > tempArray( 100 );
      return tempArray.toView();  // DANGER: tempArray destroyed, view dangles!
    }

    // WRONG - storing view of temporary
    class DataHolder
    {
      arrayView1d< real64 > m_view;
      
      void setData()
      {
        array1d< real64 > temp( 100 );
        m_view = temp.toView();  // DANGER: temp destroyed at end of scope!
      }
    };

    // CORRECT - return the array
    array1d< real64 > getSafeArray()
    {
      array1d< real64 > result( 100 );
      // populate...
      return result;  // Move semantics ensure safety
    }

    // CORRECT - store array, create view when needed
    class SafeDataHolder
    {
      array1d< real64 > m_data;
      
      arrayView1d< real64 > getView() { return m_data.toView(); }
      arrayView1d< real64 const > getViewConst() const { return m_data.toViewConst(); }
    };

GPU Memory Management
---------------------

**Rule: Minimize Host-Device Transfers**
  Keep data on the appropriate memory space and minimize transfers.

  **Rationale:** PCIe bandwidth is often the bottleneck in GPU applications.

  .. code-block:: cpp

    class GPUOptimizedSolver
    {
      // WRONG - repeated transfers
      void inefficientCompute( array1d< real64 > & data )
      {
        for( integer iter = 0; iter < 100; ++iter )
        {
          data.move( LvArray::MemorySpace::cuda );     // Transfer to GPU
          computeOnGPU( data.toView() );
          data.move( LvArray::MemorySpace::host );     // Transfer back
          real64 norm = computeNorm( data );           // On CPU
        }
      }
      
      // CORRECT - minimize transfers
      void efficientCompute( array1d< real64 > & data )
      {
        // Single transfer to GPU
        data.move( LvArray::MemorySpace::cuda );
        
        arrayView1d< real64 > dataView = data.toView();
        array1d< real64 > gpuNorms( 100 );
        arrayView1d< real64 > normsView = gpuNorms.toView();
        
        for( integer iter = 0; iter < 100; ++iter )
        {
          computeOnGPU( dataView );
          computeNormOnGPU( dataView, normsView[iter] );  // Stay on GPU
        }
        
        // Single transfer back
        gpuNorms.move( LvArray::MemorySpace::host );
        // Process norms on host
      }
    };

====================
Testing Requirements
====================

Numerical Testing
-----------------

**Rule: Test Edge Cases and Known Solutions**
  Include tests for boundary conditions, special cases, and problems with analytical solutions.

  **Rationale:** Ensures correctness and catches regressions.

  .. code-block:: cpp

    TEST( MatrixOperations, Inversion )
    {
      // Test identity matrix
      R2SymTensor identity = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
      R2SymTensor inverse;
      real64 det = identity.invert( inverse );
      
      EXPECT_NEAR( det, 1.0, 1e-14 );
      for( integer i = 0; i < 6; ++i )
      {
        EXPECT_NEAR( inverse[i], identity[i], 1e-14 );
      }
      
      // Test singular matrix detection
      R2SymTensor singular = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
      det = singular.invert( inverse );
      EXPECT_NEAR( det, 0.0, 1e-14 );
      
      // Test numerical stability
      R2SymTensor nearSingular = { 1.0, 1.0, 1.0 + 1e-10, 0.0, 0.0, 0.0 };
      det = nearSingular.invert( inverse );
      EXPECT_GT( LvArray::math::abs( det ), 0.0 );
      
      // Verify A * A^-1 = I
      R2SymTensor shouldBeIdentity;
      shouldBeIdentity.Rij_eq_AikBkj( nearSingular, inverse );
      for( integer i = 0; i < 3; ++i )
      {
        EXPECT_NEAR( shouldBeIdentity[i], 1.0, 1e-8 );  // Relaxed tolerance for ill-conditioned
      }
    }

GPU Testing  
-----------

**Rule: Test CPU and GPU Paths Identically**
  Ensure both execution paths produce identical results within tolerance.

  **Rationale:** Catches GPU-specific bugs and ensures portability.

  .. code-block:: cpp

    template< typename POLICY >
    class KernelTest : public ::testing::Test
    {
    protected:
      void testDotProduct()
      {
        constexpr localIndex size = 10000;
        array1d< real64 > a( size );
        array1d< real64 > b( size );
        
        // Initialize with deterministic values
        for( localIndex i = 0; i < size; ++i )
        {
          a[i] = std::sin( i * 0.1 );
          b[i] = std::cos( i * 0.1 );
        }
        
        // Compute on device
        a.move( LvArray::MemorySpace::cuda, false );
        b.move( LvArray::MemorySpace::cuda, false );
        
        arrayView1d< real64 const > aView = a.toViewConst();
        arrayView1d< real64 const > bView = b.toViewConst();
        
        RAJA::ReduceSum< typename POLICY::ReducePolicy, real64 > dotProduct( 0.0 );
        
        forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
        {
          dotProduct += aView[i] * bView[i];
        });
        
        real64 const result = dotProduct.get();
        
        // Verify against CPU reference
        real64 reference = 0.0;
        for( localIndex i = 0; i < size; ++i )
        {
          reference += a[i] * b[i];
        }
        
        // GPU reduction may have different rounding
        EXPECT_NEAR( result, reference, 1e-10 * LvArray::math::abs( reference ) );
      }
    };

    using TestTypes = ::testing::Types< serialPolicy,
    #ifdef GEOS_USE_CUDA
                                        cuda_policy,
    #endif
    #ifdef GEOS_USE_HIP
                                        hip_policy,
    #endif
                                        openmp_policy >;
    
    TYPED_TEST_SUITE( KernelTest, TestTypes );
    
    TYPED_TEST( KernelTest, DotProduct )
    {
      this->testDotProduct();
    }

=========================
Contribution Requirements
=========================

Documentation Standards
-----------------------

**Rule: Document All Public APIs**
  Provide complete Doxygen documentation for all public interfaces.

  **Rationale:** Good documentation is essential for maintainability and usability.

  .. code-block:: cpp

    /**
     * @brief Solve a nonlinear system using Newton's method
     * @tparam MATRIX_TYPE type of the system matrix
     * @param[in] jacobian the Jacobian matrix of the system
     * @param[in,out] solution initial guess on input, solution on output
     * @param[in] residualFunction functor to compute residual
     * @param[in] parameters solver configuration parameters
     * @return convergence information including iterations and final residual
     * 
     * @note The solution vector is modified in place
     * @warning The jacobian must be non-singular
     * 
     * This function implements a damped Newton method with line search.
     * The convergence criteria are specified in the parameters struct.
     */
    template< typename MATRIX_TYPE >
    ConvergenceInfo solveNonlinearSystem( MATRIX_TYPE const & jacobian,
                                          arrayView1d< real64 > const & solution,
                                          std::function< void( arrayView1d< real64 const >,
                                                              arrayView1d< real64 > ) > const & residualFunction,
                                          NonlinearSolverParameters const & parameters );

Unit Test Coverage
------------------

**Rule: Test All Code Paths**
  Ensure unit tests cover normal operation, edge cases, and error conditions.

  **Rationale:** Comprehensive testing prevents regressions and ensures reliability.

  .. code-block:: cpp

    class SolverTest : public ::testing::Test
    {
    protected:
      void SetUp() override
      {
        // Initialize test fixtures
        m_solver = std::make_unique< LinearSolver >();
        setupTestMatrix( m_matrix );
        setupTestRHS( m_rhs );
      }
      
      std::unique_ptr< LinearSolver > m_solver;
      CRSMatrix< real64 > m_matrix;
      array1d< real64 > m_rhs;
    };

    TEST_F( SolverTest, SolveIdentitySystem )
    {
      // Test trivial case
      m_matrix.setIdentity();
      array1d< real64 > solution = m_rhs;
      
      EXPECT_TRUE( m_solver->solve( m_matrix, solution ) );
      
      for( localIndex i = 0; i < m_rhs.size(); ++i )
      {
        EXPECT_NEAR( solution[i], m_rhs[i], 1e-12 );
      }
    }

    TEST_F( SolverTest, HandleSingularMatrix )
    {
      // Test error handling
      m_matrix.zero();  // Singular matrix
      array1d< real64 > solution = m_rhs;
      
      EXPECT_FALSE( m_solver->solve( m_matrix, solution ) );
      EXPECT_GT( m_solver->lastResidualNorm(), m_solver->tolerance() );
    }

=========================
Additional Best Practices
=========================

Const Correctness
-----------------

**Rule: Use Const Wherever Possible**
  Mark variables, parameters, and member functions const when they don't modify state.

  **Rationale:** Prevents bugs, enables optimizations, and documents intent.

  .. code-block:: cpp

    class DataProcessor
    {
    public:
      // Const member function
      real64 computeAverage() const  // Promises not to modify object
      {
        return m_sum / m_count;
      }
      
      // Const parameters for input
      void processData( arrayView1d< real64 const > const & input,  // Can't modify input
                       arrayView1d< real64 > const & output )        // Can modify output
      {
        // Const local variables
        real64 const scale = computeScale( input );
        localIndex const size = input.size();
        
        for( localIndex i = 0; i < size; ++i )
        {
          output[i] = scale * input[i];
        }
      }
      
    private:
      real64 m_sum;
      mutable real64 m_cachedAverage;  // Can be modified in const functions
      localIndex m_count;
    };

Thread Safety
-------------

**Rule: Avoid Global Mutable State**
  Minimize shared mutable state and use proper synchronization when necessary.

  **Rationale:** Prevents data races and ensures correctness in parallel execution.

  .. code-block:: cpp

    // WRONG - global mutable state
    namespace globals
    {
      real64 convergenceTolerance = 1e-6;  // DANGER: race condition!
    }

    // CORRECT - pass parameters explicitly
    class Solver
    {
      struct Parameters
      {
        real64 convergenceTolerance = 1e-6;
        integer maxIterations = 100;
      };
      
      void solve( Parameters const & params );  // Thread-safe
    };

    // CORRECT - use thread-local storage when necessary
    namespace stats
    {
      thread_local integer numKernelCalls = 0;  // Per-thread counter
      
      integer getTotalKernelCalls()
      {
        // Need synchronization to aggregate
        return MpiWrapper::sum( numKernelCalls );
      }
    }

Header Dependencies
-------------------

**Rule: Minimize Header Inclusions**
  Use forward declarations when possible and include only necessary headers.

  **Rationale:** Reduces compilation time and prevents circular dependencies.

  .. code-block:: cpp

    // header file: Solver.hpp
    #pragma once
    
    // WRONG - including full headers unnecessarily
    #include "Matrix.hpp"           // Full class definition
    #include "Vector.hpp"           // Full class definition  
    #include "Preconditioner.hpp"   // Full class definition
    
    // CORRECT - forward declarations when possible
    namespace geos
    {
      class Matrix;          // Forward declaration sufficient
      class Vector;          // Forward declaration sufficient
      class Preconditioner;  // Forward declaration sufficient
      
      class Solver
      {
      public:
        void solve( Matrix const & A, Vector & x );
        void setPreconditioner( Preconditioner * precond );
        
      private:
        Preconditioner * m_preconditioner = nullptr;  // Pointer - forward decl OK
      };
    }
    
    // Include full headers only in implementation file

