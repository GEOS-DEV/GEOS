###############################################################################
Code Rules
###############################################################################

.. note::
   This document uses collapsible code blocks. Click on "Show/Hide Code" to expand examples.

Introduction
============
This document provides comprehensive coding rules and best practices for developing with GEOS. These rules complement the :doc:`CodeStyle` guidelines and focus on typing, error handling, parallelism, performance, and general architecture principles.

All contributors should familiarize themselves with these rules to maintain code quality, consistency, and reliability across the GEOS codebase.

Typing
======

Type Correspondence Table
-------------------------

GEOS defines its own type system for **portability**, **clarity**, and **performance**.

**Always use GEOS types** instead of standard C++ types.

Basic Types
^^^^^^^^^^^

.. list-table:: Standard C++ vs GEOS Types
   :header-rows: 1
   :widths: 30 30 40

   * - Standard C++ Type
     - GEOS Type
     - Description
   * - ``int``
     - ``integer``
     - Signed integer
   * - ``std::size_t``
     - ``localIndex``
     - Local array indexing (MPI partition)
   * - N/A
     - ``globalIndex``
     - Global indexing (across MPI)
   * - ``float``
     - ``real32``
     - 32-bit floating point
   * - ``double``
     - ``real64``
     - 64-bit floating point
   * - ``std::string``
     - ``string``
     - String type
   * - ``std::string_view``
     - ``string_view``
     - Non-owning string view

.. dropdown:: Example: Using GEOS basic types
   :icon: code

   .. code-block:: c++

      // BAD - Do not use standard types
      int count = 0;
      double value = 1.5;
      
      // GOOD - Use GEOS types
      integer count = 0;
      real64 value = 1.5;
      localIndex elementIndex = 42;

Array Types
^^^^^^^^^^^

GEOS provides multi-dimensional CHAI array types for managed host-device data.

Refer to the rule of :ref:`CHAI Memory Management`.

.. list-table:: GEOS Array Types
   :header-rows: 1
   :widths: 35 65

   * - GEOS Type
     - Description
   * - ``array1d<T>``
     - 1D owning array
   * - ``arrayView1d<T>``
     - 1D non-owning view (modifiable)
   * - ``arraySlice1d<T>``
     - 1D slice (from multi-dim array)
   * - ``array2d<T>``
     - 2D owning array
   * - ``arrayView2d<T>``
     - 2D non-owning view
   * - ``arraySlice2d<T>``
     - 2D slice
   * - ``array3d<T>``, ``array4d<T>``, ``array5d<T>``
     - 3D/4D/5D owning arrays
   * - ``arrayView3d<T>``, etc.
     - 3D/4D/5D views
   * - ``SortedArray<T>``
     - Sorted array with search capabilities
   * - ``ArrayOfArrays<T>``
     - Jagged 2D array (variable row sizes)
   * - ``ArrayOfSets<T>``
     - Array of sets (unique elements per row)
   * - ``stackArray1d<T, N>``
     - Stack-allocated array (max size N)

.. dropdown:: Example: Using GEOS arrays
   :icon: code

   .. code-block:: c++

      // BAD - Do not use std::vector for data that can be packed or kernel addressable
      std::vector<real64> values;
      std::array<integer, 10> fixedData;

      // GOOD - Use CHAI managed arrays for kernels
      array1d<real64> values;
      arrayView1d<real64> valuesView = values.toView();
      arrayView1d<real64 const> constView = values.toViewConst();
      stackArray1d<integer, 10> fixedData;

Tensor Types
^^^^^^^^^^^^

Use GEOS tensor types for geometric and mechanical calculations:

.. list-table:: GEOS Tensor Types
   :header-rows: 1
   :widths: 30 70

   * - Type
     - Description
   * - ``R1Tensor``
     - 3D vector (real64)
   * - ``R1Tensor32``
     - 3D vector (real32)
   * - ``R2SymTensor``
     - Symmetric 6-component tensor (Voigt notation)

``std`` Standard Library
------------------------

It is possible to use ``std`` library components for data on host memory, for doing so, the rule is 
to **only use GEOS ``std`` container wrappers** instead of direct standard library containers.

This rule allow us to control bounds checking depending on ``GEOS_USE_BOUNDS_CHECK`` macro / ``GEOS_ENABLE_BOUNDS_CHECK`` cmake option..

   * - Standard C++ Type
     - GEOS Type
     - Description
   * - ``std::vector``
     - ``stdVector``
     - ``std`` dynamically allocated array
   * - ``std::array``
     - ``stdArray``
     - ``std`` fixed constexpr sized array
   * - ``std::map``
     - ``map``
     - ``std`` sorted dictionary
   * - ``std::unordered_map``
     - ``unordered_map``
     - ``std`` unsorted dictionary (hash-map)

.. dropdown:: Example: Container usage
   :icon: code

   .. code-block:: c++

      // BAD - Direct std containers
      std::vector<real64> data;
      std::array<integer, 10> fixedData;
      std::map<string, real64> orderedMap;
      
      // GOOD - GEOS wrappers
      stdVector<real64> data;
      stdArray<integer, 10> fixedData;
      map<string, real64> orderedMap;

The following standard library components are allowed:
- ``std::pair``, ``std::tuple`` for temporary structures, but beware of memory allocation, sometimes struct are to prefer.
- ``std::function`` for callbacks
- ``std::optional`` for optional return values
- ``std::move``, ``std::forward`` for move semantics
- Algorithms from ``<algorithm>`` where appropriate
- ``std::unique_ptr``, ``std::shared_ptr`` for memory management

Contribution
============

Documentation
-------------

Always Provide Comprehensive Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All new features must include appropriate documentation at multiple levels:

Wrapper Documentation (User-Oriented)
""""""""""""""""""""""""""""""""""""""

All data repository wrappers must be documented with ``setDescription()``.
As much as possible, valid values rules should be provided, and default values should be provided with ``setApplyDefaultValue()``.

.. dropdown:: Example: Wrapper documentation
   :icon: code

   .. code-block:: c++

      registerWrapper( "timeStep", &m_timeStep ).
        setInputFlag( InputFlags::REQUIRED ).
        setApplyDefaultValue( 1.0e-6 ).
        setDescription( "Initial time step value. This parameter controls "
                        "the starting time increment for the simulation. "
                        "Valid range: (0, 1e-3]. See User Guide Chapter 5." );

RST Documentation (User Guide)
"""""""""""""""""""""""""""""""

User-facing features must be documented in the Sphinx RST documentation (``src/docs/sphinx/``):

- Explain **what** the feature does
- Provide **usage examples**
- Document **input parameters**
- Include **expected output** or behaviour
- Add **references** to related features

Doxygen & Naming (Developer-Oriented)
""""""""""""""""""""""""""""""""""""""

Keep Doxygen comments and naming as clear as possible. Keep in mind that it is target to developers from other domains.

.. dropdown:: Example: Doxygen documentation
   :icon: code

   .. code-block:: c++
      // GOOD - The documentation gives useful information about what the function *does*, without
      //        going into theorical knowledge nor implementation details.
      /**
       * @brief Computes the residual for the nonlinear solver.
       * @param[in] domain The domain partition containing mesh and field data
       * @param[in] time Current simulation time
       * @param[in] dt Time step size
       * @param[out] residual The computed residual vector
       * @return Residual L2 norm
       * 
       * This method assembles the discrete residual for the physics equations
       * by looping over all elements in the target regions. The residual is
       * stored in the provided vector and its L2 norm is returned.
       */
      real64 computeResidual( DomainPartition & domain,
                              real64 const time,
                              real64 const dt,
                              arrayView1d<real64> const & residual );

      // BAD - Documentation does not tell anything about the point of this function.
      /**
        * @brief Compute CFL numbers
        * @param[in] time current time
        * @param[in] dt the time step size
        * @param[in] domain the domain partition
        */
      void computeCFLNumbers( real64 const time,
                              real64 const dt,
                              DomainPartition & domain ) const;


Testing
-------

Unit Tests
^^^^^^^^^^

**Always provide unit tests** to ensure code behaves according to expectations.
Unit testing is not about testing every single line of code (an unrealistic goal), but about *verifying assumptions and preventing regressions*. Don't assume your code works as intended—prove it with tests.
- Focus on **critical logic** (core functionality and algorithms), **edge cases** (extreme, erroneous, empty values).
- **Use GoogleTest framework**.
- In GEOS, the goal of unit test is never to take position on the validity of physical results (this is the point of integrated-tests).
- Place tests in ``src/coreComponents/<component>/tests/``.
- Test GPU code paths if applicable (use ``#ifdef GEOS_USE_DEVICE``).

.. dropdown:: Example: Unit test structure
   :icon: code

   .. code-block:: c++

      #include <gtest/gtest.h>
      
      // Test core functionality
      TEST( MyComponentTest, ComputesCorrectResult )
      {
        MyComponent component;
        EXPECT_DOUBLE_EQ( component.compute( 2.0 ), 4.0 );
      }
      
      // Test edge case
      TEST( MyComponentTest, HandlesZeroInput )
      {
        MyComponent component;
        EXPECT_DOUBLE_EQ( component.compute( 0.0 ), 0.0 );
      }
      
      // Test error detection
      TEST( MyComponentTest, RejectsInvalidInput )
      {
        MyComponent component;
        EXPECT_THROW( component.compute( -1.0 ), std::invalid_argument );
      }

      // entry point of the test binary
      int main( int argc, char * * argv )
      {
        ::testing::InitGoogleTest( &argc, argv );
        g_commandLineOptions = *geos::basicSetup( argc, argv );
        int const result = RUN_ALL_TESTS();
        geos::basicCleanup();
        return result;
      }

Code Coverage
^^^^^^^^^^^^^

**Code coverage should never decrease.** New code contributions must maintain or improve overall code coverage.
Use code-cov reports to identify untested code paths.

Integrated Tests & examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Every XML-usable & serializable component must have examples (at least, have an integrated-test case):**
- Each different usage pattern requires its own example,
- Examples must be ran as integrated tests to guaranty they keep valid,
- Place examples in ``examples/`` directory, and integrated tests in ``inputFiles/``,

For more info, please refer to :ref:`IntegratedTests.rst`.

Logging
=======

Rules
-----

Use GEOS Logging Macros (``GEOS_LOG*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Use of ``std::cout``, ``std::cerr``, or ``printf`` for logging must never appear** on the develop branch.
Always use GEOS logging macros defined in ``Logger.hpp``:

.. dropdown:: Example: Correct logging usage
   :icon: code

   .. code-block:: c++

      // BAD - Do not use
      std::cout << "Starting simulation" << std::endl;
      printf("Value = %f\n", value);

      // GOOD - Use GEOS macros
      GEOS_LOG("Starting simulation");
      GEOS_LOG_RANK_0("Master rank initialization complete");

**Available logging macros:**

- ``GEOS_LOG(...)`` - Log on all ranks
- ``GEOS_LOG_VAR(...)`` - Log variable name and value
- ``GEOS_LOG_IF(condition, msg)`` - Conditional logging
- ``GEOS_LOG_RANK_0(msg)`` - Log only on rank 0 (MPI master)
- ``GEOS_LOG_RANK_0_IF(condition, msg)`` - Conditional logging on rank 0
- ``GEOS_LOG_RANK(msg)`` - Log to rank-specific output stream
- ``GEOS_LOG_RANK_VAR(var)`` - Log variable on rank stream

Use Log Levels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Logging should be meaningful and controlled**:
- When possible, **use ``GEOS_LOG_LEVEL_RANK*`` macros with appropriate ``logInfo``** to allow output filtering,
- **Avoid unnecessary logging**, prefer rank 0 logging for global information to avoid redundant output,
- **Consider logs performance impact** (output flush frequency, formatting).

See :doc:`LogLevels` for using / adding log level (e.g., ``logInfo::Convergence``, ``logInfo::TimeStep``).

.. dropdown:: Example: Using log levels
   :icon: code

   .. code-block:: c++

      // physics solver cpp - effective logging
      GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                             GEOS_FMT( "Newton iteration {}: residual = {:4.2e}",
                                       iter, residual ) );

      // src/coreComponents/physicsSolvers/LogLevelsInfo.hpp - convergence logging config
      struct Convergence
      {
        static constexpr int getMinLogLevel() { return 1; }
        static constexpr std::string_view getDescription() { return "Convergence information"; }
      };

.. dropdown:: Example: Avoiding log spam
   :icon: code

   .. code-block:: c++

      // BAD - Will spam output for every element on every ranks
      forAll< parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        GEOS_LOG("Processing element " << ei); // DO NOT DO THIS
      });

      // GOOD - Log summary information
      GEOS_LOG_RANK_0(GEOS_FMT("Processed {} elements", numElems));

Logging in RAJA Kernels
^^^^^^^^^^^^^^^^^^^^^^^^

**Logging inside RAJA kernels only for DEBUG purposes** and must never appear on the develop branch.

**Rationale:** Logging on GPU can:
- Reserve cache lines and degrade performance,
- Cause race conditions in output,
- Generate massive amounts of output,
- Produce unpredictable behaviour (e.g. runtime crashes on AMD).

.. dropdown:: Example: Debug-only kernel logging
   :icon: code

   .. code-block:: c++

      forAll< parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        GEOS_LOG("processing..."); // Compiles, but remove before merge to develop!
        GEOS_LOG(GEOS_FMT("processing elem {}/{}", ei, numElems)); // Will crash on CUDA builds (unsupported GEOS_FMT)
        printf("DEBUG: ei=%d\n", ei); // Compiles, but remove before merge to develop!
      });

Structured Data Output
----------------------

For tabular or statistical output, use structured logging facilities instead of plain text logging:

- Use ``TableTextFormatter`` with a ``logLevel`` XML setting for formatted text table output,
- Use ``TableCSVFormatter`` with a ``writeCSV()`` XML setting for structured data output.

.. dropdown:: Example: Log table output
   :icon: code

   .. code-block:: c++
      // log output preparation, layout precomputation can be kept for periodic outputs
      TableLayout statsLogLayout( "", { "Flux(es)", "Region", "Element Count", massColumn, rateColumn } );
      m_logLayout = statsLogLayout;

      // filling table data (log format is for current timestep)
      m_logData.addRow( fluxName,
                             elementSetName,
                             GEOS_FMT( "{}", wrappedStats.stats().m_elementCount ),
                             GEOS_FMT( "{}", wrappedStats.stats().m_producedMass ),
                             GEOS_FMT( "{}", wrappedStats.stats().m_productionRate ) );

      // statistics CSV output (only if enabled and from rank 0)
      if( logLevelActive && MpiWrapper::commRank() == 0 )
      {
        string const title = GEOS_FMT( "{}, flux statistics for: {}",
                                      getName(), fluxesStr );
        string const multilineTitle = wrapTextToMaxLength( title, 80 ); // enabling automatic line wrapping
        m_logLayout.setTitle( multilineTitle );
        TableTextFormatter const tableStatFormatter( m_logLayout );
        GEOS_LOG_RANK( tableStatFormatter.toString( statsData ) );
      }

.. dropdown:: Example: CSV output
   :icon: code

   .. code-block:: c++

      // CSV output preparation, layout precomputation can be kept for periodic outputs
      TableLayout const statsCSVLayout( "", {"Time [s]", "Flux(es)", "Region", "Element Count", massColumn, rateColumn} );
      m_csvLayout = statsCSVLayout;

      // filling table data (CSV format is all timestep series)
      m_csvData.addRow( wrappedStats.getStatsPeriodStart(),
                             fluxName,
                             elementSetName,
                             wrappedStats.stats().m_elementCount,
                             wrappedStats.stats().m_producedMass,
                             wrappedStats.stats().m_productionRate );

      // statistics CSV output (only if enabled and from rank 0)
      if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
      {
        std::ofstream outputFile( m_csvFilename );
        TableCSVFormatter const tableStatFormatter( m_csvLayout );
        outputFile << tableStatFormatter.toString( csvData );
        outputFile.close();
        csvData.clear();
      }

Error Management
================

- **Developpers must provide runtime validations** of their code behaviour,
- for runtime validation fails that can be encountered by users, **meaningful messages and context must be provided**.

Errors (GEOS_ERROR*)
--------------------

Why Use Errors?
^^^^^^^^^^^^^^^

Use ``GEOS_ERROR*`` macros when encountering a **blocking error** where:

- The simulation result would be certifiably invalid,
- The application state is unrecoverable,
- Continuing execution would be unsafe or meaningless.

.. dropdown:: Example: Error handling
   :icon: code

   .. code-block:: c++

      // Example: not implemented feature
      GEOS_ERROR( "Rounding interpolation with derivatives not implemented" );

      // Example: invalid mesh topology
      GEOS_ERROR_IF( numFaces == 0, 
                     GEOS_FMT("Element {} has zero faces - invalid mesh", elemId) );

      // Example: comparison-based error
      GEOS_ERROR_IF_LT( dt, 0.0, "Time step must be positive" );

Exceptions (GEOS_THROW*)
------------------------

Why Use Exceptions?
^^^^^^^^^^^^^^^^^^^

Use ``GEOS_THROW*`` macros for the same reasons as ``GEOS_ERROR*`` (**unrecoverable state**), and when:

- You want to **add context** at higher call stack levels,
- You need to **propagate detailed error information** upward
- If you use custom exception class, you need to document them here.

.. list-table:: Available exception types
   :header-rows: 1
   :widths: 35 65

  * - Error class
    - Used for
  * - ``std::runtime_error``
    - General runtime errors
  * - ``geos::InputError`` 
    - Errors in user input (XML, data files)
  * - ``geos::SimulationError``
    - Errors during simulation execution
  * - ``geos::BadTypeError`` 
    - Type conversion errors
  * - ``geos::NotAnError``
    - Non critical, volontary stopping state

.. dropdown:: Example: Throwing exceptions
   :icon: code

   .. code-block:: c++

      // Example 1: throw with context
      GEOS_THROW_IF( !file.is_open(),
                     GEOS_FMT("Cannot open file {}", filename),
                     InputError );

      // Example 2: rethrow with context info (because of unclear exception)
      try
      {
        tolerance = stod( inputParams[8] );
      }
      catch( const std::invalid_argument & e )
      {
        GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), 
                    InputError, getDataContext() );
      }

Never Catch and Continue
^^^^^^^^^^^^^^^^^^^^^^^^^

**Do not catch exceptions and continue execution**. Exceptions in GEOS indicate serious problems, they
do not serve as an alternative code path to consult a behaviour state.

**Rationale:**

- **Performance:** Exception handling is expensive, especially with stack unwinding,
- **Predictability:** Catching and continuing makes behaviour unpredictable and hard to debug,
- **Safety:** The error state may have corrupted simulation data, broke the code path at an unexpected place, and/or left the system in an invalid state.

.. dropdown:: Example: Exception handling anti-pattern
   :icon: code

   .. code-block:: c++

      // BAD - Do not catch to continue
      // Here, the "solver" should have a "isFailed" flag or a returned resulting state
      try
      {
        solver.solve();
      }
      catch( std::exception const & e ) // <- Note: this exception could have been throw for unexpected reason
      {
        GEOS_LOG("Solver failed, trying alternative");
        alternativeSolver.solve();
      }

      // GOOD - Manage your code paths manually
      if(! solver.solve() )
      {
        GEOS_LOG("Solver failed, trying alternative");
        alternativeSolver.solve();
      }

Warnings (GEOS_WARNING*)
------------------------

Why Use Warnings?
^^^^^^^^^^^^^^^^^

Use ``GEOS_WARNING*`` macros when:

- An issue is detected but **simulation can continue**
- Simulation may be **affected but not completly invalid**
- User should be notified of a **potential problem**

.. dropdown:: Example: Warning usage
   :icon: code

   .. code-block:: c++

      // Example: suboptimal but valid configuration
      GEOS_WARNING_IF( dt > recommendedDt,
                       GEOS_FMT("Time step {} exceeds recommended value {}. This may affect solution accuracy.",
                                dt, recommendedDt),
                       getDataContext() );

Errors & Warning in RAJA Kernels
--------------------------------

Same rule as logging: **Error and warning output inside RAJA kernels only for DEBUG purposes** and must never appear on the develop branch.

- GPU kernel errors cause immediate kernel termination, adding potentially costly branches,
- Error handling on device has significant performance & cache impact,
- GPU error handling may be unsupported depending on the platform / device,
- Can cause deadlocks in parallel execution.

.. dropdown:: Example: Kernel error handling (debug only)
   :icon: code

   .. code-block:: c++

      // On GPU, errors will cause kernel abort
      forAll< parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        #if !defined( NDEBUG )
        // This will trap on GPU - use sparingly for debugging only
        GEOS_ERROR_IF( badCondition, "Bad condition detected" );
        #endif
      });

Contextualization with DataContext
-----------------------------------

GEOS provides a abstract contextualization mechanism using ``DataContext`` to provide detailed error
information on the source of the data, which can be:

- the file & line of the simulation XML file,
- if no source exists, ``Group`` / ``Wrapper`` path in the data repository.

**Add ``DataContext`` when** Errors occur related to a **data repository component**:

- a ``Wrapper`` when validation fails for a given user input setting,
- a ``Group`` when the error is implied by the instance state,
- multiple data-contexts when the error is due to multiple objects / settings (first = more important).

Use ``getDataContext()`` or ``getWrapperDataContext()`` as last parameters of the logging macros to give context.

.. dropdown:: Example: Using DataContext
   :icon: code

   .. code-block:: c++

      // Example 1, using the getWrapperDataContext() because the error is
      // related to a precise input setting.
      GEOS_ERROR_IF( !meshBodies.hasGroup( meshBodyName ),
                     GEOS_FMT( "MeshBody ({}) is specified in targetRegions, but does not exist.", 
                               meshBodyName ),
                     getWrapperDataContext( viewKeyStruct::targetRegionsString() ) );

      // Exemple 2, using getDataContext() because the error is related with the entire instance
      // state (a physics solver for instance)
      GEOS_ERROR_IF( convergenceFailed, 
                     "Failed to converge!",
                     getDataContext() );

      // Exemple 3, using multiple DataContext because the error is related with multiple objects
      GEOS_THROW_IF( m_isThermal && !isFluidModelThermal,
                     GEOS_FMT( "CompositionalMultiphaseBase {}: the thermal option is enabled in the solver, but the fluid model {} is incompatible with the thermal option",
                               getDataContext(), fluid.getDataContext() ),
                     InputError, getDataContext(), fluid.getDataContext() );

Parallelism
===========

RAJA Kernels / CHAI Memory
---------------------------

Do Not Launch Kernels from Within Kernels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Never launch RAJA kernels from within other RAJA kernels.** This causes portability issues and is not supported on GPU.

.. dropdown:: Example: Nested kernel anti-pattern
   :icon: code

   .. code-block:: c++

      // BAD - Do not nest kernels
      forAll<parallelDevicePolicy<>>( n, [=] GEOS_HOST_DEVICE ( localIndex i )
      {
        forAll<parallelDevicePolicy<>>( m, [=] GEOS_HOST_DEVICE ( localIndex j )
        {
          // ...
        });
      });
      
      // GOOD - Use nested loops or restructure algorithm
      forAll<parallelDevicePolicy<>>( n, [=] GEOS_HOST_DEVICE ( localIndex i )
      {
        for( localIndex j = 0; j < m; ++j )
        {
          // ...
        }
      });

CHAI Memory Management
^^^^^^^^^^^^^^^^^^^^^^

GEOS uses CHAI for **automatic kernel memory migration**. Follow these rules:

- **Use CHAI-managed arrays:** ``array1d``, ``array2d``, etc.
- **Call** ``.move(space)`` **to explicitly control memory location**
- **Never use manual CUDA calls** - RAJA has wrappers for device operations

.. dropdown:: Example: CHAI memory management
   :icon: code

   .. code-block:: c++

      array1d<real64> data( 1000 );
      
      // Explicitly move to device before kernel
      // (would have been done automatically when giving an LvArray view to a RAJA kernel)
      data.move( parallelDeviceMemorySpace );
      
      // Use in kernel
      auto dataView = data.toView();
      forAll<parallelDevicePolicy<>>( numElems, [=] GEOS_HOST_DEVICE ( localIndex i )
      {
        dataView[i] = computeValue( i );
      });
      
      // Move back to host if needed for CPU operations (costly)
      data.move( hostMemorySpace );

Use RAJA Reductions
^^^^^^^^^^^^^^^^^^^

Use RAJA reductions for summation, min, max, and other parallel operations:

.. dropdown:: Example: RAJA reductions
   :icon: code

   .. code-block:: c++

      RAJA::ReduceSum<parallelDeviceReduce, real64> totalSum( 0.0 );
      RAJA::ReduceMin<parallelDeviceReduce, real64> minValue( std::numeric_limits<real64>::max() );
      RAJA::ReduceMax<parallelDeviceReduce, real64> maxValue( std::numeric_limits<real64>::lowest() );
      
      forAll<parallelDevicePolicy<>>( n, [=] GEOS_HOST_DEVICE ( localIndex i )
      {
        totalSum += values[i];
        minValue.min( values[i] );
        maxValue.max( values[i] );
      });
      
      real64 sum = totalSum.get();
      real64 min = minValue.get();
      real64 max = maxValue.get();

MPI Communication
-----------------

Ensure All Ranks Participate in Collective Operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**All MPI ranks must participate in collective operations.** Beware of branch deadlocks.

.. dropdown:: Example: MPI collective operations
   :icon: code

   .. code-block:: c++

      // BAD - only rank 0 calls MPI_AllReduce() - DEADLOCK!
      if( MpiWrapper::commRank() == 0 )
      {
        MpiWrapper::allReduce( ... );
      }

      // BAD - some ranks may have this condition to false - not every ranks call MPI_Allreduce - DEADLOCK!
      if( localWellElementsCount > 0 )
      {
        MpiWrapper::allReduce( ... ); // WRONG - other ranks waiting!
      }

      // GOOD - All ranks participate, even if they may not have influence in the result
      real64 localSum = computeLocalSum();
      real64 globalSum = MpiWrapper::allReduce( localSum );
      
      // GOOD - Conditional execution AFTER collective operation
      real64 globalSum = MpiWrapper::allReduce( localSum );
      if( MpiWrapper::commRank() == 0 )
      {
        GEOS_LOG_RANK_0( GEOS_FMT("Global sum = {}", globalSum) );
      }

Always Use MpiWrapper
^^^^^^^^^^^^^^^^^^^^^

**Always use ``MpiWrapper`` instead of raw MPI calls.**
Refer to ``MpiWrapper.hpp`` for more.

.. dropdown:: Example: Using MpiWrapper
   :icon: code

   .. code-block:: c++

      // BAD - Raw MPI calls
      int rank;
      MPI_Comm_rank( MPI_COMM_WORLD, &rank );
      
      // GOOD - Use MpiWrapper
      int rank = MpiWrapper::commRank();
      int size = MpiWrapper::commSize();
      MpiWrapper::barrier();
      
      // Common patterns
      real64 globalMax = MpiWrapper::max( localValue );
      real64 globalSum = MpiWrapper::sum( localValue );
      MpiWrapper::allGather( sendBuf, recvBuf, count );

Communication Batching
^^^^^^^^^^^^^^^^^^^^^^^

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


Validation / Precision
======================

Computation validation
-----------------------

Use Tolerance-Based Comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Always consider proper tolerance** for floating-point numbers comparisons, taking into account rounding errors, even for extreme values.

.. dropdown:: Example: Correct float comparison
   :icon: code

   .. code-block:: c++

      // BAD - Direct comparison
      if( value == 1.0 ) { ... }
      if( a == b ) { ... }

      // GOOD - Absolute tolerance
      real64 const absTol = 1.0e-12;
      if( std::abs(value) < absTol ) { ... }

      // GOOD - Relative tolerance
      real64 const relTol = 1.0e-8;
      real64 const scale = std::max( std::abs(a), std::abs(b) );
      if( std::abs(a - b) < relTol * scale ) { ... }

Division Safety, NaN/Inf values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Always verify denominators are not zero or near-zero before division**.

- In General, **we should not make any computation that could result in NaN/Inf values**.

.. dropdown:: Example: Division Safety
   :icon: code

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
^^^^^^^^^^^^^^^^^^^^

Overflow leads to undefined behavior and memory corruption.

**Validate that operations won't exceed type limits**, especially for index calculations.

Use LvArray Math Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``LvArray::math`` namespace provides safe numerical utilities, avoid LvArray math utilities instead of platform specific math function.

.. dropdown:: Example: LvArray math utilities
   :icon: code

   .. code-block:: c++

      #include "LvArray/src/math.hpp"

      real64 const absValue = LvArray::math::abs( value );
      real64 const maxValue = LvArray::math::max( a, b );
      real64 const minValue = LvArray::math::min( a, b );
      
Data repository validation
-------------------------------

Validate in postInputInitialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Perform cross-parameter validation as soon as possible**, and when possible **only once **.
This can be done in the ``postInputInitialization()`` method, or later if needed, (i.e. in the 
``registerDataOnMesh()`` or during the simulation).

At the ``postInputInitialization()`` stage:

- All XML input has been parsed
- Default values have been applied
- Cross-references can be resolved
- User can receive clear error messages before simulation starts.

.. dropdown:: Example: Post-input validation
   :icon: code

   .. code-block:: c++

      void MyClass::postInputInitialization()
      {
        // Call parent first
        Base::postInputInitialization();
        
        // Validate individual parameters
        GEOS_THROW_IF_LT( m_timeStep, 0.0,
                          "timeStep must be positive",
                          InputError, getDataContext() );
        
        // Validate cross-parameter constraints
        GEOS_THROW_IF( m_endTime <= m_startTime,
                       "endTime must be greater than startTime",
                       InputError, getDataContext() );
        
        // Validate against other components
        GEOS_THROW_IF( !hasRequiredSolver(),
                       "Required solver not found in configuration",
                       InputError, getDataContext() );
      }

Additional Validation Guidelines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Use ``setInputFlag()`` to mark parameters as ``REQUIRED`` or ``OPTIONAL``
- Use ``setApplyDefaultValue()`` to provide sensible defaults
- Use ``setRTTypeName()`` for runtime type validation (e.g., ``rtTypes::CustomTypes::positive``)
- Document valid ranges in ``setDescription()``

Performance
===========

Profiling
---------

Profile When Affecting Kernel behaviour
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When modifying performance-critical code (kernels, assembly loops, solvers):

1. **Use Caliper profiling** on representative test cases before changes
2. **Measure performance** after changes
3. **Document performance impact** in pull request description

.. dropdown:: Example: Performance profiling workflow
   :icon: code

   .. code-block:: bash

      # Profile before changes
      geos -i test_case.xml --caliper-output=baseline.cali
      
      # Make code changes
      
      # Profile after changes
      geos -i test_case.xml --caliper-output=modified.cali
      
      # Compare results
      cali-query -q "SELECT time.duration WHERE annotation='myKernel'" baseline.cali
      cali-query -q "SELECT time.duration WHERE annotation='myKernel'" modified.cali

Caliper Integration
^^^^^^^^^^^^^^^^^^^^

Use ``GEOS_MARK_FUNCTION`` and ``Timer`` for performance tracking for the main computation functions / scopes.

.. dropdown:: Example: Performance instrumentation
   :icon: code

   .. code-block:: c++

      void PhysicsSolverBase::solverStep( ... )
      {
        GEOS_MARK_FUNCTION; // Caliper tracking for entire function

        {
          GEOS_CALIPER_MARK_SCOPE("assemble system");
          Timer timer( m_timers.get_inserted( "assembly" ) );
          assembleSystem();
        }

        {
          GEOS_CALIPER_MARK_SCOPE("linear solve");
          Timer timer( m_timers.get_inserted( "linear solver" ) );
          linearSolve();
        }
      }

Speed Optimization Rules
-----------------------------

- **Hoist Loop Invariants**, move computations that don't change during iterations outside the loop.
- When it does not critically affect the code architecture and clarity, **Fuse multiple related kernels to reduce memory traffic and launch overhead** (i.e., statistics kernels process all physics field at once).
- **Optimize Memory Access for Cache and Coalescing**. Access memory sequentially and ensure coalesced access, especially on GPUs.
- **Minimize Host-Device Transfers**. Keep data on the appropriate memory space and minimize transfers.

Memory Management Rules
---------------------------------

Generalities
^^^^^^^^^^^^^^^^^^^^

- **Minimize dynamic memory allocation** as much as possible, particularly in performance-critical code,
- For frequently allocated/deallocated objects, **consider memory pools**,
- For containers with known size, **reserve capacity upfront**.

Be Const / Constexpr
^^^^^^^^^^^^^^^^^^^^

**Use ``const`` and ``constexpr`` when possible**, enabling:

- **Compiler optimization,** enables better code generation by giving constraints to the code,
- **Code safety,** prevents accidental modification for constant contexts,
- **Show code intention,** make code more readable.

Also, **Mark methods const if the method is not designed to modify** the object state.

Value / Const Function Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Pass large and non-trivially copyable objects by const reference**:
- A "large" size can be defined as "more that 16 bytes, ", which is 2 pointers, 2 integer / index, or 2 double,
- An object is non-trivially copyable objects when it needs to perform sub-allocation when being copied,

.. dropdown:: Value / Const Reference function parameters practices
   :icon: code

   .. code-block:: c++

      // 16 bytes (string_view are equivalent to 2 constant pointers): PASS BY VALUE
      static constexpr string_view str0 = "Hello";
      string const str1 = getTimeStr(); // constant string are passed as string_view, beware of lifetime!
      void func( string_view param );

      arrayView2d< int > myDataView;

      // 16 bytes: PASS BY VALUE
      stdArray< double, 2 > vec2D;
      void func( vec2D param );

      // 12 to 16 bytes (globalIndex is 8, integer is 4 to 8): PASS BY VALUE
      struct SmallStruct
      {
        globalIndex id;
        integer count;
      };
      void func( SmallStruct param );

      // 16 to 20 bytes (globalIndex is 8, integer is 4 to 8): PASS BY REFERENCE
      struct BigStruct
      {
        globalIndex id;
        integer countA;
        integer countB;
      };
      void func( BigStruct const & param );

      // Does dynamic allocation: PASS BY REFERENCE
      map< int, long > myMap;
      stdVector< int > myList { 123 };
      void func( map< int, long > const & myMap, 
                 stdVector< int > const & myList );

      // LvArray types parameter practices depends on the intent
      array2d< int > myConstantData;
      array2d< int > myMutableData;,
      // passing as non-view means that we will potencially resize the arrays: PASS BY REFERENCE
      void arrayResizer( array2d< const int > & myConstantData
                         array2d< int > & myMutableData );
      // passing as view means that we may just affect the data, only if non-const template: PASS BY VALUE
      void arrayProcess( arrayView2d< const int > myConstantData
                         arrayView2d< int > myMutableData );

      // 16 bytes superficially... but does dynamic allocation: PASS BY REFERENCE
      struct DynAllocStruct
      {
        string const a;
      };
      void func( DynAllocStruct const & param );

      // Any Group subclass is large: PASS BY REFERENCE
      void modify( DomainPartition & domain );

Pointer / Reference Function Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When passing nullable objects as parameter, pointers can be passed.

**Pointer parameters must always be null-checked when dereferenced**.
If applying this rule produces repetitive code, consider using a const reference.

This rule imply the following:
- **A reference may never be null**,
- In a given method, **``this`` may never be null**.

Passing a pointer or reference parameter show the intention of modifying the underlying object.
The difference of the two practices is that passing a pointer means that we consider that the 
underlying object can be null.

Constexpr for Compile-Time Constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Use ``constexpr`` for values known at compile time**.

The rule is not absolute: when the impact is not significative, and the code needs is getting really unclear, the rule can be ignored.

.. dropdown:: Example: Constexpr usage
   :icon: code

   .. code-block:: c++

      // Compile-time constants
      // A rule of thumb is that oftenly any time a value is constexpr, it can also be static.
      static constexpr localIndex numDimensions = 3;
      static constexpr real64 tolerance = 1.0e-10;
      
      // Constexpr functions (evaluated at compile time when possible)
      constexpr localIndex computeArraySize( localIndex n )
      { return n * n + 2 * n + 1; }

Provide Views to Arrays
^^^^^^^^^^^^^^^^^^^^^^^^

- When possible, prefer provide views to arrays (to const data if possible).
- This **must** be done especially for RAJA kernels, and views must be **captured by value**.

The rule is generalizable to ``string_view`` for strings, but not applicable in GPU context.

Why Use Views?

- **No memory allocation:** Views are lightweight references
- **Mutability correctness:** Can provide ``const`` read-only access to inner data
- **GPU compatibility:** Views work seamlessly on device

.. dropdown:: Example: Views for arrays
   :icon: code

   .. code-block:: c++

      // BAD - Creates a copy
      array1d< real64 > copy = originalArray;
      
      // GOOD - Creates a view (lightweight)
      arrayView1d< real64 > view = originalArray.toView();
      
      // GOOD - Const view for read-only access
      arrayView1d< real64 const > constView = originalArray.toViewConst();

View Lifetime Management
------------------------

**Never Outlive Parent Arrays**
Dangling views cause segmentation faults and undefined behavior, that can be particularly hard to diagnose.
The rule is applicable to ``arrayView*`` and ``string_view``.

.. dropdown:: Example: Lifetime Management
   :icon: code

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

General Architecture
====================

Avoid Coupling
--------------

Minimize coupling between components when possible **without compromising memory efficiency or execution speed**.

Principles
^^^^^^^^^^

- **Loose coupling:** Components should depend on interfaces, not concrete implementations,
- **No circular dependancies:** Consider the existing GEOS dependancies to not make components co-dependant (headers inclusion, packages referencing in ``CMakeLists.txt``, avoid tightly coupled objects),
- **Dependency injection:** Public components should receive their dependencies from external sources. Pass required dependencies using intermediate types instead of direct implementation types, using lambda, templates, and minimal interfaces, relying on **lambda**, **templates** and **minimal interfaces** (loose coupling, testability),
- **Performance exceptions:** Tight coupling is acceptable when required for performance,
- **Minimize header inclusions and dependancies**.

.. dropdown:: Example: Reducing coupling
   :icon: code

   .. code-block:: c++

      // BAD - Tight coupling to concrete class / implementation
      class SolverA
      {
        void registerSomeMeshData( VTKMesh & mesh );
        void solveHypre( HypreInterface * solver, HypreMatrix * matrix );
      };
      
      // GOOD - Depends on interface
      class SolverB
      {
        void registerSomeMeshData( MeshLevel & mesh ); // More general interface
        void solve( LAInterface * solver, MatrixBase * matrix );
      };
      
      // Performance-critical tight coupling
      template< typename FIELD_T >
      class Kernel
      {
        // Direct access to specific data layout for performance
        void compute( FIELD_T const & data );
      };
      
      template<> Kernel::compute( arrayView1d< double > & field )
      { /* template specialization */ }

      template<> Kernel::compute( arrayView2d< double > & field );
      { /* template specialization */ }

Avoid Globals and Mutable State
--------------------------------

**Minimize mutable global states when possible.**
Prefer passing context explicitly.

Why to avoid globals:

- **Thread safety:** Globals can cause race conditions in parallel code
- **Testability:** Hard to test code that depends on global state
- **Predictability:** Global state makes code behaviour harder to understand
- **Encapsulation:** Violates modularity principles

.. dropdown:: Example: Avoiding global state
   :icon: code

   .. code-block:: c++

      // BAD - Global mutable state
      static real64 g_tolerance = 1.0e-6;
      
      void solve()
      {
        if( residual < g_tolerance ) { ... } // Depends on global
      }
      
      // GOOD - Pass as parameter
      void solve( real64 tolerance )
      {
        if( residual < tolerance ) { ... }
      }
      
      // GOOD - Member variable
      class Solver
      {
        real64 const m_tolerance;
      public:
        Solver( real64 tolerance ) : m_tolerance( 1.0e-6 ) {}
        
        void solve()
        {
          if( residual < m_tolerance ) { ... }
        }
      };

Acceptable Global Usage:

- **Library wrappers** (MPI wrapper, system-level resources, dependancies interface)
- **Read-only configuration** (immutable after initialization)
- **Performance counters** (for profiling)

Magic values, ``Group`` paths
-------------------------------

**Avoid magic values**:
- **arbitrary values should not be written more than once**, define constants, consider using or extending ``PhysicsConstants.hpp`` / ``Units.hpp``,
- Because of architecture considerations, **Avoid using full raw ``Group`` path**, especially when the path has parts unrelated to the component consulting it -> prefer as relative paths as possible,
- ``Group`` / ``Wrapper`` name in the data repository is considered as a magic value, **use the ``getCatalogName()`` / ``getName()`` getters**,
- **Prefer to let appear the calculus of constants** rather than writing its value directly without explaination (constexpr has no runtime cost).

======================
Physics-Specific Rules
======================

Unit Consistency
----------------

**Verify physical units**, and **document any non-S.I. units** (variable name, comments).
**Non-S.I. units are only for internal computation use.**

Physical Bounds Checking
------------------------

Unphysical values indicate errors and can cause solver failures.
**Validate that state variables remain within physically meaningful bounds**.

If a value is not strictly disallowed but does not seem possible for the model, **show a warning to the user that he can disable**.


  .. code-block:: cpp

    class ConstitutiveModel
    {
      void updateState( real64 const pressure,
                       real64 const temperature,
                       real64 const porosity )
      {
        // Check physical bounds
        GEOS_ERROR_IF( pressure < 1e-5, 
                       GEOS_FMT( "Pressure {:.2e} Pa below cavitation limit", pressure ) );
        
        GEOS_ERROR_IF( temperature < 0.0,
                       GEOS_FMT( "Absolute temperature {:.2f} K is negative", temperature ) );
        
        GEOS_ERROR_IF( porosity < 0.0 || porosity > 1.0,
                       GEOS_FMT( "Porosity {:.3f} outside [0,1]", porosity ) );
        
        // Warn about unusual but valid values, allow user to disable warning
        GEOS_WARNING_IF( isHighTemperatureWarningEnabled && temperature > 1000.0,
                        GEOS_FMT( "High temperature {:.0f} K may be outside model validity", 
                                  temperature ) );
      }
    };
