.. _CodingPractices:


###############################################################################
Coding Practice Rules
###############################################################################

This section covers general C++ coding practices in place for the GEOS project.

All contributors should familiarize themselves with these rules to maintain code quality, consistency, and reliability across the GEOS codebase.

.. note::
   This document uses collapsible code blocks. Click on "Show/Hide Code" to expand examples.

Memory Managment
================

Allocations
-----------

**Minimize dynamic memory allocation** as much as possible, particularly in performance-critical code,

- For frequently allocated/deallocated objects, **consider memory pools**,
- For containers with known size, **reserve capacity upfront**.

As they are intended to low-level allocator codes, **avoid** ``new`` **/** ``delete`` **and** ``malloc`` **/** ``free`` **as much as possible**:

- When applicable, **Prefer composition** over pointers,
- **Use RAII principles** as much as possible,
- **Prefer** ``std::unique_ptr`` **over** ``std::shared_ptr``.

Pointer / Reference Function Arguments
--------------------------------------

**Pass object by reference rather than passing pointer types.**

When passing nullable objects as parameter, pointers can be passed.

**Pointer parameters must always be null-checked when dereferenced**.
If applying this rule produces repetitive code, consider using a const reference.

This rule imply the following:

- **A reference may never be null**,
- In a given method, ``this`` **may never be null**.

.. dropdown:: Why?
  :icon: info

  Passing a pointer or reference parameter show the intention of modifying the underlying object.
  The difference of the two practices is that passing a pointer means that we consider that the 
  underlying object can be null.

Provide Views to Arrays
-----------------------

- When possible, **prefer provide views to arrays** (to const data if possible).
- **views must be passed / captured by value** for function / lambdas.

The rule is generalizable to ``string_view`` for strings, but not applicable in GPU context.

.. dropdown:: Why Use Views?
  :icon: info

  - **No memory allocation:** Views are lightweight references
  - **Mutability correctness:** Depending on the view type, can provide ``const`` read-only access to inner data
  - **GPU compatibility:** LvArray views work seamlessly on device

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

.. dropdown:: Why?
  :icon: info

  Dangling views cause segmentation faults and undefined behavior, that can be particularly hard to diagnose.
  The rule is applicable to ``arrayView*`` types and ``string_view``.

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
    public:
      arrayView1d< real64 > getView() { return m_data.toView(); }
      arrayView1d< real64 const > getViewConst() const { return m_data.toViewConst(); }
    private:
      array1d< real64 > m_data;
    };

Value / Const Function Arguments
--------------------------------

**Pass large and non-trivially copyable objects by const reference**:

- A "large" size can be defined as "more that 16 bytes, ", which is 2 pointers, 2 integer / index, or 2 double,
- An object is non-trivially copyable objects when it needs to perform sub-allocation when being copied,

.. dropdown:: Example: Value / Const Reference function parameters practices
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

Avoid Undesired Mutability
--------------------------

**Enforce** ``const`` **and** ``constexpr`` **when possible**

.. dropdown:: Why?
  :icon: info

  ``const`` and ``constexpr`` declaration enables:

  - **Enables compiler optimization**, improving code generation with explicit code constraints,
  - **Improves code safety,** preventing accidental modification for constant contexts,
  - **Show code intention,** making code clearer.

Also, **mark methods const if the method is not designed to modify** the object state.  

Use of ``mutable``
------------------

``mutable`` **keyword usage is to avoid, excepted for specific patterns.**

Acceptable uses:

- **Thread synchronization primitives** (``std::mutex``, ``std::atomic`` for counters),
- **Caching** when it preserves logical constness, is thread-safe and clearly documented.

Avoid ``mutable`` for:

- Working around design issues (consider redesigning instead, i.e. splitting data structures),
- GPU-accessible data (creates host/device LvArray synchronization issues).

.. dropdown:: Why mutable is discouraged?
  :icon: info

  - **Breaks const methods contract:** Makes reasoning about code harder
  - **Thread-safety confusion:** ``const`` methods may not be thread-safe
  - **GPU incompatibility:** Can cause memory coherence issues with LvArray views

Constexpr for Compile-Time Constants
------------------------------------

**Use** ``constexpr`` **for values known at compile time**.

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

Validation / Precision
======================

Use Tolerance-Based Comparisons
-------------------------------

**Always consider proper tolerance** for floating-point numbers comparisons, taking into account rounding errors, even for extreme values.

.. dropdown:: Example: Correct float comparison
  :icon: code

  .. code-block:: c++

      // BAD - Direct comparison
      if( value == 1.0 ) { ... }
      if( a == b ) { ... }

      // GOOD - Absolute tolerance
      real64 const absTol = 1.0e-12;
      if( LvArray::math::abs(value) < absTol ) { ... }

      // GOOD - Relative tolerance
      real64 const relTol = 1.0e-8;
      real64 const scale = LvArray::math::max( LvArray::math::abs(a), LvArray::math::abs(b) );
      if( LvArray::math::abs(a - b) < relTol * scale ) { ... }

Division Safety, NaN/Inf Values
-------------------------------

- **Always verify that denominators are not zero or near-zero before a division**.  
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
-------------------

Overflow leads to undefined behavior and memory corruption. **Validate that operations won't exceed type / container limits**, especially for index calculations.

Performance
===========

Speed Optimization Rules
-----------------------------

- **Hoist Loop Invariants**, move computations that don't change during iterations outside the loop.
- When it does not critically affect the code architecture and clarity, **fuse multiple related kernels to reduce memory traffic and launch overhead** (i.e., statistics kernels process all physics field at once).  
- **Optimize Memory Access for Cache and Coalescing**. Access memory sequentially and ensure coalesced access, especially on GPUs.
- **Minimize Host-Device Transfers**. Keep data on the appropriate memory space and minimize transfers.

General Architecture
====================

Avoid Coupling
--------------

Minimize coupling between components when possible **without compromising memory efficiency or execution speed**.

Principles:

- **Loose coupling:** Components should depend on interfaces, not concrete implementations,
- **No circular dependencies:** Consider the existing GEOS dependencies to not make components co-dependent (headers inclusion, packages referencing in ``CMakeLists.txt``, avoid tightly coupled objects),  
- **Dependency injection:** Public components should receive their dependencies from external sources. Pass required dependencies using intermediate types instead of direct implementation types, relying on **lambda**, **templates** and **minimal interfaces** (loose coupling, testability),
- **Performance exceptions:** Tight coupling is acceptable when required for performance,
- **Minimize header inclusions and dependencies**.

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

.. dropdown:: Why to avoid globals?
  :icon: info

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

Magic values
------------

**Avoid magic values**:

- **arbitrary values should not be written more than once**, define constants, consider using or extending ``PhysicsConstants.hpp`` / ``Units.hpp``,
- **Prefer to let appear the calculus of constants** rather than writing its value directly without explaination (constexpr has no runtime cost).
