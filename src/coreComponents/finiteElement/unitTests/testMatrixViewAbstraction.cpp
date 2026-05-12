/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMatrixViewAbstraction.cpp
 *
 * @brief Verifies that the kernel-hierarchy @c MATRIX_VIEW template parameter
 *        is genuinely generic — i.e. that any type satisfying the documented
 *        contract can be used in place of @c DefaultGlobalMatrixView, not
 *        only @c LvArray::CRSMatrixView.
 *
 * The MATRIX_VIEW contract (see @c ImplicitKernelBase doxygen) requires:
 *   - default-constructible (used by SparsityKernelBase in sparsity-only mode)
 *   - copy-constructible and safe to capture by value into device lambdas
 *   - @c numRows() , @c numColumns() , @c numNonZeros() , @c numNonZeros(row)
 *   - @c getColumns(row)  -> indexable, yields globalIndex
 *   - @c getEntries(row)  -> indexable, yields a mutable real64-like
 *   - @c addToRow< AtomicPolicy >( row, cols, vals, nCols ) const
 *   - @c addToRowBinarySearchUnsorted< AtomicPolicy >( row, cols, vals, nCols ) const
 *
 * This test pins that contract by (a) defining a minimal @c MockMatrixView
 * type whose only relationship to @c CRSMatrixView is conformance to the
 * above interface, (b) declaring a @c hasMatrixViewInterface SFINAE trait
 * and asserting both @c MockMatrixView and the production
 * @c DefaultGlobalMatrixView satisfy it, (c) forcing a production
 * @c LaplaceFEMKernel @c complete() instantiation against the mock, and
 * (d) exercising the mock functionally to confirm the contract is internally
 * consistent.
 *
 * If a future change to a kernel reaches past the documented interface
 * (e.g. into a CRSMatrix-internal layout assumption), the mock-backed
 * compile path here will fail.
 *
 * NOTE: this test deliberately does NOT run a full physics-kernel launch
 *       against MockMatrixView — that path requires mesh and constitutive
 *       fixtures and belongs in an integration test.  The compile-time
 *       LaplaceFEMKernel instantiation proves the leaf kernel stays within
 *       the documented interface; the functional smoke test confirms the mock
 *       behaves like a real sparse matrix view.
 */

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/logger/Logger.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "constitutive/NullModel.hpp"
#include "finiteElement/elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "physicsSolvers/simplePDE/LaplaceFEMKernels.hpp"
#include "gtest/gtest.h"

#include <type_traits>
#include <vector>

using namespace geos;

//==============================================================================
//  MockMatrixView
//
//  A host-only, hand-rolled CSR-backed implementation of the kernel-facing
//  MATRIX_VIEW contract.  Storage is borrowed (non-owning) — the test owns
//  the backing std::vector<>'s and constructs a view over them.  This mirrors
//  the relationship CRSMatrixView has to CRSMatrix.
//==============================================================================
class MockMatrixView
{
public:

  /// Default constructor — required by SparsityKernelBase's MATRIX_VIEW().
  MockMatrixView() = default;

  /**
   * @brief Construct a view over externally-owned CSR storage.
   * @param nRows     number of local rows
   * @param nCols     number of global columns (informational)
   * @param offsets   row offsets, length nRows+1
   * @param sizes     per-row populated counts, length nRows (sizes[i] <= offsets[i+1]-offsets[i])
   * @param cols      column-index storage, length offsets[nRows]
   * @param vals      value storage, length offsets[nRows]
   */
  MockMatrixView( localIndex const nRows,
                  localIndex const nCols,
                  localIndex const * const offsets,
                  localIndex * const sizes,
                  globalIndex * const cols,
                  real64 * const vals )
    : m_numRows( nRows )
    , m_numCols( nCols )
    , m_offsets( offsets )
    , m_sizes( sizes )
    , m_cols( cols )
    , m_vals( vals )
  {}

  // Trivially-copyable contract is enforced via static_assert below; let
  // the compiler synthesize the copy/move members.

  // --- attribute querying ---------------------------------------------------

  localIndex numRows()    const { return m_numRows; }
  localIndex numColumns() const { return m_numCols; }
  localIndex numNonZeros() const { return m_offsets == nullptr ? 0 : m_offsets[ m_numRows ]; }
  localIndex numNonZeros( localIndex const row ) const { return m_sizes[ row ]; }

  // --- row access -----------------------------------------------------------

  arraySlice1d< globalIndex const > getColumns( localIndex const row ) const
  {
    return arraySlice1d< globalIndex const >( m_cols + m_offsets[ row ],
                                              &m_sizes[ row ],
                                              nullptr );
  }

  arraySlice1d< real64 > getEntries( localIndex const row ) const
  {
    return arraySlice1d< real64 >( m_vals + m_offsets[ row ],
                                   &m_sizes[ row ],
                                   nullptr );
  }

  // --- atomic adds (host-only; AtomicPolicy is required for signature
  //     compatibility but is not honored as a real atomic in this mock) ----

  template< typename AtomicPolicy >
  void addToRow( localIndex const row,
                 globalIndex const * const cols,
                 real64 const * const vals,
                 localIndex const nCols ) const
  {
    addToRowImpl( row, cols, vals, nCols );
  }

  template< typename AtomicPolicy >
  void addToRowBinarySearchUnsorted( localIndex const row,
                                     globalIndex const * const cols,
                                     real64 const * const vals,
                                     localIndex const nCols ) const
  {
    addToRowImpl( row, cols, vals, nCols );
  }

private:

  // Naive O(nCols * nnz) linear search; correctness, not perf.
  void addToRowImpl( localIndex const row,
                     globalIndex const * const cols,
                     real64 const * const vals,
                     localIndex const nCols ) const
  {
    localIndex const nnz = m_sizes[ row ];
    globalIndex const * const rowCols = m_cols + m_offsets[ row ];
    real64 * const           rowVals = m_vals + m_offsets[ row ];
    for( localIndex i = 0; i < nCols; ++i )
    {
      for( localIndex j = 0; j < nnz; ++j )
      {
        if( rowCols[ j ] == cols[ i ] )
        {
          rowVals[ j ] += vals[ i ];
          break;
        }
      }
    }
  }

  localIndex          m_numRows = 0;
  localIndex          m_numCols = 0;
  localIndex const *  m_offsets = nullptr;
  localIndex *        m_sizes   = nullptr;
  globalIndex *       m_cols    = nullptr;
  real64 *            m_vals    = nullptr;
};

//==============================================================================
//  hasMatrixViewInterface< T >
//
//  Compile-time SFINAE trait that succeeds iff T provides every operation
//  the kernel hierarchy invokes on a MATRIX_VIEW.  This trait IS the
//  documentation of the contract.
//==============================================================================
template< typename T, typename = void >
struct hasMatrixViewInterface : std::false_type {};

template< typename T >
struct hasMatrixViewInterface< T,
                               std::void_t<
                                 // Default-constructible (required by SparsityKernelBase)
                                 decltype( T() ),
                                 // Query methods.
                                 decltype( std::declval< T const >().numRows() ),
                                 decltype( std::declval< T const >().numColumns() ),
                                 decltype( std::declval< T const >().numNonZeros() ),
                                 decltype( std::declval< T const >().numNonZeros( std::declval< localIndex >() ) ),
                                 // Row views.
                                 decltype( std::declval< T const >().getColumns( std::declval< localIndex >() ) ),
                                 decltype( std::declval< T const >().getEntries( std::declval< localIndex >() ) ),
                                 decltype( std::declval< T const >().getColumns( std::declval< localIndex >() )[ std::declval< localIndex >() ] ),
                                 decltype( std::declval< T const >().getEntries( std::declval< localIndex >() )[ std::declval< localIndex >() ] = real64{} ),
                                 // addToRow< AtomicPolicy >( row, cols, vals, nCols )
                                 decltype( std::declval< T const >().template addToRow< parallelDeviceAtomic >(
                                             std::declval< localIndex >(),
                                             static_cast< globalIndex const * >( nullptr ),
                                             static_cast< real64 const * >( nullptr ),
                                             std::declval< localIndex >() ) ),
                                 // addToRowBinarySearchUnsorted< AtomicPolicy >( row, cols, vals, nCols )
                                 decltype( std::declval< T const >().template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                                             std::declval< localIndex >(),
                                             static_cast< globalIndex const * >( nullptr ),
                                             static_cast< real64 const * >( nullptr ),
                                             std::declval< localIndex >() ) )
                                 >
                               > : std::integral_constant< bool,
                                                           std::is_default_constructible< T >::value &&
                                                           std::is_copy_constructible< T >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().numRows() ), localIndex >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().numColumns() ), localIndex >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().numNonZeros() ), localIndex >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().numNonZeros( std::declval< localIndex >() ) ), localIndex >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().getColumns( std::declval< localIndex >() ) ), arraySlice1d< globalIndex const > >::value &&
                                                           std::is_convertible< decltype( std::declval< T const >().getEntries( std::declval< localIndex >() ) ), arraySlice1d< real64 > >::value
                                                           > {};

//==============================================================================
//  Pin the contract at compile time.
//
//  These static_asserts fire if either:
//    (a) MockMatrixView's interface drifts from the documented contract, or
//    (b) DefaultGlobalMatrixView (a.k.a. CRSMatrixView<real64,globalIndex const>)
//        ever loses a member the kernels depend on.
//
//  Any future MATRIX_VIEW type (e.g. HypreParCSRMatrixView) is required to
//  satisfy hasMatrixViewInterface<> in order to be usable in the kernel
//  templates.  Add a static_assert for it next to these.
//==============================================================================

static_assert( std::is_trivially_copyable< MockMatrixView >::value,
               "MockMatrixView must be trivially copyable so it can be captured "
               "by value into device lambdas without dynamic allocation or "
               "refcount side effects." );

static_assert( std::is_copy_constructible< DefaultGlobalMatrixView >::value,
               "DefaultGlobalMatrixView must remain copy-constructible because "
               "kernels store matrix views by value." );

static_assert( hasMatrixViewInterface< MockMatrixView >::value,
               "MockMatrixView does not satisfy the MATRIX_VIEW contract — "
               "update either the mock or the trait." );

static_assert( hasMatrixViewInterface< DefaultGlobalMatrixView >::value,
               "DefaultGlobalMatrixView (CRSMatrixView<real64,globalIndex const>) "
               "no longer satisfies the MATRIX_VIEW contract — kernels will not "
               "compile against it." );

using MockLaplaceKernel = LaplaceFEMKernel< CellElementSubRegion,
                                           constitutive::NullModel,
                                           finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2_impl,
                                           MockMatrixView >;

static_assert( std::is_constructible< MockLaplaceKernel,
                                      NodeManager const &,
                                      EdgeManager const &,
                                      FaceManager const &,
                                      localIndex const,
                                      CellElementSubRegion const &,
                                      finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2_impl const &,
                                      constitutive::NullModel &,
                                      arrayView1d< globalIndex const > const,
                                      globalIndex const,
                                      MockMatrixView const,
                                      arrayView1d< real64 > const,
                                      real64 const,
                                      string const >::value,
               "LaplaceFEMKernel must be constructible with MockMatrixView in "
               "the MATRIX_VIEW slot." );

static_assert( std::is_same<
                 decltype( std::declval< MockLaplaceKernel const >().complete(
                             std::declval< localIndex >(),
                             std::declval< MockLaplaceKernel::StackVariables & >() ) ),
                 real64 >::value,
               "LaplaceFEMKernel::complete must compile against MockMatrixView. "
               "A failure here means production kernel code reached outside the "
               "documented MATRIX_VIEW contract." );

// Force instantiation of LaplaceFEMKernel< ..., MockMatrixView >::complete().
//
// The two static_asserts above only verify signatures (decltype is unevaluated
// context, so they do not force function template body instantiation).  The
// body of complete() is where the kernel actually calls into MATRIX_VIEW —
// `m_matrix.numRows()` and `m_matrix.template addToRowBinarySearchUnsorted<
// parallelDeviceAtomic >( ... )`.  To catch a future kernel edit that reaches
// past the documented contract, we need the body to type-check against the
// mock.
//
// Taking the address of `instantiateComplete< MockLaplaceKernel >` at
// namespace scope is an odr-use in evaluated context and therefore forces the
// function template to be instantiated, which in turn instantiates
// MockLaplaceKernel::complete().  The variable itself is otherwise unused;
// [[maybe_unused]] suppresses the warning.
template< typename KERNEL >
real64 instantiateComplete( KERNEL const & kernel,
                            typename KERNEL::StackVariables & stack )
{
  return kernel.complete( 0, stack );
}

using MockLaplaceCompleteFunction = real64 (*)( MockLaplaceKernel const &,
                                                MockLaplaceKernel::StackVariables & );

[[maybe_unused]] MockLaplaceCompleteFunction const mockLaplaceCompleteFunction
  = &instantiateComplete< MockLaplaceKernel >;

//==============================================================================
//  Functional tests: confirm MockMatrixView behaves like a real sparse matrix
//  view when exercised through the contract interface.  These are smoke tests
//  for the mock itself — if they fail, the mock is broken, which would
//  invalidate the static_asserts as a proof.
//==============================================================================

namespace
{

/**
 * Build a small 3x3 sparse matrix with the following structure:
 *
 *      col:   0   1   2
 *   row 0: [  .   1   2 ]      (2 entries)
 *   row 1: [  3   4   . ]      (2 entries)
 *   row 2: [  5   .   6 ]      (2 entries)
 *
 * Returns a MockMatrixView and stashes the backing storage in the out-params.
 */
MockMatrixView build3x3( std::vector< localIndex >  & offsets,
                         std::vector< localIndex >  & sizes,
                         std::vector< globalIndex > & cols,
                         std::vector< real64 >      & vals )
{
  offsets = { 0, 2, 4, 6 };
  sizes   = { 2, 2, 2 };
  cols    = { 1, 2,   0, 1,   0, 2 };
  vals    = { 1, 2,   3, 4,   5, 6 };
  return MockMatrixView( 3, 3,
                         offsets.data(),
                         sizes.data(),
                         cols.data(),
                         vals.data() );
}

} // namespace

TEST( MatrixViewAbstraction, mock_attribute_queries )
{
  std::vector< localIndex >  offsets;
  std::vector< localIndex >  sizes;
  std::vector< globalIndex > cols;
  std::vector< real64 >      vals;
  MockMatrixView const m = build3x3( offsets, sizes, cols, vals );

  EXPECT_EQ( m.numRows(),    3 );
  EXPECT_EQ( m.numColumns(), 3 );
  EXPECT_EQ( m.numNonZeros(), 6 );
  EXPECT_EQ( m.numNonZeros( 0 ), 2 );
  EXPECT_EQ( m.numNonZeros( 1 ), 2 );
  EXPECT_EQ( m.numNonZeros( 2 ), 2 );
}

TEST( MatrixViewAbstraction, mock_row_access )
{
  std::vector< localIndex >  offsets;
  std::vector< localIndex >  sizes;
  std::vector< globalIndex > cols;
  std::vector< real64 >      vals;
  MockMatrixView const m = build3x3( offsets, sizes, cols, vals );

  // row 0 columns: { 1, 2 }, entries: { 1, 2 }
  auto const c0 = m.getColumns( 0 );
  auto const e0 = m.getEntries( 0 );
  EXPECT_EQ( c0[ 0 ], 1 );
  EXPECT_EQ( c0[ 1 ], 2 );
  EXPECT_DOUBLE_EQ( e0[ 0 ], 1.0 );
  EXPECT_DOUBLE_EQ( e0[ 1 ], 2.0 );

  // row 2 columns: { 0, 2 }, entries: { 5, 6 }
  auto const c2 = m.getColumns( 2 );
  auto const e2 = m.getEntries( 2 );
  EXPECT_EQ( c2[ 0 ], 0 );
  EXPECT_EQ( c2[ 1 ], 2 );
  EXPECT_DOUBLE_EQ( e2[ 0 ], 5.0 );
  EXPECT_DOUBLE_EQ( e2[ 1 ], 6.0 );
}

TEST( MatrixViewAbstraction, mock_addToRow )
{
  std::vector< localIndex >  offsets;
  std::vector< localIndex >  sizes;
  std::vector< globalIndex > cols;
  std::vector< real64 >      vals;
  MockMatrixView const m = build3x3( offsets, sizes, cols, vals );

  // Add { col 0: +10, col 1: +20 } into row 1.
  // Existing row 1: cols {0,1}, vals {3,4}.  After: vals {13,24}.
  globalIndex const addCols[] = { 0, 1 };
  real64      const addVals[] = { 10.0, 20.0 };
  m.template addToRow< serialAtomic >( 1, addCols, addVals, 2 );

  auto const e1 = m.getEntries( 1 );
  EXPECT_DOUBLE_EQ( e1[ 0 ], 13.0 );
  EXPECT_DOUBLE_EQ( e1[ 1 ], 24.0 );
}

TEST( MatrixViewAbstraction, mock_addToRowBinarySearchUnsorted )
{
  std::vector< localIndex >  offsets;
  std::vector< localIndex >  sizes;
  std::vector< globalIndex > cols;
  std::vector< real64 >      vals;
  MockMatrixView const m = build3x3( offsets, sizes, cols, vals );

  // Add unsorted { col 2: +100, col 0: +200 } into row 2.
  // Existing row 2: cols {0,2}, vals {5,6}.  After: vals {205, 106}.
  globalIndex const addCols[] = { 2, 0 };
  real64      const addVals[] = { 100.0, 200.0 };
  m.template addToRowBinarySearchUnsorted< serialAtomic >( 2, addCols, addVals, 2 );

  auto const e2 = m.getEntries( 2 );
  EXPECT_DOUBLE_EQ( e2[ 0 ], 205.0 );
  EXPECT_DOUBLE_EQ( e2[ 1 ], 106.0 );
}

TEST( MatrixViewAbstraction, mock_dirichlet_style_row_mutation )
{
  // Exercise the real FieldSpecificationOps Dirichlet path. This is the
  // hardest part of the contract to satisfy on non-CSR backings because it
  // requires mutable row-entry access paired with global row-column indices.
  std::vector< localIndex >  offsets;
  std::vector< localIndex >  sizes;
  std::vector< globalIndex > cols;
  std::vector< real64 >      vals;
  MockMatrixView const m = build3x3( offsets, sizes, cols, vals );

  globalIndex const targetDof = 1;
  localIndex  const localRow  = 1;   // row 1 has columns {0,1}; diag-of-interest is col 1
  real64 rhs = 0.0;
  FieldSpecificationEqual::SpecifyFieldValue< MockMatrixView >( targetDof,
                                                                0,
                                                                m,
                                                                rhs,
                                                                3.0,
                                                                2.0 );

  // Re-read and check.
  auto const eAfter = m.getEntries( localRow );
  EXPECT_DOUBLE_EQ( eAfter[ 0 ], 0.0 );    // col 0 -> zeroed
  EXPECT_DOUBLE_EQ( eAfter[ 1 ], 4.0 );    // col 1 -> diagonal
  EXPECT_DOUBLE_EQ( rhs, -4.0 );           // -diagonal * (bcValue - fieldValue)
}

TEST( MatrixViewAbstraction, default_constructed_view_is_safe_to_hold )
{
  // SparsityKernelBase passes a default-constructed MATRIX_VIEW() to the
  // ImplicitKernelBase base in sparsity-only mode.  No method is invoked on
  // it; the requirement is purely that the type be default-constructible
  // and copyable.  Verify both hold for the mock.
  MockMatrixView const m;
  MockMatrixView const copy = m;
  EXPECT_EQ( m.numRows(),    0 );
  EXPECT_EQ( copy.numRows(), 0 );
}

//==============================================================================
//  Negative compile-time check: a type missing required members must NOT
//  satisfy hasMatrixViewInterface<>.  This protects the trait from silently
//  matching anything.
//==============================================================================

namespace
{
struct IncompleteMatrixView
{
  IncompleteMatrixView() = default;
  localIndex numRows() const { return 0; }
  // Deliberately missing: numColumns, numNonZeros, getColumns, getEntries,
  // addToRow, addToRowBinarySearchUnsorted.
};
}

static_assert( !hasMatrixViewInterface< IncompleteMatrixView >::value,
               "hasMatrixViewInterface accepted a type that is missing required "
               "members — the trait is too permissive." );
