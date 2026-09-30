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

/** @file VTKRefinementFields.hpp */
#ifndef GEOS_VTK_REFINEMENT_FIELDS_HPP
#define GEOS_VTK_REFINEMENT_FIELDS_HPP

#include "VTKRefinementTopology.hpp"
#include <vtkSmartPointer.h>

#include <set>
#include <string>

class vtkCellData;
class vtkPointData;
class vtkDataArray;
class vtkAbstractArray;

namespace geos::vtk::refinement
{

enum class PointTransferPolicy
{
  continuous,
  equal,
  nodeSet
};

/** Initialization semantics belong to GEOS, independently of array storage.
 * Floating point fields default to continuous interpolation; integral fields
 * require equal support values. Cell values are intensive unless declared
 * extensive. Relational arrays are supplied by their specialized handlers.
 */
struct TransferPolicies
{
  std::map< std::string, PointTransferPolicy > pointArrays;
  std::set< std::string > extensiveCellArrays;
  std::set< std::string > excludedPointArrays{ "collocated_nodes" };
  std::set< std::string > excludedCellArrays;
};

vtkSmartPointer< vtkPointData > transferPointData( vtkPointData & input, PointRegistry const & points, TransferPolicies const & policies );
vtkSmartPointer< vtkCellData > transferCellData( vtkCellData & input, vtkIdType parentCount, Connectivity const & parents,
                                                 std::vector< double > const & fractions, TransferPolicies const & policies );

/** Canonical typed tuples for authoritative shared-point records.
 * Full array schema, component names, and active roles are checked on install.
 * Integral values never pass through VTK's double-valued tuple interface.
 */
class PointFieldLayout
{
public:
  explicit PointFieldLayout( vtkPointData & data );
  std::vector< unsigned char > pack( vtkIdType point ) const;
  void install( vtkIdType point, std::vector< unsigned char > const & tuple ) const;

private:
  vtkSmartPointer< vtkPointData > m_data;
  std::vector< vtkDataArray * > m_arrays;
  std::vector< unsigned char > m_schema;
};

/** Typed surface-cell tuples, including string and bit arrays.
 * Authoritative replica transfer does not convert integral values to double.
 */
class CellFieldLayout
{
public:
  explicit CellFieldLayout( vtkCellData & data );
  std::vector< unsigned char > pack( vtkIdType cell ) const;
  void install( vtkIdType cell, std::vector< unsigned char > const & tuple ) const;
private:
  vtkSmartPointer< vtkCellData > m_data;
  std::vector< vtkAbstractArray * > m_arrays;
  std::vector< unsigned char > m_schema;
};

} // namespace geos::vtk::refinement
#endif
