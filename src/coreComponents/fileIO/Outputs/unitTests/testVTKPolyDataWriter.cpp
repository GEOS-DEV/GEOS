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

#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"

#include "mesh/CellElementRegion.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ParticleRegion.hpp"
#include "mesh/ParticleSubRegion.hpp"

#include <conduit.hpp>
#include <gtest/gtest.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>
#include <set>

using namespace geos;
using namespace geos::dataRepository;

namespace
{

/**
 * Keep this list synchronized with the layouts in types::StandardArrays. Expanding
 * the operations directly in the test body avoids adding a templated test helper
 * instantiation for every supported array type to coverage's function denominator.
 */
#define GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, VALUE_TYPE ) \
  OP( VALUE_TYPE, 1, RAJA::PERM_I )                         \
  OP( VALUE_TYPE, 2, RAJA::PERM_IJ )                        \
  OP( VALUE_TYPE, 2, RAJA::PERM_JI )                        \
  OP( VALUE_TYPE, 3, RAJA::PERM_IJK )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_IKJ )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_JIK )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_JKI )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_KIJ )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_KJI )                       \
  OP( VALUE_TYPE, 4, RAJA::PERM_IJKL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IJLK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IKJL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IKLJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_ILJK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_ILKJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JIKL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JILK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JKIL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JKLI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JLIK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JLKI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KIJL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KILJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KJIL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KJLI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KLIJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KLJI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LIJK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LIKJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LJIK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LJKI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LKIJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LKJI )

#define GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, VALUE_TYPE ) \
  OP( VALUE_TYPE, 1, RAJA::PERM_I )                             \
  OP( VALUE_TYPE, 2, RAJA::PERM_IJ )                            \
  OP( VALUE_TYPE, 2, RAJA::PERM_JI )

#define GEOS_FOR_EACH_STANDARD_ARRAY( OP )                \
  GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, real32 )         \
  GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, real64 )         \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, integer )    \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, localIndex ) \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, globalIndex )

class TestVTKWriter final : public geos::vtk::VTKPolyDataWriterInterface
{
  using Base = geos::vtk::VTKPolyDataWriterInterface;

public:
  TestVTKWriter()
    : Base( "unused" )
  {
    setOnlyPlotSpecifiedFieldNamesFlag( 1 );
  }

  using Base::writeElementFields;
  using Base::writeNodeFields;
  using Base::writeParticleFields;
};

struct FieldExpectation
{
  string name;
  integer numDims;
  localIndex numComponents;
  bool labeled;
  real64 nodeValue;
  real64 elementValue;
  real64 particleValue;
};

void checkField( vtkDataSetAttributes * const attributes,
                 FieldExpectation const & expected,
                 stdVector< real64 > const & expectedTupleValues )
{
  vtkDataArray * const data = attributes->GetArray( expected.name.c_str() );
  ASSERT_NE( data, nullptr );
  ASSERT_EQ( data->GetNumberOfComponents(), expected.numComponents );
  ASSERT_EQ( data->GetNumberOfTuples(), expectedTupleValues.size() );

  for( vtkIdType tuple = 0; tuple < data->GetNumberOfTuples(); ++tuple )
  {
    for( localIndex component = 0; component < expected.numComponents; ++component )
    {
      EXPECT_DOUBLE_EQ( data->GetComponent( tuple, component ), expectedTupleValues[tuple] );
    }
  }

  if( expected.numDims == 1 )
  {
    EXPECT_EQ( data->GetComponentName( 0 ), nullptr );
    return;
  }

  if( expected.numDims == 2 && !expected.labeled )
  {
    for( localIndex component = 0; component < expected.numComponents; ++component )
    {
      EXPECT_EQ( data->GetComponentName( component ), nullptr );
    }
    return;
  }

  std::set< string > expectedNames;
  for( localIndex component = 0; component < expected.numComponents; ++component )
  {
    string componentName;
    for( integer dim = 1; dim < expected.numDims; ++dim )
    {
      if( dim > 1 )
      {
        componentName += "/";
      }
      bool const high = ( component & ( 1 << ( dim - 1 ) ) ) != 0;
      componentName += expected.labeled ? ( high ? "high" : "low" ) : ( high ? "1" : "0" );
    }
    expectedNames.insert( std::move( componentName ) );
  }

  std::set< string > actualNames;
  for( localIndex component = 0; component < expected.numComponents; ++component )
  {
    char const * const componentName = data->GetComponentName( component );
    ASSERT_NE( componentName, nullptr );
    actualNames.insert( componentName );
  }
  EXPECT_EQ( actualNames, expectedNames );
}

TEST( VTKPolyDataWriter, StandardArrayFieldsAndMetadata )
{
  conduit::Node rootNode;
  Group root( "root", rootNode );

  NodeManager nodeManager( "nodes", &root );
  nodeManager.resize( 4 );
  array1d< localIndex > nodeIndices( 3 );
  nodeIndices[0] = 2;
  nodeIndices[1] = 0;
  nodeIndices[2] = 3;

  CellElementRegion elementRegion( "elements", &root );
  CellElementSubRegion & elementSubRegion0 =
    elementRegion.createElementSubRegion< CellElementSubRegion >( "elementSubRegion0" );
  CellElementSubRegion & elementSubRegion1 =
    elementRegion.createElementSubRegion< CellElementSubRegion >( "elementSubRegion1" );
  elementSubRegion0.resize( 2 );
  elementSubRegion1.resize( 1 );
  std::array< ObjectManagerBase *, 2 > const elementSubRegions =
  { { &elementSubRegion0, &elementSubRegion1 } };

  ParticleRegion particleRegion( "particles", &root );
  Group & particleSubRegionGroup =
    particleRegion.getGroup( ParticleRegionBase::viewKeyStruct::particleSubRegions() );
  ParticleSubRegion & particleSubRegion0 =
    particleSubRegionGroup.registerGroup< ParticleSubRegion >( "particleSubRegion0" );
  ParticleSubRegion & particleSubRegion1 =
    particleSubRegionGroup.registerGroup< ParticleSubRegion >( "particleSubRegion1" );
  particleSubRegion0.resize( 1 );
  particleSubRegion1.resize( 2 );
  std::array< ObjectManagerBase *, 2 > const particleSubRegions =
  { { &particleSubRegion0, &particleSubRegion1 } };

  stdVector< string > const componentLabels{ "low", "high" };
  string_array fieldNames;
  stdVector< FieldExpectation > expectations;
  localIndex fieldIndex = 0;

  #define GEOS_REGISTER_VTK_FIELD( VALUE_TYPE, NDIM, PERMUTATION )                              \
    do                                                                                          \
    {                                                                                           \
      using ArrayType = Array< VALUE_TYPE, NDIM, PERMUTATION >;                                 \
      string const fieldName = "standardArray_" + std::to_string( fieldIndex );                 \
      bool const labeled = fieldIndex % 2 == 0;                                                 \
      real64 const nodeValue = 100 + fieldIndex;                                                \
      real64 const elementValue = 1000 + 10 * fieldIndex;                                       \
      real64 const particleValue = 2000 + 10 * fieldIndex;                                      \
      fieldNames.emplace_back( fieldName );                                                     \
      expectations.push_back( { fieldName, NDIM, 1 << ( NDIM - 1 ), labeled,                    \
                                nodeValue, elementValue, particleValue } );                      \
                                                                                                \
      Wrapper< ArrayType > & nodeWrapper = nodeManager.registerWrapper< ArrayType >( fieldName ); \
      std::array< localIndex, NDIM > nodeDims;                                                   \
      nodeDims.fill( 2 );                                                                        \
      nodeDims[0] = nodeManager.size();                                                          \
      nodeWrapper.reference().resize( NDIM, nodeDims.data() );                                  \
      std::fill( nodeWrapper.reference().data(),                                                 \
                 nodeWrapper.reference().data() + nodeWrapper.reference().size(),               \
                 static_cast< VALUE_TYPE >( nodeValue ) );                                      \
      if( labeled )                                                                              \
      {                                                                                          \
        for( integer dim = 1; dim < NDIM; ++dim )                                               \
        {                                                                                        \
          nodeWrapper.setDimLabels( dim, { componentLabels.begin(), componentLabels.end() } );  \
        }                                                                                        \
      }                                                                                          \
                                                                                                \
      for( std::size_t subRegionIndex = 0; subRegionIndex < elementSubRegions.size();            \
           ++subRegionIndex )                                                                    \
      {                                                                                          \
        ObjectManagerBase & subRegion = *elementSubRegions[subRegionIndex];                      \
        Wrapper< ArrayType > & wrapper = subRegion.registerWrapper< ArrayType >( fieldName );    \
        std::array< localIndex, NDIM > dims;                                                     \
        dims.fill( 2 );                                                                           \
        dims[0] = subRegion.size();                                                               \
        wrapper.reference().resize( NDIM, dims.data() );                                         \
        std::fill( wrapper.reference().data(),                                                    \
                   wrapper.reference().data() + wrapper.reference().size(),                      \
                   static_cast< VALUE_TYPE >( elementValue + subRegionIndex ) );                 \
        if( labeled )                                                                             \
        {                                                                                         \
          for( integer dim = 1; dim < NDIM; ++dim )                                              \
          {                                                                                       \
            wrapper.setDimLabels( dim, { componentLabels.begin(), componentLabels.end() } );     \
          }                                                                                       \
        }                                                                                         \
      }                                                                                           \
                                                                                                \
      for( std::size_t subRegionIndex = 0; subRegionIndex < particleSubRegions.size();            \
           ++subRegionIndex )                                                                     \
      {                                                                                           \
        ObjectManagerBase & subRegion = *particleSubRegions[subRegionIndex];                      \
        Wrapper< ArrayType > & wrapper = subRegion.registerWrapper< ArrayType >( fieldName );     \
        std::array< localIndex, NDIM > dims;                                                      \
        dims.fill( 2 );                                                                            \
        dims[0] = subRegion.size();                                                                \
        wrapper.reference().resize( NDIM, dims.data() );                                          \
        std::fill( wrapper.reference().data(),                                                     \
                   wrapper.reference().data() + wrapper.reference().size(),                       \
                   static_cast< VALUE_TYPE >( particleValue + subRegionIndex ) );                 \
        if( labeled )                                                                              \
        {                                                                                          \
          for( integer dim = 1; dim < NDIM; ++dim )                                               \
          {                                                                                        \
            wrapper.setDimLabels( dim, { componentLabels.begin(), componentLabels.end() } );      \
          }                                                                                        \
        }                                                                                          \
      }                                                                                           \
      ++fieldIndex;                                                                               \
    } while( false );

  GEOS_FOR_EACH_STANDARD_ARRAY( GEOS_REGISTER_VTK_FIELD )
#undef GEOS_REGISTER_VTK_FIELD

  ASSERT_EQ( fieldIndex, 75 );

  string const indexedFieldName = "indexedNodeField";
  array1d< real64 > & indexedField =
    nodeManager.registerWrapper< array1d< real64 > >( indexedFieldName ).reference();
  indexedField.resize( nodeManager.size() );
  for( localIndex node = 0; node < nodeManager.size(); ++node )
  {
    indexedField[node] = 10 + node;
  }
  fieldNames.emplace_back( indexedFieldName );

  TestVTKWriter writer;
  writer.setFieldNames( fieldNames );

  vtkNew< vtkUnstructuredGrid > nodeOutput;
  writer.writeNodeFields( nodeManager, nodeIndices.toViewConst(), nodeOutput->GetPointData() );

  vtkNew< vtkUnstructuredGrid > elementOutput;
  writer.writeElementFields( elementRegion, elementOutput->GetCellData() );

  vtkNew< vtkUnstructuredGrid > particleOutput;
  writer.writeParticleFields( particleRegion, particleOutput->GetCellData() );

  for( FieldExpectation const & expected : expectations )
  {
    SCOPED_TRACE( expected.name );
    checkField( nodeOutput->GetPointData(), expected,
                { expected.nodeValue, expected.nodeValue, expected.nodeValue } );
    checkField( elementOutput->GetCellData(), expected,
                { expected.elementValue, expected.elementValue, expected.elementValue + 1 } );
    checkField( particleOutput->GetCellData(), expected,
                { expected.particleValue, expected.particleValue + 1, expected.particleValue + 1 } );
  }

  FieldExpectation const indexedExpectation{ indexedFieldName, 1, 1, false, 0, 0, 0 };
  checkField( nodeOutput->GetPointData(), indexedExpectation, { 12, 10, 13 } );
}

#undef GEOS_FOR_EACH_STANDARD_ARRAY
#undef GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY
#undef GEOS_FOR_EACH_REAL_STANDARD_ARRAY

} // namespace
