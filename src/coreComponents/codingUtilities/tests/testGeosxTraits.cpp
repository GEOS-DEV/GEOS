/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/traits.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::traits;

TEST( testGeosxTraits, Pointer )
{
  static_assert( std::is_same< Pointer< int >, int * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< R1Tensor >, R1Tensor * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< std::vector< double > >, double * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< std::string >, char * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< array3d< std::string > >, std::string * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< SortedArray< float > >, float const * >::value, "Should be true." );

  static_assert( std::is_same< ConstPointer< int >, int const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< R1Tensor >, R1Tensor const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< std::vector< double > >, double  const* >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< std::string >, char const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< array3d< std::string > >, std::string const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< SortedArray< float > >, float const * >::value, "Should be true." );
}

TEST( testGeosxTraits, has_alias_value_type )
{
  static_assert( has_alias_value_type< array1d< double > >, "Should be true." );
  static_assert( has_alias_value_type< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_alias_value_type< SortedArray< string > >, "Should be true." );
  static_assert( has_alias_value_type< std::vector< int > >, "Should be true." );
  static_assert( has_alias_value_type< std::map< string, string > >, "Should be true." );

  static_assert( !has_alias_value_type< int >, "Should be false." );
  static_assert( !has_alias_value_type< double >, "Should be false." );
  static_assert( !has_alias_value_type< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_data_method )
{
  static_assert( has_data_method< array1d< double > >, "Should be true." );
  static_assert( has_data_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_data_method< std::vector< int > >, "Should be true." );

  static_assert( !has_data_method< std::map< string, string > >, "Should be false." );
  static_assert( !has_data_method< SortedArray< string > >, "Should be false." );
  static_assert( !has_data_method< int >, "Should be false." );
  static_assert( !has_data_method< double >, "Should be false." );
  static_assert( !has_data_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_chai_move_method )
{
  static_assert( has_chai_move_method< array1d< double > >, "Should be true." );
  static_assert( has_chai_move_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_chai_move_method< SortedArray< string > >, "Should be true." );
  static_assert( has_chai_move_method< ArrayOfArrays< int > >, "Should be true." );
  
  static_assert( !has_chai_move_method< std::vector< int > >, "Should be true." );
  static_assert( !has_chai_move_method< std::map< string, string > >, "Should be true." );
  static_assert( !has_chai_move_method< int >, "Should be false." );
  static_assert( !has_chai_move_method< double >, "Should be false." );
  static_assert( !has_chai_move_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_empty_method )
{
  static_assert( has_empty_method< array1d< double > >, "Should be true." );
  static_assert( has_empty_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_empty_method< SortedArray< string > >, "Should be true." );
  static_assert( has_empty_method< std::vector< int > >, "Should be true." );
  static_assert( has_empty_method< std::map< string, string > >, "Should be true." );

  static_assert( !has_empty_method< int >, "Should be false." );
  static_assert( !has_empty_method< double >, "Should be false." );
  static_assert( !has_empty_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_size_method )
{
  static_assert( has_size_method< array1d< double > >, "Should be true." );
  static_assert( has_size_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_size_method< SortedArray< string > >, "Should be true." );
  static_assert( has_size_method< std::vector< int > >, "Should be true." );
  static_assert( has_size_method< std::map< string, string > >, "Should be true." );

  static_assert( !has_size_method< int >, "Should be false." );
  static_assert( !has_size_method< double >, "Should be false." );
  static_assert( !has_size_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_dimension_size_method )
{
  static_assert( has_dimension_size_method< array1d< double > >, "Should be true." );
  static_assert( has_dimension_size_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  
  static_assert( !has_dimension_size_method< SortedArray< string > >, "Should be false." );
  static_assert( !has_dimension_size_method< std::vector< int > >, "Should be false." );
  static_assert( !has_dimension_size_method< std::map< string, string > >, "Should be false." );
  static_assert( !has_dimension_size_method< int >, "Should be false." );
  static_assert( !has_dimension_size_method< double >, "Should be false." );
  static_assert( !has_dimension_size_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, has_resize_method )
{
  static_assert( has_resize_method< array1d< double > >, "Should be true." );
  static_assert( has_resize_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( has_resize_method< std::vector< int > >, "Should be true." );
  static_assert( has_resize_method< ArrayOfArrays< int > >, "Should be true." );
  static_assert( has_resize_method< ArrayOfSets< int > >, "Should be true." );
  
  static_assert( !has_resize_method< SortedArray< string > >, "Should be false." );
  static_assert( !has_resize_method< std::map< string, string > >, "Should be false." );
  static_assert( !has_resize_method< int >, "Should be false." );
  static_assert( !has_resize_method< double >, "Should be false." );
  static_assert( !has_resize_method< R2SymTensor >, "Should be false." );

}

TEST( testGeosxTraits, has_resize_default_method )
{
  static_assert( has_resize_default_method< array1d< double >, double >, "Should be true." );
  static_assert( has_resize_default_method< array5d< array1d< R1Tensor > >, array1d< R1Tensor > >, "Should be true." );
  static_assert( has_resize_default_method< array2d< std::string >, std::string >, "Should be true." );

  static_assert( !has_resize_default_method< array1d< double >, int >, "Should be false." );
  static_assert( !has_resize_default_method< array5d< array1d< R1Tensor > >, array1d< float > >, "Should be false." );
  static_assert( !has_resize_default_method< array2d< std::string >, char >, "Should be false." );

  static_assert( !has_resize_default_method< std::vector< int >, int >, "Should be false." );
  static_assert( !has_resize_default_method< SortedArray< std::string >, std::string >, "Should be false." );
  static_assert( !has_resize_default_method< std::map< string, string >, float >, "Should be false." );
  static_assert( !has_resize_default_method< int, int >, "Should be false." );
  static_assert( !has_resize_default_method< double, std::string >, "Should be false." );
  static_assert( !has_resize_default_method< R2SymTensor, char >, "Should be false." );
}

TEST( testGeosxTraits, has_resize_dimensions_method )
{
  static_assert( has_resize_dimensions_method< array1d< double > >, "Should be true." );
  static_assert( has_resize_dimensions_method< array5d< array1d< R1Tensor > > >, "Should be true." );
  
  static_assert( !has_resize_dimensions_method< SortedArray< string > >, "Should be false." );
  static_assert( !has_resize_dimensions_method< std::vector< int > >, "Should be false." );
  static_assert( !has_resize_dimensions_method< std::map< string, string > >, "Should be false." );
  static_assert( !has_resize_dimensions_method< int >, "Should be false." );
  static_assert( !has_resize_dimensions_method< double >, "Should be false." );
  static_assert( !has_resize_dimensions_method< R2SymTensor >, "Should be false." );
}

TEST( testGeosxTraits, is_tensorT )
{
  static_assert( is_tensorT< R1Tensor >, "Should be true." );
  static_assert( is_tensorT< R2Tensor >, "Should be true." );
  static_assert( is_tensorT< R2SymTensor >, "Should be true." );

  static_assert( !is_tensorT< int >, "Should be false." );
  static_assert( !is_tensorT< double >, "Should be false." );
  static_assert( !is_tensorT< void >, "Should be false." );
}


TEST( testGeosxTraits, is_string )
{
  static_assert( is_string< string >, "Should be true." );

  static_assert( !is_string< int >, "Should be false." );
  static_assert( !is_string< double >, "Should be false." );
  static_assert( !is_string< void >, "Should be false." );
}

TEST( testGeosxTraits, is_array )
{
  static_assert( is_array< LvArray::Array< int, 1, camp::make_idx_seq< 0 >::type, int > >, "Should be true." );
  
  static_assert( !is_array< int >, "Should be false." );
  static_assert( !is_array< double >, "Should be false." );
  static_assert( !is_array< void >, "Should be false." );
}

TEST( testGeosxTraits, is_map )
{
  static_assert( is_map< map< string, int > >, "Should be true." );
  static_assert( is_map< unordered_map< string, int > >, "Should be true." );

  static_assert(!is_map< int >, "Should be false." );
  static_assert(!is_map< double >, "Should be false." );
  static_assert(!is_map< void >, "Should be false." );
  SUCCEED();
}

TEST( testGeosxTraits, is_set )
{
  static_assert( is_set< SortedArray< string > >, "Should be true." );

  static_assert( !is_set< int >, "Should be false." );
  static_assert( !is_set< double >, "Should be false." );
  static_assert( !is_set< void >, "Should be false." );
  SUCCEED();
}

TEST( testGeosxTraits, is_pair )
{
  static_assert( is_pair< std::pair< string, int > >, "Should be true." );

  static_assert( !is_pair< int >, "Should be false" );
  static_assert( !is_pair< double >, "Should be false" );
  static_assert( !is_pair< void >, "Should be false" );
}

