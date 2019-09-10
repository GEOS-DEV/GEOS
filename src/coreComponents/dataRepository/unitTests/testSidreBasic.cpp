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

#include <gtest/gtest.h>

#include <mpi.h>

#include "dataRepository/Group.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "dataRepository/Wrapper.hpp"
#include "common/DataTypes.hpp"



namespace geosx
{
namespace dataRepository
{

#ifdef GEOSX_USE_ATK
TEST( testSidreBasic, testSidreBasic )
{
  MPI_Init( nullptr, nullptr );
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  const string path = "test_sidre_basic";
  const string protocol = "sidre_hdf5";
  const int group_size = 44;
  const int sized_from_parent = 55;
  const int num_items = 10;
  const uint expected_size = num_items * sizeof(globalIndex);
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  Group * root = new Group( std::string( "data" ), nullptr );
  root->resize( group_size );

  /* Create a Wrapper which creates the associated sidre::View */
  Wrapper< globalIndex_array > * data_view = root->registerWrapper< globalIndex_array >( "globalIndex_data" );
  data_view->setSizedFromParent( sized_from_parent );

  /* Resize the array */
  data_view->resize( num_items );

  /* Check that the Wrapper size and byteSize functions return the proper values */
  EXPECT_EQ( data_view->size(), num_items );
  EXPECT_EQ( data_view->byteSize(), expected_size );

  /* Set the data */
  globalIndex_array & data = data_view->reference();
  for( int i = 0 ; i < num_items ; i++ )
  {
    data[i] = i;
  }

  /* Check that the Wrapper dataPtr points to the right thing */
  globalIndex * dataPtr = data_view->dataPtr();
  EXPECT_EQ( dataPtr, &(data[0]));
  for( int i = 0 ; i < data_view->size() ; i++ )
  {
    EXPECT_EQ( dataPtr[i], data[i] );
  }

  /* Save the sidre tree */
  root->prepareToWrite();
  SidreWrapper::writeTree( 1, path, protocol, MPI_COMM_GEOSX );
  root->finishWriting();

  /* Delete geos tree and reset sidre tree. */
  delete root;
  ds.destroyAllAttributes();
  ds.destroyAllBuffers();
  ds.getRoot()->destroyGroups();

  /* Restore the sidre tree */
  SidreWrapper::reconstructTree( path + ".root", protocol, MPI_COMM_GEOSX );
  root = new Group( std::string( "data" ), nullptr );

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  Wrapper< globalIndex_array > * data_view_new = root->registerWrapper< globalIndex_array >( "globalIndex_data" );

  /* Load the data */
  root->prepareToRead();
  SidreWrapper::loadExternalData( path + ".root", MPI_COMM_GEOSX );
  root->finishReading();

  /* Should be the same as stored. */
  EXPECT_EQ( data_view_new->size(), num_items );

  globalIndex_array & data_new = data_view_new->reference();
  for( int i = 0 ; i < data_view_new->size() ; i++ )
  {
    EXPECT_EQ( data_new[i], i );
  }

  EXPECT_EQ( data_view_new->sizedFromParent(), sized_from_parent );
  EXPECT_EQ( root->size(), group_size );

  delete root;
  MPI_Finalize();
}

#endif /* GEOSX_USE_ATK */


int main( int argc, char * argv[] ) {
  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();
  return result;
}


} /* end namespace dataRepository */
} /* end namespace goesx */
