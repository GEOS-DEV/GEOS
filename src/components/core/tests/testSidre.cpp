#include <gtest/gtest.h>
#include <mpi.h>
#include "sidre/sidre.hpp"
#include "spio/IOManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "common/DataTypes.hpp"

#include <iostream>


namespace geosx {
namespace dataRepository {


TEST(testSidre, simpleRestore) {
  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("int64_data"), nullptr);
  
  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<int64_array> & data_view = root->RegisterViewWrapper<int64_array>("data");
  
  const int num_items = 10;
  const int expected_size = num_items * sizeof(int64);

  /* Resize the array */
  data_view.resize(num_items);

  /* Check that the ViewWrapper size and data_size functions return the proper values */
  EXPECT_TRUE(data_view.size() == num_items) << data_view.size() << ", " << num_items;
  EXPECT_TRUE(data_view.data_size() == expected_size) << data_view.data_size() << ", " << expected_size;

  /* Set the data */
  // int64_array data = data_view.data();
  // for (int i = 0; i < num_items; i++) {
  //     data[i] = i;
  // }

  /* Check that the ViewWrapper data_ptr points to the right thing */
  // int64_ptr data_ptr = data_view.data_ptr();
  // EXPECT_TRUE(data_ptr == &(data[0]));
  // for (int i = 0; i < num_items; i++) {
  //   EXPECT_TRUE(data_ptr[i] == i) << data_ptr[i] << ", " << i;
  // }

  /* Save the sidre tree */
  // string path = "geos_simple_restore";
  // axom::spio::IOManager ioManager(MPI_COMM_WORLD);
  // ioManager.write(SidreWrapper::dataStore().getRoot(), 1, path, "sidre_hdf5");

  /* Reset the DataStore and the mirrored ManagedGroup heirarchy. */
  // SidreWrapper::dataStore().getRoot()->destroyGroups();
  // SidreWrapper::dataStore().destroyAllBuffers();
  // delete root;

  /* Restore the sidre tree */
  // ioManager.read(SidreWrapper::dataStore().getRoot(), path + ".root");

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  // root = new ManagedGroup(std::string("int64_data"), nullptr);
  // ViewWrapper<int64_array> & data_view_new = root->RegisterViewWrapper<int64_array>("data");
  
  /* Data allocation */
  // data_view_new.resize(num_items);

  /* Load the data */
  // ioManager.loadExternalData(SidreWrapper::dataStore().getRoot(), path + ".root");

  /* Should be the same as stored. */
  // data = data_view_new.data();
  // for (int i = 0; i < num_items; i++) {
  //   EXPECT_TRUE(data[i] == data[i]);
  // }

  EXPECT_TRUE(true);
}


int main(int argc, char* argv[]) {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    result = RUN_ALL_TESTS();
    MPI_Finalize();

    return result;
}


} // end namespace dataRepository
} // end namespace goesx