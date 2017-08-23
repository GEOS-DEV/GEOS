#include <gtest/gtest.h>
#include <mpi.h>
#include "sidre/sidre.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "common/DataTypes.hpp"

#include <iostream>


namespace geosx {
namespace dataRepository {


TEST(testSidreWrite, simpleWrite) {
  MPI_Init(0, nullptr);
  const string path = "geos_simple_write";
  const string protocol = "sidre_hdf5";
  const int num_items = 10;
  const int expected_size = num_items * sizeof(int64);
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);

  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<int64_array> & data_view = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Resize the array */
  data_view.resize(num_items);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(data_view.size() == num_items) << data_view.size() << ", " << num_items << std::endl;
  EXPECT_TRUE(data_view.dataSize() == expected_size) << data_view.dataSize() << ", " << expected_size << std::endl;

  /* Set the data */
  int64_array & data = data_view.data();
  for (int i = 0; i < num_items; i++) {
      data[i] = i;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  int64_ptr dataPtr = data_view.dataPtr();
  EXPECT_TRUE(dataPtr == &(data[0])) << dataPtr << ", " << &(data[0]);
  for (int i = 0; i < data_view.size(); i++) {
    EXPECT_TRUE(dataPtr[i] == data[i]) << dataPtr[i] << ", " << data[i] << std::endl;
  }

  /* Check that the size ViewWrapper is self consistent. */
  ViewWrapper<int32> & size_wrapper = root->getWrapper<int32>(string("size"));
  int & size_int = *size_wrapper.data();
  int * size_ptr = size_wrapper.dataPtr();
  EXPECT_TRUE(&size_int == size_ptr) << &size_int << ", " << size_ptr << std::endl;
  EXPECT_TRUE(size_int == *size_ptr) << size_int << ", " << *size_ptr << std::endl;

  /* Save the sidre tree */
  root->writeRestart(1, path, protocol, MPI_COMM_WORLD);

  /* Delete geos tree and reset sidre tree. */
  delete root;
  ds.destroyAllAttributes();
  ds.destroyAllBuffers();
  ds.getRoot()->destroyGroups();

  /* Restore the sidre tree */
  root = new ManagedGroup(std::string("data"), nullptr);
  root->reconstructSidreTree(path + ".root", protocol, MPI_COMM_WORLD);

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  ViewWrapper<int64_array> & data_view_new = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Load the data */
  root->loadSidreExternalData(path + ".root", MPI_COMM_WORLD);

  /* Should be the same as stored. */
  data = data_view_new.data();
  for (int i = 0; i < data_view.size(); i++) {
    EXPECT_TRUE(data[i] == data[i]);
  }

  delete root;
  MPI_Finalize();
}


int main(int argc, char* argv[]) {
  int result = 0;

  testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}


} // end namespace dataRepository
} // end namespace goesx
