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


TEST(testSidre, simpleRestore) {
  MPI_Init(0, nullptr);
  const string path = "geos_simple_restore";
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

  /* Check that the name ViewWrapper is self consistent. */
  ViewWrapper<string> & name_wrapper = root->getWrapper<string>(string("name"));
  string & name_str = name_wrapper.data();
  const char * name_ptr = name_wrapper.dataPtr();
  EXPECT_TRUE(name_str.c_str() == name_ptr) << name_str.c_str() << ", " << name_ptr << std::endl;
  EXPECT_TRUE(name_str.size() == (uint32) name_wrapper.dataSize()) << name_str.size() << ", " << name_wrapper.dataSize() << std::endl;

  for (uint32 i = 0; i < name_str.size(); i++) {
    EXPECT_TRUE(name_str[i] == name_ptr[i]) << name_str[i] << ", " << name_ptr[i] << std::endl;
  }

  /* Check that the path ViewWrapper is self consistent. */
  ViewWrapper<string> & path_wrapper = root->getWrapper<string>(string("path"));
  string & path_str = path_wrapper.data();
  const char * path_ptr = path_wrapper.dataPtr();
  EXPECT_TRUE(path_str.c_str() == path_ptr) << path_str.c_str() << ", " << path_ptr << std::endl;
  EXPECT_TRUE(path_str.size() == (uint32) path_wrapper.dataSize()) << path_str.size() << ", " << path_wrapper.dataSize() << std::endl;

  for (uint32 i = 0; i < path_str.size(); i++) {
    EXPECT_TRUE(path_str[i] == path_ptr[i]) << path_str[i] << ", " << path_ptr[i] << std::endl;
  }

  /* Check that the size ViewWrapper is self consistent. */
  ViewWrapper<int32> & size_wrapper = root->getWrapper<int32>(string("size"));
  int & size_int = *size_wrapper.data();
  int * size_ptr = size_wrapper.dataPtr();
  EXPECT_TRUE(&size_int == size_ptr) << &size_int << ", " << size_ptr << std::endl;
  EXPECT_TRUE(size_int == *size_ptr) << size_int << ", " << *size_ptr << std::endl;

  /* Save the sidre tree */
  root->writeRestart(1, path, protocol, MPI_COMM_WORLD);

#if 1
  /* Reset the DataStore and the mirrored ManagedGroup heirarchy. */
  SidreWrapper::dataStore().getRoot()->destroyGroups();
  SidreWrapper::dataStore().destroyAllBuffers();
  delete root;

  /* Restore the sidre tree */
  root = new ManagedGroup(std::string("data"), nullptr);
  root->reconstructSidreTree(path + ".root", protocol, MPI_COMM_WORLD);

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  ViewWrapper<int64_array> & data_view_new = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Load the data */
  root->resizeSubViews();
  root->registerSubViews();
  root->loadSidreExternalData(path + ".root", MPI_COMM_WORLD);

  /* Should be the same as stored. */
  data = data_view_new.data();
  for (int i = 0; i < num_items; i++) {
    EXPECT_TRUE(data[i] == data[i]);
  }
#endif

  MPI_Finalize();
  EXPECT_TRUE(true);
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
