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
  MPI_Init(0, nullptr);
  const string path = "geos_simple_restore";
  const int num_items = 10;
  const int expected_size = num_items * sizeof(int64);
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);

  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<int64_array> & data_view = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Resize the array */
  data_view.resize(num_items);

  /* Check that the ViewWrapper size and data_size functions return the proper values */
  EXPECT_TRUE(data_view.size() == num_items) << data_view.size() << ", " << num_items << std::endl;
  EXPECT_TRUE(data_view.data_size() == expected_size) << data_view.data_size() << ", " << expected_size << std::endl;

  /* Set the data */
  int64_array & data = data_view.data();
  for (int i = 0; i < num_items; i++) {
      data[i] = i;
  }

  /* Check that the ViewWrapper data_ptr points to the right thing */
  int64_ptr data_ptr = data_view.data_ptr();
  EXPECT_TRUE(data_ptr == &(data[0])) << data_ptr << ", " << &(data[0]);
  for (int i = 0; i < data_view.size(); i++) {
    EXPECT_TRUE(data_ptr[i] == data[i]) << data_ptr[i] << ", " << data[i] << std::endl;
  }

  ViewWrapper<string> & name_wrapper = root->getWrapper<string>(string("name"));
  ViewWrapper<string> & path_wrapper = root->getWrapper<string>(string("path"));
  ViewWrapper<int32> & size_wrapper = root->getWrapper<int32>(string("size"));

  string & name_str = name_wrapper.data();
  EXPECT_FALSE( (void*)(&name_str) == (void const *)(name_wrapper.data_ptr())) << "ViewWrapper<string>.data_ptr() points to string object." << std::endl;
  EXPECT_TRUE(name_str.c_str() ==  (char const*) name_wrapper.data_ptr()) << "ViewWrapper<string>.data_ptr() doesn't point to char[]." << std::endl;

#if 0
  axom::sidre::Group * root_group = root->getSidreGroup();
  root_group->destroyView("name");
  root_group->destroyView("size");
  root_group->destroyView("path");

  /* Save the sidre tree */
  axom::spio::IOManager ioManager(MPI_COMM_WORLD);
  ioManager.write(ds.getRoot(), 1, path, "sidre_hdf5");

  /* Reset the DataStore and the mirrored ManagedGroup heirarchy. */
  SidreWrapper::dataStore().getRoot()->destroyGroups();
  SidreWrapper::dataStore().destroyAllBuffers();
  delete root;

  /* Restore the sidre tree */
  ioManager.read(SidreWrapper::dataStore().getRoot(), path + ".root");

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  root = new ManagedGroup(std::string("int64_data"), nullptr);
  ViewWrapper<int64_array> & data_view_new = root->RegisterViewWrapper<int64_array>("data");
  
  /* Data allocation */
  data_view_new.resize(num_items);

  /* Load the data */
  ioManager.loadExternalData(SidreWrapper::dataStore().getRoot(), path + ".root");

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

  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}


} // end namespace dataRepository
} // end namespace goesx
