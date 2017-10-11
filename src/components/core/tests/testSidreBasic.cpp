#include <gtest/gtest.h>
#include <mpi.h>
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "common/DataTypes.hpp"


namespace geosx {
namespace dataRepository {

#ifdef USE_ATK
TEST(testSidreBasic, testSidreBasic) {
  MPI_Init(0, nullptr);
  const string path = "test_sidre_basic";
  const string protocol = "sidre_hdf5";
  const int group_size = 44;
  const int sized_from_parent = 55;
  const int num_items = 10;
  const uint expected_size = num_items * sizeof(int64);
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);
  root->resize(group_size);

  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<int64_array> * data_view = root->RegisterViewWrapper<int64_array>("int64_data");
  data_view->setSizedFromParent(sized_from_parent);

  /* Resize the array */
  data_view->resize(num_items);

  /* Check that the ViewWrapper size and byteSize functions return the proper values */
  EXPECT_EQ(data_view->size(), num_items);
  EXPECT_EQ(data_view->byteSize(), expected_size);

  /* Set the data */
  ViewWrapper<int64_array>::rtype data = data_view->data();
  for (int i = 0; i < num_items; i++) {
      data[i] = i;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  int64_ptr dataPtr = data_view->dataPtr();
  EXPECT_EQ(dataPtr, &(data[0]));
  for (int i = 0; i < data_view->size(); i++) {
    EXPECT_EQ(dataPtr[i], data[i]);
  }

  /* Save the sidre tree */
  root->prepareToWriteRestart();
  SidreWrapper::writeTree(1, path, protocol, MPI_COMM_WORLD);

  /* Delete geos tree and reset sidre tree. */
  delete root;
  ds.destroyAllAttributes();
  ds.destroyAllBuffers();
  ds.getRoot()->destroyGroups();

  /* Restore the sidre tree */
  SidreWrapper::reconstructTree(path + ".root", protocol, MPI_COMM_WORLD);
  root = new ManagedGroup(std::string("data"), nullptr);

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  ViewWrapper<int64_array> * data_view_new = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Load the data */
  root->prepareToLoadExternalData();
  SidreWrapper::loadExternalData(path + ".root", MPI_COMM_WORLD);
  root->finishLoadingExternalData();

  /* Should be the same as stored. */
  EXPECT_EQ(data_view_new->size(), num_items);
  ViewWrapper<int64_array>::rtype_const data_new = data_view_new->data();
  for (int i = 0; i < data_view_new->size(); i++) {
    EXPECT_EQ(data_new[i], i);
  }

  EXPECT_EQ(data_view_new->sizedFromParent(), sized_from_parent);
  EXPECT_EQ(root->size(), group_size);

  delete root;
  MPI_Finalize();
}
#endif /* ATK_FOUND */


int main(int argc, char* argv[]) {
  int result = 0;
  testing::InitGoogleTest(&argc, argv);
  result = RUN_ALL_TESTS();
  return result;
}


} /* end namespace dataRepository */
} /* end namespace goesx */
