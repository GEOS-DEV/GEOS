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

  /* Restore the sidre tree */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);
  root->reconstructSidreTree(path + ".root", protocol, MPI_COMM_WORLD);

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  ViewWrapper<int64_array> & data_view = root->RegisterViewWrapper<int64_array>("int64_data");

  /* Load the data */
  root->resizeSubViews();
  root->registerSubViews();
  root->loadSidreExternalData(path + ".root", MPI_COMM_WORLD);

  /* Should be the same as stored. */
  int64_array & data = data_view.data();
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
