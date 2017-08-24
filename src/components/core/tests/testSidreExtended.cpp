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


TEST(testSidreExtended, testSidreExtended) {
  MPI_Init(0, nullptr);
  const string path = "test_sidre_extended";
  const string protocol = "sidre_hdf5";
  const int group_size = 44;
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);
  root->resize(group_size);

  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<int64_array> & view_int64 = root->RegisterViewWrapper<int64_array>("int64");

  /* Resize the array */
  int view_int64_size = 10;
  uint32 expected_size = view_int64_size * sizeof(int64);
  view_int64.resize(view_int64_size);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_int64.size() == view_int64_size);
  EXPECT_TRUE(view_int64.dataSize() == expected_size);

  /* Set the data */
  for (int i = 0; i < view_int64.size(); i++) {
      view_int64.data()[i] = i;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_int64.dataPtr() == &(view_int64.data()[0]));



  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<string> & view_hope = root->RegisterViewWrapper<string>("hope");
  
  string hope = "I hope this works";
  int view_hope_size = hope.size();
  expected_size = view_hope_size * sizeof(char);

  /* Set the data */
  view_hope.data() = hope;

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_hope.size() == view_hope_size);
  EXPECT_TRUE(view_hope.dataSize() == expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_hope.dataPtr() == view_hope.data().c_str());





  ManagedGroup & strings_group = root->RegisterGroup("strings");
  strings_group.resize(group_size);


  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<string> & view_hello = strings_group.RegisterViewWrapper<string>("hello");
  
  string hello = "Hello, how are you doing today?";
  int view_hello_size = hello.size();
  expected_size = view_hello_size * sizeof(char);

  /* Set the data */
  view_hello.data() = hello;

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_hello.size() == view_hello_size);
  EXPECT_TRUE(view_hello.dataSize() == expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_hello.dataPtr() == view_hello.data().c_str());



  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<string> & view_goodbye = strings_group.RegisterViewWrapper<string>("goodbye");
  
  string goodbye = "I am not doing very well so I'll bid you goodbye.";
  int view_goodbye_size = goodbye.size();
  expected_size = view_goodbye_size * sizeof(char);

  /* Set the data */
  view_goodbye.data() = goodbye;

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_goodbye.size() == view_goodbye_size);
  EXPECT_TRUE(view_goodbye.dataSize() == expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_goodbye.dataPtr() == view_goodbye.data().c_str());




  ManagedGroup & doubles_group = root->RegisterGroup("doubles");
  doubles_group.resize(group_size);


  ViewWrapper<real64_array> & view_real64_1 = doubles_group.RegisterViewWrapper<real64_array>("doubles1");

  /* Resize the array */
  int view_real64_1_size = 100;
  expected_size = view_real64_1_size * sizeof(real64);
  view_real64_1.resize(view_real64_1_size);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_real64_1.size() == view_real64_1_size);
  EXPECT_TRUE(view_real64_1.dataSize() == expected_size);

  /* Set the data */
  for (int i = 0; i < view_real64_1.size(); i++) {
      view_real64_1.data()[i] = double(i) / 100.0;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_real64_1.dataPtr() == &(view_real64_1.data()[0]));


  ViewWrapper<real64_array> & view_real64_2 = doubles_group.RegisterViewWrapper<real64_array>("doubles2");

  /* Resize the array */
  int view_real64_2_size = 100;
  expected_size = view_real64_2_size * sizeof(real64);
  view_real64_2.resize(view_real64_2_size);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_real64_2.size() == view_real64_2_size);
  EXPECT_TRUE(view_real64_2.dataSize() == expected_size);

  /* Set the data */
  for (int i = 0; i < view_real64_2.size(); i++) {
      view_real64_2.data()[i] = double(i) * double(i) / 100.0;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_real64_2.dataPtr() == &(view_real64_2.data()[0]));




  ManagedGroup & mixed_group = doubles_group.RegisterGroup("mixed");
  mixed_group.resize(group_size);


  ViewWrapper<int32_array> & view_int32 = mixed_group.RegisterViewWrapper<int32_array>("int32");

  /* Resize the array */
  int view_int32_size = 1000;
  expected_size = view_int32_size * sizeof(int32);
  view_int32.resize(view_int32_size);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_int32.size() == view_int32_size);
  EXPECT_TRUE(view_int32.dataSize() == expected_size);

  /* Set the data */
  for (int i = 0; i < view_int32.size(); i++) {
      view_int32.data()[i] = i * i;
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_int32.dataPtr() == &(view_int32.data()[0]));


  ViewWrapper<real32_array> & view_real32 = mixed_group.RegisterViewWrapper<real32_array>("real32");

  /* Resize the array */
  int view_real32_size = 1000;
  expected_size = view_real32_size * sizeof(int32);
  view_real32.resize(view_real32_size);

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_real32.size() == view_real32_size);
  EXPECT_TRUE(view_real32.dataSize() == expected_size);

  /* Set the data */
  for (int i = 0; i < view_real32.size(); i++) {
      view_real32.data()[i] = float(i) / float(i + 1.0 - i * i);
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_real32.dataPtr() == &(view_real32.data()[0]));


  /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<string> & view_whatup = mixed_group.RegisterViewWrapper<string>("whatup");
  
  string whatup = "Whatup man? I like you're hat.";
  int view_whatup_size = whatup.size();
  expected_size = view_whatup_size * sizeof(char);

  /* Set the data */
  view_whatup.data() = whatup;

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_whatup.size() == view_whatup_size);
  EXPECT_TRUE(view_whatup.dataSize() == expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_whatup.dataPtr() == view_whatup.data().c_str());


   /* Create a ViewWrapper which creates the associated sidre::View */
  ViewWrapper<real64> & view_pi = mixed_group.RegisterViewWrapper<real64>("pi");

  int view_pi_size = 1;
  expected_size = sizeof(real64);

  /* Set the data */
  *(view_pi.data()) = 3.14159;

  /* Check that the ViewWrapper size and dataSize functions return the proper values */
  EXPECT_TRUE(view_pi.size() == view_pi_size);
  EXPECT_TRUE(view_pi.dataSize() == expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view_pi.dataPtr() == view_pi.data());





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
  ViewWrapper<int64_array> & view_int64_new= root->RegisterViewWrapper<int64_array>("int64");
  ViewWrapper<string> & view_hope_new = root->RegisterViewWrapper<string>("hope");

  ManagedGroup & strings_group_new = root->RegisterGroup("strings");
  ViewWrapper<string> & view_hello_new = strings_group_new.RegisterViewWrapper<string>("hello");
  ViewWrapper<string> & view_goodbye_new = strings_group_new.RegisterViewWrapper<string>("goodbye");
  
  ManagedGroup & doubles_group_new = root->RegisterGroup("doubles");
  ViewWrapper<real64_array> & view_real64_1_new = doubles_group_new.RegisterViewWrapper<real64_array>("doubles1");
  ViewWrapper<real64_array> & view_real64_2_new = doubles_group_new.RegisterViewWrapper<real64_array>("doubles2");

  ManagedGroup & mixed_group_new = doubles_group_new.RegisterGroup("mixed");
  ViewWrapper<int32_array> & view_int32_new = mixed_group_new.RegisterViewWrapper<int32_array>("int32");
  ViewWrapper<real32_array> & view_real32_new = mixed_group_new.RegisterViewWrapper<real32_array>("real32");
  ViewWrapper<string> & view_whatup_new = mixed_group_new.RegisterViewWrapper<string>("whatup");
  ViewWrapper<real64> & view_pi_new = mixed_group_new.RegisterViewWrapper<real64>("pi");


  /* Load the data */
  root->loadSidreExternalData(path + ".root", MPI_COMM_WORLD);


  /* Group sizes should have carried over. */
  EXPECT_TRUE(root->size() == group_size);
  EXPECT_TRUE(strings_group_new.size() == group_size);
  EXPECT_TRUE(doubles_group_new.size() == group_size);
  EXPECT_TRUE(mixed_group_new.size() == group_size);

  

  /* Should be the same as stored. */
  EXPECT_TRUE(view_int64_new.size() == view_int64_size);
  for (int i = 0; i < view_int64_new.size(); i++) {
    EXPECT_TRUE(view_int64_new.data()[i] == i);
  }

  EXPECT_TRUE(view_hope_new.data().compare(hope) == 0);

  EXPECT_TRUE(view_hello_new.data().compare(hello) == 0);

  EXPECT_TRUE(view_goodbye_new.data().compare(goodbye) == 0);

  EXPECT_TRUE(view_real64_1_new.size() == view_real64_1_size);
  for (int i = 0; i < view_real64_1_new.size(); i++) {
    EXPECT_TRUE(view_real64_1_new.data()[i] == double(i) / 100.0);
  }

  EXPECT_TRUE(view_real64_2_new.size() == view_real64_2_size);
  for (int i = 0; i < view_real64_2_new.size(); i++) {
    EXPECT_TRUE(view_real64_2_new.data()[i] == double(i) * double(i) / 100.0);
  }

  EXPECT_TRUE(view_int32_new.size() == view_int32_size);
  for (int i = 0; i < view_int32_new.size(); i++) {
    EXPECT_TRUE(view_int32_new.data()[i] == i * i);
  }

  EXPECT_TRUE(view_real32_new.size() == view_real32_size);
  for (int i = 0; i < view_real32_new.size(); i++) {
    EXPECT_TRUE(view_real32_new.data()[i] == float(i) / float(i + 1.0 - i * i));
  }

  EXPECT_TRUE(view_whatup_new.data().compare(whatup) == 0);

  // EXPECT_TRUE(*(view_pi_new.data()) == 3.14159);

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
