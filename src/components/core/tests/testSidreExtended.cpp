#include <gtest/gtest.h>
#include <mpi.h>
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "common/DataTypes.hpp"


namespace geosx {
namespace dataRepository {

#ifdef USE_ATK
template<typename T> 
ViewWrapper<T> & createArrayView(ManagedGroup * parent, const string name,
                                 int sfp, const T & data)
{
  ViewWrapper<T> & view = *parent->RegisterViewWrapper<T>(name);
  view.setSizedFromParent(sfp);

  /* Resize the array */
  localIndex expected_size = data.size() * sizeof(typename T::value_type);
  view.resize(data.size());

  /* Check that the ViewWrapper size and dataSize return the proper values */
  EXPECT_EQ(view.size(), data.size());
  EXPECT_EQ(view.dataSize(), expected_size);

  /* Set the data */
  view_rtype<T> view_data = view.data();
  for (int i = 0; i < view.size(); i++) {
      view_data[i] = data[i];
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view.dataPtr(), &(view.data()[0]));

  return view;
}


template<typename T> 
void checkArrayView(ViewWrapper<T> & view, int sfp, const T & data) 
{
  EXPECT_EQ(view.sizedFromParent(), sfp);
  EXPECT_EQ(view.size(), data.size());
  view_rtype<T> view_data = view.data();
  for (int i = 0; i < view.size(); i++) {
    EXPECT_DOUBLE_EQ(view_data[i], data[i]);
  }
}


ViewWrapper<string> & createStringview(ManagedGroup * parent, const string name,
                                       int sfp, const string str) 
{
  ViewWrapper<string> & view = *parent->RegisterViewWrapper<string>(name);
  view.setSizedFromParent(sfp);
  
  localIndex expected_size = str.size() * sizeof(char);

  /* Set the data */
  view.data() = str;

  /* Check that the ViewWrapper size and dataSize return the proper values */
  EXPECT_EQ(view.size(), str.size());
  EXPECT_EQ(view.dataSize(), expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view.dataPtr(), view.data().c_str());

  return view;
}


void checkStringView(ViewWrapper<string> & view, const int sfp, const string str) {
  EXPECT_EQ(view.sizedFromParent(), sfp);
  EXPECT_EQ(view.data().compare(str), 0);
}


template<typename T>
ViewWrapper<T> & createScalarView(ManagedGroup * parent, const string name, 
                                  int sfp, const T value) {
  ViewWrapper<T> & view = *parent->RegisterViewWrapper<T>(name);
  view.setSizedFromParent(sfp);

  /* Set the data */
  *(view.data()) = value;

  /* Check that the ViewWrapper size and dataSize return the proper values */
  EXPECT_EQ(view.size(), 1);
  EXPECT_EQ(view.dataSize(), sizeof(T));

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view.dataPtr(), view.data());

  return view;
}


template<typename T>
void checkScalarView(ViewWrapper<T> & view, int sfp, const T value) {
  EXPECT_EQ(view.sizedFromParent(), sfp);
  EXPECT_EQ(*(view.data()), value);
}


//
//TEST(testSidreExtended, testSidreExtended) {
//  MPI_Init(0, nullptr);
//  const string path = "test_sidre_extended";
//  const string protocol = "sidre_hdf5";
//  const int group_size = 44;
//  int sfp = 55;
//  axom::sidre::DataStore & ds = SidreWrapper::dataStore();
//
//  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
//  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);
//  root->resize(group_size);
//
//  /* Create a new int64_array ViewWrapper. */
//  string view_int64_name = "int64";
//  int view_int64_sfp = sfp++;
//  int view_int64_size = 100;
//  int64_array view_int64_data(view_int64_size);
//  for (int i = 0; i < view_int64_size; i++) {
//    view_int64_data[i] = i * i * i;
//  }
//  createArrayView(root, view_int64_name, view_int64_sfp, view_int64_data);
//
//  /* Create a new string ViewWrapper. */
//  string view_hope_name = "hope";
//  int view_hope_sfp = sfp++;
//  string view_hope_str = "I sure hope these tests pass.";
//  createStringview(root, view_hope_name, view_hope_sfp, view_hope_str);
//
//  /* Create a new group. */
//  ManagedGroup * strings_group = root->RegisterGroup("strings");
//  strings_group->resize(group_size + 1);
//
//  /* Create a new string ViewWrapper. */
//  string view_hello_name = "hello";
//  int view_hello_sfp = sfp++;
//  string view_hello_str = "Hello, how are you doing on this fine day?";
//  createStringview(strings_group, view_hello_name, view_hello_sfp,
//                   view_hello_str);
//
//  /* Create a new string ViewWrapper. */
//  string view_goodbye_name = "goodbye";
//  int view_goodbye_sfp = sfp++;
//  string view_goodbye_str = "I hate this weather so I'm heading inside. Goodbye.";
//  createStringview(strings_group, view_goodbye_name, view_goodbye_sfp,
//                   view_goodbye_str);
//
//  /* Create a new group. */
//  ManagedGroup * real64_group = root->RegisterGroup("real64");
//  real64_group->resize(group_size + 2);
//
//  /* Create a new real64_array ViewWrapper. */
//  string view_real641_name = "real641";
//  int view_real641_sfp = sfp++;
//  int view_real641_size = 1000;
//  real64_array view_real641_data(view_real641_size);
//  for (real64 i = 0; i < view_real641_size; i++) {
//    view_real641_data[i] = i * i / (i + 5);
//  }
//  createArrayView(real64_group, view_real641_name, view_real641_sfp,
//                  view_real641_data);
//
//  /* Create a new real64_array ViewWrapper. */
//  string view_real642_name = "real642";
//  int view_real642_sfp = sfp++;
//  int view_real642_size = 1000;
//  real64_array view_real642_data(view_real642_size);
//  for (real64 i = 0; i < view_real642_size; i++) {
//    view_real642_data[i] = i * i / (5 + 5 * i + i * i);
//  }
//  createArrayView(real64_group, view_real642_name, view_real642_sfp,
//                  view_real642_data);
//
//  /* Create a new group. */
//  ManagedGroup * mixed_group = real64_group->RegisterGroup("mixed");
//  mixed_group->resize(group_size + 3);
//
//  /* Create a new int32_array ViewWrapper. */
//  string view_int32_name = "int32";
//  int view_int32_sfp = sfp++;
//  int view_int32_size = 953;
//  int32_array view_int32_data(view_int32_size);
//  for (int32 i = 0; i < view_int32_size; i++) {
//    view_int32_data[i] = i * i - 100 * i + 3;
//  }
//  createArrayView(mixed_group, view_int32_name, view_int32_sfp, view_int32_data);
//
//  /* Create a new real32_array ViewWrapper. */
//  string view_real32_name = "real32";
//  int view_real32_sfp = sfp++;
//  int view_real32_size = 782;
//  real32_array view_real32_data(view_real32_size);
//  for (real32 i = 0; i < view_real32_size; i++) {
//    view_real32_data[i] = (i * i - 100 * i + 3) / (i + 3);
//  }
//  createArrayView(mixed_group, view_real32_name, view_real32_sfp,
//                  view_real32_data);
//
//  /* Create a new string ViewWrapper. */
//  string view_what_name = "what";
//  int view_what_sfp = sfp++;
//  string view_what_str = "What are you talking about? Who doesn't like storms?";
//  createStringview(mixed_group, view_what_name, view_what_sfp, view_what_str);
//
//  /* Create a new real64 ViewWrapper. */
//  string view_pi_name = "pi";
//  int view_pi_sfp = sfp++;
//  real64 view_pi_value = 3.14159;
//  createScalarView(mixed_group, view_pi_name, view_pi_sfp, view_pi_value);
//
//
//  /* Save the sidre tree */
//  root->writeRestart(1, path, protocol, MPI_COMM_WORLD);
//
//  /* Delete geos tree and reset sidre tree. */
//  delete root;
//  ds.destroyAllAttributes();
//  ds.destroyAllBuffers();
//  ds.getRoot()->destroyGroups();
//
//  /* Restore the sidre tree */
//  root = new ManagedGroup(std::string("data"), nullptr);
//  root->reconstructSidreTree(path + ".root", protocol, MPI_COMM_WORLD);
//
//  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
//  ViewWrapper<int64_array> & view_int64_new = *root->RegisterViewWrapper<int64_array>(view_int64_name);
//  ViewWrapper<string> & view_hope_new = *root->RegisterViewWrapper<string>(view_hope_name);
//
//  ManagedGroup * strings_group_new = root->RegisterGroup("strings");
//  ViewWrapper<string> & view_hello_new = *strings_group_new->RegisterViewWrapper<string>(view_hello_name);
//  ViewWrapper<string> & view_goodbye_new = *strings_group_new->RegisterViewWrapper<string>(view_goodbye_name);
//
//  ManagedGroup * real64_group_new = root->RegisterGroup("real64");
//  ViewWrapper<real64_array> & view_real641_new = *real64_group_new->RegisterViewWrapper<real64_array>(view_real641_name);
//  ViewWrapper<real64_array> & view_real642_new = *real64_group_new->RegisterViewWrapper<real64_array>(view_real642_name);
//
//  ManagedGroup * mixed_group_new = real64_group_new->RegisterGroup("mixed");
//  ViewWrapper<int32_array> & view_int32_new = *mixed_group_new->RegisterViewWrapper<int32_array>(view_int32_name);
//  ViewWrapper<real32_array> & view_real32_new = *mixed_group_new->RegisterViewWrapper<real32_array>(view_real32_name);
//  ViewWrapper<string> & view_what_new = *mixed_group_new->RegisterViewWrapper<string>(view_what_name);
//  ViewWrapper<real64> & view_pi_new = *mixed_group_new->RegisterViewWrapper<real64>(view_pi_name);
//
//  /* Load the data */
//  root->loadSidreExternalData(path + ".root", MPI_COMM_WORLD);
//
//  /* Group sizes should have carried over. */
//  EXPECT_EQ(root->size(), group_size);
//  EXPECT_EQ(strings_group_new->size(), group_size + 1);
//  EXPECT_EQ(real64_group_new->size(), group_size + 2);
//  EXPECT_EQ(mixed_group_new->size(), group_size + 3);
//
//  /* Check that ViewWrapper values were restored. */
//  checkArrayView(view_int64_new, view_int64_sfp, view_int64_data);
//  checkStringView(view_hope_new, view_hope_sfp, view_hope_str);
//  checkStringView(view_hello_new, view_hello_sfp, view_hello_str);
//  checkStringView(view_goodbye_new, view_goodbye_sfp, view_goodbye_str);
//  checkArrayView(view_real641_new, view_real641_sfp, view_real641_data);
//  checkArrayView(view_real642_new, view_real642_sfp, view_real642_data);
//  checkArrayView(view_int32_new, view_int32_sfp, view_int32_data);
//  checkArrayView(view_real32_new, view_real32_sfp, view_real32_data);
//  checkStringView(view_what_new, view_what_sfp, view_what_str);
//  checkScalarView(view_pi_new, view_pi_sfp, view_pi_value);
//
//  delete root;
//  MPI_Finalize();
//}
#endif /* ATK_FOUND */

int main(int argc, char* argv[]) {
  int result = 0;
  testing::InitGoogleTest(&argc, argv);
  result = RUN_ALL_TESTS();
  return result;
}


} /* end namespace dataRepository */
} /* end namespace goesx */
