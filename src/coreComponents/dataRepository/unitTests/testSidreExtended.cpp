/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include <gtest/gtest.h>

#ifdef __clang__
#pragma clang diagnostic push
#endif

#include <mpi.h>
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "common/DataTypes.hpp"


namespace geosx {
namespace dataRepository {

#ifdef GEOSX_USE_ATK
template<typename T> 
ViewWrapper<array1d<T>> * createArrayView(ManagedGroup * parent, const string & name,
                                        int sfp, const array1d<T> & data)
{
  ViewWrapper<array1d<T>> * view = parent->RegisterViewWrapper<array1d<T>>(name);
  view->setSizedFromParent(sfp);

  /* Resize the array */
  localIndex expected_size = data.size() * sizeof(T);
  view->resize(data.size());

  /* Check that the ViewWrapper size and byteSize return the proper values */
  EXPECT_EQ(view->size(), data.size());
  EXPECT_EQ(view->byteSize(), expected_size);

  /* Set the data */
  view_rtype<T> view_data = view->data();
  for (int i = 0; i < view->size(); i++) {
    view_data[i] = data[i];
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view->dataPtr(), &(view->data()[0]));

  return view;
}


template<typename T> 
void checkArrayView(const ViewWrapper<array1d<T>> * view, int sfp, const array1d<T> & data) 
{
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(view->size(), data.size());
  view_rtype_const<T> view_data = view->data();
  for (int i = 0; i < view->size(); i++) {
    EXPECT_EQ(view_data[i], data[i]);
  }
}

template<typename T> 
ViewWrapper<array2d<T>> * createArray2dView(ManagedGroup * parent, const string & name,
                                          int sfp, const array2d<T> & data)
{
  ViewWrapper<array2d<T>> * view = parent->RegisterViewWrapper<array2d<T>>(name);
  view->setSizedFromParent(sfp);

  /* Resize the array */
  localIndex expected_size = data.size() * sizeof(T);
  long long dims[2];
  dims[0] = data.size(0);
  dims[1] = data.size(1);
  view->resize(2, dims);

  /* Check that the ViewWrapper size and byteSize return the proper values */
  EXPECT_EQ(view->size(0), data.size(0));
  EXPECT_EQ(view->size(1), data.size(1));
  EXPECT_EQ(view->size(), data.size());
  EXPECT_EQ(view->byteSize(), expected_size);

  /* Set the data */
  array2d<T> & view_data = view->reference();
  for (int i = 0; i < dims[0]; i++) 
  {
    for (int j = 0; j < dims[1]; j++) 
    {
      view_data[i][j] = data[i][j];
    }
  }

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_TRUE(view->dataPtr() == &view_data[0][0]);

  return view;
}


template<typename T> 
void checkArray2dView(const ViewWrapper<array2d<T>> * view, int sfp, const array2d<T> & data) 
{
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(view->size(), data.size());
  EXPECT_EQ(view->size(0), data.size(0));
  EXPECT_EQ(view->size(1), data.size(1));
  
  const array2d<T> & view_data = view->reference();
  for (int i = 0; i < data.size(0); i++) 
  {
    for (int j = 0; j < data.size(1); j++) 
    {
      EXPECT_EQ(view_data[i][j], data[i][j]);
    }
  }

}


template<typename T> 
ViewWrapper<set<T>> * createSetView(ManagedGroup * parent, const string & name,
                                    localIndex sfp, const set<T> & data)
{
  ViewWrapper<set<T>> * view = parent->RegisterViewWrapper<set<T>>(name);
  view->setSizedFromParent(int(sfp));

  /* Resize the array */
  localIndex expected_size = data.size() * sizeof(T);

  /* Set the data */
  view->reference().insert(data.begin(), data.end());

  /* Check that the ViewWrapper size and byteSize return the proper values */
  EXPECT_EQ(view->size(), data.size());
  EXPECT_EQ(view->byteSize(), expected_size);

  /* Check that the set is sorted */
  EXPECT_TRUE(view->reference().isSorted());

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view->dataPtr(), &(view->data()[0]));

  return view;
}


template<typename T> 
void checkSetView(const ViewWrapper<set<T>> * view, localIndex sfp, const set<T> & data) 
{
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(view->size(), data.size());
  const set<T> & view_data = view->reference();
  for (int i = 0; i < view->size(); i++) {
    EXPECT_EQ(view_data[i], data[i]);
  }

  /* Check that the set is sorted */
  EXPECT_TRUE(view->reference().isSorted());
}


ViewWrapper<string> * createStringView(ManagedGroup * parent, const string & name,
                                       int sfp, const string & str) 
{
  ViewWrapper<string> * view = parent->RegisterViewWrapper<string>(name);
  view->setSizedFromParent(sfp);
  
  localIndex expected_size = static_cast<localIndex>(str.size()) * sizeof(char);

  /* Set the data */
  view->data() = str;

  /* Check that the ViewWrapper size and byteSize return the proper values */
  EXPECT_EQ(static_cast<uint>(view->size()), str.size());
  EXPECT_EQ(view->byteSize(), expected_size);

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view->dataPtr(), view->data().c_str());

  return view;
}


void checkStringView(const ViewWrapper<string> * view, const int sfp, const string & str) {
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(view->reference(), str);
}


ViewWrapper<string_array> * createStringArrayView(ManagedGroup * parent, const string & name,
                                                int sfp, const string_array & arr)
{
  ViewWrapper<string_array> * view = parent->RegisterViewWrapper<string_array>(name);
  view->setSizedFromParent(sfp);

  unsigned long expected_size = arr.size() * sizeof(string);
  view->resize( static_cast<localIndex>(arr.size()));

  EXPECT_EQ(static_cast<uint>(view->size()), arr.size());
  EXPECT_EQ(view->byteSize(), expected_size);

  view_rtype<string_array> view_data = view->data();
  for (string_array::size_type i = 0; i < arr.size(); ++i)
  {
    view_data[i] = arr[i];
  }
  
  EXPECT_EQ(view->dataPtr(), &(view->data()[0]));
  return view;
}


void checkStringArrayView(const ViewWrapper<string_array> * view, const int sfp, const string_array & arr)
{
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(view->size(), arr.size());
  view_rtype_const<string_array> view_data = view->data();
  for (int i = 0; i < view->size(); i++) {
    EXPECT_EQ(view_data[i], arr[i]);
  }
}


template<typename T>
ViewWrapper<T> * createScalarView(ManagedGroup * parent, const string & name, 
                                  int sfp, const T value) {
  ViewWrapper<T> * view = parent->RegisterViewWrapper<T>(name);
  view->setSizedFromParent(sfp);

  /* Set the data */
  *(view->data()) = value;

  /* Check that the ViewWrapper size and byteSize return the proper values */
  EXPECT_EQ(view->size(), 1);
  EXPECT_EQ(view->byteSize(), sizeof(T));

  /* Check that the ViewWrapper dataPtr points to the right thing */
  EXPECT_EQ(view->dataPtr(), view->data());

  return view;
}


template<typename T>
void checkScalarView(const ViewWrapper<T> * view, int sfp, const T value) {
  EXPECT_EQ(view->sizedFromParent(), sfp);
  EXPECT_EQ(*(view->data()), value);
}


TEST(testSidreExtended, testSidreExtended)
{
  MPI_Init( nullptr, nullptr);
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  const string path = "test_sidre_extended";
  const string protocol = "sidre_hdf5";
  const int group_size = 44;
  int sfp = 55;
  axom::sidre::DataStore & ds = SidreWrapper::dataStore();

  /* Create a new ManagedGroup directly below the sidre::DataStore root. */
  ManagedGroup * root = new ManagedGroup(std::string("data"), nullptr);
  root->resize(group_size);

  /* Create a new globalIndex_array ViewWrapper. */
  string view_globalIndex_name = "globalIndex";
  int view_globalIndex_sfp = sfp++;
  int view_globalIndex_size = 100;
  globalIndex_array view_globalIndex_data(view_globalIndex_size);
  for (int i = 0; i < view_globalIndex_size; i++) {
    view_globalIndex_data[i] = i * i * i;
  }
  createArrayView(root, view_globalIndex_name, view_globalIndex_sfp, view_globalIndex_data);

  /* Create a new string ViewWrapper. */
  string view_hope_name = "hope";
  int view_hope_sfp = sfp++;
  string view_hope_str = "I sure hope these tests pass.";
  createStringView(root, view_hope_name, view_hope_sfp, view_hope_str);

  /* Create a new string array ViewWrapper. */
  string view_restart_name = "restart";
  int view_restart_sfp = sfp++;
  string_array view_restart_arr(6);
  view_restart_arr[0] = "Man ";
  view_restart_arr[1] = "this ";
  view_restart_arr[2] = "restart ";
  view_restart_arr[3] = "stuff ";
  view_restart_arr[4] = "better ";
  view_restart_arr[5] = "work.";
  createStringArrayView(root, view_restart_name, view_restart_sfp, view_restart_arr);

  /* Create a new group. */
  ManagedGroup * strings_group = root->RegisterGroup("strings");
  strings_group->resize(group_size + 1);

  /* Create a new string ViewWrapper. */
  string view_hello_name = "hello";
  int view_hello_sfp = sfp++;
  string view_hello_str = "Hello, how are you doing on this fine day?";
  createStringView(strings_group, view_hello_name, view_hello_sfp,
                   view_hello_str);

  /* Create a new string ViewWrapper. */
  string view_goodbye_name = "goodbye";
  int view_goodbye_sfp = sfp++;
  string view_goodbye_str = "I hate this weather so I'm heading inside. Goodbye.";
  createStringView(strings_group, view_goodbye_name, view_goodbye_sfp, 
                   view_goodbye_str);

  /* Create a new group. */
  ManagedGroup * real64_group = root->RegisterGroup("real64");
  real64_group->resize(group_size + 2);

  /* Create a new real64_array ViewWrapper. */
  string view_real641_name = "real641";
  int view_real641_sfp = sfp++;
  int view_real641_size = 1000;
  real64_array view_real641_data(view_real641_size);
  for (int i = 0; i < view_real641_size; i++) {
    view_real641_data[i] = i * i / (i + 5.0);
  }
  createArrayView(real64_group, view_real641_name, view_real641_sfp, 
                  view_real641_data);

  /* Create a new real64_array ViewWrapper. */
  string view_real642_name = "real642";
  int view_real642_sfp = sfp++;
  int view_real642_size = 1000;
  real64_array view_real642_data(view_real642_size);
  for (int i = 0; i < view_real642_size; i++) {
    view_real642_data[i] = i * i / (5.0 + 5.0 * i + i * i);
  }
  createArrayView(real64_group, view_real642_name, view_real642_sfp, 
                  view_real642_data);

  /* Create a new group. */
  ManagedGroup * mixed_group = real64_group->RegisterGroup("mixed");
  mixed_group->resize(group_size + 3);

  /* Create a new localIndex_array ViewWrapper. */
  string view_localIndex_name = "localIndex";
  int view_localIndex_sfp = sfp++;
  int view_localIndex_size = 953;
  localIndex_array view_localIndex_data(view_localIndex_size);
  for (localIndex i = 0; i < view_localIndex_size; i++) {
    view_localIndex_data[i] = i * i - 100 * i + 3;
  }
  createArrayView(mixed_group, view_localIndex_name, view_localIndex_sfp, view_localIndex_data);

  /* Create a new real32_array ViewWrapper. */
  string view_real32_name = "real32";
  int view_real32_sfp = sfp++;
  int view_real32_size = 782;
  real32_array view_real32_data(view_real32_size);
  for (int i = 0; i < view_real32_size; i++) {
    view_real32_data[i] = (i * i - 100.0f * i + 3.0f) / (i + 3.0f);
  }
  createArrayView(mixed_group, view_real32_name, view_real32_sfp, 
                  view_real32_data);

  /* Create a new string ViewWrapper. */
  string view_what_name = "what";
  int view_what_sfp = sfp++;
  string view_what_str = "What are you talking about? Who doesn't like storms?";
  createStringView(mixed_group, view_what_name, view_what_sfp, view_what_str);

  /* Create a new real64 ViewWrapper. */
  string view_pi_name = "pi";
  int view_pi_sfp = sfp++;
  real64 view_pi_value = 3.14159;
  createScalarView(mixed_group, view_pi_name, view_pi_sfp, view_pi_value);


  /* Create a new set<int> ViewWrapper. */
  string view_setlocalIndex_name = "view_setlocalIndex";
  localIndex view_setlocalIndex_sfp = sfp++;
  set<localIndex> view_setlocalIndex_set;
  view_setlocalIndex_set.insert(4);
  view_setlocalIndex_set.insert(3);
  view_setlocalIndex_set.insert(10);
  view_setlocalIndex_set.insert(0);
  view_setlocalIndex_set.insert(1000);
  view_setlocalIndex_set.insert(-1000);
  createSetView(mixed_group, view_setlocalIndex_name, view_setlocalIndex_sfp, view_setlocalIndex_set);


  /* Create a new set<string> ViewWrapper. */
  string view_setString_name = "view_setString";
  int view_setString_sfp = sfp++;
  set<string> view_setString_set;
  view_setString_set.insert("zaa");
  view_setString_set.insert("aab");
  view_setString_set.insert("QD");
  view_setString_set.insert("az");
  view_setString_set.insert("g");
  view_setString_set.insert("this better be sorted");
  createSetView(mixed_group, view_setString_name, view_setString_sfp, view_setString_set);


  string view_real642d_name = "view_real642d";
  int view_real642d_sfp = sfp++;
  localIndex dim0 = 10;
  localIndex dim1 = 20;
  array2d<real64> view_real642d_arr(dim0, dim1);
  for (localIndex i = 0; i < dim0; i++)
  {
    for (localIndex j = 0; j < dim1; j++)
    {
      view_real642d_arr[i][j] = i * i / 3.0 + j;
    }
  }
  createArray2dView(mixed_group, view_real642d_name, view_real642d_sfp, view_real642d_arr);


  string view_r1t2d_name = "view_r1t2d";
  int view_r1t2d_sfp = sfp++;
  dim0 = 10;
  dim1 = 20;
  array2d<R1Tensor> view_r1t2d_arr(dim0, dim1);
  for (localIndex i = 0; i < dim0; i++)
  {
    for (localIndex j = 0; j < dim1; j++)
    {
      for (localIndex k = 0; k < 3; k++) 
      {
        view_r1t2d_arr[i][j] = i * i / 3.0 + j + i * j * k / 7.0;
      }
    }
  }
  createArray2dView(mixed_group, view_r1t2d_name, view_r1t2d_sfp, view_r1t2d_arr);



  /* Save the sidre tree */
  root->prepareToWrite();
  SidreWrapper::writeTree(1, path, protocol, MPI_COMM_GEOSX);
  root->finishWriting();

  /* Delete geos tree and reset sidre tree. */
  delete root;
  ds.destroyAllAttributes();
  ds.destroyAllBuffers();
  ds.getRoot()->destroyGroups();

  /* Restore the sidre tree */
  SidreWrapper::reconstructTree(path + ".root", protocol, MPI_COMM_GEOSX);
  root = new ManagedGroup(std::string("data"), nullptr);

  /* Create dual GEOS tree. ManagedGroups automatically register with the associated sidre::View. */
  ViewWrapper<globalIndex_array> * view_globalIndex_new = root->RegisterViewWrapper<globalIndex_array>(view_globalIndex_name);
  ViewWrapper<string> * view_hope_new = root->RegisterViewWrapper<string>(view_hope_name);
  ViewWrapper<string_array> * view_restart_new = root->RegisterViewWrapper<string_array>(view_restart_name);

  ManagedGroup * strings_group_new = root->RegisterGroup("strings");
  ViewWrapper<string> * view_hello_new = strings_group_new->RegisterViewWrapper<string>(view_hello_name);
  ViewWrapper<string> * view_goodbye_new = strings_group_new->RegisterViewWrapper<string>(view_goodbye_name);
  
  ManagedGroup * real64_group_new = root->RegisterGroup("real64");
  ViewWrapper<real64_array> * view_real641_new = real64_group_new->RegisterViewWrapper<real64_array>(view_real641_name);
  ViewWrapper<real64_array> * view_real642_new = real64_group_new->RegisterViewWrapper<real64_array>(view_real642_name);

  ManagedGroup * mixed_group_new = real64_group_new->RegisterGroup("mixed");
  ViewWrapper<localIndex_array> * view_localIndex_new = mixed_group_new->RegisterViewWrapper<localIndex_array>(view_localIndex_name);
  ViewWrapper<real32_array> * view_real32_new = mixed_group_new->RegisterViewWrapper<real32_array>(view_real32_name);
  ViewWrapper<string> * view_what_new = mixed_group_new->RegisterViewWrapper<string>(view_what_name);
  ViewWrapper<real64> * view_pi_new = mixed_group_new->RegisterViewWrapper<real64>(view_pi_name);
  ViewWrapper<set<localIndex>> * view_setlocalIndex_new = mixed_group_new->RegisterViewWrapper<set<localIndex>>(view_setlocalIndex_name);
  ViewWrapper<set<string>> * view_setString_new = mixed_group_new->RegisterViewWrapper<set<string>>(view_setString_name);
  ViewWrapper<array2d<real64>> * view_real642d_new = mixed_group_new->RegisterViewWrapper<array2d<real64>>(view_real642d_name);
  ViewWrapper<array2d<R1Tensor>> * view_r1t2d_new = mixed_group_new->RegisterViewWrapper<array2d<R1Tensor>>(view_r1t2d_name);


  /* Load the data */
  root->prepareToRead();
  SidreWrapper::loadExternalData(path + ".root", MPI_COMM_GEOSX);
  root->finishReading();

  /* Group sizes should have carried over. */
  EXPECT_EQ(root->size(), group_size);
  EXPECT_EQ(strings_group_new->size(), group_size + 1);
  EXPECT_EQ(real64_group_new->size(), group_size + 2);
  EXPECT_EQ(mixed_group_new->size(), group_size + 3);

  /* Check that ViewWrapper values were restored. */
  checkArrayView(view_globalIndex_new, view_globalIndex_sfp, view_globalIndex_data);  
  checkStringView(view_hope_new, view_hope_sfp, view_hope_str);
  checkStringArrayView(view_restart_new, view_restart_sfp, view_restart_arr);
  checkStringView(view_hello_new, view_hello_sfp, view_hello_str);
  checkStringView(view_goodbye_new, view_goodbye_sfp, view_goodbye_str);
  checkArrayView(view_real641_new, view_real641_sfp, view_real641_data);  
  checkArrayView(view_real642_new, view_real642_sfp, view_real642_data);  
  checkArrayView(view_localIndex_new, view_localIndex_sfp, view_localIndex_data);  
  checkArrayView(view_real32_new, view_real32_sfp, view_real32_data);  
  checkStringView(view_what_new, view_what_sfp, view_what_str);
  checkScalarView(view_pi_new, view_pi_sfp, view_pi_value);
  checkSetView(view_setlocalIndex_new, view_setlocalIndex_sfp, view_setlocalIndex_set);
  checkSetView(view_setString_new, view_setString_sfp, view_setString_set);
  checkArray2dView(view_real642d_new, view_real642d_sfp, view_real642d_arr);
  checkArray2dView(view_r1t2d_new, view_r1t2d_sfp, view_r1t2d_arr);


  delete root;
  MPI_Finalize();
}

#endif /* GEOSX_USE_ATK */

int main(int argc, char* argv[]) {
  int result = 0;
  testing::InitGoogleTest(&argc, argv);
  result = RUN_ALL_TESTS();
  return result;
}


} /* end namespace dataRepository */
} /* end namespace goesx */
