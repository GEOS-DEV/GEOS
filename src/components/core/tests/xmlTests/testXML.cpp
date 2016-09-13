/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#endif


std::string filename;

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  for( int i=0 ; i<argc ; ++i )
  {
    std::cout<<argv[i]<<std::endl;
  }
  filename = argv[1];
  return RUN_ALL_TESTS();
}

TEST(testXML,testXML)
{
  std::cout<<filename<<std::endl;

}
