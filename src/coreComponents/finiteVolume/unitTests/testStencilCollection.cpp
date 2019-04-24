/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include <gtest/gtest.h>

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "finiteVolume/FluxStencil.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include <chrono>
#include <iostream>

#define TEST_SIZE 1000000

using namespace geosx;

struct Cell
{
  Cell(): er(0), esr(0), ei(0) {}
  Cell(localIndex i): er(0), esr(0), ei(i) {}

  localIndex er;
  localIndex esr;
  localIndex ei;
};

template<typename INDEX, typename T>
void makeStencilTPFA(localIndex size, FluxStencil<INDEX, T> & stencil)
{
  stencil.reserve(size, 2);

  stackArray1d<INDEX, 2> cells(2);
  stackArray1d<real64, 2> weights(2);

  for (localIndex kf = 0; kf < size; ++kf)
  {
    cells[0] = INDEX(kf);
    cells[1] = INDEX(kf+1);
    weights[0] = 1.0/2.0;
    weights[1] = 1.0/2.0;

    stencil.add(2, cells.data(), weights.data(), kf);
  }
  stencil.compress();
}

template<typename LAMBDA>
void testStencilLoop( LAMBDA && compute )
{
  FluxStencil<Cell, double> stencil;
  makeStencilTPFA(TEST_SIZE, stencil);

  array1d<double> in(TEST_SIZE + 1);
  array1d<double> out(TEST_SIZE);
  for (localIndex i = 0; i < TEST_SIZE; ++i)
  {
    in[i] = double(i);
  }
  out = 0.0;

  auto start = std::chrono::high_resolution_clock::now();
  compute( in, out, stencil );
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff = end - start;
  std::cout << "Time elapsed: " << diff.count() << " s" << std::endl;

  for (localIndex kf = 0; kf < TEST_SIZE; ++kf)
  {
    EXPECT_FLOAT_EQ( out[kf], (in[kf]/2 + in[kf+1]/2) );
  }
}

TEST(testStencilCollection, noAcessorIterationTPFA)
{
  testStencilLoop( [&] ( arrayView1d<double const> const & dataIn,
                         arrayView1d<double> const & dataOut,
                         FluxStencil<Cell, double> const & stencil )
  {
    csArrayView2d<FluxStencil<Cell, double>::Entry const> const & connections = stencil.getConnections();

    forall_in_range<stencilPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      for (localIndex i = 0; i < connections.size(iconn); ++i)
      {
        FluxStencil<Cell, double>::Entry const & entry = connections(iconn, i);
        dataOut[iconn] += entry.weight * dataIn[entry.index.ei];
      }
    } );
  } );
}

int main( int argc, char* argv[] )
{
  logger::InitializeLogger();

  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();

  logger::FinalizeLogger();
  return result;
}
