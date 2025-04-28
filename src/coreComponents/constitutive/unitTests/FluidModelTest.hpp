/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluidModelTest.hpp
 */

#ifndef GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_
#define GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "TestFluid.hpp"

#include <conduit.hpp>
#include <gtest/gtest.h>

namespace geos
{
namespace testing
{

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE = 2 >
class FluidModelTest : public ::testing::Test
{
public:
  using FluidModel = FLUID_TYPE;
  using FluidWrapper = typename FLUID_TYPE::KernelWrapper;
  static constexpr integer numComp = NUM_COMP;
  static constexpr integer numDof = numComp + 2;
  static constexpr integer numPhase = NUM_PHASE;

public:
  using TestPoint = std::tuple<
    real64 const,           // pressure
    real64 const,           // temperature
    Feed< NUM_COMP > const  // composition
    >;

public:
  FluidModelTest();
  ~FluidModelTest() override = default;

  FluidModel * getFluid( string const & name );

  template< typename LAMBDA >
  FluidModel * createFluid( string const & name, LAMBDA && function );

  void testNumericalDerivatives( FluidModel * fluid, TestPoint const & data );

protected:
  static void populateLinearScale( array1d< real64 > & array, real64 const x0, real64 const x1, integer const & n );
  static void populateLogScale( array1d< real64 > & array, real64 const x0, real64 const x1, integer const & n );

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
};

} // namespace testing

} //namespace geos

#include "FluidModelTest_impl.hpp"

#endif // GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_
