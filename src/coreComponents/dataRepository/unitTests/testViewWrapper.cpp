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

#include <gtest/gtest.h>


#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/ManagedGroup.hpp"
using namespace geosx;
using namespace dataRepository;


TEST(testViewWrapper, testSetters)
{
  ManagedGroup group("group",nullptr);
  ViewWrapper<int> wrapper( "wrapper", &group );
  ViewWrapperBase * wrapperBasePtr = &wrapper;

  {
    {
      auto rval = wrapper.setSizedFromParent(true);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.sizedFromParent() );
    }
    {
      auto rval = wrapper.setSizedFromParent(false);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_FALSE( wrapper.sizedFromParent() );
    }

    {
      auto rval = wrapperBasePtr->setSizedFromParent(true);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->sizedFromParent() );
    }
    {
      auto rval = wrapperBasePtr->setSizedFromParent(false);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_FALSE( wrapperBasePtr->sizedFromParent() );
    }
  }

  {
    {
      auto rval = wrapper.setRestartFlags(RestartFlags::NO_WRITE);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getRestartFlags()==RestartFlags::NO_WRITE );
    }
    {
      auto rval = wrapper.setRestartFlags(RestartFlags::WRITE_AND_READ);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getRestartFlags()==RestartFlags::WRITE_AND_READ );
    }

    {
      auto rval = wrapperBasePtr->setRestartFlags(RestartFlags::NO_WRITE);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getRestartFlags()==RestartFlags::NO_WRITE );
    }
    {
      auto rval = wrapperBasePtr->setRestartFlags(RestartFlags::WRITE_AND_READ);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getRestartFlags()==RestartFlags::WRITE_AND_READ );
    }
  }

  {
    {
      auto rval = wrapper.setPlotLevel(PlotLevel::LEVEL_0);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getPlotLevel()==PlotLevel::LEVEL_0 );
    }
    {
      auto rval = wrapper.setPlotLevel(PlotLevel::LEVEL_1);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getPlotLevel()==PlotLevel::LEVEL_1 );
    }

    {
      auto rval = wrapperBasePtr->setPlotLevel(PlotLevel::LEVEL_0);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getPlotLevel()==PlotLevel::LEVEL_0 );
    }
    {
      auto rval = wrapperBasePtr->setPlotLevel(PlotLevel::LEVEL_1);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getPlotLevel()==PlotLevel::LEVEL_1 );
    }
  }

  {
    {
      auto rval = wrapper.setInputFlag(InputFlags::OPTIONAL);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getInputFlag()==InputFlags::OPTIONAL );
    }
    {
      auto rval = wrapper.setInputFlag(InputFlags::REQUIRED);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getInputFlag()==InputFlags::REQUIRED );
    }

    {
      auto rval = wrapperBasePtr->setInputFlag(InputFlags::OPTIONAL);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getInputFlag()==InputFlags::OPTIONAL );
    }
    {
      auto rval = wrapperBasePtr->setInputFlag(InputFlags::REQUIRED);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getInputFlag()==InputFlags::REQUIRED );
    }
  }

  {
    {
      string description("Description of wrapped object 1");
      auto rval = wrapper.setDescription(description);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getDescription()==description );
    }
    {
      string description("Description of wrapped object 2");
      auto rval = wrapper.setDescription(description);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapper<int> * >::value) );
      EXPECT_TRUE( wrapper.getDescription()==description );
    }

    {
      string description("Description of wrapped object 3");
      auto rval = wrapperBasePtr->setDescription(description);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getDescription()==description );
    }
    {
      string description("Description of wrapped object 4");
      auto rval = wrapperBasePtr->setDescription(description);
      EXPECT_TRUE( (std::is_same< decltype(rval),ViewWrapperBase * >::value) );
      EXPECT_TRUE( wrapperBasePtr->getDescription()==description );
    }
  }
}
