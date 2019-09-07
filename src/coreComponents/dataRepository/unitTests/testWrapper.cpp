/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2019 Total, S.A
* Copyright (c) 2019-     GEOSX Contributors
* All right reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
*/

#include <gtest/gtest.h>
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"

using namespace geosx;
using namespace dataRepository;


TEST( testWrapper, testSetters )
{
Group group( "group", nullptr );
Wrapper< int > wrapper( "wrapper", &group );
WrapperBase * wrapperBasePtr = &wrapper;

{
{
auto rval = wrapper.setSizedFromParent( true );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.sizedFromParent() );
}
{
auto rval = wrapper.setSizedFromParent( false );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_FALSE( wrapper.sizedFromParent() );
}

{
auto rval = wrapperBasePtr->setSizedFromParent( true );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->sizedFromParent() );
}
{
auto rval = wrapperBasePtr->setSizedFromParent( false );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_FALSE( wrapperBasePtr->sizedFromParent() );
}
}

{
{
auto rval = wrapper.setRestartFlags( RestartFlags::NO_WRITE );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getRestartFlags()==RestartFlags::NO_WRITE );
}
{
auto rval = wrapper.setRestartFlags( RestartFlags::WRITE_AND_READ );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getRestartFlags()==RestartFlags::WRITE_AND_READ );
}

{
auto rval = wrapperBasePtr->setRestartFlags( RestartFlags::NO_WRITE );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getRestartFlags()==RestartFlags::NO_WRITE );
}
{
auto rval = wrapperBasePtr->setRestartFlags( RestartFlags::WRITE_AND_READ );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getRestartFlags()==RestartFlags::WRITE_AND_READ );
}
}

{
{
auto rval = wrapper.setPlotLevel( PlotLevel::LEVEL_0 );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getPlotLevel()==PlotLevel::LEVEL_0 );
}
{
auto rval = wrapper.setPlotLevel( PlotLevel::LEVEL_1 );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getPlotLevel()==PlotLevel::LEVEL_1 );
}

{
auto rval = wrapperBasePtr->setPlotLevel( PlotLevel::LEVEL_0 );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getPlotLevel()==PlotLevel::LEVEL_0 );
}
{
auto rval = wrapperBasePtr->setPlotLevel( PlotLevel::LEVEL_1 );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getPlotLevel()==PlotLevel::LEVEL_1 );
}
}

{
{
auto rval = wrapper.setInputFlag( InputFlags::OPTIONAL );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getInputFlag()==InputFlags::OPTIONAL );
}
{
auto rval = wrapper.setInputFlag( InputFlags::REQUIRED );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getInputFlag()==InputFlags::REQUIRED );
}

{
auto rval = wrapperBasePtr->setInputFlag( InputFlags::OPTIONAL );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getInputFlag()==InputFlags::OPTIONAL );
}
{
auto rval = wrapperBasePtr->setInputFlag( InputFlags::REQUIRED );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getInputFlag()==InputFlags::REQUIRED );
}
}

{
{
string description( "Description of wrapped object 1" );
auto rval = wrapper.setDescription( description );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getDescription()==description );
}
{
string description( "Description of wrapped object 2" );
auto rval = wrapper.setDescription( description );
EXPECT_TRUE( (std::is_same< decltype(rval), Wrapper< int > * >::value) );
EXPECT_TRUE( wrapper.getDescription()==description );
}

{
string description( "Description of wrapped object 3" );
auto rval = wrapperBasePtr->setDescription( description );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getDescription()==description );
}
{
string description( "Description of wrapped object 4" );
auto rval = wrapperBasePtr->setDescription( description );
EXPECT_TRUE( (std::is_same< decltype(rval), WrapperBase * >::value) );
EXPECT_TRUE( wrapperBasePtr->getDescription()==description );
}
}
}
