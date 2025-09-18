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

#include "fileIO/timeHistory/ScatterDataProvider.hpp"
#include <gtest/gtest.h>

using namespace geos;

class MinimalScatterDataProvider : public ScatterDataProvider
{
public:
  MinimalScatterDataProvider()
  {
    m_scatterData.resize( 0 );
    m_coordinates.resize( 0, 3 );
    m_metadata.resize( 0 );
  }
  localIndex getNumScatterPoints() const override { return 0; }
  array1d< real64 > const & getScatterData() const override { return m_scatterData; }
  array2d< real64 > const & getScatterCoordinates() const override { return m_coordinates; }
  string_array const & getScatterMetadata() const override { return m_metadata; }
private:
  array1d< real64 > m_scatterData;
  array2d< real64 > m_coordinates;
  string_array m_metadata;
};

TEST( ScatterDataProviderTest, ZeroScatterPoints )
{
  MinimalScatterDataProvider provider;
  EXPECT_EQ( provider.getNumScatterPoints(), 0 );
  EXPECT_EQ( provider.getScatterData().size(), 0 );
  EXPECT_EQ( provider.getScatterCoordinates().size( 0 ), 0 );
  EXPECT_EQ( provider.getScatterCoordinates().size( 1 ), 3 );
  EXPECT_EQ( provider.getScatterMetadata().size(), 0 );
}

class EdgeCaseScatterDataProvider : public ScatterDataProvider
{
public:
  EdgeCaseScatterDataProvider()
  {
    m_scatterData.resize( 2 );
    m_scatterData[0] = -1e10; m_scatterData[1] = 1e10;
    m_coordinates.resize( 2, 3 );
    m_coordinates( 0, 0 ) = -1e5; m_coordinates( 0, 1 ) = 0.0; m_coordinates( 0, 2 ) = 1e5;
    m_coordinates( 1, 0 ) = 1e-5; m_coordinates( 1, 1 ) = -1e-5; m_coordinates( 1, 2 ) = 0.0;
    m_metadata.resize( 2 );
    m_metadata[0] = "";
    m_metadata[1] = "Special_Station_#1";
  }
  localIndex getNumScatterPoints() const override { return 2; }
  array1d< real64 > const & getScatterData() const override { return m_scatterData; }
  array2d< real64 > const & getScatterCoordinates() const override { return m_coordinates; }
  string_array const & getScatterMetadata() const override { return m_metadata; }
private:
  array1d< real64 > m_scatterData;
  array2d< real64 > m_coordinates;
  string_array m_metadata;
};

TEST( ScatterDataProviderTest, EdgeCases )
{
  EdgeCaseScatterDataProvider provider;
  EXPECT_EQ( provider.getNumScatterPoints(), 2 );
  EXPECT_DOUBLE_EQ( provider.getScatterData()[0], -1e10 );
  EXPECT_DOUBLE_EQ( provider.getScatterData()[1], 1e10 );
  EXPECT_DOUBLE_EQ( provider.getScatterCoordinates()( 0, 0 ), -1e5 );
  EXPECT_DOUBLE_EQ( provider.getScatterCoordinates()( 0, 2 ), 1e5 );
  EXPECT_EQ( provider.getScatterMetadata()[0], "" );
  EXPECT_EQ( provider.getScatterMetadata()[1], "Special_Station_#1" );
}

int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
