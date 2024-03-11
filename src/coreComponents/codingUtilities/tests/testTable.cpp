/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/TableData.hpp"
#include "codingUtilities/TableFormatter.hpp"
#include "codingUtilities/TableLayout.hpp"
#include "dataRepository/Group.hpp"
// TPL includes
#include <gtest/gtest.h>

using namespace geos;


TEST( testTable, tableClass )
{

  {
    TableLayout tableLayout( {"Well\nelement no.\nPV weighted\nbar",
                              "CordX",
                              "CoordZ",
                              "Prev\nelement",
                              "Next\nelement"}
                             );
    tableLayout.setTitle( "InternalWellGenerator well_injector1" );

    TableData tableData;
    tableData.addRow( "value1", "[30.21543]", "3.0", 54, 0 );
    tableData.addRow( "", "", "", "", "" );
    tableData.addRow( "Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante", "[30.21543]", "30.45465142",
                      787442, 10 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString(
                 tableData ),
               "\n+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+\n"
               "|                                                                                InternalWellGenerator well_injector1                                                                                 |\n"
               "+----------------------------------------------------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n"
               "|                                                                     Well                                                                     |    CordX     |    CoordZ     |   Prev    |   Next    |\n"
               "|                                                                 element no.                                                                  |              |               |  element  |  element  |\n"
               "|                                                                 PV weighted                                                                  |              |               |           |           |\n"
               "|                                                                     bar                                                                      |              |               |           |           |\n"
               "+----------------------------------------------------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n"
               "|                                                                    value1                                                                    |  [30.21543]  |      3.0      |    54     |     0     |\n"
               "|                                                                                                                                              |              |               |           |           |\n"
               "|  Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |  [30.21543]  |  30.45465142  |  787442   |    10     |\n"
               "+----------------------------------------------------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n\n"
               );
  };

  {
    TableLayout tableLayout( {"Duis fringilla, ligula sed porta fringilla,\nligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante",
                              "CordX",
                              "CoordZ",
                              "Prev\nelement",
                              "Next\nelement"}
                             );
    tableLayout.setTitle( "InternalWellGenerator well_injector1" );

    TableData tableData;
    tableData.addRow( "value1", "[30.21543]", "3.0", 54, 0 );
    tableData.addRow( "", "", "", "", "" );
    tableData.addRow( "value23", "[30.21543]", "30.45465142", 787442, 10 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "\n+---------------------------------------------------------------------------------------------------------------------------------------------------------+\n"
               "|                                                          InternalWellGenerator well_injector1                                                           |\n"
               "+--------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n"
               "|                           Duis fringilla, ligula sed porta fringilla,                            |    CordX     |    CoordZ     |   Prev    |   Next    |\n"
               "|  ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |              |               |  element  |  element  |\n"
               "+--------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n"
               "|                                              value1                                              |  [30.21543]  |      3.0      |    54     |     0     |\n"
               "|                                                                                                  |              |               |           |           |\n"
               "|                                             value23                                              |  [30.21543]  |  30.45465142  |  787442   |    10     |\n"
               "+--------------------------------------------------------------------------------------------------+--------------+---------------+-----------+-----------+\n\n"
               );
  }

  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::right}
    }
                             );
    tableLayout.setTitle( "InternalWellGenerator well_injector1" );

    TableData tableData;
    tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
    tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "\n+-----------------------------------------------------------------------------------------+\n"
               "|                          InternalWellGenerator well_injector1                           |\n"
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
               "|  Cras egestas  |  CoordX  |  C                    |  CoordZ     |     Prev  |     Next  |\n"
               "|                |          |                       |             |  element  |  element  |\n"
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
               "|     value1     |          |  3.0                  |  3.0129877  |        2  |        1  |\n"
               "|      val1      |  v       |  [3.045,42.02,89.25]  |  3          |       10  |        3  |\n"
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n\n"
               );
  }

  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left, false},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::middle, false},
    }
                             );
    tableLayout.setTitle( "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

    TableData tableData;
    tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
    tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "\n+------------------------------------------------------------------------------------------------------------------+\n"
               "|    Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis     |\n"
               "+---------------------------+---------------------+----------------------------------+-----------------------------+\n"
               "|       Cras egestas        |             CoordX  |                C                 |  CoordZ                     |\n"
               "+---------------------------+---------------------+----------------------------------+-----------------------------+\n"
               "|          value1           |                     |               3.0                |  3.0129877                  |\n"
               "|           val1            |                  v  |       [3.045,42.02,89.25]        |  3                          |\n"
               "+---------------------------+---------------------+----------------------------------+-----------------------------+\n\n"
               );
  }

  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::middle},
    }
                             );
    tableLayout.setTitle( "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

    TableData tableData;
    tableData.addRow( "value1" );
    tableData.addRow( "val1" );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "\n+-------------------------------------------------------------------------------------------------------------------+\n"
               "|     Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis      |\n"
               "+-------------------------------------------------------------------------------------------------------------------+\n"
               "|                                                  Cras egestas                                                     |\n"
               "+-------------------------------------------------------------------------------------------------------------------+\n"
               "|                                                     value1                                                        |\n"
               "|                                                      val1                                                         |\n"
               "+-------------------------------------------------------------------------------------------------------------------+\n\n"
               );
  }

  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::middle},
    }
                             );

    TableData tableData;
    tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
    tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
               "|  Cras egestas  |  CoordX  |           C           |  CoordZ     |  Prev     |   Next    |\n"
               "|                |          |                       |             |  element  |  element  |\n"
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
               "|     value1     |          |          3.0          |  3.0129877  |  2        |     1     |\n"
               "|      val1      |       v  |  [3.045,42.02,89.25]  |  3          |  10       |     3     |\n"
               "+----------------+----------+-----------------------+-------------+-----------+-----------+\n\n"
               );
  };

  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::middle},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::middle},
    } );
    tableLayout.setTitle( "InternalWellGenerator well_injector1" );
    tableLayout.setMargin( TableLayout::MarginValue::tiny );

    TableData tableData;
    tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
    tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.ToString( tableData ),
               "\n+-----------------------------------------------------------------+\n"
               "|              InternalWellGenerator well_injector1               |\n"
               "+------------+------+-------------------+---------+-------+-------+\n"
               "|Cras egestas|CoordX|         C         |CoordZ   |Prev   | Next  |\n"
               "|            |      |                   |         |element|element|\n"
               "+------------+------+-------------------+---------+-------+-------+\n"
               "|   value1   |      |        3.0        |3.0129877|2      |   1   |\n"
               "|    val1    |     v|[3.045,42.02,89.25]|3        |10     |   3   |\n"
               "+------------+------+-------------------+---------+-------+-------+\n\n"
               );
  };
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
