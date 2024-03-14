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

  //table with empty row
  {
    TableLayout tableLayout( {"Well\nelement no.\nPV weighted\nbar",
                              "CordX",
                              "CoordZ",
                              "Prev\nelement",
                              "Next\nelement"},
                             "InternalWellGenerator well_injector1"
                             );

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

  //same but with different values
  {
    TableLayout tableLayout( {"Duis fringilla, ligula sed porta fringilla,\nligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante",
                              "CordX",
                              "CoordZ",
                              "Prev\nelement",
                              "Next\nelement"}, "InternalWellGenerator well_injector1" );

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

  //table with TableLayout::ColumnParam
  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::right}
    }, "InternalWellGenerator well_injector1" );

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

  //test with hidden column
  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left, false},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::center, false},
    }, "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

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
  
  //test with 1 column
  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
    }, "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

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

  //test without title
  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::center},
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

  //test with tiny margin
  {
    TableLayout tableLayout( {
      TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
      TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::center},
      TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left},
      TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::center},
    }, "InternalWellGenerator well_injector1" );

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

  //test 2D table
  {
    TableLayout tableLayout( {"FakePressure", "Value1", "Value2"} );
    TableData2D tableData;

    for( real64 p = 10000; p<20000; p+=5000 )
    {
      for( real64 t = 400; t>=270; t+=-50.0 )
      {
        real64 value = t/p;
        tableData.addCell( t, p, value );
      }
    }

    TableTextFormatter tableLog( tableLayout );
    EXPECT_EQ( tableLog.ToString( tableData.buildTableData()),
               "+----------------+----------+------------------------+\n"
               "|  FakePressure  |  Value1  |         Value2         |\n"
               "+----------------+----------+------------------------+\n"
               "|      300       |   0.03   |          0.02          |\n"
               "|      350       |  0.035   |  0.023333333333333334  |\n"
               "|      400       |   0.04   |  0.02666666666666667   |\n"
               "+----------------+----------+------------------------+\n\n"
               );

  }
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
