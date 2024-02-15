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
#include "../Table.hpp"
#include "../../dataRepository/Group.hpp"
// TPL includes
#include <gtest/gtest.h>

using namespace geos;


TEST( testTable, tableClass )
{

  std::vector< string > tableTestsOutput;

  std::ostringstream oss;

  Table tableTest( {"Well\nelement no.\nPV weighted\nbar",
                    "CordX",
                    "CoordZ",
                    "Prev\nelement",
                    "Next\nelement"}
                   );
  tableTest.setTitle( "InternalWellGenerator well_injector1" );
  tableTest.addRow< 5 >( "value1", "[30.21543]", "3.0", 54, 0 );
  tableTest.addRow< 5 >( "", "", "", "", "" );
  tableTest.addRow< 5 >( "Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante", "[30.21543]", "30.45465142",
                         787442, 10 );
  tableTest.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[0],
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

  Table tableTest2( {"Duis fringilla, ligula sed porta fringilla,\nligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante",
                     "CordX",
                     "CoordZ",
                     "Prev\nelement",
                     "Next\nelement"}
                    );
  tableTest2.setTitle( "InternalWellGenerator well_injector1" );
  tableTest2.addRow< 5 >( "value1", "[30.21543]", "3.0", 54, 0 );
  tableTest2.addRow< 5 >( "", "", "", "", "" );
  tableTest2.addRow< 5 >( "value23", "[30.21543]", "30.45465142", 787442, 10 );
  tableTest2.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[1],
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

  Table tableTest3 (
  {
    "Cras egestas ipsum a nisl. Vivamus variu\ndolor utsisicdis parturient montes,\nnascetur ridiculus mus. Duis fringilla, ligula sed porta fringilla, ligula wisicommodo felis,\nut adi\npiscing felis dui in enim. Suspendisse malesuada ultrices ante",
    "CoordX", "C", "CoordZ",
  } );
  tableTest3.setTitle( "InternalWellGenerator well_injector1" );
  tableTest3.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ ( tableTestsOutput[2],
              "\n+-----------------------------------------------------------------------------------------------------------------------------+\n"
              "|                                            InternalWellGenerator well_injector1                                             |\n"
              "+-------------------------------------------------------------------------------------------------+----------+-----+----------+\n"
              "|                            Cras egestas ipsum a nisl. Vivamus variu                             |  CoordX  |  C  |  CoordZ  |\n"
              "|                               dolor utsisicdis parturient montes,                               |          |     |          |\n"
              "|  nascetur ridiculus mus. Duis fringilla, ligula sed porta fringilla, ligula wisicommodo felis,  |          |     |          |\n"
              "|                                             ut adi                                              |          |     |          |\n"
              "|                 piscing felis dui in enim. Suspendisse malesuada ultrices ante                  |          |     |          |\n"
              "+-------------------------------------------------------------------------------------------------+----------+-----+----------+\n\n"
              );

  Table tableTest4( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::left},
    Table::ColumnParam{{"C"}, Table::Alignment::left},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::right},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::right}
  } );
  tableTest4.setTitle( "InternalWellGenerator well_injector1" );
  tableTest4.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest4.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest4.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[3],
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

  Table tableTest5( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::left},
    Table::ColumnParam{{"C"}, Table::Alignment::left},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::right, false},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::right, false},

  } );
  tableTest5.setTitle( "InternalWellGenerator well_injector1" );
  tableTest5.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest5.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest5.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[4],
             "\n+-----------------------------------------------------------------+\n"
             "|              InternalWellGenerator well_injector1               |\n"
             "+----------------+----------+-----------------------+-------------+\n"
             "|  Cras egestas  |  CoordX  |  C                    |  CoordZ     |\n"
             "+----------------+----------+-----------------------+-------------+\n"
             "|     value1     |          |  3.0                  |  3.0129877  |\n"
             "|      val1      |  v       |  [3.045,42.02,89.25]  |  3          |\n"
             "+----------------+----------+-----------------------+-------------+\n\n"
             );

  Table tableTest6( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::left},
    Table::ColumnParam{{"C"}, Table::Alignment::left},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::middle},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::right},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::middle},
  } );
  tableTest6.setTitle( "InternalWellGenerator well_injector1" );
  tableTest6.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest6.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest6.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[5],
             "\n+-----------------------------------------------------------------------------------------+\n"
             "|                          InternalWellGenerator well_injector1                           |\n"
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
             "|  Cras egestas  |  CoordX  |  C                    |   CoordZ    |     Prev  |   Next    |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
             "|     value1     |          |  3.0                  |  3.0129877  |        2  |     1     |\n"
             "|      val1      |  v       |  [3.045,42.02,89.25]  |      3      |       10  |     3     |\n"
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n\n"
             );

  Table tableTest7( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::right},
    Table::ColumnParam{{"C"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::left, false},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::middle, false},
  } );
  tableTest7.setTitle( "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );
  tableTest7.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest7.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest7.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[6],
             "\n+------------------------------------------------------------------------------------------------------------------+\n"
             "|    Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis     |\n"
             "+---------------------------+---------------------+----------------------------------+-----------------------------+\n"
             "|       Cras egestas        |             CoordX  |                C                 |  CoordZ                     |\n"
             "+---------------------------+---------------------+----------------------------------+-----------------------------+\n"
             "|          value1           |                     |               3.0                |  3.0129877                  |\n"
             "|           val1            |                  v  |       [3.045,42.02,89.25]        |  3                          |\n"
             "+---------------------------+---------------------+----------------------------------+-----------------------------+\n\n"
             );

  Table tableTest8( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::right},
    Table::ColumnParam{{"C"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::left},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::middle},
  } );
  tableTest8.setTitle( "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );
  tableTest8.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest8.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest8.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[7],
             "\n+----------------------------------------------------------------------------------------------------------------+\n"
             "|    Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis    |\n"
             "+-------------------+-------------+--------------------------+----------------+--------------+-------------------+\n"
             "|   Cras egestas    |     CoordX  |            C             |  CoordZ        |  Prev        |       Next        |\n"
             "|                   |             |                          |                |  element     |      element      |\n"
             "+-------------------+-------------+--------------------------+----------------+--------------+-------------------+\n"
             "|      value1       |             |           3.0            |  3.0129877     |  2           |         1         |\n"
             "|       val1        |          v  |   [3.045,42.02,89.25]    |  3             |  10          |         3         |\n"
             "+-------------------+-------------+--------------------------+----------------+--------------+-------------------+\n\n"
             );

  Table tableTest9( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},

  } );
  tableTest9.setTitle( "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );
  tableTest9.addRow< 1 >( "value1" );
  tableTest9.addRow< 1 >( "val1" );
  tableTest9.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[8],
             "\n+-------------------------------------------------------------------------------------------------------------------+\n"
             "|     Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis      |\n"
             "+-------------------------------------------------------------------------------------------------------------------+\n"
             "|                                                  Cras egestas                                                     |\n"
             "+-------------------------------------------------------------------------------------------------------------------+\n"
             "|                                                     value1                                                        |\n"
             "|                                                      val1                                                         |\n"
             "+-------------------------------------------------------------------------------------------------------------------+\n\n"
             );

  Table tableTest10( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
  } );
  tableTest10.setTitle( "title1" );
  tableTest10.addRow< 1 >( "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );
  tableTest10.addRow< 1 >( "val1" );
  tableTest10.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[9],
             "\n+--------------------------------------------------------------------------------------------------------------+\n"
             "|                                                    title1                                                    |\n"
             "+--------------------------------------------------------------------------------------------------------------+\n"
             "|                                                Cras egestas                                                  |\n"
             "+--------------------------------------------------------------------------------------------------------------+\n"
             "|  Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis    |\n"
             "|                                                    val1                                                      |\n"
             "+--------------------------------------------------------------------------------------------------------------+\n\n"
             );

  Table tableTest11( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
  } );
  tableTest11.setTitle( "title1" );
  tableTest11.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[10],
             "\n+------------------+\n"
             "|      title1      |\n"
             "+------------------+\n"
             "|  Cras egestas    |\n"
             "+------------------+\n\n"
             );

  Table tableTest12( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::right},
    Table::ColumnParam{{"C"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::left},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::middle},
  } );
  tableTest12.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest12.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest12.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[11],
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
             "|  Cras egestas  |  CoordX  |           C           |  CoordZ     |  Prev     |   Next    |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n"
             "|     value1     |          |          3.0          |  3.0129877  |  2        |     1     |\n"
             "|      val1      |       v  |  [3.045,42.02,89.25]  |  3          |  10       |     3     |\n"
             "+----------------+----------+-----------------------+-------------+-----------+-----------+\n\n"
             );

  Table tableTest13( {
    Table::ColumnParam{{"Cras egestas"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordX"}, Table::Alignment::right},
    Table::ColumnParam{{"C"}, Table::Alignment::middle},
    Table::ColumnParam{{"CoordZ"}, Table::Alignment::left},
    Table::ColumnParam{{"Prev\nelement"}, Table::Alignment::left},
    Table::ColumnParam{{"Next\nelement"}, Table::Alignment::middle},
  } );
  tableTest13.setTitle( "InternalWellGenerator well_injector1" );
  tableTest13.setMargin( Table::MarginValue::tiny );
  tableTest13.addRow< 6 >( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableTest13.addRow< 6 >( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );
  tableTest13.draw( oss );
  tableTestsOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( tableTestsOutput[12],
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
}


int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
