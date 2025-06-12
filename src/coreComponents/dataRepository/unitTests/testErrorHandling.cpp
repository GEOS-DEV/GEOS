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
#include "common/logger/ErrorHandling.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace dataRepository;

TEST( ErrorHandling, testYaml )
{
  g_errorLogger.setFilename( "errorsOutput.yaml" );
  g_errorLogger.setWriteValue( true );
  double minPrecision = 1e-6;
  double maxPrecision = 1e-3;
  int x = 5;

  DataFileContext const context = DataFileContext( "Base Test Class", __FILE__, __LINE__ );

  if( g_errorLogger.writeFile() )
  {
    g_errorLogger.createFile();
  }

  GEOS_WARNING( "Conflicting pressure boundary conditions" );
  GEOS_WARNING_IF( x == 5, "Pressure value is too small." );
  GEOS_WARNING_CTX_IF( x == 5,
                       GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                 context.toString(), minPrecision, maxPrecision, minPrecision ),
                       context );
  try
  {
    GEOS_THROW_CTX_IF( x == 5,
                       "Group " << context.toString() << " has no wrapper named" << std::endl,
                       std::domain_error,
                       context );
  }
  catch( std::domain_error const & ex )
  {
    string const errorMsg = "Table input error.\n";
    g_errorLogger.currentErrorMsg()
      .addToMsg( errorMsg )
      .addContextInfo( context.getContextInfo().setPriority( 2 ) );
  }

  if( g_errorLogger.writeFile() )
  {
    g_errorLogger.write( g_errorLogger.currentErrorMsg() );
  }
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();

  return result;
}
