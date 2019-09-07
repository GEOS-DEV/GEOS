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

#include "common/initialization.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include "managers/ProblemManager.hpp"

using namespace geosx;


int main( int argc, char *argv[] )
{
basicSetup( argc, argv );

printTypeSummary();

{
timeval tim;
gettimeofday(&tim, nullptr);
real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);

std::string restartFileName;
bool restart = ProblemManager::ParseRestart( argc, argv, restartFileName );
if (restart) {
GEOS_LOG_RANK_0("Loading restart file " << restartFileName);
dataRepository::SidreWrapper::reconstructTree( restartFileName, "sidre_hdf5", MPI_COMM_GEOSX );
}

ProblemManager problemManager( "Problem", nullptr );

problemManager.InitializePythonInterpreter();
problemManager.ParseCommandLineInput( argc, argv );

if ( !problemManager.getSchemaFileName().empty() )
{
problemManager.GenerateDocumentation();
}
else
{
problemManager.ParseInputFile();

problemManager.ProblemSetup();

if (restart) {
problemManager.ReadRestartOverwrite( restartFileName );
}

MPI_Barrier(MPI_COMM_GEOSX);
GEOS_LOG_RANK_0("Running simulation");

gettimeofday(&tim, nullptr);
const real64 t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

problemManager.RunSimulation();

gettimeofday(&tim, nullptr);
const real64 t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

GEOS_LOG_RANK_0("\ninit time = " << std::setprecision(5) << t_initialize-t_start <<
"s, run time = " << t_run-t_initialize << "s");
}


problemManager.ClosePythonInterpreter();
}

basicCleanup();


return 0;
}

