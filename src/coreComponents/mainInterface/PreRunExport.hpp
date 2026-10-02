/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MAININTERFACE_PRERUNEXPORT_HPP_
#define GEOS_MAININTERFACE_PRERUNEXPORT_HPP_

#include <memory>

namespace geos
{
struct CommandLineOptions;

/** Handle the standalone JSON capability query before initializing MPI or runtime libraries. */
bool handleCapabilitiesCommand( int argc, char * argv[], int & exitCode );

/** Emit a JSON diagnostic for a failing opt-in export command; ordinary runs are unchanged. */
void reportPreRunExportError( int argc, char * argv[], char const * message );

/** Export global input metadata without entering the event/time loop. */
void runPreRunExport( std::unique_ptr< CommandLineOptions > options );
}
#endif
