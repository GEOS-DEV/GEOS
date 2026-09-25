/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_BOUNDARY_MORTAR_COMPONENTS_HPP
#define GEOS_MPM_BOUNDARY_MORTAR_COMPONENTS_HPP

#include "MPMBoundaryMortar.hpp"
#include <climits>
#include <sstream>

namespace geos
{
namespace mpm
{
namespace boundaryMortar
{

/** Coarse-grained scheduling only: each component is still a complete solve.
 * A future algebraic-interface-block scheduler can reuse this separation of
 * local computation, collective error agreement, and result distribution.
 */
inline int componentWorker( int component, int dimensions, int ranks, int requestedWorkers )
{
  require( dimensions >= 1 && dimensions <= 3 && ranks >= 1 &&
           requestedWorkers >= 1 && requestedWorkers <= 3 &&
           component >= 0 && component < dimensions, "Invalid mortar component schedule." );
  return component % std::min( dimensions, std::min( ranks, requestedWorkers ) );
}

/** Run independent x/y/z projections on up to three ranks.
 *
 * Communicator supplies rank(), size(), allMinimum(int), and typed
 * broadcast(T*,int,root). solve(component) MUST NOT call collectives or mutate
 * the shared trial field. All local solves finish before the first broadcast:
 * moving a root inside a blocking solve/broadcast loop would remain serial.
 *
 * All solve failures are distributed before any results are returned. Nonworker
 * ranks, including ranks without local particles or interface, participate in
 * the same collective sequence. Dense factors stay only on their worker and
 * are destroyed by solve(); all ranks receive the small Result arrays.
 */
template< class Communicator, class Solve >
std::array< Result, 3 > solveComponents( Communicator const & comm, int dimensions,
                                       int requestedWorkers, Index dofs, Index rows,
                                       Solve solve )
{
  require( dofs <= Index( INT_MAX ) && rows <= Index( INT_MAX ),
           "Boundary mortar result exceeds the MPI count limit." );
  (void)componentWorker( 0, dimensions, comm.size(), requestedWorkers );
  std::array< Result, 3 > results;
  std::array< int, 3 > success{{ 1, 1, 1 }};
  std::array< std::string, 3 > errors;
  int const rank = comm.rank();
  for( int c = 0; c < dimensions; ++c )
  {
    int const owner = componentWorker( c, dimensions, comm.size(), requestedWorkers );
    if( rank != owner ) continue;
    try
    {
      results[c] = solve( c );
      Result const & r = results[c];
      require( r.velocity.size() == dofs && r.mortarImpulse.size() == dofs &&
               r.reactionImpulse.size() == dofs && r.fallbackDof.size() == dofs &&
               r.preJump.size() == rows && r.postJump.size() == rows &&
               r.retainedFraction.size() == rows && r.accepted.size() == rows,
               "Boundary mortar worker returned inconsistent result dimensions." );
      require( r.rank <= Index( INT_MAX ) && r.unknowns <= Index( INT_MAX ) &&
               r.assembledRows <= Index( INT_MAX ) && r.rejectedRows <= Index( INT_MAX ),
               "Boundary mortar result metadata exceeds the MPI count limit." );
    }
    catch( std::exception const & e ) { success[c] = 0; errors[c] = e.what(); }
    catch( ... ) { success[c] = 0; errors[c] = "Unknown exception in component projection."; }
  }

  // Bounded messages avoid unbounded receives following a failed factorization.
  // Exchange every component's status before throwing, so no rank is left in a
  // later component's collective while another rank has already terminated.
  std::string firstError;
  for( int c = 0; c < dimensions; ++c )
  {
    int const owner = componentWorker( c, dimensions, comm.size(), requestedWorkers );
    comm.broadcast( &success[c], 1, owner );
    if( !success[c] )
    {
      std::array< char, 2048 > message{};
      if( rank == owner ) std::copy_n( errors[c].data(), std::min( errors[c].size(), message.size()-1 ), message.data() );
      comm.broadcast( message.data(), int( message.size() ), owner );
      if( firstError.empty() )
      {
        std::ostringstream text;
        text << "Boundary mortar component " << c << " on worker " << owner << ": " << message.data();
        firstError = text.str();
      }
    }
  }
  if( !firstError.empty() ) throw std::runtime_error( firstError );

  // Allocate receive storage before data broadcasts. A recoverable allocation
  // exception must be agreed collectively, not strand another rank in bcast.
  int ready = 1;
  try
  {
    for( int c = 0; c < dimensions; ++c )
    {
      if( rank == componentWorker( c, dimensions, comm.size(), requestedWorkers ) ) continue;
      Result & r = results[c];
      r.velocity.resize( dofs ); r.mortarImpulse.resize( dofs ); r.reactionImpulse.resize( dofs );
      r.fallbackDof.resize( dofs ); r.preJump.resize( rows ); r.postJump.resize( rows );
      r.retainedFraction.resize( rows ); r.accepted.resize( rows );
    }
  }
  catch( ... ) { ready = 0; }
  require( comm.allMinimum( ready ) != 0, "A rank could not allocate boundary mortar result buffers." );

  for( int c = 0; c < dimensions; ++c )
  {
    int const owner = componentWorker( c, dimensions, comm.size(), requestedWorkers );
    Result & r = results[c];
    std::array< int, 4 > counts{};
    std::array< double, 7 > scalars{};
    if( rank == owner )
    {
      counts = {{ int(r.rank), int(r.unknowns), int(r.assembledRows), int(r.rejectedRows) }};
      scalars = {{ r.maximumJump, r.relativeResidual, r.boundaryResidual, r.smallestRelativePivot,
                   r.setupSeconds, r.qrSeconds, r.recoverySeconds }};
    }
    comm.broadcast( counts.data(), int(counts.size()), owner );
    comm.broadcast( scalars.data(), int(scalars.size()), owner );
    r.rank = counts[0]; r.unknowns = counts[1]; r.assembledRows = counts[2]; r.rejectedRows = counts[3];
    r.maximumJump = scalars[0]; r.relativeResidual = scalars[1]; r.boundaryResidual = scalars[2];
    r.smallestRelativePivot = scalars[3]; r.setupSeconds = scalars[4];
    r.qrSeconds = scalars[5]; r.recoverySeconds = scalars[6];
    if( dofs > 0 )
    {
      comm.broadcast( r.velocity.data(), int(dofs), owner );
      comm.broadcast( r.mortarImpulse.data(), int(dofs), owner );
      comm.broadcast( r.reactionImpulse.data(), int(dofs), owner );
      comm.broadcast( r.fallbackDof.data(), int(dofs), owner );
    }
    if( rows > 0 )
    {
      comm.broadcast( r.preJump.data(), int(rows), owner );
      comm.broadcast( r.postJump.data(), int(rows), owner );
      comm.broadcast( r.retainedFraction.data(), int(rows), owner );
      comm.broadcast( r.accepted.data(), int(rows), owner );
    }
  }
  return results;
}

} // namespace boundaryMortar
} // namespace mpm
} // namespace geos
#endif
