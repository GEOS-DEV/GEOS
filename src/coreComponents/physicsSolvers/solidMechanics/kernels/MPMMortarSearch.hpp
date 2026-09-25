/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_MORTAR_SEARCH_HPP
#define GEOS_MPM_MORTAR_SEARCH_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>

namespace geos
{
namespace mpm
{
namespace mortarSearch
{

/** Gather bins into the ORIGINAL sorted/unique marker order, with cheap fast paths.
 *
 * The marker builder appends each index to one bin in increasing order. First
 * deduplicate the at-most-27 bin keys: periodic wrapping can revisit a bin.
 * One bin (or already ordered disjoint index ranges) needs no marker sort at all.
 * Keep std::sort for interleaved bins: benchmarks showed a heap-based k-way merge
 * was slower on that case, despite its better asymptotic comparison count.
 *
 * The complete candidate sequence is unchanged, not just the candidate set.
 * That preserves the later partial-sort tie behavior and floating-point fit
 * accumulation order. No search radius, sample count, or surface fit changes.
 */
template< class Key, class Index, class BinMap >
void gatherSortedMarkers( std::array< Key, 27 > keys, std::size_t count,
                          BinMap const & bins, std::vector< Index > & markers )
{
  if( count > keys.size() ) throw std::runtime_error( "Too many mortar neighbor bins." );
  std::sort( keys.begin(), keys.begin()+count );
  auto const equivalent = []( Key const & a, Key const & b ) { return !(a < b) && !(b < a); };
  count = std::size_t( std::unique( keys.begin(), keys.begin()+count, equivalent )-keys.begin() );
  std::array< std::vector< Index > const *, 27 > streams{};
  std::size_t numberOfStreams = 0, total = 0;
  for( std::size_t k = 0; k < count; ++k )
  {
    auto const bin = bins.find( keys[k] );
    if( bin == bins.end() || bin->second.empty() ) continue;
    streams[numberOfStreams++] = &bin->second;
    total += bin->second.size();
  }
  markers.clear();
  markers.reserve( total );
  bool sorted = true;
  for( std::size_t k = 0; k < numberOfStreams; ++k )
  {
    auto const & values = *streams[k];
    if( !markers.empty() && values.front() < markers.back() ) sorted = false;
    markers.insert( markers.end(), values.begin(), values.end() );
  }
  if( !sorted ) std::sort( markers.begin(), markers.end() );
  // Production bins have unique indices. Keep this defensive pass so the helper
  // also matches the legacy contract for repeated entries in synthetic fixtures.
  markers.erase( std::unique( markers.begin(), markers.end() ), markers.end() );
}

} // namespace mortarSearch
} // namespace mpm
} // namespace geos
#endif
