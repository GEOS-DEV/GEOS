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

#ifndef GEOS_INDICES_HPP
#define GEOS_INDICES_HPP

#include "common/GeosxMacros.hpp"

#include "LvArray/src/limits.hpp"

// TODO remove in the end.
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// TODO make conditional in the end.
#include <NamedType/named_type.hpp>

namespace geos
{

template< typename OUTPUT, typename INPUT >
inline GEOS_HOST_DEVICE
OUTPUT intConv( INPUT input )
{
  return LvArray::integerConversion< OUTPUT >( input );
}

using NodeLocIdx = fluent::NamedType< localIndex, struct NodeLocIdxTag, fluent::Comparable, fluent::Printable >;
using NodeGlbIdx = fluent::NamedType< globalIndex, struct NodeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using CellLocIdx = fluent::NamedType< localIndex, struct CellLocIdxTag, fluent::Comparable, fluent::Printable >;
using CellGlbIdx = fluent::NamedType< globalIndex, struct CellGlbIdxTag, fluent::Comparable, fluent::Printable >;

using MpiRank = fluent::NamedType< int, struct MpiRankTag, fluent::Comparable, fluent::Printable, fluent::Addable >;

EdgeGlbIdx operator "" _egi( unsigned long long int i )
{
  return EdgeGlbIdx{ EdgeGlbIdx::UnderlyingType( i ) };
}

FaceGlbIdx operator "" _fgi( unsigned long long int i )
{
  return FaceGlbIdx{ FaceGlbIdx::UnderlyingType( i ) };
}

MpiRank operator "" _mpi( unsigned long long int i )
{
  return MpiRank{ MpiRank::UnderlyingType( i ) };
}

void to_json( json & j,
              const MpiRank & v )
{
  j = v.get();
}

void from_json( const json & j,
                MpiRank & v )
{
  v = MpiRank{ j.get< MpiRank::UnderlyingType >() };  // TODO use a `traits` instead
}

void from_json( const json & j,
                NodeGlbIdx & v )
{
  v = NodeGlbIdx{ j.get< NodeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const NodeGlbIdx & v )
{
  j = v.get();
}

void to_json( json & j,
              const EdgeGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                EdgeGlbIdx & v )
{
  v = EdgeGlbIdx{ j.get< EdgeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const FaceGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                FaceGlbIdx & v )
{
  v = FaceGlbIdx{ j.get< FaceGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const CellGlbIdx & v )
{
  j = v.get();
}

}

#endif //GEOS_INDICES_HPP
