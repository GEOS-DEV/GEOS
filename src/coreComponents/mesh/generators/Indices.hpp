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
#include "common/DataTypes.hpp"

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

using NodeLocIdx = fluent::NamedType< localIndex, struct NodeLocIdxTag, fluent::Comparable, fluent::Printable, fluent::PreIncrementable >;
using NodeGlbIdx = fluent::NamedType< globalIndex, struct NodeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable, fluent::PreIncrementable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable, fluent::PreIncrementable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using CellLocIdx = fluent::NamedType< localIndex, struct CellLocIdxTag, fluent::Comparable, fluent::Printable, fluent::PreIncrementable >;
using CellGlbIdx = fluent::NamedType< globalIndex, struct CellGlbIdxTag, fluent::Comparable, fluent::Printable >;

using MpiRank = fluent::NamedType< int, struct MpiRankTag, fluent::Comparable, fluent::Printable, fluent::Addable >;

inline NodeLocIdx operator "" _nli( unsigned long long int i )
{
  return NodeLocIdx{ NodeLocIdx::UnderlyingType( i ) };
}

inline EdgeLocIdx operator "" _eli( unsigned long long int i )
{
  return EdgeLocIdx{ EdgeLocIdx::UnderlyingType( i ) };
}

inline FaceLocIdx operator "" _fli( unsigned long long int i )
{
  return FaceLocIdx{ FaceLocIdx::UnderlyingType( i ) };
}

inline CellLocIdx operator "" _cli( unsigned long long int i )
{
  return CellLocIdx{ CellLocIdx::UnderlyingType( i ) };
}

inline NodeGlbIdx operator "" _ngi( unsigned long long int i )
{
  return NodeGlbIdx{ NodeGlbIdx::UnderlyingType( i ) };
}

inline EdgeGlbIdx operator "" _egi( unsigned long long int i )
{
  return EdgeGlbIdx{ EdgeGlbIdx::UnderlyingType( i ) };
}

inline FaceGlbIdx operator "" _fgi( unsigned long long int i )
{
  return FaceGlbIdx{ FaceGlbIdx::UnderlyingType( i ) };
}

inline CellGlbIdx operator "" _cgi( unsigned long long int i )
{
  return CellGlbIdx{ CellGlbIdx::UnderlyingType( i ) };
}

inline MpiRank operator "" _mpi( unsigned long long int i )
{
  return MpiRank{ MpiRank::UnderlyingType( i ) };
}

inline void to_json( json & j,
                     const MpiRank & v )
{
  j = v.get();
}

inline void from_json( const json & j,
                       MpiRank & v )
{
  v = MpiRank{ j.get< MpiRank::UnderlyingType >() };  // TODO use a `traits` instead
}

inline void from_json( const json & j,
                       NodeGlbIdx & v )
{
  v = NodeGlbIdx{ j.get< NodeGlbIdx::UnderlyingType >() };
}

inline void to_json( json & j,
                     const NodeGlbIdx & v )
{
  j = v.get();
}

inline void to_json( json & j,
                     const NodeLocIdx & v )
{
  j = v.get();
}

inline void to_json( json & j,
                     const EdgeGlbIdx & v )
{
  j = v.get();
}
inline void to_json( json & j,
                     const EdgeLocIdx & v )
{
  j = v.get();
}

inline void from_json( const json & j,
                       EdgeGlbIdx & v )
{
  v = EdgeGlbIdx{ j.get< EdgeGlbIdx::UnderlyingType >() };
}

inline void to_json( json & j,
                     const FaceGlbIdx & v )
{
  j = v.get();
}

inline void to_json( json & j,
                     const FaceLocIdx & v )
{
  j = v.get();
}

inline void from_json( const json & j,
                       FaceGlbIdx & v )
{
  v = FaceGlbIdx{ j.get< FaceGlbIdx::UnderlyingType >() };
}

inline void to_json( json & j,
                     const CellGlbIdx & v )
{
  j = v.get();
}

inline void to_json( json & j,
                     const CellLocIdx & v )
{
  j = v.get();
}

}

#endif //GEOS_INDICES_HPP
