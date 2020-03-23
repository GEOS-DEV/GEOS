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

/**
 * @file ToElementRelation.hpp
 */

#ifndef GEOSX_MESH_TOELEMENTRELATION_HPP_
#define GEOSX_MESH_TOELEMENTRELATION_HPP_

#include "InterObjectRelation.hpp"

namespace geosx
{

class ElementRegionManager;

template< typename BASETYPE >
class ToElementRelation
{
public:

  using base_type = BASETYPE;

  ToElementRelation();
  ~ToElementRelation();

  template< typename ... DIMS >
  void resize( DIMS... newdims );

  localIndex size() const
  {
    return m_toElementRegion.size();
  }

  localIndex size( int const dim ) const
  {
    return m_toElementRegion.size( dim );
  }

  void setElementRegionManager( ElementRegionManager const * const input )
  {
    m_elemRegionManager = input;
  }

  ElementRegionManager const * getElementRegionManager() const
  {
    return m_elemRegionManager;
  }

  BASETYPE m_toElementRegion;
  BASETYPE m_toElementSubRegion;
  BASETYPE m_toElementIndex;

  ElementRegionManager const * m_elemRegionManager;
};

template< typename BASETYPE >
ToElementRelation< BASETYPE >::ToElementRelation():
  m_toElementRegion(),
  m_toElementSubRegion(),
  m_toElementIndex(),
  m_elemRegionManager( nullptr )
{}

template< typename BASETYPE >
ToElementRelation< BASETYPE >::~ToElementRelation()
{}


template< typename BASETYPE >
template< typename ... DIMS >
void ToElementRelation< BASETYPE >::resize( DIMS... newdims )
{
  m_toElementRegion.resize( newdims ... );
  m_toElementSubRegion.resize( newdims ... );
  m_toElementIndex.resize( newdims ... );
}

typedef ToElementRelation< array2d< localIndex > > FixedToManyElementRelation;
typedef ToElementRelation< ArrayOfArrays< localIndex > > OrderedVariableToManyElementRelation;

void erase( OrderedVariableToManyElementRelation & relation,
            localIndex const firstIndex,
            localIndex const er,
            localIndex const esr,
            localIndex const ei );

void insert( OrderedVariableToManyElementRelation & relation,
             localIndex const firstIndex,
             localIndex const er,
             localIndex const esr,
             localIndex const ei );



} /* namespace geosx */

#endif /* GEOSX_MESH_TOELEMENTRELATION_HPP_ */
