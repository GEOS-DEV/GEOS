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
 * @file InterObjectRelationship.hpp
 */

#ifndef GEOSX_MESH_INTEROBJECTRELATION_HPP_
#define GEOSX_MESH_INTEROBJECTRELATION_HPP_


//#include "Common/typedefs.h"
#include <map>

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{
template< typename BASETYPE >
class InterObjectRelation : public BASETYPE
{
public:

  using base_type = BASETYPE;

  template< typename ... ARGS >
  InterObjectRelation( ARGS && ... args ):
    BASETYPE( std::forward< ARGS >( args )... )
  {}

  /// equals operator that sets *this to a single value of any type
  template< typename rTYPE > InterObjectRelation & operator=( const rTYPE & rhs )
  {
    BASETYPE::operator=( rhs );
    return (*this);
  }

  const base_type & Base() const { return static_cast< const BASETYPE & >(*this); }
  base_type & Base() { return dynamic_cast< BASETYPE & >(*this); }

  void SetRelatedObject( ObjectManagerBase const * const relatedObject )
  { m_relatedObject = relatedObject; }

  const ObjectManagerBase * RelatedObject() const
  { return m_relatedObject; }

  arrayView1d< globalIndex const > const & RelatedObjectLocalToGlobal() const
  { return this->m_relatedObject->localToGlobalMap(); }

  unordered_map< globalIndex, localIndex > const & RelatedObjectGlobalToLocal() const
  { return this->m_relatedObject->globalToLocalMap(); }

private:
  ObjectManagerBase const * m_relatedObject = nullptr;
};

typedef InterObjectRelation< array2d< localIndex > >                FixedOneToManyRelation;
}

#endif /* GEOSX_MESH_INTEROBJECTRELATION_HPP_ */
