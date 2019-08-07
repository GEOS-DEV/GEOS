/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * InterObjectRelationship.h
 *
 *  Created on: May 4, 2012
 *      Author: settgast1
 */

#ifndef INTEROBJECTRELATION_H_
#define INTEROBJECTRELATION_H_


//#include "Common/typedefs.h"
#include <map>

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{
template < typename BASETYPE >
class InterObjectRelation : public BASETYPE
{
public:


  /// equals operator that sets *this to a single value of any type
  template<typename rTYPE> InterObjectRelation& operator=( const rTYPE& rhs )
  {
    BASETYPE::operator=(rhs);
    return (*this);
  }

  const BASETYPE& Base() const { return static_cast<const BASETYPE&>(*this); }
  BASETYPE& Base() { return dynamic_cast<BASETYPE&>(*this); }

  void SetRelatedObject( ObjectManagerBase const * const relatedObject )
  { m_relatedObject = relatedObject; }

  const ObjectDataStructureBaseT* RelatedObject() const
  { return m_relatedObject; }

  globalIndex_array const & RelatedObjectLocalToGlobal() const
  { return this->m_relatedObject->m_localToGlobalMap; }

  const std::map<globalIndex,localIndex>& RelatedObjectGlobalToLocal() const
  { return this->m_relatedObject->m_globalToLocalMap; }

private:
  ObjectDataStructureBaseT const * m_relatedObject = nullptr;
};

typedef InterObjectRelation<array1d<localIndex>>                OneToOneRelation;
typedef InterObjectRelation<array1d<localIndex const>>          OneToOneConstRelation;

typedef InterObjectRelation<array2d<localIndex>>                FixedOneToManyRelation;
typedef InterObjectRelation<array2d<localIndex const>>          FixedOneToManyConstRelation;

typedef InterObjectRelation<array1d<array1d<localIndex>>>       OrderedVariableOneToManyRelation;
typedef InterObjectRelation<array1d<array1d<localIndex const>>> OrderedVariableOneToManyConstRelation;

typedef InterObjectRelation<array1d<set<localIndex>>>           UnorderedVariableOneToManyRelation;
typedef InterObjectRelation<array1d<set<localIndex const>>>     UnorderedVariableOneToManyConstRelation;

}

#endif /* INTEROBJECTRELATION_H_ */
