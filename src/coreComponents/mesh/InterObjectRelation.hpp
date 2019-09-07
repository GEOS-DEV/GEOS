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

using base_type = BASETYPE;

/// equals operator that sets *this to a single value of any type
template<typename rTYPE> InterObjectRelation& operator=( const rTYPE& rhs )
{
BASETYPE::operator=(rhs);
return (*this);
}

const base_type & Base() const { return static_cast<const BASETYPE&>(*this); }
base_type & Base() { return dynamic_cast<BASETYPE&>(*this); }

void SetRelatedObject( ObjectManagerBase const * const relatedObject )
{ m_relatedObject = relatedObject; }

const ObjectDataStructureBaseT* RelatedObject() const
{ return m_relatedObject; }

globalIndex_array const & RelatedObjectLocalToGlobal() const
{ return this->m_relatedObject->m_localToGlobalMap; }

const unordered_map<globalIndex,localIndex>& RelatedObjectGlobalToLocal() const
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
