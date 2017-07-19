/*
 * InterObjectRelationship.cpp
 *
 *  Created on: May 4, 2012
 *      Author: settgast1
 */

#include "InterObjectRelation.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"


template < typename BASETYPE >
const gArray1d& InterObjectRelation<BASETYPE>::RelatedObjectLocalToGlobal() const
{
  return this->m_relatedObject->m_localToGlobalMap;
}
template const gArray1d& InterObjectRelation<lArray1d>::RelatedObjectLocalToGlobal() const;
template const gArray1d& InterObjectRelation<lArray2d>::RelatedObjectLocalToGlobal() const;
template const gArray1d& InterObjectRelation<Array1dT<lArray1d> >::RelatedObjectLocalToGlobal() const;
template const gArray1d& InterObjectRelation<Array1dT<lSet> >::RelatedObjectLocalToGlobal() const;
template const gArray1d& InterObjectRelation<Array1dT<pArray1d> >::RelatedObjectLocalToGlobal() const;
template const gArray1d& InterObjectRelation<Array1dT<pSet> >::RelatedObjectLocalToGlobal() const;


template < typename BASETYPE >
const std::map<globalIndex,localIndex>& InterObjectRelation<BASETYPE>::RelatedObjectGlobalToLocal() const
{
  return this->m_relatedObject->m_globalToLocalMap;
}
template const std::map<globalIndex,localIndex>& InterObjectRelation<lArray1d>::RelatedObjectGlobalToLocal() const;
template const std::map<globalIndex,localIndex>& InterObjectRelation<lArray2d>::RelatedObjectGlobalToLocal() const;
template const std::map<globalIndex,localIndex>& InterObjectRelation<Array1dT<lArray1d> >::RelatedObjectGlobalToLocal() const;
template const std::map<globalIndex,localIndex>& InterObjectRelation<Array1dT<lSet> >::RelatedObjectGlobalToLocal() const;
