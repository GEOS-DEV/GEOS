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
  InterObjectRelation();

  //InterObjectRelation( const ObjectDataStructureBaseT& relatedObject );

  ~InterObjectRelation() {}

  InterObjectRelation( const InterObjectRelation& copiedRelationship );

  // equals operator should not be called
  InterObjectRelation& operator=( const InterObjectRelation& rhs )
  {
    BASETYPE::operator=( static_cast<BASETYPE>(rhs));
    m_relatedObject = rhs.m_relatedObject;
    return *this;
  }

  InterObjectRelation& operator=( const BASETYPE& rhs )
  {
    BASETYPE::operator=( rhs );
    return *this;
  }


  /// equals operator that sets *this to a single value of any type
  template<typename rTYPE> InterObjectRelation& operator=( const rTYPE& rhs )
  {
    BASETYPE::operator=(rhs);
    return (*this);
  }


  const BASETYPE& Base() const { return static_cast<const BASETYPE&>(*this); }
  BASETYPE& Base() { return dynamic_cast<BASETYPE&>(*this); }

  void SetRelatedObject( const ObjectDataStructureBaseT* const relatedObject )
  { m_relatedObject = relatedObject; }

  const ObjectDataStructureBaseT* RelatedObject() const
  { return m_relatedObject; }

  globalIndex_const_array& RelatedObjectLocalToGlobal() const
  {
//    return this->m_relatedObject->m_localToGlobalMap;
  }

  const std::map<globalIndex,localIndex>& RelatedObjectGlobalToLocal() const
  {
    return this->m_relatedObject->m_globalToLocalMap;
  }

private:
  const ObjectDataStructureBaseT* m_relatedObject;


};


template < typename BASETYPE >
InterObjectRelation<BASETYPE>::InterObjectRelation( ):
m_relatedObject(NULL)
{}


template < typename BASETYPE >
InterObjectRelation<BASETYPE>::InterObjectRelation( const InterObjectRelation& copiedRelationship ):
BASETYPE( static_cast<BASETYPE const&>(copiedRelationship)),
m_relatedObject(copiedRelationship.m_relatedObject)
{

}



typedef InterObjectRelation<localIndex_array> OneToOneRelation;
typedef InterObjectRelation<lArray2d> FixedOneToManyRelation;
typedef InterObjectRelation<array<localIndex_array> > OrderedVariableOneToManyRelation;
typedef InterObjectRelation<array<lSet> > UnorderedVariableOneToManyRelation;

typedef InterObjectRelation<array<array<localIndex_array> > > OrderedVariableOneToManyToManyRelation;

typedef InterObjectRelation<array< pArray1d > > OrderedVariableOneToManyPairRelation;
typedef InterObjectRelation<array< pSet > > UnorderedVariableOneToManyPairRelation;
}


#endif /* INTEROBJECTRELATION_H_ */
