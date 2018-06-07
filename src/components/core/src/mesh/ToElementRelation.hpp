/*
 * ToElementRelation.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESH_TOELEMENTRELATION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESH_TOELEMENTRELATION_HPP_

#include "InterObjectRelation.hpp"

namespace geosx
{

class ElementRegionManager;

template< typename BASETYPE >
class ToElementRelation
{
public:
  ToElementRelation();
  ~ToElementRelation();

  template< typename... DIMS >
  void resize( DIMS... newdims );

  localIndex size() const
  {
    return m_toElementRegion.size();
  }

//  template< bool DOPACK >
//  friend localIndex Pack( char *& buffer,
//                   localIndex_array const & packList,
//                   ElementRegionManager const * const elemManager );

//private:
  BASETYPE m_toElementRegion;
  BASETYPE m_toElementSubRegion;
  BASETYPE m_toElementIndex;

};

template< typename BASETYPE >
ToElementRelation<BASETYPE>::ToElementRelation():
  m_toElementRegion(),
  m_toElementSubRegion(),
  m_toElementIndex()
{

}

template< typename BASETYPE >
ToElementRelation<BASETYPE>::~ToElementRelation()
{
}


template< typename BASETYPE >
template< typename... DIMS >
void ToElementRelation<BASETYPE>::resize( DIMS... newdims )
{
  m_toElementRegion.resize(newdims...);
  m_toElementSubRegion.resize(newdims...);
  m_toElementIndex.resize(newdims...);
}



//typedef ToElementRelation<localIndex_array> OneToOneRelation;
typedef ToElementRelation<lArray2d> FixedToManyElementRelation;
typedef ToElementRelation<array<localIndex_array> > OrderedVariableToManyElementRelation;
typedef ToElementRelation<array<lSet> > UnorderedVariableToManyElementRelation;


class ElementRegionManager;
template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 UnorderedVariableToManyElementRelation const & var,
                 ElementRegionManager const * const elementRegionManager );


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESH_TOELEMENTRELATION_HPP_ */
