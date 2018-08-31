/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

  void setElementRegionManager( ElementRegionManager const * const input )
  {
    m_elemRegionManager = input;
  }

  ElementRegionManager const * getElementRegionManager() const
  {
    return m_elemRegionManager;
  }


//  template< bool DOPACK >
//  friend localIndex Pack( char *& buffer,
//                   localIndex_array const & packList,
//                   ElementRegionManager const * const elemManager );

//private:
  BASETYPE m_toElementRegion;
  BASETYPE m_toElementSubRegion;
  BASETYPE m_toElementIndex;

  ElementRegionManager const * m_elemRegionManager;
};

template< typename BASETYPE >
ToElementRelation<BASETYPE>::ToElementRelation():
  m_toElementRegion(),
  m_toElementSubRegion(),
  m_toElementIndex(),
  m_elemRegionManager(nullptr)
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
typedef ToElementRelation<array2d<localIndex>> FixedToManyElementRelation;
typedef ToElementRelation<array1d<localIndex_array> > OrderedVariableToManyElementRelation;
typedef ToElementRelation<array1d<set<localIndex>> > UnorderedVariableToManyElementRelation;

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESH_TOELEMENTRELATION_HPP_ */
