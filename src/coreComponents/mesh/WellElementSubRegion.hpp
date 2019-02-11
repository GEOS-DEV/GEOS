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

#ifndef WELLELEMENTSUBREGION_HPP_
#define WELLELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto wellElementSubRegion = "wellElementSubRegion";
}
}

class WellElementSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType=OneToOneRelation;
  
  static const string CatalogName()
  { return "WellCell"; }

  virtual const string getCatalogName() const override
  {
    return WellElementSubRegion::CatalogName();
  }

  WellElementSubRegion( string const & name, ManagedGroup * const parent );
  virtual ~WellElementSubRegion() override;
  
  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                                   const NodeManager& nodeManager,
                                                   const bool useReferencePos = true) const override
  {
    return m_elementCenter[k]; // TODO
  }

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
    // TODO
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override
  {
    // TODO
  }

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
  };

  /*!
   * @brief returns the element to node relations.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const override
  { 
    return arraySlice1dRval<localIndex const>(0,0,0); // TODO
  }

  /*!
   * @brief returns the element to node relations.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) override
  {
    return arraySlice1dRval<localIndex>(0,0,0); // TODO
  }

private:
  /// The elements to nodes relation is one to one relation.
  NodeMapType  m_toNodesRelation;
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MESH_WELLELEMENTSUBREGION_HPP_ */

