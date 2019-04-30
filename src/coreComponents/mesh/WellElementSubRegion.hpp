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

#include "../wells/WellElement.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto wellElementData = "wellElementData";
}
}

/**
 * @class WellElementSubRegion
 *
 * This class describes a collection of local well elements
 */  
class WellElementSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType=FixedOneToManyRelation;
  
  static const string CatalogName() { return dataRepository::keys::wellElementData; }

  virtual const string getCatalogName() const override { return WellElementSubRegion::CatalogName(); }

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */  
  WellElementSubRegion( string const & name, 
                        ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  virtual ~WellElementSubRegion() override;

  /**
   * @brief Getter for the number of well elements on this rank
   * @return the number of well elememts on this rank
   */
  localIndex numWellElementsLocal()  const
  { return integer_conversion<localIndex>(size()); }

  
  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                                   const NodeManager& nodeManager,
                                                   const bool useReferencePos = true) const override
  {
    return m_elementCenter[k]; 
  }

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override
  {
  }

  /*!
   * @brief returns the element to node relations.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const override
  { 
    return m_toNodesRelation[k];
  }

  /*!
   * @brief returns the element to node relations.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) override
  {
    return m_toNodesRelation[k];
  }

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nextWellElementIndexString = "nextWellElementIndex";
    static constexpr auto gravityDepthString = "gravityDepth";
    
    dataRepository::ViewKey nextWellElementIndex = { nextWellElementIndexString };
    dataRepository::ViewKey gravityDepth = { gravityDepthString };

  } viewKeysWellElementData;
  
protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:
  /// The elements to nodes relation is one to one relation.
  NodeMapType  m_toNodesRelation;

  // indices of the next well element 
  array1d<localIndex> m_nextWellElementIndex;

  // depth of the current well element
  array1d<real64> m_gravityDepth;
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MESH_WELLELEMENTSUBREGION_HPP_ */

