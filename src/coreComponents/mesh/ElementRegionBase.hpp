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

#ifndef ELEMENT_REGION_BASE_H
#define ELEMENT_REGION_BASE_H

#include "CellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "FaceElementSubRegion.hpp"

namespace geosx
{

class StableTimeStep;

class FaceManager;

/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

//  static const string CatalogName()
//  { return "ElementRegionBase"; }
//
//  virtual const string getCatalogName() const override
//  { return ElementRegionBase::CatalogName(); }


  ///@}


  ElementRegionBase() = delete;

  ElementRegionBase( string const & name, ManagedGroup * const parent );


  ElementRegionBase(const ElementRegionBase& init);

  virtual ~ElementRegionBase() override;

  virtual void GenerateMesh( ManagedGroup const * const cellBlocks )
  {
    GEOS_ERROR( "ElementRegionBase::GenerateMesh() should be overriden if called.");
  }

//  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const NodeManager )
//  {
//    GEOS_ERROR( "ElementRegionBase::GenerateAggregates() should be overriden if called.");
//  }

  subGroupMap & GetSubRegions()
  {
    return GetGroup(viewKeyStruct::elementSubRegions)->GetSubGroups();
  }

  subGroupMap const & GetSubRegions() const
  {
    return GetGroup(viewKeyStruct::elementSubRegions)->GetSubGroups();
  }

  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(viewKeyStruct::elementSubRegions)->GetGroup<SUBREGIONTYPE>(regionName);
  }
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(viewKeyStruct::elementSubRegions)->GetGroup<SUBREGIONTYPE>(regionName);
  }

  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup(viewKeyStruct::elementSubRegions)->GetGroup<SUBREGIONTYPE>(index);
  }
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup(viewKeyStruct::elementSubRegions)->GetGroup<SUBREGIONTYPE>(index);
  }

  localIndex numSubRegions() const
  {
    return this->GetGroup(viewKeyStruct::elementSubRegions)->GetSubGroups().size();
  }


  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions<CellElementSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions<CellElementSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    ManagedGroup const * const elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
    elementSubRegions->forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES...>( std::forward<LAMBDA>(lambda) );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    ManagedGroup * const elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
    elementSubRegions->forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES...>( std::forward<LAMBDA>(lambda) );
  }


  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    forElementSubRegionsIndex<CellElementSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    forElementSubRegionsIndex<CellElementSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase const * const subRegion = this->GetSubRegion(esr);
      applyLambdaToContainer<ElementSubRegionBase, SUBREGIONTYPE, SUBREGIONTYPES...>( subRegion, [&]( auto const * const castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      });
    }
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase * const subRegion = this->GetSubRegion(esr);
      applyLambdaToContainer<ElementSubRegionBase, SUBREGIONTYPE,SUBREGIONTYPES...>( subRegion, [&]( auto * const castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      });
    }
  }


  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto materialListString = "materialList";
    static constexpr auto elementSubRegions = "elementSubRegions";
  };

  string_array & getMaterialList() {return m_materialList;}
  string_array const & getMaterialList() const {return m_materialList;}

protected:

private:

  ElementRegionBase& operator=(const ElementRegionBase& rhs);

  string_array m_materialList;
  string m_numericalMethod;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENT_REGION_BASE_H_ */
