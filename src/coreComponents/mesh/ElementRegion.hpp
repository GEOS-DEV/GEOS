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

/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTREGION_H
#define ELEMENTREGION_H

#include "CellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "FaceManager.hpp"
#include "FaceElementSubRegion.hpp"

namespace geosx
{

//template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//
//  SUBREGIONTYPE * const subRegion = dynamic_cast<SUBREGIONTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval = true;
//  }
//  else
//  {
//    rval = applyLambdaToCellBlocks< SUBREGIONTYPES..., LAMBDA >( cellSubRegion, std::forward<LAMBDA>(lambda) );
//  }
//
//  return rval;
//}
//
//template< typename SUBREGIONTYPE, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//  SUBREGIONTYPE * const subRegion = dynamic_cast<SUBREGIONTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval=true;
//  }
//
//  return rval;
//}
//
//template< typename... SUBREGIONTYPES, typename LAMBDA >
//void forSomeCellBlocks( LAMBDA && lambda ) const
//{
//  ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);
//
//  for( auto const & subGroupIter : cellBlockSubRegions->GetSubGroups() )
//  {
//    bool isNull =
//    !applyLambdaToCellBlocks<SUBREGIONTYPES...>( subGroupIter.second, [&]( auto * const subRegion )
//    {
//      lambda( subRegion );
//    });
//    GEOS_ERROR_IF( isNull, "subRegion "<<subGroupIter.second->getName()<<" is can not be casted to any "
//                   "types specified in parameter pack.");
//  }
//}






class StableTimeStep;


/**
 * Class to manage the data stored at the element level.
 */
class ElementRegion : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ElementRegion"; }

  virtual const string getCatalogName() const override final
  { return ElementRegion::CatalogName(); }


  ///@}


  ElementRegion() = delete;

  ElementRegion( string const & name, ManagedGroup * const parent );


  ElementRegion(const ElementRegion& init);

  virtual ~ElementRegion() override;

  void GenerateMesh( ManagedGroup const * const cellBlocks );

  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const NodeManager ); 

  void GenerateFractureMesh( FaceManager const * const faceManager );

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

  void AddCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.push_back( cellBlockName );
  }


  template< typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup const * const cellSubRegion, LAMBDA&& lambda )
  { return false; }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup const * const cellSubRegion, LAMBDA&& lambda )
  {
    bool rval = false;
    SUBREGIONTYPE const * const subRegion = dynamic_cast<SUBREGIONTYPE const *>( cellSubRegion );
    if( subRegion!= nullptr )
    {
      lambda( subRegion );
      rval = true;
    }
    else
    {
      rval = applyLambdaToCellBlocks< SUBREGIONTYPES... >( cellSubRegion, std::forward<LAMBDA>(lambda) );
    }
    return rval;
  }

  template< typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
  { return false; }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
  {
    bool rval = false;
    SUBREGIONTYPE * const subRegion = dynamic_cast<SUBREGIONTYPE *>( cellSubRegion );
    if( subRegion!= nullptr )
    {
      lambda( subRegion );
      rval = true;
    }
    else
    {
      rval = applyLambdaToCellBlocks< SUBREGIONTYPES... >( cellSubRegion, std::forward<LAMBDA>(lambda) );
    }
    return rval;
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
    ManagedGroup const * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);

    for( auto const & subGroupIter : elementSubRegions->GetSubGroups() )
    {
      applyLambdaToCellBlocks<SUBREGIONTYPE, SUBREGIONTYPES...>( subGroupIter.second, [&]( auto const * const subRegion )
      {
        lambda( subRegion );
      });
    }
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);

    for( auto & subGroupIter : elementSubRegions->GetSubGroups() )
    {
      applyLambdaToCellBlocks<SUBREGIONTYPE, SUBREGIONTYPES...>( subGroupIter.second, [&]( auto * const subRegion )
      {
        lambda( subRegion );
      });
    }
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
      applyLambdaToCellBlocks<SUBREGIONTYPE,SUBREGIONTYPES...>( subRegion, [&]( auto const * const castedSubRegion )
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
      applyLambdaToCellBlocks<SUBREGIONTYPE,SUBREGIONTYPES...>( subRegion, [&]( auto * const castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      });
    }
  }


  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto materialListString = "materialList";
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto coarseningRatioString = "coarseningRatio"; 
    static constexpr auto elementSubRegions = "elementSubRegions";
    static constexpr auto sourceCellBlockNames = "cellBlocks";

    static constexpr auto fractureConnectorIndicesString = "fractureConnectorIndices";
    static constexpr auto fractureElementConnectorString = "fractureElementConnectors";
    static constexpr auto fractureToCellConnectorString = "fractureCellConnectors";
    static constexpr auto fractureCellConnectorIndicesString = "fractureCellConnectorIndices";


  } m_regionViewKeys;

  string_array & getMaterialList() {return m_materialList;}
  string_array const & getMaterialList() const {return m_materialList;}

protected:
  virtual void PostProcessInput() override;

private:

  ElementRegion& operator=(const ElementRegion& rhs);

  string_array m_cellBlockNames;
  string_array m_fractureSetNames;
  string_array m_materialList;
  string m_numericalMethod;
  real64 m_coarseningRatio;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
