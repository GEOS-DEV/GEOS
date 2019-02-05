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

#include "managers/ObjectManagerBase.hpp"
#include "FaceManager.hpp"
#include "CellBlockSubRegion.hpp"
#include "FaceElementSubRegion.hpp"

namespace geosx
{

//template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//
//  CELLTYPE * const subRegion = dynamic_cast<CELLTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval = true;
//  }
//  else
//  {
//    rval = applyLambdaToCellBlocks< CELLTYPES..., LAMBDA >( cellSubRegion, std::forward<LAMBDA>(lambda) );
//  }
//
//  return rval;
//}
//
//template< typename CELLTYPE, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//  CELLTYPE * const subRegion = dynamic_cast<CELLTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval=true;
//  }
//
//  return rval;
//}
//
//template< typename... CELLTYPES, typename LAMBDA >
//void forSomeCellBlocks( LAMBDA && lambda ) const
//{
//  ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);
//
//  for( auto const & subGroupIter : cellBlockSubRegions->GetSubGroups() )
//  {
//    bool isNull =
//    !applyLambdaToCellBlocks<CELLTYPES...>( subGroupIter.second, [&]( auto * const subRegion )
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

  void GenerateFractureMesh( FaceManager const * const faceManager );

  subGroupMap & GetSubRegions()
  {
    return GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups();
  }

  subGroupMap const & GetSubRegions() const
  {
    return GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups();
  }

  template< typename CELLTYPE=ElementSubRegionBase >
  CELLTYPE const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(regionName);
  }
  template< typename CELLTYPE=ElementSubRegionBase >
  CELLTYPE * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(regionName);
  }

  template< typename CELLTYPE=ElementSubRegionBase >
  CELLTYPE const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(index);
  }
  template< typename CELLTYPE=ElementSubRegionBase >
  CELLTYPE * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(index);
  }

  localIndex numSubRegions() const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups().size();
  }


  template< typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup const * const cellSubRegion, LAMBDA&& lambda )
  { return false; }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup const * const cellSubRegion, LAMBDA&& lambda )
  {
    bool rval = false;
    CELLTYPE const * const subRegion = dynamic_cast<CELLTYPE const *>( cellSubRegion );
    if( subRegion!= nullptr )
    {
      lambda( subRegion );
      rval = true;
    }
    else
    {
      rval = applyLambdaToCellBlocks< CELLTYPES... >( cellSubRegion, std::forward<LAMBDA>(lambda) );
    }
    return rval;
  }

  template< typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
  { return false; }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
  {
    bool rval = false;
    CELLTYPE * const subRegion = dynamic_cast<CELLTYPE *>( cellSubRegion );
    if( subRegion!= nullptr )
    {
      lambda( subRegion );
      rval = true;
    }
    else
    {
      rval = applyLambdaToCellBlocks< CELLTYPES... >( cellSubRegion, std::forward<LAMBDA>(lambda) );
    }
    return rval;
  }

  template< typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda ) const
  {
    forCellBlocks<CellBlockSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda )
  {
    forCellBlocks<CellBlockSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda ) const
  {
    ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    for( auto const & subGroupIter : cellBlockSubRegions->GetSubGroups() )
    {
      applyLambdaToCellBlocks<CELLTYPE, CELLTYPES...>( subGroupIter.second, [&]( auto const * const subRegion )
      {
        lambda( subRegion );
      });
    }
  }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda )
  {
    ManagedGroup * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    for( auto & subGroupIter : cellBlockSubRegions->GetSubGroups() )
    {
      applyLambdaToCellBlocks<CELLTYPE, CELLTYPES...>( subGroupIter.second, [&]( auto * const subRegion )
      {
        lambda( subRegion );
      });
    }
  }


  template< typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda ) const
  {
    forCellBlocksIndex<CellBlockSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda )
  {
    forCellBlocksIndex<CellBlockSubRegion, FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase const * const subRegion = this->GetSubRegion(esr);
      applyLambdaToCellBlocks<CELLTYPE,CELLTYPES...>( subRegion, [&]( auto const * const castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      });
    }
  }

  template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase * const subRegion = this->GetSubRegion(esr);
      applyLambdaToCellBlocks<CELLTYPE,CELLTYPES...>( subRegion, [&]( auto * const castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      });
    }
  }


  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto materialListString = "materialList";
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto cellBlockSubRegions = "cellBlockSubRegions";
    static constexpr auto cellBlockSubRegionNames = "cellBlocks";

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

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
