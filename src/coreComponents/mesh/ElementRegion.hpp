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



class StableTimeStep;

namespace geosx
{

namespace dataRepository
{
namespace keys
{
}
}

class CellBlockSubRegion;



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

  CellBlockSubRegion const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }
  CellBlockSubRegion * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }

  CellBlockSubRegion const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }
  CellBlockSubRegion * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }

  localIndex numSubRegions() const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups().size();
  }


  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda )
  {
    ManagedGroup * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion * subRegion ) -> void
      {
        lambda( subRegion );
      });
  }


  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda ) const
  {
    ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion const * subRegion ) -> void
      {
        lambda( subRegion );
      });
  }

  template< typename LAMBDA >
  void forCellBlocksIndex( LAMBDA lambda ) const
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * cellBlock = this->GetSubRegion(esr);
      lambda( esr, cellBlock );
    }
  }

  template< typename LAMBDA >
  void forCellBlocksIndex( LAMBDA lambda )
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * cellBlock = this->GetSubRegion(esr);
      lambda( esr, cellBlock );
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
