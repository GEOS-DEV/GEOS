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
 * CellBlockSubRegion.hpp
 *
 *  Created on: May 11, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_

#include "CellBlock.hpp"

namespace geosx
{

class CellBlockSubRegion : public CellBlock
{
public:
  CellBlockSubRegion( string const & name, ManagedGroup * const parent );
  virtual ~CellBlockSubRegion() override;

  void ProcessInputFile_PostProcess() override;

  void InitializePreSubGroups( ManagedGroup * const ) override final;

  void InitializePostSubGroups( ManagedGroup * const ) override final;

  void CopyFromCellBlock( CellBlock const * source );

  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  void MaterialPassThru( string const & matName,
                         string const & setName,
                         set<localIndex> & materialSet,
                         ManagedGroup * material );


  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) override final;

  struct viewKeyStruct : public CellBlock::viewKeyStruct
  {
    static constexpr auto constitutivePointVolumeFraction = "ConstitutivePointVolumeFraction";
    static constexpr auto dNdXString = "dNdX";

    static constexpr auto constitutiveGroupingString = "ConstitutiveGrouping";
    static constexpr auto constitutiveMapString = "ConstitutiveMap";



    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString };
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString };
    dataRepository::ViewKey dNdX                  = { dNdXString };

  } m_CellBlockSubRegionViewKeys;

  struct groupKeyStruct : public CellBlock::groupKeyStruct
  {
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
//    constitutiveModelsString = constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString;


  } m_CellBlockSubRegionGroupKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

//  virtual groupKeyStruct & groupKeys() { return m_ObjectManagerBaseGroupKeys; }
//  virtual groupKeyStruct const & groupKeys() const { return m_ObjectManagerBaseGroupKeys; }


  map< string, localIndex_array > m_constitutiveGrouping;

  array3d< real64 > m_constitutivePointVolumeFraction;

  // TODO this needs to be stored by the FiniteElementManager!!
  std::pair< array2d< localIndex >, array2d< localIndex > > m_constitutiveMapView;

  array3d< R1Tensor > m_dNdX;

  dataRepository::ManagedGroup const * GetConstitutiveModels() const
  { return &m_constitutiveModels; }

  dataRepository::ManagedGroup * GetConstitutiveModels()
  { return &m_constitutiveModels; }

private:
  dataRepository::ManagedGroup m_constitutiveModels;

  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInNodelist;
  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInFacelist;

  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_ */
