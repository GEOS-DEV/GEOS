/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOSX_MESH_CELLELEMENTSUBREGION_HPP_
#define GEOSX_MESH_CELLELEMENTSUBREGION_HPP_

#include "CellBlock.hpp"

namespace geosx
{

class CellElementSubRegion : public CellBlock
{
public:
  CellElementSubRegion( string const & name, Group * const parent );
  virtual ~CellElementSubRegion() override;

  void CopyFromCellBlock( CellBlock const * source );

  void ConstructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                      string const & setName );

  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  virtual void ViewPackingExclusionList( SortedArray<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) final override;

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


  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

//  virtual groupKeyStruct & groupKeys() { return m_ObjectManagerBaseGroupKeys; }
//  virtual groupKeyStruct const & groupKeys() const { return m_ObjectManagerBaseGroupKeys; }


  map< string, localIndex_array > m_constitutiveGrouping;

  array3d< real64 > m_constitutivePointVolumeFraction;

  // TODO this needs to be stored by the FiniteElementManager!!
  std::pair< array2d< localIndex >, array2d< localIndex > > m_constitutiveMapView;

  array3d< R1Tensor > m_dNdX;


private:

  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInNodelist;
  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInFacelist;

  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;


};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
