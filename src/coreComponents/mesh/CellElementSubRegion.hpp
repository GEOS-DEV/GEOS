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


#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_

#include "CellBlock.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#if STANDARD_ELEMENT_DNDX_LAYOUT
#define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX(k, q, n, i)
#else
#if CALC_SHAPE_FUNCTION_DERIVATIVES
#define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX[i][n]
#else
#define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX(i, n, q, k)
#endif
#endif

#if STANDARD_ELEMENT_DETJ_LAYOUT
#define DETJ_ACCESSOR(detJ, k, q) detJ(k, q)
#else
#define DETJ_ACCESSOR(detJ, k, q) detJ(q, k)
#endif

#if STANDARD_ELEMENT_MEANSTRESS_LAYOUT
#define MEANSTRESS_ACCESSOR(meanStress, k, q) meanStress(k, q)
#else
#define MEANSTRESS_ACCESSOR(meanStress, k, q) meanStress(q, k)
#endif

#if STANDARD_ELEMENT_DEVIATORSTRESS_LAYOUT
#define DEVIATORSTRESS_ACCESSOR(deviatorStress, k, q, i) deviatorStress(k, q, i)
#else
#define DEVIATORSTRESS_ACCESSOR(deviatorStress, k, q, i) deviatorStress(i, q, k)
#endif

#if STANDARD_ELEMENT_TONODESRELATION_LAYOUT
#define TONODESRELATION_ACCESSOR(toNodes, k, i) toNodes(k, i)
#else
#define TONODESRELATION_ACCESSOR(toNodes, k, i) toNodes(i, k)
#endif

namespace geosx
{

class CellElementSubRegion : public CellBlock
{
public:
  CellElementSubRegion( string const & name, ManagedGroup * const parent );
  virtual ~CellElementSubRegion() override;

  void initializeDNDXReordered();

  arrayView4d< double const > const & getDNDXReordered() const
  { return m_dNdX_reordered; }


  void initializeDetJReordered();

  arrayView2d< double const > const & getDetJReordered() const
  { return m_detJ_reordered; }

  void initializeMeanStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex);

  void outputMeanStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex);

  arrayView2d< double > const & getMeanStressReordered() const
  { return m_meanStress_reordered; }

  void initializeDeviatorStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex);

  void outputDeviatorStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex);

  arrayView3d< double > const & getDeviatorStressReordered() const
  { return m_deviatorStress_reordered; }

  void initializeToNodesRelationReordered();

  arrayView2d< localIndex const > const & getToNodesRelationReordered() const
  { return m_toNodesRelation_reordered; }

#if SSLE_USE_PATCH_KERNEL
  arrayView2d< localIndex const > const & getPatchToNodesRelationReordered() const
  { return m_patchToNodesRelation_reordered; }
#endif

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

  void MaterialPassThru( string const & matName,
                         string const & setName,
                         set<localIndex> & materialSet,
                         ManagedGroup * material );


  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

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

  array4d< double > m_dNdX_reordered;
  array2d< double > m_detJ_reordered;
  array2d< double > m_meanStress_reordered;
  array3d< double > m_deviatorStress_reordered;
  array2d< localIndex > m_toNodesRelation_reordered;
#if SSLE_USE_PATCH_KERNEL
  array2d< localIndex > m_patchToNodesRelation_reordered;
#endif

private:

  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInNodelist;
  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInFacelist;

  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_ */
