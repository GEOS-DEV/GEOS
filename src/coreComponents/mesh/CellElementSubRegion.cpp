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


#include "CellElementSubRegion.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, ManagedGroup * const parent ):
  CellBlock( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping, 0)->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::constitutiveMapString,
                       &m_constitutiveMapView, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0);


  RegisterViewWrapper( viewKeyStruct::constitutivePointVolumeFraction,
                       &m_constitutivePointVolumeFraction, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0)->setSizedFromParent(1);

}

CellElementSubRegion::~CellElementSubRegion()
{
  // TODO Auto-generated destructor stub
}

void CellElementSubRegion::initializeDNDXReordered()
{
  static bool initialized = false;
  if ( initialized ) return;

#if STANDARD_ELEMENT_DNDX_LAYOUT
  GEOS_LOG("Using the standard layout for dNdX!!!!");
  m_dNdX_reordered.resize(m_dNdX.size(0), m_dNdX.size(1), m_dNdX.size(2), 3);
#else
  GEOS_LOG("Using the alternate layout for dNdX!!!!");
  m_dNdX_reordered.resize(3, m_dNdX.size(2), m_dNdX.size(1), m_dNdX.size(0));
#endif

  m_dNdX_reordered.setUserCallBack("dNdX_reordered");

  for (localIndex k = 0; k < m_dNdX.size(0); ++k)
  {
    for (localIndex q = 0; q < m_dNdX.size(1); ++q)
    {
      for (localIndex n = 0; n < m_dNdX.size(2); ++n)
      {
        for (int i = 0; i < 3; ++i)
        {
          DNDX_ACCESSOR(m_dNdX_reordered, k, q, n, i) = m_dNdX(k, q, n)[i];
        }
      }
    }
  }

  initialized = true;
}

void CellElementSubRegion::initializeDetJReordered()
{
  static bool initialized = false;
  if ( initialized ) return;

  arrayView2d< double const > const & detJ = getReference< array2d<real64> >("detJ");

#if STANDARD_ELEMENT_DETJ_LAYOUT
  GEOS_LOG("Using the standard layout for detJ!!!!");
  m_detJ_reordered.resize(detJ.size(0), detJ.size(1));
#else
  GEOS_LOG("Using the alternate layout for detJ!!!!");
  m_detJ_reordered.resize(detJ.size(1), detJ.size(0));
#endif

  m_detJ_reordered.setUserCallBack("detJ_reordered");

  for (localIndex k = 0; k < detJ.size(0); ++k)
  {
    for (localIndex q = 0; q < detJ.size(1); ++q)
    {
      DETJ_ACCESSOR(m_detJ_reordered, k, q) = detJ(k, q);
    }
  }

  initialized = true;
}

void CellElementSubRegion::initializeMeanStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex)
{
  static bool initialized = false;
  if ( initialized ) return;

  arrayView2d<real64 const> const & meanStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("MeanStress", constitutiveManager)[er][esr][matIndex];

#if STANDARD_ELEMENT_MEANSTRESS_LAYOUT
  GEOS_LOG("Using the standard layout for meanStress!!!!");
  m_meanStress_reordered.resize(meanStress.size(0), meanStress.size(1));
#else
  GEOS_LOG("Using the alternate layout for meanStress!!!!");
  m_meanStress_reordered.resize(meanStress.size(1), meanStress.size(0));
#endif

  m_meanStress_reordered.setUserCallBack("meanStress_reordered");

  for (localIndex k = 0; k < meanStress.size(0); ++k)
  {
    for (localIndex q = 0; q < meanStress.size(1); ++q)
    {
      MEANSTRESS_ACCESSOR(m_meanStress_reordered, k, q) = meanStress(k, q);
    }
  }

  initialized = true;
}

void CellElementSubRegion::outputMeanStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex)
{
  m_meanStress_reordered.move(chai::CPU, false);

  arrayView2d<real64> const & meanStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("MeanStress", constitutiveManager)[er][esr][matIndex];

  for (localIndex k = 0; k < meanStress.size(0); ++k)
  {
    for (localIndex q = 0; q < meanStress.size(1); ++q)
    {
      meanStress(k, q) = MEANSTRESS_ACCESSOR(m_meanStress_reordered, k, q);
    }
  }
}


void CellElementSubRegion::initializeDeviatorStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex)
{
  static bool initialized = false;
  if ( initialized ) return;

  arrayView2d<R2SymTensor const> const & deviatorStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<R2SymTensor>, arrayView2d<R2SymTensor> >("DeviatorStress", constitutiveManager)[er][esr][matIndex];

#if STANDARD_ELEMENT_DEVIATORSTRESS_LAYOUT
  GEOS_LOG("Using the standard layout for deviatorStress!!!!");
  m_deviatorStress_reordered.resize(deviatorStress.size(0), deviatorStress.size(1), 6);
#else
  GEOS_LOG("Using the alternate layout for deviatorStress!!!!");
  m_deviatorStress_reordered.resize(6, deviatorStress.size(1), deviatorStress.size(0));
#endif

  m_deviatorStress_reordered.setUserCallBack("deviatorStress_reordered");

  for (localIndex k = 0; k < deviatorStress.size(0); ++k)
  {
    for (localIndex q = 0; q < deviatorStress.size(1); ++q)
    {
      for (localIndex i = 0; i < 6; ++i)
      {
        DEVIATORSTRESS_ACCESSOR(m_deviatorStress_reordered, k, q, i) = deviatorStress(k, q).Data()[i];
      }
    }
  }

  initialized = true;
}

void CellElementSubRegion::outputDeviatorStressReordered(ElementRegionManager * elemManager, constitutive::ConstitutiveManager * constitutiveManager, localIndex const er, localIndex const esr, localIndex const matIndex)
{
  m_deviatorStress_reordered.move(chai::CPU, false);

  arrayView2d<R2SymTensor> const & deviatorStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<R2SymTensor>, arrayView2d<R2SymTensor> >("DeviatorStress", constitutiveManager)[er][esr][matIndex];

  for (localIndex k = 0; k < deviatorStress.size(0); ++k)
  {
    for (localIndex q = 0; q < deviatorStress.size(1); ++q)
    {
      for (localIndex i = 0; i < 6; ++i)
      {
        deviatorStress(k, q).Data()[i] = DEVIATORSTRESS_ACCESSOR(m_deviatorStress_reordered, k, q, i);
      }
    }
  }
}


void CellElementSubRegion::CopyFromCellBlock( CellBlock const * source )
{
  this->SetElementType(source->GetElementTypeString());
  this->numNodesPerElement() = source->numNodesPerElement();
  this->numFacesPerElement() = source->numFacesPerElement();
  this->resize(source->size());
  this->nodeList() = source->nodeList();
  this->m_localToGlobalMap = source->m_localToGlobalMap;
  this->ConstructGlobalToLocalMap();
  source->forExternalProperties([&]( const dataRepository::ViewWrapperBase * vw )->void
  {
    std::type_index typeIndex = std::type_index( vw->get_typeid());
    rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                    true,
                                    [&]( auto type, auto baseType ) -> void
    {
      using fieldType = decltype(type);
      const dataRepository::ViewWrapper<fieldType> & field = dataRepository::ViewWrapper< fieldType >::cast( *vw );
      const fieldType & fieldref = field.reference();
      this->RegisterViewWrapper( vw->getName(), &const_cast< fieldType & >( fieldref ), 0 ); //TODO remove const_cast
    });
  });
}

void CellElementSubRegion::ConstructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                                        string const & setName )
{
  set<localIndex> const & targetSet = faceManager->sets()->getReference<set<localIndex> >(setName);
  m_toFacesRelation.resize(0,2);
  this->resize( targetSet.size() );


}

void CellElementSubRegion::MaterialPassThru( string const & matName,
                                           string const & setName,
                                           set<localIndex> & materialSet,
                                           ManagedGroup * material )
{}





void CellElementSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
}


localIndex CellElementSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}


localIndex CellElementSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex CellElementSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                      arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         nodeList().Base(),
                                         m_unmappedGlobalIndicesInNodelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         nodeList().RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         faceList().Base(),
                                         m_unmappedGlobalIndicesInFacelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         faceList().RelatedObjectLocalToGlobal() );

  return packedSize;
}


localIndex CellElementSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                 localIndex_array & packList,
                                                 bool const overwriteUpMaps,
                                                 bool const overwriteDownMaps )
{
  localIndex unPackedSize = 0;



  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->m_globalToLocalMap,
                                     nodeList().RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->m_globalToLocalMap,
                                     faceList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void CellElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}


} /* namespace geosx */
