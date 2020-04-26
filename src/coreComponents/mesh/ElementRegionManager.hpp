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

/**
 * @file ElementRegionManager.hpp
 */

#ifndef GEOSX_MESH_ELEMENTREGIONMANAGER_HPP
#define GEOSX_MESH_ELEMENTREGIONMANAGER_HPP

#include "CellBlock.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "FaceElementRegion.hpp"
#include "fileIO/schema/schemaUtilities.hpp"
#include "wells/WellElementRegion.hpp"
#include "EmbeddedSurfaceRegion.hpp"

namespace geosx
{

class MeshManager;

/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionManager : public ObjectManagerBase
{
public:

  constexpr static int maxNumNodesPerElem = 8;

  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< VIEWTYPE > >;

  template< typename VIEWTYPE >
  using ElementReferenceAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > >;

  /**
   * The MaterialViewAccessor at the ElementRegionManager level is an 3d array that contains a
   * ReferenceWrapper around the VIEWTYPE. The dimensions are denoted as follows:
   * var[elementRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d< VIEWTYPE > > >;

  template< typename CONSTITUTIVE_TYPE >
  using ConstitutiveRelationAccessor = array1d< array1d< array1d< CONSTITUTIVE_TYPE * > > >;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ZoneManager"; }

  virtual const string getCatalogName() const override final
  { return ElementRegionManager::CatalogName(); }



  ///@}

  ElementRegionManager( string const &, Group * const parent );
  virtual ~ElementRegionManager() override;

  /**
   * @brief Get the number of elements within all ElementSubRegions of type T
   */
  template< typename T = ElementSubRegionBase >
  localIndex getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forElementSubRegions< T >( [&]( ElementSubRegionBase const & cellBlock )
    {
      numElem += cellBlock.size();
    } );
    return numElem;
  }

//  void Initialize(  ){}

  void GenerateMesh( Group * const cellBlockManager );

  void GenerateCellToEdgeMaps( FaceManager const * const faceManager );

  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const nodeManager );

  void GenerateWells( MeshManager * const meshManager, MeshLevel * const meshLevel );

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;
//  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;


  virtual void ExpandObjectCatalogs() override;

  virtual void SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType ) override;

  using Group::resize;

  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

  void SetMaxGlobalIndex();

  subGroupMap const & GetRegions() const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups();
  }
  subGroupMap & GetRegions()
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups();
  }

  template< typename T=ElementRegionBase >
  T const * GetRegion( string const & regionName ) const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( regionName );
  }

  template< typename T=ElementRegionBase >
  T * GetRegion( string const & regionName )
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( regionName );
  }

  template< typename T=ElementRegionBase >
  T const * GetRegion( localIndex const index ) const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( index );
  }

  template< typename T=ElementRegionBase >
  T * GetRegion( localIndex const index )
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( index );
  }

  localIndex numRegions() const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups().size();
  }

  localIndex numCellBlocks() const;


  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegions( LAMBDA && lambda )
  {
    Group * const elementRegions = this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegions( LAMBDA && lambda ) const
  {
    Group const * const elementRegions = this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    Group * const elementRegions = this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    Group const * const elementRegions = this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }



  template< typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion, FaceElementRegion, EmbeddedSurfaceRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion, FaceElementRegion, EmbeddedSurfaceRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }



  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = *this->GetRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( elementRegion, [&]( auto & castedRegion )
      {
        lambda( er, castedRegion );
      } );
    }
  }

  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = *this->GetRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( elementRegion, [&]( auto const & castedRegion )
      {
        lambda( er, castedRegion );
      } );
    }
  }


  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion, FaceElementRegion, EmbeddedSurfaceRegion,
                               WellElementRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion, FaceElementRegion, EmbeddedSurfaceRegion,
                               WellElementRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto & elementRegion )
    {
      lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
    } );
  }

  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto const & elementRegion )
    {
      lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
    } );
  }

  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase &,
                                                   auto & subRegion )
    {
      lambda( subRegion );
    }
      );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase const &,
                                                   auto const & subRegion )
    {
      lambda( subRegion );
    } );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ElementRegionBase &,
                                                                                                                   auto & subRegion )
    {
      lambda( targetIndex, subRegion );
    } );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ElementRegionBase const &,
                                                                                                                   auto const & subRegion )
    {
      lambda( targetIndex, subRegion );
    } );
  }

  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion, WellElementSubRegion >( targetRegions,
                                                                                                                                std::forward< LAMBDA >( lambda ) );
  }

  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion, WellElementSubRegion >( targetRegions,
                                                                                                                                std::forward< LAMBDA >( lambda ) );
  }


  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = *this->GetRegion( er );

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        } );
      }
    }
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = *this->GetRegion( er );

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase const & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        } );
      }
    }
  }


  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementRegions( targetRegions, [&] ( localIndex const targetIndex, ElementRegionBase & elementRegion )
    {
      localIndex const er = elementRegion.getIndexInParent();

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
        {
          lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
        } );
      }
    } );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementRegions( targetRegions, [&] ( localIndex const targetIndex, ElementRegionBase const & elementRegion )
    {
      localIndex const er = elementRegion.getIndexInParent();

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase const & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
        {
          lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
        } );
      }
    } );
  }


  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS > ConstructViewAccessor( string const & name,
                                                    string const & neighborName = string() ) const;
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS > ConstructViewAccessor( string const & name,
                                                    string const & neighborName = string() );

  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  ConstructReferenceAccessor( string const & viewName, string const & neighborName = string() ) const;

  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  ConstructReferenceAccessor( string const & viewName, string const & neighborName = string() );

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  MaterialViewAccessor< LHS >
  ConstructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const * const cm ) const;

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  MaterialViewAccessor< LHS >
  ConstructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const * const cm );

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructMaterialViewAccessor( string const & viewName,
                                 arrayView1d< string const > const & regionNames,
                                 arrayView1d< string const > const & materialNames,
                                 bool const allowMissingViews = false ) const;

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructMaterialViewAccessor( string const & viewName,
                                 arrayView1d< string const > const & regionNames,
                                 arrayView1d< string const > const & materialNames,
                                 bool const allowMissingViews = false );

  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  ConstructFullConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm );

  using Group::PackSize;
  using Group::Pack;
  using ObjectManagerBase::PackGlobalMapsSize;
  using ObjectManagerBase::PackGlobalMaps;
  using ObjectManagerBase::UnpackGlobalMaps;
  using ObjectManagerBase::PackUpDownMapsSize;
  using ObjectManagerBase::PackUpDownMaps;
  using ObjectManagerBase::UnpackUpDownMaps;



  int PackSize( string_array const & wrapperNames,
                ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  int Pack( buffer_unit_type * & buffer,
            string_array const & wrapperNames,
            ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  using ObjectManagerBase::Unpack;
  int Unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor< arrayView1d< localIndex > > & packList );

  int Unpack( buffer_unit_type const * & buffer,
              ElementReferenceAccessor< array1d< localIndex > > & packList );



  int PackGlobalMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  int PackGlobalMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;


  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                        ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList );

  int PackUpDownMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;
  int PackUpDownMapsSize( ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;
  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;


  int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementReferenceAccessor< localIndex_array > & packList,
                        bool const overwriteMap );


  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto elementRegionsGroup = "elementRegionsGroup";
  } m_ElementRegionManagerKeys;


private:
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   string_array const & wrapperNames,
                   ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  template< bool DOPACK, typename T >
  int
  PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                         T const & packList ) const;

  template< typename T >
  int UnpackPrivate( buffer_unit_type const * & buffer,
                     T & packList );

  ElementRegionManager( const ElementRegionManager & );
  ElementRegionManager & operator=( const ElementRegionManager & );
};


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructViewAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      Group const * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference< VIEWTYPE >( viewName );
      }
    }
  }
  return viewAccessor;
}


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::
  ConstructViewAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      Group * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference< VIEWTYPE >( viewName );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::
  ConstructReferenceAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      Group const * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ));
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::
  ConstructReferenceAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      Group * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ));
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::
  ConstructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const * const cm ) const
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase const * const subRegion = elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group const * const constitutiveGroup = subRegion->GetConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm->numSubGroups(); ++matIndex )
      {
        string constitutiveName = cm->GetGroup( matIndex )->getName();
        dataRepository::Group const * const constitutiveRelation = constitutiveGroup->GetGroup( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > const * const
          wrapper = constitutiveRelation->getWrapper< VIEWTYPE >( viewName );

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::
  ConstructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const * const cm )
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group * const constitutiveGroup = subRegion->GetConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm->numSubGroups(); ++matIndex )
      {
        string constitutiveName = cm->GetGroup( matIndex )->getName();
        dataRepository::Group * const constitutiveRelation = constitutiveGroup->GetGroup( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > * const
          wrapper = constitutiveRelation->getWrapper< VIEWTYPE >( viewName );

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructMaterialViewAccessor( string const & viewName,
                                                     arrayView1d< string const > const & regionNames,
                                                     arrayView1d< string const > const & materialNames,
                                                     bool const allowMissingViews ) const
{
  GEOSX_ASSERT_EQ( regionNames.size(), materialNames.size() );
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( GetRegion( kReg )->numSubRegions() );
  }

  subGroupMap const & regionMap = GetRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    GEOSX_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
    ElementRegionBase const & region = *GetRegion( er );

    region.forElementSubRegionsIndex( [&]( localIndex const esr,
                                           ElementSubRegionBase const & subRegion )
    {
      dataRepository::Group const & constitutiveGroup = *subRegion.GetConstitutiveModels();
      dataRepository::Group const * const constitutiveRelation = constitutiveGroup.GetGroup( materialNames[k] );
      GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                      "Material " << materialNames[k] << " not found in " << regionNames[k] << '/' << subRegion.getName() );
      dataRepository::Wrapper< VIEWTYPE > const * const wrapper = constitutiveRelation->getWrapper< VIEWTYPE >( viewName );
      GEOSX_ERROR_IF( !allowMissingViews && wrapper == nullptr, "Material " << materialNames[k] << " does not contain " << viewName );
      if( wrapper != nullptr )
      {
        accessor[er][esr] = wrapper->reference();
      }
    } );
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructMaterialViewAccessor( string const & viewName,
                                                     arrayView1d< string const > const & regionNames,
                                                     arrayView1d< string const > const & materialNames,
                                                     bool const allowMissingViews )
{
  GEOSX_ASSERT_EQ( regionNames.size(), materialNames.size() );
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( GetRegion( kReg )->numSubRegions() );
  }

  subGroupMap const & regionMap = GetRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    GEOSX_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
    ElementRegionBase & region = *GetRegion( er );

    region.forElementSubRegionsIndex( [&]( localIndex const esr, ElementSubRegionBase & subRegion )
    {
      dataRepository::Group & constitutiveGroup = *subRegion.GetConstitutiveModels();
      dataRepository::Group * const constitutiveRelation = constitutiveGroup.GetGroup( materialNames[k] );
      GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                      "Material " << materialNames[k] << " not found in " << regionNames[k] << '/' << subRegion.getName() );
      dataRepository::Wrapper< VIEWTYPE > * const wrapper = constitutiveRelation->getWrapper< VIEWTYPE >( viewName );
      GEOSX_ERROR_IF( !allowMissingViews && wrapper == nullptr, "Material " << materialNames[k] << " does not contain " << viewName );
      if( wrapper != nullptr )
      {
        accessor[er][esr] = wrapper->reference();
      }
    } );
  }
  return accessor;
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ElementRegionManager::ConstructFullConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm )
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group * const
      constitutiveGroup = subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm->numSubGroups(); ++matIndex )
      {
        string const constitutiveName = cm->GetGroup( matIndex )->getName();

        CONSTITUTIVE_TYPE * const
        constitutiveRelation = constitutiveGroup->GetGroup< CONSTITUTIVE_TYPE >( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

}
#endif /* GEOSX_MESH_ELEMENTREGIONMANAGER_HPP */
