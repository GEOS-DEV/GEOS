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

/**
 * @file DofManager.cpp
 */

#include "DofManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationOps.hpp"
#include "mesh/MeshLevel.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

namespace geosx
{

using namespace dataRepository;

DofManager::DofManager( string name, localIndex const verbosity )
  : m_name( std::move( name ) ),
    m_verbosity( verbosity ),
    m_domain( nullptr ),
    m_mesh( nullptr ),
    m_closed( false )
{
  initializeDataStructure();
}

void DofManager::initializeDataStructure()
{
  // we pre-allocate an oversized array to store connectivity type
  // instead of resizing it dynamically as fields are added.
  m_connectivity.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );
  m_connectivity = Connectivity::None;

  // we pre-allocate an oversized array to store sparsity pattern type
  // instead of resizing it dynamically as fields are added.
  m_sparsityPattern.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );
}

void DofManager::clear()
{
  // deallocate index arrays from the mesh
  for( FieldDescription const & field : m_fields )
  {
    removeIndexArray( field );
  }

  // delete internal data
  m_fields.clear();
  m_connectivity.clear();
  m_sparsityPattern.clear();

  initializeDataStructure();

  m_closed = false;
}

void DofManager::setMesh( DomainPartition * const domain,
                          localIndex const meshLevelIndex,
                          localIndex const meshBodyIndex )
{
  // TODO: this should be m_domain != domain
  if( m_domain != nullptr )
  {
    // Domain is changed! Delete old data structure and create new
    clear();
  }
  m_domain = domain;
  m_mesh = m_domain->getMeshBody( meshBodyIndex )->getMeshLevel( meshLevelIndex );
}

// .... DOF MANAGER :: FIELD INDEX
localIndex DofManager::getFieldIndex( string const & key ) const
{
  for( localIndex i = 0; i < m_fields.size(); ++i )
  {
    if( m_fields[i].name == key )
    {
      return i;
    }
  }
  GEOS_ERROR( "Field's string key not found in list of active fields." );
  return -1;
}

bool DofManager::keyInUse( string const & key ) const
{
  for( localIndex i = 0; i < m_fields.size(); ++i )
  {
    if( m_fields[i].name == key )
    {
      return true;
    }
  }
  return false;
}

string DofManager::getKey( string const & fieldName ) const
{
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( fieldName ), "getKey: requested field name must be already existing." );

  // get field index
  localIndex fieldIndex = getFieldIndex( fieldName );

  return m_fields[fieldIndex].key;
}

globalIndex DofManager::numGlobalDofs( string const & fieldName ) const
{
  if( fieldName.length() > 0 )
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( fieldName ), "numGlobalDofs: requested field name must be already existing." );

    // get field index
    localIndex fieldIndex = getFieldIndex( fieldName );

    return m_fields[fieldIndex].numGlobalRows;
  }
  else
  {
    globalIndex sumGlobalDofs = 0;
    for( FieldDescription const & field : m_fields )
    {
      sumGlobalDofs += field.numGlobalRows;
    }
    return sumGlobalDofs;
  }
}

localIndex DofManager::numLocalDofs( string const & fieldName ) const
{
  if( fieldName.length() > 0 )
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( fieldName ), "numLocalDofs: requested field name must be already existing." );

    // get field index
    localIndex fieldIndex = getFieldIndex( fieldName );

    return m_fields[fieldIndex].numLocalRows;
  }
  else
  {
    localIndex sumLocalDofs = 0;
    for( FieldDescription const & field : m_fields )
    {
      sumLocalDofs += field.numLocalRows;
    }
    return sumLocalDofs;
  }
}

localIndex DofManager::offsetLocalDofs( string const & fieldName ) const
{
  if( fieldName.length() > 0 )
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( fieldName ), "offsetLocalDofs: requested field name must be already existing." );

    // get field index
    localIndex fieldIndex = getFieldIndex( fieldName );

    return m_fields[fieldIndex].firstLocalRow;
  }
  else
  {
    localIndex sumOffsetLocalDofs = 0;
    for( FieldDescription const & field : m_fields )
    {
      sumOffsetLocalDofs += field.firstLocalRow;
    }
    return sumOffsetLocalDofs;
  }
}

/* ================================================================================== */

namespace
{

/**
 * @brief A struct to abstract away some details of mesh access
 * @tparam LOC type of mesh location
 */
template< DofManager::Location LOC >
struct MeshHelper
{
};

template<>
struct MeshHelper<DofManager::Location::Node>
{
  using ManagerType = NodeManager;

  static constexpr auto managerGroupName = MeshLevel::groupStructKeys::nodeManagerString;
  static constexpr auto mapViewKey = ElementSubRegionBase::viewKeyStruct::nodeListString;
  static constexpr auto syncObjName = "node";

  template< typename MANAGER >
  using MapType = typename MANAGER::NodeMapType;

  HAS_ALIAS( NodeMapType )

  template< typename MANAGER >
  static bool constexpr hasMapTypeAlias()
  {
    return has_alias_NodeMapType<MANAGER>::value;
  }
};

template<>
struct MeshHelper<DofManager::Location::Edge>
{
  using ManagerType = EdgeManager;

  static constexpr auto managerGroupName = MeshLevel::groupStructKeys::edgeManagerString;
  static constexpr auto mapViewKey = ElementSubRegionBase::viewKeyStruct::edgeListString;
  static constexpr auto syncObjName = "edge";

  template< typename MANAGER >
  using MapType = typename MANAGER::EdgeMapType;

  HAS_ALIAS( EdgeMapType )

  template< typename MANAGER >
  static bool constexpr hasMapTypeAlias()
  {
    return has_alias_EdgeMapType<MANAGER>::value;
  }
};

template<>
struct MeshHelper<DofManager::Location::Face>
{
  using ManagerType = FaceManager;

  static constexpr auto managerGroupName = MeshLevel::groupStructKeys::faceManagerString;
  static constexpr auto mapViewKey = ElementSubRegionBase::viewKeyStruct::faceListString;
  static constexpr auto syncObjName = "face";

  template< typename MANAGER >
  using MapType = typename MANAGER::FaceMapType;

  HAS_ALIAS( FaceMapType )

  template< typename MANAGER >
  static bool constexpr hasMapTypeAlias()
  {
    return has_alias_FaceMapType<MANAGER>::value;
  }
};

template<>
struct MeshHelper<DofManager::Location::Elem>
{
  using ManagerType = ElementSubRegionBase;

  static constexpr auto managerGroupName = MeshLevel::groupStructKeys::elemManagerString;
  static constexpr auto mapViewKey = FaceManager::viewKeyStruct::elementListString; // TODO which key
  static constexpr auto syncObjName = "elems";

  template< typename MANAGER >
  using MapType = typename MANAGER::ElemMapType;

  HAS_ALIAS( ElemMapType )

  template< typename MANAGER >
  static bool constexpr hasMapTypeAlias()
  {
    return has_alias_ElemMapType<MANAGER>::value;
  }
};

template< DofManager::Location LOC, typename MANAGER, bool >
struct MapTypeHelper
{
  using type = FixedOneToManyRelation; // dummy type
};

template< DofManager::Location LOC, typename MANAGER >
struct MapTypeHelper<LOC, MANAGER, true>
{
  using type = typename MeshHelper<LOC>::template MapType<MANAGER>;
};

// return dummy type if target manager type does not declare a type alias to map to LOC objects
// this allows all switchyards to compile, but one shouldn't attempt to access a non-existent map
template< DofManager::Location LOC, typename MANAGER >
using MapType = typename MapTypeHelper< LOC, MANAGER, MeshHelper<LOC>::template hasMapTypeAlias<MANAGER>() >::type;

// some helper crust to extract underlying type from InterObjectRelation and the likes
template< typename T, bool >
struct BaseTypeHelper
{
  using type = T;
};

template< typename T >
struct BaseTypeHelper<T, true>
{
  using type = typename T::base_type;
};

HAS_ALIAS( base_type )

template< typename MAP >
using BaseType = typename BaseTypeHelper<MAP, has_alias_base_type<MAP>::value>::type;

/**
 * @brief Helper struct that specializes access to various map types
 * @tparam MAP type of the map
 */
template< typename MAP >
struct MapHelperImpl
{
};

template< typename T >
struct MapHelperImpl< array2d<T> >
{
  static localIndex size0( array2d<T> const & map )
  {
    return map.size( 0 );
  }

  static localIndex size1( array2d<T> const & map,
                           localIndex const i0 )
  {
    return map.size( 1 );
  }

  static T const & value( array2d<T> const & map,
                          localIndex const i0,
                          localIndex const i1 )
  {
    return map( i0, i1 );
  }
};

template< typename T >
struct MapHelperImpl< array1d< array1d<T> > >
{
  static localIndex size0( array1d<array1d<T> > const & map )
  {
    return map.size();
  }

  static localIndex size1( array1d< array1d<T> > const & map,
                           localIndex const i0 )
  {
    return map[i0].size();
  }

  static T const & value( array1d< array1d<T> > const & map,
                          localIndex const i0,
                          localIndex const i1 )
  {
    return map[i0][i1];
  }
};

template< typename T >
struct MapHelperImpl< ArrayOfArrays<T> >
{
  static localIndex size0( ArrayOfArrays<T> const & map )
  {
    return map.size();
  }

  static localIndex size1( ArrayOfArrays<T> const & map,
                           localIndex const i0 )
  {
    return map.sizeOfArray( i0 );
  }

  static T const & value( ArrayOfArrays<T> const & map,
                          localIndex const i0,
                          localIndex const i1 )
  {
    return map( i0, i1 );
  }
};

template< typename T >
struct MapHelperImpl< array1d< set<T> > >
{
  static localIndex size0( array1d< set<T> > const & map )
  {
    return map.size();
  }

  static localIndex size1( array1d< set<T> > const & map,
                           localIndex const i0 )
  {
    return map[i0].size();
  }

  static T const & value( array1d< set<T> > const & map,
                          localIndex const i0,
                          localIndex const i1 )
  {
    return map[i0][i1];
  }
};

template< typename BASETYPE >
struct MapHelperImpl< ToElementRelation<BASETYPE> >
{
  static localIndex size0( ToElementRelation<BASETYPE> const & map )
  {
    return MapHelperImpl<BASETYPE>::size0( map.m_toElementIndex );
  }

  static localIndex size1( ToElementRelation<BASETYPE> const & map,
                           localIndex const i0 )
  {
    return MapHelperImpl<BASETYPE>::size1( map.m_toElementIndex, i0 );
  }

  static auto value( ToElementRelation<BASETYPE> const & map,
                     localIndex const i0,
                     localIndex const i1 )
  {
    return std::make_tuple( MapHelperImpl<BASETYPE>::value( map.m_toElementRegion, i0, i1 ),
                            MapHelperImpl<BASETYPE>::value( map.m_toElementSubRegion, i0, i1 ),
                            MapHelperImpl<BASETYPE>::value( map.m_toElementIndex, i0, i1 ) );
  }
};

/**
 * @brief Helper struct that specializes access to various map types
 * @tparam MAP type of the map
 *
 * @note We may need to strip off InterObjectRelation and get the underlying map type, hence the extra layer
 */
template< typename MAP >
using MapHelper = MapHelperImpl< BaseType<MAP> >;

/**
 * @brief Switchyard for Location (Node/Edge/Face/Elem)
 * @tparam LAMBDA generic functor type
 * @param loc location type
 * @param lambda functor to be called
 */
template< typename LAMBDA >
bool LocationSwitch( DofManager::Location const loc,
                     LAMBDA lambda )
{
  switch( loc )
  {
    case DofManager::Location::Node:
      lambda( std::integral_constant<DofManager::Location, DofManager::Location::Node>() );
      return true;
    case DofManager::Location::Edge:
      lambda( std::integral_constant<DofManager::Location, DofManager::Location::Edge>() );
      return true;
    case DofManager::Location::Face:
      lambda( std::integral_constant<DofManager::Location, DofManager::Location::Face>() );
      return true;
    case DofManager::Location::Elem:
      lambda( std::integral_constant<DofManager::Location, DofManager::Location::Elem>() );
      return true;
    default:
      return false;
  }
}

template< typename LAMBDA >
bool LocationSwitch( DofManager::Location const loc1,
                     DofManager::Location const loc2,
                     LAMBDA lambda )
{
  bool ret2;
  bool const ret1 =
  LocationSwitch( loc1, [&]( auto const loc_type1 )
  {
    ret2 =
    LocationSwitch( loc2, [&]( auto const loc_type2 )
    {
      lambda( loc_type1, loc_type2 );
    } );
  } );
  return ret1 && ret2;
}

template< DofManager::Location LOC >
typename MeshHelper<LOC>::ManagerType const * getObjectManager( MeshLevel const * const meshLevel )
{
  using ObjectManager = typename MeshHelper<LOC>::ManagerType;
  GEOS_ASSERT( meshLevel != nullptr );
  ObjectManager const * manager = meshLevel->GetGroup<ObjectManager>( MeshHelper<LOC>::managerGroupName );
  GEOS_ASSERT( manager != nullptr );
  return manager;
}

template< DofManager::Location LOC >
typename MeshHelper<LOC>::ManagerType * getObjectManager( MeshLevel * const meshLevel )
{
  using ObjectManager = typename MeshHelper<LOC>::ManagerType;
  return const_cast<ObjectManager *>( getObjectManager<LOC>( const_cast<MeshLevel const *>( meshLevel ) ) );
}

template< DofManager::Location LOC, DofManager::Location CONN_LOC >
struct MeshLoopHelper;

template< DofManager::Location LOC >
struct MeshLoopHelper<LOC, LOC>
{
  template< typename ... SUBREGIONTYPES, typename LAMBDA >
  static void visit( MeshLevel * const meshLevel,
                     array1d<string> const & regions,
                     LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper<LOC>::ManagerType;

    // get access to location ghost rank (we don't want to visit ghosted locations
    ObjectManagerLoc const * const objectManager = getObjectManager<LOC>( meshLevel );
    array1d<integer> const & ghostRank = objectManager->GhostRank();

    // create an array to track previously visited locations (to avoid multiple visits)
    array1d<integer> locVisited( objectManager->size() );
    locVisited = 0;

    meshLevel->getElemManager()->
      forElementSubRegionsComplete<SUBREGIONTYPES...>( regions, [&]( localIndex const er,
                                                                     localIndex const esr,
                                                                     ElementRegionBase *,
                                                                     auto * subRegion )
    {
      // derive some more useful, subregion-dependent type aliases
      using ElementSubRegionType = std::remove_pointer_t<decltype( subRegion )>;
      using ElemToLocMapType = MapType<LOC, ElementSubRegionType>;

      // get access to element-to-location map
      ElemToLocMapType const & elemToLocMap =
        subRegion->template getReference<ElemToLocMapType>( MeshHelper<LOC>::mapViewKey );

      // loop over all elements (including ghosts, which may be necessary to access some locally owned locations)
      for( localIndex ei = 0; ei < subRegion->size(); ++ei )
      {
        // loop over all locations incident on an element
        for( localIndex a = 0; a < MapHelper<ElemToLocMapType>::size1( elemToLocMap, ei ); ++a )
        {
          localIndex const locIdx = MapHelper<ElemToLocMapType>::value( elemToLocMap, ei, a );

          // check if we should visit this location
          if( ghostRank[locIdx] < 0 && !std::exchange( locVisited[locIdx], 1 ) )
          {
            lambda( locIdx, locIdx, 0 );
          }
        }
      }
    } );
  }
};

template< DofManager::Location LOC, DofManager::Location CONN_LOC >
struct MeshLoopHelper
{
  template< typename ... SUBREGIONTYPES, typename LAMBDA >
  static void visit( MeshLevel * const meshLevel,
                     array1d<string> const & regions,
                     LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper<LOC>::ManagerType;
    using LocToConnMapType = MapType<CONN_LOC, ObjectManagerLoc>;

    // get access to location ghost rank (we don't want to visit ghosted locations
    ObjectManagerLoc const * const objectManager = getObjectManager<LOC>( meshLevel );

    // get access to location-to-connected map
    LocToConnMapType const & locToConnMap =
      objectManager->template getReference<LocToConnMapType>( MeshHelper<CONN_LOC>::mapViewKey );

    // call the specialized version first, then add an extra loop over connected objects
    MeshLoopHelper<LOC, LOC>::template visit<SUBREGIONTYPES...>( meshLevel, regions,
                                                                 [&]( localIndex const locIdx,
                                                                      localIndex const,
                                                                      localIndex const )
    {
      // loop over all connected locations
      for( localIndex b = 0; b < MapHelper<LocToConnMapType>::size1( locToConnMap, locIdx ); ++b )
      {
        auto const connIdx = MapHelper<LocToConnMapType>::value( locToConnMap, locIdx, b );
        lambda( locIdx, connIdx, b );
      }
    } );
  }
};

template< DofManager::Location LOC >
struct MeshLoopHelper<LOC, DofManager::Location::Elem>
{
  template< typename ... SUBREGIONTYPES, typename LAMBDA >
  static void visit( MeshLevel * const meshLevel,
                     array1d<string> const & regions,
                     LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper<LOC>::ManagerType;
    using LocToConnMapType = BaseType< MapType< DofManager::Location::Elem, ObjectManagerLoc> >;

    // get access to location ghost rank (we don't want to visit ghosted locations
    ObjectManagerLoc const * const objectManager = getObjectManager<LOC>( meshLevel );

    // get access to location-to-connected map
    LocToConnMapType const & elemRegionList =
      objectManager->template getReference<LocToConnMapType>( ObjectManagerLoc::viewKeyStruct::elementRegionListString );
    LocToConnMapType const & elemSubRegionList =
      objectManager->template getReference<LocToConnMapType>( ObjectManagerLoc::viewKeyStruct::elementSubRegionListString );
    LocToConnMapType const & elemIndexList =
      objectManager->template getReference<LocToConnMapType>( ObjectManagerLoc::viewKeyStruct::elementListString );

    // call the specialized version first, then add an extra loop over connected objects
    MeshLoopHelper<LOC, LOC>::template visit<SUBREGIONTYPES...>( meshLevel, regions,
                                                                 [&]( localIndex const locIdx,
                                                                      localIndex const,
                                                                      localIndex const )
    {
      // loop over all connected locations
      for( localIndex b = 0; b < MapHelper<LocToConnMapType>::size1( elemIndexList, locIdx ); ++b )
      {
        localIndex const er  = MapHelper<LocToConnMapType>::value( elemRegionList, locIdx, b );
        localIndex const esr = MapHelper<LocToConnMapType>::value( elemSubRegionList, locIdx, b );
        localIndex const ei  = MapHelper<LocToConnMapType>::value( elemIndexList, locIdx, b );

        if( er >= 0 && esr >= 0 && ei >= 0 )
        {
          lambda( locIdx, std::make_tuple( er, esr, ei ), b );
        }
      }
    } );
  }
};

template< DofManager::Location CONN_LOC >
struct MeshLoopHelper<DofManager::Location::Elem, CONN_LOC>
{
  template< typename ... SUBREGIONTYPES, typename LAMBDA >
  static void visit( MeshLevel * const meshLevel,
                     array1d<string> const & regions,
                     LAMBDA lambda )
  {
    meshLevel->getElemManager()->
      forElementSubRegionsComplete<SUBREGIONTYPES...>( regions, [&]( localIndex const er,
                                                                     localIndex const esr,
                                                                     ElementRegionBase *,
                                                                     auto * subRegion )
    {
      // derive some more useful, subregion-dependent type aliases
      using ElementSubRegionType = std::remove_pointer_t<decltype( subRegion )>;
      using ElemToConnMapType = MapType<CONN_LOC, ElementSubRegionType>;

      // get access to element-to-location map
      ElemToConnMapType const & elemToConnMap =
        subRegion->template getReference<ElemToConnMapType>( MeshHelper<CONN_LOC>::mapViewKey );

      arrayView1d<integer const> const & elemGhostRank = subRegion->GhostRank();

      for( localIndex ei = 0; ei < subRegion->size(); ++ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          auto const elemIdx = std::make_tuple( er, esr, ei );

          for( localIndex a = 0; a < MapHelper<ElemToConnMapType>::size1( elemToConnMap, ei ); ++a )
          {
            localIndex const locIdx = MapHelper<ElemToConnMapType>::value( elemToConnMap, ei, a );
            lambda( elemIdx, locIdx, a );
          }
        }
      }
    } );
  }
};

template<>
struct MeshLoopHelper<DofManager::Location::Elem, DofManager::Location::Elem>
{
  template< typename ... SUBREGIONTYPES, typename LAMBDA >
  static void visit( MeshLevel * const meshLevel,
                     array1d<string> const & regions,
                     LAMBDA lambda )
  {
    meshLevel->getElemManager()->
      forElementSubRegionsComplete<SUBREGIONTYPES...>( regions, [&]( localIndex const er,
                                                                     localIndex const esr,
                                                                     ElementRegionBase *,
                                                                     auto * subRegion )
    {
      arrayView1d<integer const> const & elemGhostRank = subRegion->GhostRank();

      for( localIndex ei = 0; ei < subRegion->size(); ++ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          auto const elemIdx = std::make_tuple( er, esr, ei );
          lambda( elemIdx, elemIdx, 0 );
        }
      }
    } );
  }
};

/**
 * @brief Visit mesh locations within active regions and their connected locations (e.g. nodes of an element)
 *
 * @tparam LOC type of location (Node/Edge/Face/Elem)
 * @tparam CONN_LOC type of connected location (Node/Edge/Face/Elem)
 * @tparam SUBREGIONTYPES variadic pack of subregion types to visit;
 *         if omitted, uses ElementRegionManager's default behavior
 * @tparam LAMBDA type of user-provided functor or lambda to call
 * @param mesh the mesh to loop over
 * @param regions list of region names to visit
 * @param lambda functor or lambda to call
 *
 * Calls user-provided lambda with 3 parameters:
 * - mesh local index of primary location (one index for node/edge/face, reg/subreg/index tuple for elems)
 * - mesh local index of connected location
 * - index of connected location inside primary location's map (e.g. which node of an element)
 * If CONN_LOC == LOC, lambda is called only once per location, as lambda(loc, loc, 0)
 *
 * @note This function will exclude ghosted locations, but will not perform a similar check for
 *       connected locations - these may include ghosts, it is up to the user to filter as needed.
 *       Similarly, while primary loop is limited to @p regions, adjacent locations may not belong.
 */
template< DofManager::Location LOC, DofManager::Location CONN_LOC, typename ... SUBREGIONTYPES, typename LAMBDA >
void forMeshLocation( MeshLevel * const mesh,
                      array1d<string> const & regions,
                      LAMBDA && lambda )
{
  MeshLoopHelper<LOC, CONN_LOC>::template visit<SUBREGIONTYPES...>( mesh,
                                                                    regions,
                                                                    std::forward<LAMBDA>( lambda ) );
}

/**
 * @brief A shortcut for previous function with CONN_LOC == LOC
 */
template< DofManager::Location LOC, typename ... SUBREGIONTYPES, typename LAMBDA >
void forMeshLocation( MeshLevel * const mesh,
                      array1d<string> const & regions,
                      LAMBDA && lambda )
{
  forMeshLocation<LOC, LOC, SUBREGIONTYPES...>( mesh,
                                                regions,
                                                std::forward<LAMBDA>( lambda ) );
}

/**
 * @brief This template abstracts some differences between index arrays that live
 *        on elements vs other mesh objects: (de)registration, value access, etc.
 * @tparam LOC target location
 */
template< typename INDEX, DofManager::Location LOC >
struct IndexArrayHelper
{
  using ArrayType = array1d< std::remove_const_t<INDEX> >;
  using ViewType = arrayView1d<INDEX>;
  using Accessor = ViewType const;
  using Mesh = add_const_if_t< MeshLevel, std::is_const<INDEX>::value >;

  template< typename ... SUBREGIONTYPES >
  static void
  create( Mesh * const mesh, DofManager::FieldDescription & field )
  {
    GEOS_ASSERT( field.location == LOC );

    ObjectManagerBase * baseManager = getObjectManager<LOC>( mesh );
    baseManager->RegisterViewWrapper<ArrayType>( field.key )->
      setApplyDefaultValue( -1 )->
      setPlotLevel( PlotLevel::LEVEL_1 )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setDescription( field.docstring );
  }

  static Accessor get( Mesh * const mesh, DofManager::FieldDescription const & field )
  {
    auto * baseManager = getObjectManager<LOC>( mesh );
    return baseManager->template getReference<ArrayType>( field.key );
  }

  static inline INDEX value( Accessor & indexArray, localIndex const i )
  {
    return indexArray[i];
  }

  static inline INDEX & reference( Accessor & indexArray, localIndex const i )
  {
    return indexArray[i];
  }

  template< typename ... SUBREGIONTYPES >
  static void
  remove( Mesh * const mesh, DofManager::FieldDescription const & field )
  {
    getObjectManager<LOC>( mesh )->DeregisterViewWrapper( field.key );
  }
};

template< typename INDEX >
struct IndexArrayHelper< INDEX, DofManager::Location::Elem >
{
  using ArrayType = array1d< std::remove_const_t<INDEX> >;
  using ViewType = arrayView1d<INDEX>;
  using Accessor = ElementRegionManager::ElementViewAccessor< ViewType > const;
  using Mesh = add_const_if_t< MeshLevel, std::is_const<INDEX>::value >;

  template< typename ... SUBREGIONTYPES >
  static void
  create( Mesh * const mesh, DofManager::FieldDescription & field )
  {
    GEOS_ASSERT( field.location == DofManager::Location::Elem );

    mesh->getElemManager()->
      template forElementSubRegionsComplete<SUBREGIONTYPES...>( field.regionNames,
                                                                [&]( localIndex const er,
                                                                     localIndex const esr,
                                                                     auto * const region,
                                                                     auto * const subRegion )
    {
      subRegion->template RegisterViewWrapper<ArrayType>( field.key )->
        setApplyDefaultValue( -1 )->
        setPlotLevel( PlotLevel::LEVEL_1 )->
        setRestartFlags( RestartFlags::NO_WRITE )->
        setDescription( field.docstring );
    } );
  }

  static Accessor get( Mesh * const mesh,
                       DofManager::FieldDescription const & field )
  {
    return mesh->getElemManager()->template ConstructViewAccessor< ArrayType, ViewType >( field.key );
  }

  static inline INDEX value( Accessor & indexArray,
                             std::tuple<localIndex, localIndex, localIndex> const & e )
  {
    if( indexArray[std::get<0>( e )].empty() || indexArray[std::get<0>( e )][std::get<1>( e )].empty() )
    {
      return -1;
    }
    return indexArray[std::get<0>( e )][std::get<1>( e )][std::get<2>( e )];
  }

  static inline INDEX & reference( Accessor & indexArray,
                                   std::tuple<localIndex, localIndex, localIndex> const & e )
  {
    return indexArray[std::get<0>( e )][std::get<1>( e )][std::get<2>( e )];
  }

  template< typename ... SUBREGIONTYPES >
  static void
  remove( Mesh * const mesh, DofManager::FieldDescription const & field )
  {
    GEOS_ASSERT( field.location == DofManager::Location::Elem );

    mesh->getElemManager()->
      template forElementSubRegionsComplete<SUBREGIONTYPES...>( field.regionNames,
                                                               [&]( localIndex const er,
                                                                    localIndex const esr,
                                                                    auto * const,
                                                                    auto * const subRegion )
    {
      subRegion->DeregisterViewWrapper( field.key );
    } );
  }
};

template< DofManager::Location LOC, typename ... SUBREGIONTYPES >
void createIndexArrayImpl( DomainPartition * const domain,
                           MeshLevel * const mesh,
                           DofManager::FieldDescription & field )
{
  using helper = IndexArrayHelper<globalIndex, LOC>;

  // 0. register index arrays
  helper::template create<SUBREGIONTYPES ...>( mesh, field );
  typename helper::Accessor & indexArray = helper::get( mesh, field );

  // step 1. loop over all active regions, determine number of local mesh objects
  field.numLocalNodes = 0;
  forMeshLocation<LOC, SUBREGIONTYPES...>( mesh, field.regionNames,
                                           [&]( auto const locIdx,
                                                auto const,
                                                localIndex const )
  {
    helper::reference( indexArray, locIdx ) = field.numComponents * field.numLocalNodes++;
  } );
  field.numLocalRows = field.numComponents * field.numLocalNodes;

  // step 2. gather row counts across ranks
  std::tie( field.firstLocalRow, field.numGlobalRows ) =
    CommunicationTools::PrefixSum<globalIndex>( field.numLocalRows );

  // step 3. adjust local dof offsets to reflect processor offset
  forMeshLocation<LOC, SUBREGIONTYPES...>( mesh, field.regionNames,
                                           [&]( auto const locIdx,
                                                auto const,
                                                localIndex const )
  {
    helper::reference( indexArray, locIdx ) += field.firstLocalRow;
  } );

  // step 4. synchronize across ranks
  std::map<string, string_array> fieldNames;
  fieldNames[ MeshHelper<LOC>::syncObjName ].push_back( field.key );

  CommunicationTools::
  SynchronizeFields( fieldNames, mesh,
                     domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
}

} // namespace

/* ================================================================================== */

template< typename ... SUBREGIONTYPES >
void DofManager::createIndexArray( FieldDescription & field )
{
  bool const success =
  LocationSwitch( field.location, [&]( auto const loc )
  {
    Location constexpr LOC = decltype(loc)::value;
    createIndexArrayImpl<LOC, SUBREGIONTYPES...>( m_domain, m_mesh, field );
  } );
  GEOS_ERROR_IF( !success, "createIndexArray: invalid location type" );
}

template< typename ... SUBREGIONTYPES >
void DofManager::removeIndexArray( FieldDescription const & field )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    Location constexpr LOC = decltype(loc)::value;
    IndexArrayHelper<globalIndex, LOC>::template remove<SUBREGIONTYPES...>( m_mesh, field );
  } );
}

// Just an interface to allow only three parameters
void DofManager::addField( string const & fieldName,
                           Location const location,
                           Connectivity const connectivity )
{
  addField( fieldName, location, connectivity, 1, array1d<string>() );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addField( string const & fieldName,
                           Location const location,
                           Connectivity const connectivity,
                           localIndex const components )
{
  addField( fieldName, location, connectivity, components, array1d<string>() );
}

// Just another interface to allow four parameters (no components)
void DofManager::addField( string const & fieldName,
                           Location const location,
                           Connectivity const connectivity,
                           string_array const & regions )
{
  addField( fieldName, location, connectivity, 1, regions );
}

// The real function, allowing the creation of self-connected blocks
void DofManager::addField( string const & fieldName,
                           Location const location,
                           Connectivity const connectivity,
                           localIndex const components,
                           string_array const & regions )
{
  GEOS_ERROR_IF( m_closed, "addField: cannot add fields after DofManager has been closed." );
  GEOS_ERROR_IF( keyInUse( fieldName ), "addField: requested field name matches an existing field in the DofManager." );
  GEOS_ERROR_IF( m_fields.size() >= MAX_NUM_FIELDS, "addField: limit on DofManager's MAX_NUM_FIELDS exceeded." );

  localIndex fieldIndex = m_fields.size();
  m_fields.resize( fieldIndex + 1 );

  string suffix;
  for( string const & regionName : regions )
  {
    suffix.append( "_" + regionName );
  }

  FieldDescription & field = m_fields.back();

  field.name = fieldName;
  field.location = location;
  field.numComponents = components;
  field.key = m_name + '_' + fieldName + "_dofIndex" + suffix;
  field.docstring = fieldName + " DoF indices";

  if( components > 1 )
  {
    field.docstring += " (with " + std::to_string( components ) + "-component blocks)";
  }

  // save pointers to "active" element regions
  ElementRegionManager * const elemManager = m_mesh->getElemManager();

  // retrieve full list of regions
  if( regions.empty() )
  {
    elemManager->forElementRegions( [&]( ElementRegionBase const * const region )
    {
      field.regionNames.push_back( region->getName() );
    } );
  }
  else
  {
    field.regionNames = regions;
  }

  // sort and remove duplicates
  std::sort( field.regionNames.begin(), field.regionNames.end() );
  auto end_it = std::unique( field.regionNames.begin(), field.regionNames.end() );
  localIndex const numActiveRegions = std::distance( field.regionNames.begin(), end_it );
  field.regionNames.resize( numActiveRegions );

  // check region existence
  for( string const & regionName : field.regionNames )
  {
    GEOS_ERROR_IF( elemManager->GetRegion( regionName ) == nullptr,
                   "addField: specified element region not found: " << regionName );
  }

  // based on location, allocate an index array for this field
  createIndexArray( field );

  // determine field's global offset
  if( fieldIndex > 0 )
  {
    FieldDescription & prev = m_fields[fieldIndex - 1];
    field.fieldOffset = prev.fieldOffset + prev.numGlobalRows;
  }
  else
  {
    field.fieldOffset = 0;
  }

  // save field's connectivity type (self-to-self)
  m_connectivity[fieldIndex][fieldIndex] = connectivity;

  // add sparsity pattern (LC matrix)
  std::unique_ptr<ParallelMatrix> & connLocPattern = m_sparsityPattern( fieldIndex, fieldIndex ).first;
  connLocPattern = std::make_unique<ParallelMatrix>();
  makeConnLocPattern( field, connectivity, field.regionNames, *connLocPattern );

  // log some basic info
  if( m_verbosity > 0 )
  {
    GEOS_LOG_RANK_0( "DofManager :: Added field .... " << field.docstring );
    GEOS_LOG_RANK_0( "DofManager :: Global dofs .... " << field.numGlobalRows );
    GEOS_LOG_RANK_0( "DofManager :: Field offset ... " << field.fieldOffset );
  }
}

// addField: allow the usage of a predefine location-connection pattern (user-defined)
// Interface to allow only two parameters
void DofManager::addField( string const & fieldName,
                           ParallelMatrix const & connLocInput )
{
  DofManager::addField( fieldName,
                        connLocInput,
                        1,
                        Connectivity::USER_DEFINED );
}

// Just another interface to allow three parameters (no connectivity)
void DofManager::addField( string const & fieldName,
                           ParallelMatrix const & connLocInput,
                           localIndex const components )
{
  DofManager::addField( fieldName,
                        connLocInput,
                        components,
                        Connectivity::USER_DEFINED );
}

// Just another interface to allow three parameters (no components)
void DofManager::addField( string const & fieldName,
                           ParallelMatrix const & connLocInput,
                           Connectivity const connectivity )
{
  DofManager::addField( fieldName,
                        connLocInput,
                        1,
                        connectivity );
}

// The real function
void DofManager::addField( string const & fieldName,
                           ParallelMatrix const & connLocInput,
                           localIndex const components,
                           Connectivity const connectivity )
{
  GEOS_ERROR_IF( m_closed, "addField: cannot add fields after DofManager has been closed." );
  GEOS_ERROR_IF( keyInUse( fieldName ), "addField: requested field name matches an existing field in the DofManager." );
  GEOS_ERROR_IF( m_fields.size() > MAX_NUM_FIELDS, "addField: limit on DofManager's MAX_NUM_FIELDS exceeded." );

  localIndex const fieldIndex = m_fields.size();
  m_fields.resize( fieldIndex + 1 );

  FieldDescription & field = m_fields.back();

  field.name = fieldName;
  field.location = Location::USER_DEFINED;
  field.numComponents = components;
  field.key = m_name + '_' + fieldName + "_dofIndex";
  field.docstring = fieldName + " DoF indices";

  if( components > 1 )
  {
    field.docstring += " (with " + std::to_string( components ) + "-component blocks)";
  }

  // determine field's global offset
  if( fieldIndex > 0 )
  {
    FieldDescription & prev = m_fields[fieldIndex - 1];
    field.fieldOffset = prev.fieldOffset + prev.numGlobalRows;
  }
  else
  {
    field.fieldOffset = 0;
  }

  // save field's connectivity type (self-to-self)
  m_connectivity[fieldIndex][fieldIndex] = connectivity;

  // compute useful values (number of local and global rows)
  field.numLocalNodes = connLocInput.localCols();
  field.numLocalRows = field.numLocalNodes * components;

  std::tie( field.firstLocalRow, field.numGlobalRows ) =
    CommunicationTools::PrefixSum<globalIndex>( field.numLocalRows );

  // create the pattern from the user-provided location-connectivity matrix
  localIndex const nrows = connLocInput.localRows();
  localIndex const entriesPerRow = ( nrows > 0 ) ? connLocInput.localNonzeros() / nrows * components : 0;

  std::unique_ptr<ParallelMatrix> & connLocPattern = m_sparsityPattern( fieldIndex, fieldIndex ).first;
  connLocPattern = std::make_unique<ParallelMatrix>();
  connLocPattern->createWithLocalSize( nrows, field.numLocalRows, entriesPerRow, MPI_COMM_GEOSX );

  array1d<globalIndex> colsInput, cols;
  array1d<real64> valsInput, vals;

  // User-provided matrix has a single nnz per loc/conn pair, expand it into multiple (per-component)
  for( globalIndex irow = connLocInput.ilower(); irow < connLocInput.iupper(); ++irow )
  {
    connLocInput.getRowCopy( irow, colsInput, valsInput );
    cols.resize( colsInput.size() * components );
    vals.resize( colsInput.size() * components );
    localIndex k = 0;
    for( localIndex i = 0; i < colsInput.size(); ++i )
    {
      for( localIndex c = 0; c < components; ++c )
      {
        cols[k] = colsInput[i] * components + c;
        vals[k] = (valsInput[i] - 1) * components + c + 1;
        ++k;
      }
    }
    connLocPattern->insert( irow, cols, vals );
  }

  connLocPattern->close();

  // log some basic info
  if( m_verbosity > 0 )
  {
    GEOS_LOG_RANK_0( "DofManager :: Added field .... " << field.docstring );
    GEOS_LOG_RANK_0( "DofManager :: Global dofs .... " << field.numGlobalRows );
    GEOS_LOG_RANK_0( "DofManager :: Field offset ... " << field.fieldOffset );
  }
}

// Create the sparsity pattern (location-location). High level interface
void DofManager::setSparsityPattern( ParallelMatrix & locLocDistr,
                                     string const & rowFieldName,
                                     string const & colFieldName ) const
{
  localIndex rowFieldIndex, colFieldIndex;

  if( rowFieldName.length() > 0 )
  {
    // check if the row field name is already added
    GEOS_ERROR_IF( !keyInUse( rowFieldName ), "setSparsityPattern: requested field name must be already existing." );

    // get row field index
    rowFieldIndex = getFieldIndex( rowFieldName );
  }
  else
  {
    rowFieldIndex = -1;
  }

  if( colFieldName.length() > 0 )
  {
    // check if the col field name is already added
    GEOS_ERROR_IF( !keyInUse( colFieldName ), "setSparsityPattern: requested field name must be already existing." );

    // get col field index
    colFieldIndex = getFieldIndex( colFieldName );
  }
  else
  {
    colFieldIndex = -1;
  }

  if( rowFieldIndex * colFieldIndex < 0 )
  {
    GEOS_ERROR( "setSparsityPattern accepts both two field names and none, instead just one is provided." );
  }

  // Call the low level routine
  setSparsityPattern( locLocDistr, rowFieldIndex, colFieldIndex );
}


void DofManager::setSparsityPatternOneBlock( ParallelMatrix & locLocDistr,
                                             localIndex const rowFieldIndex,
                                             localIndex const colFieldIndex ) const
{
  GEOS_ASSERT( rowFieldIndex >= 0 );
  GEOS_ASSERT( colFieldIndex >= 0 );

  locLocDistr.createWithLocalSize( m_fields[rowFieldIndex].numLocalRows,
                                   m_fields[colFieldIndex].numLocalRows,
                                   1, MPI_COMM_GEOSX );

  if( colFieldIndex == rowFieldIndex )
  {
    // Diagonal block
    ParallelMatrix const * const connLocPattDistr = m_sparsityPattern( rowFieldIndex, rowFieldIndex ).first.get();

    connLocPattDistr->MatrixMatrixMultiply( true,
                                            *connLocPattDistr,
                                            false,
                                            locLocDistr );
  }
  else
  {
    // ExtraDiagonal (coupling) block
    if( m_connectivity[rowFieldIndex][colFieldIndex] != Connectivity::None )
    {
      ParallelMatrix const * CL1;
      ParallelMatrix const * CL2;

      if( m_sparsityPattern( rowFieldIndex, colFieldIndex ).first != nullptr )
      {
        CL1 = m_sparsityPattern( rowFieldIndex, colFieldIndex ).first.get();
        CL2 = m_sparsityPattern( rowFieldIndex, colFieldIndex ).second.get();
      }
      else
      {
        CL1 = m_sparsityPattern( colFieldIndex, rowFieldIndex ).second.get();
        CL2 = m_sparsityPattern( colFieldIndex, rowFieldIndex ).first.get();
      }

      CL1->MatrixMatrixMultiply( true, *CL2, false, locLocDistr );
    }
    else
    {
      locLocDistr.close(); // empty matrix, but still needs to be closed
    }
  }
}

// Create the sparsity pattern (location-location). Low level interface
void DofManager::setSparsityPattern( ParallelMatrix & matrix,
                                     localIndex const rowFieldIndex,
                                     localIndex const colFieldIndex ) const
{
  GEOS_ASSERT( rowFieldIndex < m_fields.size() );
  GEOS_ASSERT( colFieldIndex < m_fields.size() );
  GEOS_ERROR_IF( (rowFieldIndex >= 0 && colFieldIndex < 0) || (rowFieldIndex < 0 && colFieldIndex >= 0),
                 "setSparsityPattern accepts either two non-negative values (existing field indices) or "
                 "two negative values (entire Jacobian rows/columns), instead just one index is non-negative." );

  if( rowFieldIndex >= 0 ) // both nonnegative => single row/col field
  {
    setSparsityPatternOneBlock( matrix, rowFieldIndex, colFieldIndex );
  }
  else if( m_fields.empty() ) // both negative, no fields present
  {
    matrix.createWithLocalSize( 0, 0, 0, MPI_COMM_GEOSX );
    matrix.close();
  }
  else if( m_fields.size() == 1 ) // both negative, single field present
  {
    setSparsityPatternOneBlock( matrix, 0, 0 );
  }
  else // both negative, multiple fields
  {
    GEOS_ERROR_IF( !m_closed, "setSparsityPattern: DofManager needs to be closed first" );

    // Create the matrix
    globalIndex const sumLocalDofs = numLocalDofs();

    ParallelMatrix sparsity, colPerm, localPattern;

    matrix.createWithLocalSize( sumLocalDofs, sumLocalDofs, 1, MPI_COMM_GEOSX );
    sparsity.createWithLocalSize( sumLocalDofs, sumLocalDofs, 1, MPI_COMM_GEOSX );
    colPerm.createWithLocalSize( sumLocalDofs, sumLocalDofs, 1, MPI_COMM_GEOSX  );

    array1d<globalIndex> indices;
    array1d<real64> values;

    globalIndex fieldOffset = 0;

    // Loop over all fields
    for( localIndex iGlo = 0; iGlo < m_fields.size(); ++iGlo )
    {
      FieldDescription const & field = m_fields[iGlo];
      globalIndex const row_adjustment = field.fieldOffset - field.firstLocalRow;

      // Loop over all fields
      for( localIndex jGlo = 0; jGlo < m_fields.size(); ++jGlo )
      {
        // compute single coupling block pattern
        setSparsityPatternOneBlock( localPattern, iGlo, jGlo );

        // Assemble into global pattern (with indices adjusted for field offsets)
        for( globalIndex i = localPattern.ilower(); i < localPattern.iupper(); ++i )
        {
          localPattern.getRowCopy( i, indices, values );
          if( !indices.empty() )
          {
            for( globalIndex j = 0; j < indices.size(); ++j )
            {
              indices[j] += fieldOffset;
            }
            sparsity.insert( i + row_adjustment, indices, values );
          }
        }
      }

      for( globalIndex i = 0; i < field.numLocalRows; ++i )
      {
        colPerm.insert( fieldOffset + field.firstLocalRow + i, field.fieldOffset + i, 1.0 );
      }

      fieldOffset += field.numGlobalRows;
    }
    sparsity.close();
    colPerm.close();

    // Permute the columns to adjust for rank-based ordering
    sparsity.MatrixMatrixMultiply( false, colPerm, false, matrix );
  }
}

// Allocate a vector (location-location). High level interface
void DofManager::setVector( ParallelVector & vector,
                            string const & fieldName ) const
{
  localIndex fieldIndex;

  if( !fieldName.empty() )
  {
    GEOS_ERROR_IF( !keyInUse( fieldName ), "setVector: requested field name must be already existing." );
    fieldIndex = getFieldIndex( fieldName );
  }
  else
  {
    fieldIndex = -1;
  }

  // Call the low level routine
  setVector( vector, fieldIndex );
}

// Allocate a vector (location-location). Low level interface
void DofManager::setVector( ParallelVector & vector,
                            localIndex const fieldIndex ) const
{
  if( fieldIndex >= 0 )
  {
    vector.createWithLocalSize( m_fields[fieldIndex].numLocalRows, MPI_COMM_GEOSX );
  }
  else
  {
    GEOS_ERROR_IF( !m_closed, "setVector: DofManager needs to be closed first" );

    // Create the global vector
    globalIndex const sumLocalDofs = numLocalDofs();
    vector.createWithLocalSize( sumLocalDofs, MPI_COMM_GEOSX );
  }
  vector.close();
}

template< typename FIELD_OP, typename POLICY >
void DofManager::vectorToField( ParallelVector const & vector,
                                string const & srcFieldName,
                                real64 const scalingFactor,
                                ObjectManagerBase * const manager,
                                string const & dstFieldName,
                                localIndex const loCompIndex,
                                localIndex const hiCompIndex ) const
{
  GEOS_ERROR_IF( !keyInUse( srcFieldName ), "copyVectorToField: requested field does not exist: " << srcFieldName );

  FieldDescription const & fieldDesc = m_fields[ getFieldIndex( srcFieldName ) ];

  localIndex const loComp = loCompIndex;
  localIndex hiComp = hiCompIndex;
  if( hiComp < 0 )
  {
    hiComp = fieldDesc.numComponents;
  }
  GEOS_ASSERT( loComp >= 0 && hiComp <= fieldDesc.numComponents && loComp < hiComp );

  arrayView1d<globalIndex const> const & indexArray = manager->getReference< array1d<globalIndex> >( fieldDesc.key );
  arrayView1d<integer const> const & ghostRank = manager->GhostRank();

  globalIndex const rankOffset = offsetLocalDofs();

  real64 * localVector = nullptr;
  vector.extractLocalVector( &localVector );

  ViewWrapperBase * const vw = manager->getWrapperBase( dstFieldName );
  GEOS_ASSERT( vw != nullptr );
  std::type_index typeIndex = std::type_index( vw->get_typeid() );

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                  false,
                                  [&]( auto arrayInstance, auto dataTypeInstance )
  {
    using ArrayType = decltype(arrayInstance);
    ViewWrapper<ArrayType> & view = ViewWrapper<ArrayType>::cast( *vw );
    typename ViewWrapper<ArrayType>::ViewType const & field = view.referenceAsView();

    forall_in_range<POLICY>( 0, indexArray.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex const i )
    {
      if( ghostRank[i] < 0 )
      {
        localIndex const lid = indexArray[i] - rankOffset;
        GEOS_ASSERT( lid >= 0 ); // since vectors are partitioned same as the mesh
        for( localIndex c = loComp; c < hiComp; ++c )
        {
          FIELD_OP::template SpecifyFieldValue( field,
                                                i,
                                                integer_conversion< integer >( c - loComp ),
                                                scalingFactor * localVector[lid + c] );
        }
      }
    } );
  } );
}

// Copy values from DOFs to nodes
void DofManager::copyVectorToField( ParallelVector const & vector,
                                    string const & srcFieldName,
                                    real64 const scalingFactor,
                                    ObjectManagerBase * const manager,
                                    string const & dstFieldName,
                                    localIndex const loCompIndex,
                                    localIndex const hiCompIndex ) const
{
  vectorToField< FieldSpecificationEqual, parallelHostPolicy >( vector,
                                                                srcFieldName,
                                                                scalingFactor,
                                                                manager,
                                                                dstFieldName,
                                                                loCompIndex,
                                                                hiCompIndex );
}

// Copy values from DOFs to nodes
void DofManager::addVectorToField( ParallelVector const & vector,
                                   string const & srcFieldName,
                                   real64 const scalingFactor,
                                   ObjectManagerBase * const manager,
                                   string const & dstFieldName,
                                   localIndex const loCompIndex,
                                   localIndex const hiCompIndex ) const
{
  vectorToField< FieldSpecificationAdd, parallelHostPolicy >( vector,
                                                              srcFieldName,
                                                              scalingFactor,
                                                              manager,
                                                              dstFieldName,
                                                              loCompIndex,
                                                              hiCompIndex );
}

template< typename FIELD_OP, typename POLICY >
void DofManager::fieldToVector( ObjectManagerBase const * const manager,
                                string const & srcFieldName,
                                real64 const scalingFactor,
                                ParallelVector & vector,
                                string const & dstFieldName,
                                localIndex const loCompIndex,
                                localIndex const hiCompIndex ) const
{
  GEOS_ERROR_IF( !keyInUse( dstFieldName ), "copyVectorToField: requested field does not exist: " << dstFieldName );

  FieldDescription const & fieldDesc = m_fields[ getFieldIndex( dstFieldName ) ];

  localIndex const loComp = loCompIndex;
  localIndex hiComp = hiCompIndex;
  if( hiComp < 0 )
  {
    hiComp = fieldDesc.numComponents;
  }
  GEOS_ASSERT( loComp >= 0 && hiComp <= fieldDesc.numComponents && loComp < hiComp );

  arrayView1d<globalIndex const> const & indexArray = manager->getReference< array1d<globalIndex> >( fieldDesc.key );
  arrayView1d<integer const> const & ghostRank = manager->GhostRank();

  globalIndex const rankOffset = offsetLocalDofs();

  real64 * localVector = nullptr;
  vector.extractLocalVector( &localVector );

  ViewWrapperBase const * const vw = manager->getWrapperBase( srcFieldName );
  GEOS_ASSERT( vw != nullptr );
  std::type_index typeIndex = std::type_index( vw->get_typeid() );

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                  false,
                                  [&]( auto arrayInstance, auto dataTypeInstance )
  {
    using ArrayType = decltype(arrayInstance);
    ViewWrapper<ArrayType> const & view = ViewWrapper<ArrayType>::cast( *vw );
    typename ViewWrapper<ArrayType>::ViewTypeConst const & field = view.referenceAsView();

    forall_in_range<POLICY>( 0, indexArray.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex const i )
    {
      if( ghostRank[i] < 0 )
      {
        localIndex const lid = indexArray[i] - rankOffset;
        GEOS_ASSERT( lid >= 0 ); // since vectors are partitioned same as the mesh
        for( localIndex c = loComp; c < hiComp; ++c )
        {
          FIELD_OP::template ReadFieldValue( field,
                                             i,
                                             integer_conversion< int >( c - loComp ),
                                             localVector[lid + c] );
        }
      }
    } );
  } );
}

// Copy values from nodes to DOFs
void DofManager::copyFieldToVector( ObjectManagerBase const * const manager,
                                    string const & srcFieldName,
                                    real64 const scalingFactor,
                                    ParallelVector & vector,
                                    string const & dstFieldName,
                                    localIndex const loCompIndex,
                                    localIndex const hiCompIndex ) const
{
  fieldToVector< FieldSpecificationEqual, parallelHostPolicy >( manager,
                                                                srcFieldName,
                                                                scalingFactor,
                                                                vector,
                                                                dstFieldName,
                                                                loCompIndex,
                                                                hiCompIndex );
}

// Copy values from nodes to DOFs
void DofManager::addFieldToVector( ObjectManagerBase const * const manager,
                                   string const & srcFieldName,
                                   real64 const scalingFactor,
                                   ParallelVector & vector,
                                   string const & dstFieldName,
                                   localIndex const loCompIndex,
                                   localIndex const hiCompIndex ) const
{
  fieldToVector< FieldSpecificationAdd, parallelHostPolicy >( manager,
                                                              srcFieldName,
                                                              scalingFactor,
                                                              vector,
                                                              dstFieldName,
                                                              loCompIndex,
                                                              hiCompIndex );
}

// Just an interface to allow only three parameters
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connectivity const connectivity )
{
  addCoupling( rowFieldName, colFieldName, connectivity, string_array(), true );
}

// Just another interface to allow four parameters (no symmetry)
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connectivity const connectivity,
                              string_array const & regions = string_array() )
{
  addCoupling( rowFieldName, colFieldName, connectivity, regions, true );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connectivity const connectivity,
                              bool const symmetric )
{
  addCoupling( rowFieldName, colFieldName, connectivity, string_array(), symmetric );
}

// The real function, allowing the creation of coupling blocks
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connectivity const connectivity,
                              string_array const & regions,
                              bool const symmetric )
{
  GEOS_ERROR_IF( m_closed, "addCoupling: cannot add coupling after DofManager has been closed." );
  GEOS_ERROR_IF( !keyInUse( rowFieldName ), "addCoupling: field does not exist: " << rowFieldName );
  GEOS_ERROR_IF( !keyInUse( colFieldName ), "addCoupling: field does not exist: " << colFieldName );

  // get row/col field index
  localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
  localIndex const colFieldIndex = getFieldIndex( colFieldName );

  // Check if already defined
  if( m_connectivity[rowFieldIndex][colFieldIndex] != Connectivity::None )
  {
    GEOS_ERROR( "addCoupling: coupling already defined with another connectivity" );
    return;
  }

  // get row/col field regionName
  array1d<string> const & rowRegions = m_fields[rowFieldIndex].regionNames;
  array1d<string> const & colRegions = m_fields[colFieldIndex].regionNames;

  string_array regionList( regions );
  if( regionList.empty() )
  {
    // Resize regionsList
    regionList.resize( std::min( rowRegions.size(), colRegions.size() ) );

    // Find common regions
    string_array::iterator it = std::set_intersection( rowRegions.begin(), rowRegions.end(),
                                                       colRegions.begin(), colRegions.end(),
                                                       regionList.begin() );
    regionList.resize( std::distance( regionList.begin(), it ) );
  }
  else
  {
    std::sort( regionList.begin(), regionList.end() );
    regionList.resize( std::distance( regionList.begin(), std::unique( regionList.begin(), regionList.end() ) ) );

    // Check that both fields are defined on all regions in the list
    for( string const & regionName : regionList )
    {
      GEOS_ERROR_IF( std::find( rowRegions.begin(), rowRegions.end(), regionName ) == rowRegions.end(),
                     "addCoupling: region " << regionName << " does not belong to the domain of field " << rowFieldName );
      GEOS_ERROR_IF( std::find( colRegions.begin(), colRegions.end(), regionName ) == colRegions.end(),
                     "addCoupling: region " << regionName << " does not belong to the domain of field " << colFieldName );
    }
  }

  // save field's connectivity type (rowField to colField)
  m_connectivity[rowFieldIndex][colFieldIndex] = connectivity;

  // Set connectivity with active symmetry flag
  if( symmetric )
  {
    m_connectivity[colFieldIndex][rowFieldIndex] = connectivity;
  }

  // Compute the CL matrices for row and col fields
  if( connectivity != Connectivity::None )
  {
    std::unique_ptr<ParallelMatrix> & rowConnLocPattern = m_sparsityPattern( rowFieldIndex, colFieldIndex ).first;
    rowConnLocPattern = std::make_unique<ParallelMatrix>();
    makeConnLocPattern( m_fields[rowFieldIndex], connectivity, regions, *rowConnLocPattern );

    std::unique_ptr<ParallelMatrix> & colConnLocPattern = m_sparsityPattern( rowFieldIndex, colFieldIndex ).second;
    colConnLocPattern = std::make_unique<ParallelMatrix>();
    makeConnLocPattern( m_fields[colFieldIndex], connectivity, regions, *colConnLocPattern );
  }
}

namespace
{

/**
 * Definition for entries of sparse matrix in COO format
 */
template< typename INDEX, typename VALUE >
struct MatEntry
{
  MatEntry()
    : row( -1 ),
      col( -1 ),
      val( 0 )
  {}

  MatEntry( INDEX r,
            INDEX c,
            VALUE v )
    : row( r ),
      col( c ),
      val( v )
  {}

  bool operator==( MatEntry<INDEX, VALUE> const & rhs ) const
  {
    return row == rhs.row && col == rhs.col && val == rhs.val;
  }

  bool operator<( MatEntry<INDEX, VALUE> const & rhs ) const
  {
    return row < rhs.row || (!(rhs.row < row) && col < rhs.col);
  }

  INDEX row;
  INDEX col;
  VALUE val;
};

/**
 * A temporary representation of a local portion of conn-loc sparsity
 */
template< typename INDEX, typename T >
struct LocalSparsityPattern
{
  array1d<INDEX> rowNumbers; //!< row numbers, size numLocalRows
  ArrayOfArrays<INDEX> colIndices; //!< packed column indices, size numLocalNonZeros
  ArrayOfArrays<T> nnzEntries; //!< packed values (of type localIndex), size numLocalNonZeros

  void appendRow( INDEX row,
                  array1d<INDEX> const & cols,
                  array1d<T> const & vals )
  {
    rowNumbers.push_back( row );
    colIndices.appendArray( cols.data(), cols.size() );
    nnzEntries.appendArray( vals.data(), vals.size() );
  }
};

/**
 * Convert a list of COO entries in CSR format.
 * Assumes pairs is sorted by first element of each entry
 */
template< typename INDEX, typename T, typename U >
void vectorOfPairsToCSR( std::vector<MatEntry<INDEX, T>> const & entries,
                         LocalSparsityPattern<INDEX, U> & pattern )
{
  pattern.rowNumbers.clear();
  pattern.colIndices.resize( 0 );
  pattern.nnzEntries.resize( 0 );

  if( entries.empty() )
  {
    return;
  }

  localIndex const nNz = entries.size();
  localIndex nRows = 0;
  globalIndex currRow = -1;

  // count number of distinct rows
  for( MatEntry<INDEX, T> const & entry : entries )
  {
    if( entry.row != currRow )
    {
      currRow = entry.row;
      ++nRows;
    }
  }

  pattern.rowNumbers.reserve( nRows );
  pattern.colIndices.reserve( nRows );
  pattern.nnzEntries.reserve( nRows );
  pattern.colIndices.reserveValues( nNz );
  pattern.nnzEntries.reserveValues( nNz );

  array1d<INDEX> cols;
  array1d<U> vals;
  currRow = entries.front().row;

  for( MatEntry<INDEX, T> const & entry : entries )
  {
    if( entry.row != currRow )
    {
      pattern.appendRow( currRow, cols, vals );
      cols.clear();
      vals.clear();
      currRow = entry.row;
    }
    cols.push_back( entry.col );
    vals.push_back( static_cast<U>( entry.val ) );
  }
  pattern.appendRow( currRow, cols, vals );
}

}

// Compute the sparsity pattern of the matrix connectivity - location for the specified
// field (diagonal entry in the m_connectivity collection)
void DofManager::makeConnLocPattern( FieldDescription const & fieldDesc,
                                     Connectivity const connectivity,
                                     array1d <string> const & regions,
                                     ParallelMatrix & connLocPattern )
{
  localIndex const NC = fieldDesc.numComponents;

  if( connectivity == Connectivity::None )
  {
    // Case of no connectivity (self connection)
    Connectivity selfConnectivity = static_cast<Connectivity>( fieldDesc.location );
    makeConnLocPattern( fieldDesc, selfConnectivity, regions, connLocPattern );
    return;
  }

  // Create a name decoration with all region names
  string suffix;
  for( string const & regionName : regions )
  {
    suffix.append( "_" + regionName );
  }

  // array to store the matrix in COO format
  std::vector< MatEntry< globalIndex, localIndex > > entries;
  //entries.reserve( fieldConn.numLocalNodes );
  localIndex numLocalConns;

  LocationSwitch( fieldDesc.location,
                  static_cast<Location>( connectivity ),
                  [&] ( auto const loc, auto const conn )
  {
    Location constexpr LOC  = decltype(loc)::value;
    Location constexpr CONN = decltype(conn)::value;

    // Create a unique global indexing for connectors in active regions
    FieldDescription fieldConn;
    fieldConn.location = CONN;
    fieldConn.key = "connIndex" + suffix;
    fieldConn.regionNames = regions;

    createIndexArrayImpl<CONN>( m_domain, m_mesh, fieldConn );
    numLocalConns = fieldConn.numLocalNodes;

    using ArrayHelperLoc = IndexArrayHelper<globalIndex const, LOC>;
    using ArrayHelperCon = IndexArrayHelper<globalIndex, CONN>;

    typename ArrayHelperLoc::Accessor indexArrayLoc = ArrayHelperLoc::get( m_mesh, fieldDesc );
    typename ArrayHelperCon::Accessor indexArrayCon = ArrayHelperCon::get( m_mesh, fieldConn );

    // XXX: special treatment for TPFA-style connectivity for fractures
    // We need this because edges are not added as face connectors when inserting a fracture element
    if( LOC == Location::Elem && CONN == Location::Face )
    {
      Location constexpr EDGE = Location::Edge; // for brevity
      using ArrayHelperEdge = IndexArrayHelper<globalIndex, EDGE>;

      FieldDescription fieldEdge;
      fieldEdge.location = EDGE;
      fieldEdge.key = "connectorIndices" + suffix;
      fieldEdge.regionNames = regions;

      createIndexArrayImpl<EDGE, FaceElementSubRegion>( m_domain, m_mesh, fieldEdge );
      typename ArrayHelperEdge::Accessor & indexArrayEdge = ArrayHelperEdge::get( m_mesh, fieldEdge );

      // backup old values for face connectors
      localIndex  const numLocalFaces  = fieldConn.numLocalNodes;
      globalIndex const firstLocalFace = fieldConn.firstLocalRow;

      // adjust face connector indexing to account for added edge connectors on every rank
      // the reason is we need contiguous numbering of connectors within a rank in order to create CL matrix
      numLocalConns += fieldEdge.numLocalNodes;
      localIndex const firstLocalConn = CommunicationTools::PrefixSum<globalIndex>( numLocalConns ).first;
      localIndex const adjustment = firstLocalConn - firstLocalFace;

      forMeshLocation<CONN>( m_mesh, regions,
                             [&]( auto const & faceIdx,
                                  auto const,
                                  localIndex const )
      {
        GEOS_ASSERT( ArrayHelperCon::value( indexArrayCon, faceIdx ) >= 0 );
        ArrayHelperCon::reference( indexArrayCon, faceIdx ) += adjustment;
      } );

      // also adjust edge connector indexing in a similar way
      forMeshLocation<EDGE, FaceElementSubRegion>( m_mesh, regions,
                                                   [&]( auto const & edgeIdx,
                                                        auto const,
                                                        localIndex const )
      {
        GEOS_ASSERT( ArrayHelperEdge::value( indexArrayEdge, edgeIdx ) >= 0 );
        ArrayHelperEdge::reference( indexArrayEdge, edgeIdx ) += adjustment + numLocalFaces;
      } );

      // sync index arrays after adjustment
      std::map<string, string_array> fieldNames;
      fieldNames[ MeshHelper<CONN>::syncObjName ].push_back( fieldConn.key );
      fieldNames[ MeshHelper<EDGE>::syncObjName ].push_back( fieldEdge.key );

      CommunicationTools::
      SynchronizeFields( fieldNames, m_mesh,
                         m_domain->getReference<array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

      // add links between edge connectors and fracture elements
      forMeshLocation<LOC, EDGE, FaceElementSubRegion>( m_mesh, regions,
                                                        [&]( auto const & elemIdx,
                                                             auto const & edgeIdx,
                                                             localIndex const k )
      {
        globalIndex const indexElem  = ArrayHelperLoc::value( indexArrayLoc, elemIdx );
        globalIndex const indexEdge = ArrayHelperEdge::value( indexArrayEdge, edgeIdx );
        GEOS_ASSERT( indexEdge >= 0 );
        GEOS_ASSERT( indexElem >= 0 );

        for( localIndex c = 0; c < NC; ++c )
        {
          entries.emplace_back( indexEdge, indexElem + c, NC * k + c + 1 );
        }
      } );

      // here we also have to loop through fracture elements, because face-to-element maps
      // are not updated after fracture creation (this may need to be fixed?), so the general
      // connector-based loop below misses these links. This involves assembling into remote
      // rows of sparsity pattern, which is not great, but will have to do for now.
      forMeshLocation<LOC, CONN, FaceElementSubRegion>( m_mesh, regions,
                                                        [&]( auto const & elemIdx,
                                                             auto const & faceIdx,
                                                             localIndex const k )
      {
        globalIndex const indexElem = ArrayHelperLoc::value( indexArrayLoc, elemIdx );
        globalIndex const indexConn = ArrayHelperCon::value( indexArrayCon, faceIdx );
        GEOS_ASSERT( indexElem >= 0 );
        GEOS_ASSERT( indexConn >= 0 );

        for( localIndex c = 0; c < NC; ++c )
        {
          entries.emplace_back( indexConn, indexElem + c, NC * k + c + 1 );
        }
      } );
    } //XXX: end special treatment for TPFA-style connectivity for fractures

    // loop over locally owned connectors and adjacent locations
    forMeshLocation<CONN, LOC>( m_mesh, regions,
                                [&]( auto const & connIdx,
                                     auto const & locIdx,
                                     localIndex const k )
    {
      globalIndex const indexConn = ArrayHelperCon::value( indexArrayCon, connIdx );
      globalIndex const indexLoc  = ArrayHelperLoc::value( indexArrayLoc, locIdx );
      GEOS_ASSERT( indexConn >= 0 );

      if( indexLoc >= 0 )
      {
        for( localIndex c = 0; c < NC; ++c )
        {
          entries.emplace_back( indexConn, indexLoc + c, NC * k + c + 1 );
        }
      }
    } );

    removeIndexArray( fieldConn );
  } );


  // Sort the entries lexicographically and remove duplicates (if any)
  sort( entries.begin(), entries.end() );
  entries.resize( std::unique( entries.begin(), entries.end() ) - entries.begin() );

  LocalSparsityPattern<globalIndex, real64> localPattern;
  vectorOfPairsToCSR( entries, localPattern );

  localIndex const entriesPerRow = (numLocalConns > 0 ) ? (entries.size() / numLocalConns ) : 0;
  connLocPattern.createWithLocalSize( numLocalConns, fieldDesc.numLocalRows, entriesPerRow, MPI_COMM_GEOSX );

  for( localIndex i = 0; i < localPattern.rowNumbers.size(); ++i )
  {
    connLocPattern.insert( localPattern.rowNumbers[i],
                           localPattern.colIndices[i],
                           localPattern.nnzEntries[i],
                           localPattern.colIndices.sizeOfArray( i ) );
  }
  connLocPattern.close();
}

void DofManager::close()
{
  if( m_fields.size() > 1 )
  {
    // count total number of DoFs on lower ranks
    globalIndex dofOffset = 0;
    for( FieldDescription const & field : m_fields )
    {
      dofOffset += field.firstLocalRow;
    }

    // update field offsets to account for renumbering
    for( FieldDescription & field : m_fields )
    {
      field.fieldOffset = dofOffset;
      dofOffset += field.numLocalRows;
    }

    std::map< string, string_array > fieldNames;

    // adjust index arrays for owned locations
    for( FieldDescription const & field : m_fields )
    {
      globalIndex const adjustment = field.fieldOffset - field.firstLocalRow;

      LocationSwitch( field.location, [&]( auto const loc )
      {
        Location constexpr LOC = decltype(loc)::value;
        using ArrayHelper = IndexArrayHelper< globalIndex, LOC >;

        typename ArrayHelper::Accessor indexArray = ArrayHelper::get( m_mesh, field );

        forMeshLocation< LOC >( m_mesh, field.regionNames,
                                [&]( auto const locIdx,
                                     auto const,
                                     localIndex const )
        {
          ArrayHelper::reference( indexArray, locIdx ) += adjustment;
        } );

        fieldNames[MeshHelper< LOC >::syncObjName].push_back( field.key );
      } );
    }

    // synchronize index arrays for all fields across ranks
    CommunicationTools::
    SynchronizeFields( fieldNames, m_mesh,
                       m_domain->getReference< array1d< NeighborCommunicator > >( m_domain->viewKeys.neighbors ) );
  }

  m_closed = true;
}

void DofManager::printConnectivityLocationPattern( string const & fieldName, string const & fileName ) const
{
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( fieldName ),
                 "printConnectivityLocationPattern: requested field name must be already existing." );

  string const name = !fileName.empty() ? fileName : "pattern_" + fieldName + ".mtx" ;
  localIndex const fieldIndex = getFieldIndex( fieldName );
  m_sparsityPattern( fieldIndex, fieldIndex ).first->printParallelMatrix( name );
}

// Print the coupling table on screen
void DofManager::printConnectivityMatrix( std::ostream & os ) const
{
  if( CommunicationTools::MPI_Rank(MPI_COMM_GEOSX) == 0 )
  {
    localIndex numFields = m_fields.size();

    os << std::endl;
    for( localIndex i = 0 ; i < numFields ; ++i )
    {
      for( localIndex j = 0 ; j < numFields ; ++j )
      {
        switch( m_connectivity[i][j] )
        {
          case Connectivity::Elem:
            os << " E ";
            break;
          case Connectivity::Face:
            os << " F ";
            break;
          case Connectivity::Edge:
            os << " G ";
            break;
          case Connectivity::Node:
            os << " N ";
            break;
          case Connectivity::USER_DEFINED:
            os << " U ";
            break;
          case Connectivity::None:
            os << "   ";
            break;
        }
        if( j < numFields - 1 )
        {
          os << "|";
        }
      }
      os << std::endl;
      if( i < numFields - 1 )
      {
        for( localIndex j = 0 ; j < numFields - 1 ; ++j )
        {
          os << "---|";
        }
        os << "---";
      }
      os << std::endl;
    }
  }
}

}
