/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CohesiveZoneManager.hpp
 */

#ifndef GEOS_MESH_COHESIVEZONEMANAGER_HPP
#define GEOS_MESH_COHESIVEZONEMANAGER_HPP

#include "mesh/CohesiveZoneRegionBase.hpp"
#include "mesh/CohesiveZoneRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"

#include "mesh/NodeManager.hpp"
#include "mesh/ParticleManager.hpp"

namespace geos
{

class MeshManager;
class DomainPartition;

/**
 * @class CohesiveZoneManager
 * @brief The CohesiveZoneManager class provides an interface to ObjectManagerBase in order to manage Cohesive Zone Regions
 * data
 */
class CohesiveZoneManager : public ObjectManagerBase
{
public:

    /**
    * @brief Group key associated with czRegionsGroup.
    */
    struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
    {
        /// @return cohesive zone regions group string key
        static constexpr auto cohesiveZoneRegionsGroup() { return "cohesiveZoneRegionsGroup"; }
    };

    /**
    * @brief The function is to return the name of the CohesiveZoneManager in the object catalog
    * @return string that contains the catalog name used to register/lookup this class in  the object catalog
    */
    static string catalogName()
    { return "CohesiveZoneManager"; }

    /**
    * @brief Virtual access to catalogName()
    * @return string that contains the catalog name used to register/lookup this class in the object catalog
    */
    virtual string getCatalogName() const override final
    { return CohesiveZoneManager::catalogName(); }

    /**
    * @brief Constructor.
    * @param [in] name the name of this ObjectManager
    * @param [in] parent the parent Group
    */
    CohesiveZoneManager( string const & name, Group * const parent );

    /**
    * @brief Destructor
    */
    virtual ~CohesiveZoneManager() override;

    /**
    * @brief Determines if a CohesiveZoneRegion with the input name exists.
    * @tparam T The type of CohesiveZoneRegion. May be a specific derived type of CohesiveZoneRegion.
    * @param name The name/key of the CohesiveZoneRegion
    * @return true if the region exists, false if not.
    */
    template< typename T=CohesiveZoneRegionBase >
    bool hasRegion( string const & name ) const
    {
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).hasGroup< T >( name );
    }

    /**
    * @brief Get the number of regions.
    * @return number of regions
    */
    localIndex numRegions() const
    {
        return this->getRegions().size();
    }

    /**
    * @brief Get the number of nodes in all cohesive zone regions of type T.
    * @return number of nodes
    */
    template< typename T = CohesiveZoneRegionBase >
    localIndex getNumberOfNodes() const
    {
        localIndex numNode = 0;
        this->forCohesiveZoneRegions< T >( [&]( CohesiveZoneRegionBase const & cohesiveZoneRegion )
        {
            numNode += cohesiveZoneRegion.size();
        } );
        return numNode;
    }

    /**
    * @brief Get a cohesive zone region.
    * @param key The key of cohesive zone region, either name or number.
    * @return Reference to const T.
    */
    template< typename T=CohesiveZoneRegionBase, typename KEY_TYPE=void >
    T const & getRegion( KEY_TYPE const & key ) const
    {
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).getGroup< T >( key );
    }

    /**
    * @brief Get a cohesive zone region.
    * @param key The key of the cohesive zone region, either name or number.
    * @return Reference to T.
    */
    template< typename T=CohesiveZoneRegionBase, typename KEY_TYPE=void >
    T & getRegion( KEY_TYPE const & key )
    {
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).getGroup< T >( key );
    }

    /**
    * @brief Get a collection of cohesive zone regions
    * @return reference to immutable subGroupMap
    */
    subGroupMap const & getRegions() const
    {
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).getSubGroups();
    }

    /**
    * @brief Get a collection of cohesive zone regions.
    * @return reference to mutable subGroupMap
    */
    subGroupMap & getRegions()
    {
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).getSubGroups();
    }

    /**
    * @brief Create a new CohesiveZoneRegion object as a child of this group.
    * @param childKey catalog key of the new CohesiveZoneRegion derived type to create
    * @param childName name of the new CohesiveZone object
    * @return pointer to the created CohesiveZone object
    */
    virtual Group * createChild( string const & childKey, string const & childName ) override;
    
    /**
    * @brief Add a cohesive zone region.
    * @param regionName The name of the cohesive zone region, either name or number, to create.
    * @return Reference to the new cohesive zone region.
    */
    template< typename T=CohesiveZoneRegionBase >
    T & addRegion( string const & regionName )
    {
        (void) createChild( T::catalogName(), regionName );
        return this->getGroup( groupKeyStruct::cohesiveZoneRegionsGroup() ).getGroup< T >( regionName );
    }

    /**
    * @brief Remove a cohesive zone region.
    * @param regionName The name of the cohesive zone region, either name or number, to remove.
    */
    template< typename T=CohesiveZoneRegionBase >
    void removeRegion( string const & regionName )
    {
        GEOS_ERROR_IF( !hasRegion< T >( regionName ), "Attempted to remove nonexistant cohesive zone region : " << regionName << "!" );
        Group & cohesiveZoneRegions = this->getGroup( CohesiveZoneManager::groupKeyStruct::cohesiveZoneRegionsGroup() );
        cohesiveZoneRegions.deregisterGroup( regionName );
    }

    /**
    * @brief Expand any catalogs in the data structure
    */
    virtual void expandObjectCatalogs() override;

    /**
    * @brief Inform the schema generator of any deviations between the xml and GEOS data structures.
    * @param schemaRoot        XML node corresponding to the root
    * @param schemaParent      XML node for the parent node
    * @param documentationType type of XML schema generated
    */
    virtual void setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                        xmlWrapper::xmlNode schemaParent,
                                        integer documentationType ) override;

    using Group::resize;

    /**
    * @brief Set the number of nodes for a set of cohesive zone regions.
    * @param numNodes list of the new node numbers
    * @param regionNames list of the node region names
    */
    void resize( integer_array const & numNodes,
                string_array const & regionNames );

    /**
    * @brief This function is used to launch kernel function over the cohesive zone regions of all types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename LAMBDA >
    void forCohesiveZoneRegions( LAMBDA && lambda )
    {
        forCohesiveZoneRegions< CohesiveZoneRegion >( std::forward< LAMBDA >( lambda ) );
    }

    /**
    * @brief This const function is used to launch kernel function over the cohesive zone regions of all types.
    * types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename LAMBDA >
    void forCohesiveZoneRegions( LAMBDA && lambda ) const
    {
        forCohesiveZoneRegions< CohesiveZoneRegion >( std::forward< LAMBDA >( lambda ) );
    }

    /**
    * @brief This function is used to launch kernel function over the cohesive zone bregions of the specified region
    * types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
    void forCohesiveZoneRegions( LAMBDA && lambda )
    {
        forCohesiveZoneRegionsComplete< REGIONTYPE, REGIONTYPES... >(
        [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                     auto & region )
        {
            lambda( region );
        }
        );
    }

    /**
    * @brief This const function is used to launch kernel function over the cohesive zone bregions of the specified region
    * types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
    void forCohesiveZoneRegions( LAMBDA && lambda ) const
    {
        forCohesiveZoneRegionsComplete< REGIONTYPE, REGIONTYPES... >(
        [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                     auto const & region )
        {
            lambda( region );
        } );
    }

    /**
    * @brief This function is used to launch kernel function over the cohesive zone regions of all region types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename LAMBDA >
    void forCohesiveZoneRegionsComplete( LAMBDA && lambda )
    {
        forCohesiveZoneRegionsComplete< CohesiveZoneRegion >( std::forward< LAMBDA >( lambda ) );
    }

    /**
    * @brief This const function is used to launch kernel function over the cohesive zone regions of all region types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename LAMBDA >
    void forCohesiveZoneRegionsComplete( LAMBDA && lambda ) const
    {
        forCohesiveZoneRegionsComplete< CohesiveZoneRegion >( std::forward< LAMBDA >( lambda ) );
    }

    /**
    * @brief This function is used to launch kernel function over all the cohesive zone regions that can be casted to one of
    * the specified region types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
    void forCohesiveZoneRegionsComplete( LAMBDA && lambda )
    {
        for( localIndex er=0; er<this->numRegions(); ++er )
        {
            CohesiveZoneRegionBase & region = this->getRegion( er );
            
            Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( region, [&]( auto & castedRegion )
            {
                lambda( er, castedRegion );
            } );
        }
    }

    /**
    * @brief This const function is used to launch kernel function over all the cohesive zone regions that can be casted to one of
    * the specified region types.
    * @tparam LAMBDA type of the user-provided function
    * @param lambda kernel function
    */
    template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
    void forCohesiveZoneRegionsComplete( LAMBDA && lambda ) const
    {
        for( localIndex er=0; er<this->numRegions(); ++er )
        {
            CohesiveZoneRegionBase const & region = this->getRegion( er );

            Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( region, [&]( auto const & castedRegion )
            {
                lambda( er, castedRegion );
            } );
        }
    }

private:

  /**
   * @brief Copy constructor.
   */
  CohesiveZoneManager( const CohesiveZoneManager & );

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  CohesiveZoneManager & operator=( const CohesiveZoneManager & );
};


} // End geos namespace
#endif /* GEOS_MESH_COHESIVEZONEMANAGER_HPP */
