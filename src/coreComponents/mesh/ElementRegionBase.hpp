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

#ifndef GEOSX_MESH_ELEMENTREGIONBASE_HPP
#define GEOSX_MESH_ELEMENTREGIONBASE_HPP

#include "CellElementSubRegion.hpp"
#include "FaceElementSubRegion.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "EmbeddedSurfaceSubRegion.hpp"

namespace geosx
{

class StableTimeStep;

class FaceManager;

/**
 * @class ElementRegionBase
 * @brief The ElementRegionBase is the base class to manage the data stored at the element level.
 *
 * The ElementRegion base is the base class for classes such as CellElementRegion, FaceElementRegion, 
 * WellElementRegion, and EmbeddedSurfaceRegion.
 */
class ElementRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

//  static const string CatalogName()
//  { return "ElementRegionBase"; }
//
//  virtual const string getCatalogName() const override
//  { return ElementRegionBase::CatalogName(); }


  ///@}

/**
 * @brief deleted default constructor
 */
  ElementRegionBase() = delete;

/**
 * @brief Main constructor.
 * @param name the name of the element region
 * @param parent the pointer to the parent group
 */
  ElementRegionBase( string const & name, Group * const parent );


/**
 * @brief Copy constructor.
 * @param init the element region to be copied
 */
  ElementRegionBase( const ElementRegionBase & init );

/**
 * @brief Default destructor.
 */
  virtual ~ElementRegionBase() override;

/**
 * @brief Generate mesh
 * @param cellBlocks cell blocks where the mesh is generated
 */
  virtual void GenerateMesh( Group * const cellBlocks )
  {
    GEOSX_UNUSED_VAR( cellBlocks );
    GEOSX_ERROR( "ElementRegionBase::GenerateMesh() should be overriden if called." );
  }

//  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const NodeManager )
//  {
//    GEOSX_ERROR( "ElementRegionBase::GenerateAggregates() should be overriden if called.");
//  }

/**
 * @copydoc GetSubRegions() const
 */
  subGroupMap & GetSubRegions()
  {
    return GetGroup( viewKeyStruct::elementSubRegions )->GetSubGroups();
  }

/**
 * @brief Get a collection of the subregions.
 * @return a collection of the subregions
 */
  subGroupMap const & GetSubRegions() const
  {
    return GetGroup( viewKeyStruct::elementSubRegions )->GetSubGroups();
  }

/**
 * @brief Get a pointer to a subregion by specifying its name.
 * @tparam SUBREGIONTYPE the type of the subregion
 * @param regionName the name of the subregion
 * @return a pointer to the subregion
 * @note
 */
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup( viewKeyStruct::elementSubRegions )->GetGroup< SUBREGIONTYPE >( regionName );
  }

/**
 * @copydoc GetSubRegion( string const & regionName ) const
 */
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE * GetSubRegion( string const & regionName )
  {
    return this->GetGroup( viewKeyStruct::elementSubRegions )->GetGroup< SUBREGIONTYPE >( regionName );
  }

/**
 * @brief Get a pointer to a subregion by specifying its index.
 * @tparam SUBREGIONTYPE the type of the subregion
 * @param index the index of the subregion
 * @return a pointer to the subregion with the specified subregion type
 */
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup( viewKeyStruct::elementSubRegions )->GetGroup< SUBREGIONTYPE >( index );
  }

/**
 * @copydoc GetSubRegion( localIndex const & index ) const
 */
  template< typename SUBREGIONTYPE=ElementSubRegionBase >
  SUBREGIONTYPE * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup( viewKeyStruct::elementSubRegions )->GetGroup< SUBREGIONTYPE >( index );
  }

/**
 * @brief Get the number of subregions in the region.
 * @return the number of subregions  in the region
 */
  localIndex numSubRegions() const
  {
    return this->GetGroup( viewKeyStruct::elementSubRegions )->GetSubGroups().size();
  }

/**
 * @brief Get the number of subregions  in the region
 *        for specific subregion types provided as template arguments.
 * @tparam SUBREGIONTYPES
 * @return the number of elements contained in the element region
 */
  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES >
  localIndex getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forElementSubRegions< SUBREGIONTYPE, SUBREGIONTYPES... >( [&]( Group const * cellBlock ) -> void
    {
      numElem += cellBlock->size();
    } );
    return numElem;
  }

/**
 * @brief Apply a lambda to all subregions.
 * @param lambda the lambda to be applied
 */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegions( LAMBDA && lambda ) const
 */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions with the specific subregion types
 *        listed in the template.
 * @param lambda the lambda to be applied
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    Group const * const elementSubRegions = this->GetGroup( viewKeyStruct::elementSubRegions );
    elementSubRegions->forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegions( LAMBDA && lambda ) const
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    Group * const elementSubRegions = this->GetGroup( viewKeyStruct::elementSubRegions );
    elementSubRegions->forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions, loop using subregion indices.
 * @param lambda the lambda to be applied
 */
  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    forElementSubRegionsIndex< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegionsIndex( LAMBDA && lambda ) const
 */
  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    forElementSubRegionsIndex< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions with the specific subregion types
 *        listed in the template, loop using subregion indices.
 * @param lambda the lambda to be applied
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ElementSubRegionBase const & subRegion = *this->GetSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }

/**
 * @copydoc forElementSubRegionsIndex( LAMBDA && lambda ) const
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ElementSubRegionBase & subRegion = *this->GetSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }


/**
 * @brief Struct to serve as a container for variable strings and keys.
 * @struct viewKeyStruct 
*/
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// the key for the material list
    static constexpr auto materialListString = "materialList";
    /// the key for the element subregions
    static constexpr auto elementSubRegions = "elementSubRegions";
  };

/**
 * @copydoc getMaterialList() const
 */
  string_array & getMaterialList() {return m_materialList;}

/**
 * @brief Get the material list in the element region.
 * @return the material list
 */
  string_array const & getMaterialList() const {return m_materialList;}

/**
 * @brief Get the name of the constiutive in the element region.
 * @return the string_array of the constitutive names
 */
  template< typename CONSTITUTIVE_TYPE >
  string_array getConstitutiveNames() const;

protected:

private:

  ElementRegionBase & operator=( const ElementRegionBase & rhs );

  string_array m_materialList;
  string m_numericalMethod;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template< typename CONSTITUTIVE_TYPE >
string_array ElementRegionBase::getConstitutiveNames() const
{
  string_array rval;
  for( string const & matName : m_materialList )
  {
    Group const * const matModel = this->GetSubRegion( 0 )->GetConstitutiveModels()->GetGroup( matName );
    if( dynamic_cast< CONSTITUTIVE_TYPE const * >( matModel ) != nullptr )
    {
      rval.push_back( matName );
    }
  }
  return rval;
}

}



#endif /* GEOSX_MESH_ELEMENTREGIONBASE_HPP */
