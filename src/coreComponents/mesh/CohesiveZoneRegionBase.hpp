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

#ifndef GEOS_MESH_COHESIVEZONEREGIONBASE_HPP
#define GEOS_MESH_COHESIVEZONEREGIONBASE_HPP

#include "mesh/ObjectManagerBase.hpp"

namespace geos
{

/**
 * @class CohesiveZoneRegionBase
 * @brief The CohesiveZoneRegionBase is the base class to manage the cohesive zone data stored at the node level.
 *
 * The CohesiveZoneRegionBase is the base class for the CohesiveZoneRegion class. It may be depreciated at
 * some point since no other classes are currently derived from CohesiveZoneRegionBase.
 */
class CohesiveZoneRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  CohesiveZoneRegionBase() = delete;

  /**
   * @brief Main constructor.
   * @param name the name of the particle region
   * @param parent the pointer to the parent group
   */
  CohesiveZoneRegionBase( string const & name, Group * const parent );


  /**
   * @brief Copy constructor.
   * @param init the particle region to be copied
   */
  CohesiveZoneRegionBase( const CohesiveZoneRegionBase & init );

  /**
   * @brief Default destructor.
   */
  virtual ~CohesiveZoneRegionBase() override;

  ///@}

  /**
   * @name Generation of the mesh region
   */
  ///@{

  /**
   * @brief Generate mesh.
   * @param blocks blocks where the mesh is generated
   */
  virtual void generateMesh( Group & blocks ) final
  {
    GEOS_UNUSED_VAR( blocks );
    GEOS_ERROR( "CohesiveZoneRegionBase::GenerateMesh() is non-op." );
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

//   /**
//    * @brief Get the name of the constitutive in the particle region.
//    * @tparam CONSTITUTIVE_TYPE the type of the constitutive model
//    * @return the string_array of the constitutive names
//    */
//   template< typename CONSTITUTIVE_TYPE >
//   string_array getConstitutiveNames() const;

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the particle subregions
    // static constexpr char const * particleSubRegions() { return "particleSubRegions"; }
  };

private:

  CohesiveZoneRegionBase & operator=( const CohesiveZoneRegionBase & rhs );

  /// List of materials for the particle region
  string_array m_materialList;

  /// Name of the mesh body that contains this region
  string m_meshBody;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


// /**
//  * @brief Get the names of all constitutive models of a specific type
//  *
//  * @tparam CONSTITUTIVE_TYPE type of constitutive model
//  * return string array with the names of the constitutive models
//  */
// template< typename CONSTITUTIVE_TYPE >
// string_array CohesiveZoneRegionBase::getConstitutiveNames() const
// {
//   string_array rval;
//   for( string const & matName : m_materialList )
//   {
//     if( this->getRegion( 0 ).getConstitutiveModels().hasGroup< CONSTITUTIVE_TYPE >( matName ) )
//     {
//       rval.emplace_back( matName );
//     }
//   }
//   return rval;
// }

}



#endif /* GEOS_MESH_COHESIVEZONEREGIONBASE_HPP */
