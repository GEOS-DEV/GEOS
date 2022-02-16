/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleBlock.hpp
 */

#ifndef GEOSX_MESH_PARTICLEBLOCK_HPP_
#define GEOSX_MESH_PARTICLEBLOCK_HPP_

#include "dataRepository/Group.hpp"
#include "ParticleBlockABC.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ParticleType.hpp"

namespace geosx
{

/**
 * @class ParticleBlock
 * Class deriving from ParticleBlockABC specializing the particle subregion
 * for a particle.
 */
class ParticleBlock : public ParticleBlockABC
{
public:

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static const string catalogName()
  { return "ParticleBlock"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleBlock::catalogName(); }

  ///@}

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  ParticleBlock() = delete;

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleBlock( string const & name, Group * const parent );

  /**
   * @brief Copy constructor.
   * @param[in] init the source to copy
   */
  ParticleBlock( const ParticleBlock & init ) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleBlock() override;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  virtual void setParticleType( ParticleType const particleType ) override;

  ///@}

  /**
   * @name Properties
   */
  ///@{

  /**
   * @brief Add a property to the ParticleBlock.
   * @tparam T type of the property
   * @param[in] propertyName the name of the property
   * @return a non-const reference to the property
   */
  template< typename T >
  T & addProperty( string const & propertyName )
  {
    m_externalPropertyNames.emplace_back( propertyName );
    return this->registerWrapper< T >( propertyName ).reference();
  }

  /**
   * @brief Helper function to apply a lambda function over all the external properties of the subregion
   * @tparam LAMBDA the type of the lambda function
   * @param lambda lambda function that is applied to the wrappers of external properties
   */
  template< typename LAMBDA >
  void forExternalProperties( LAMBDA && lambda )
  {
    for( auto & externalPropertyName : m_externalPropertyNames )
    {
      lambda( this->getWrapperBase( externalPropertyName ) );
    }
  }

  ///@}

protected:

private:
  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

};

}

#endif /* GEOSX_MESH_PARTICLEBLOCK_HPP_ */
