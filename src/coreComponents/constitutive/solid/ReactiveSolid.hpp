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
 * @file ReactiveSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_REACTIVESOLID_HPP_
#define GEOS_CONSTITUTIVE_SOLID_REACTIVESOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/ReactivePorosity.hpp"
#include "constitutive/NullModel.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam PORO_TYPE type of the porosity model
 * @tparam PERM_TYPE type of the permeability model
 */
template< typename PORO_TYPE,
          typename PERM_TYPE >
class ReactiveSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   */
  ReactiveSolidUpdates( NullModel const & solidModel,
                            PORO_TYPE const & porosityModel,
                            PERM_TYPE const & permModel ):
    CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >( solidModel, porosityModel, permModel )
  {}

  GEOS_HOST_DEVICE
  void updateStateFromPressureAndReactions( localIndex const k,
                                            localIndex const q,
                                            real64 const & pressure,
                                            arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & kineticReactionMolarIncrements ) const
  {
    m_porosityUpdate.updateFromReactions( k, q, kineticReactionMolarIncrements );
    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    m_permUpdate.updateFromPressureAndPorosity( k, q, pressure, porosity );
  }
  
private:
  using CoupledSolidUpdates< NullModel, ReactivePorosity, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, ReactivePorosity, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, ReactivePorosity, PERM_TYPE >::m_permUpdate;

};


/**
 * @brief ReactiveSolidBase class used for dispatch of all Reactive solids.
 */
class ReactiveSolidBase
{};


/**
 * @brief Class to represent a porous material for flow simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material
 * for flow only simulations.
 *
 * @tparam PORO_TYPE type of porosity model
 * @tparam PERM_TYPE type of the permeability model
 */

template< typename PORO_TYPE,
          typename PERM_TYPE >
class ReactiveSolid : public CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >
{
public:


  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = ReactiveSolidUpdates< PORO_TYPE, PERM_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  ReactiveSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~ReactiveSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "ReactiveSolid" ) + PERM_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


  /**
   * @brief Create a instantiation of the ReactiveSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of ReactiveSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {

    return ReactiveSolidUpdates< PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                         getPorosityModel(),
                                                         getPermModel() );
  }
private:
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getSolidModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPermModel;

};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_REACTIVESOLID_HPP_ */
