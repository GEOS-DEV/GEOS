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
#include "constitutive/solid/porosity/ReactivePorosityBase.hpp"
#include "constitutive/surfaceArea/ConstantSurfaceArea.hpp"
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
 * @tparam SURFACE_AREA_TYPE type of the reactive surface area model
 */
template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
class ReactiveSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   */
  ReactiveSolidUpdates( NullModel const & solidModel,
                        PORO_TYPE const & porosityModel,
                        PERM_TYPE const & permModel,
                        SURFACE_AREA_TYPE const & surfaceAreaModel ):
    CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >( solidModel, porosityModel, permModel ),
    m_surfaceAreaUpdate( surfaceAreaModel.createKernelWrapper() )
  {}

  GEOS_HOST_DEVICE
  virtual void updateStateFromPressureTemperatureAndReactions( localIndex const k,
                                                               localIndex const q,
                                                               real64 const & pressure,
                                                               real64 const & temperature,
                                                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements ) const override final
  {
    GEOS_UNUSED_VAR( temperature );

    m_porosityUpdate.updateFromReactions( k, q, kineticReactionMolarIncrements );
    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    m_permUpdate.updateFromPressureAndPorosity( k, q, pressure, porosity );
  }

  /**
   * @brief Update the reactive surface area of all the minerals in one quadrature point.
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   *
   * The surface areas are owned and computed by the surface area model, this only feeds it with
   * the porosity and the mineral volume fractions it needs.
   */
  GEOS_HOST_DEVICE
  void updateSurfaceArea( localIndex const k,
                          localIndex const q ) const
  {
    m_surfaceAreaUpdate.updateFromPorosityAndVolumeFractions( k, q,
                                                              m_porosityUpdate.getPorosity( k, q ),
                                                              m_porosityUpdate.getInitialPorosity( k, q ),
                                                              m_porosityUpdate.getVolumeFractions( k, q ),
                                                              m_porosityUpdate.getInitialVolumeFractions( k, q ) );
  }

private:
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_permUpdate;

  typename SURFACE_AREA_TYPE::KernelWrapper const m_surfaceAreaUpdate;

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
 * @tparam SURFACE_AREA_TYPE type of the reactive surface area model
 */

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
class ReactiveSolid : public CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >
{
public:


  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = ReactiveSolidUpdates< PORO_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >;

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
   *
   * The surface area model only shows up in the catalog name when it is not the default one,
   * so that the models relying on a constant surface area keep their historical name.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< SURFACE_AREA_TYPE, ConstantSurfaceArea > )   // default case
    {
      return string( "ReactiveSolid" ) + PERM_TYPE::catalogName();
    }
    else   // special cases
    {
      return string( "ReactiveSolid" ) + PERM_TYPE::catalogName() + SURFACE_AREA_TYPE::catalogName();
    }
  }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  /*
   * @brief get the volume fractions.
   * return a constant arrayView3d to the new volume fractions
   */
  arrayView3d< real64 const > const getVolumeFractions() const
  {
    return getPorosityModel().getVolumeFractions();
  }

  /*
   * @brief get the initial volume fractions.
   * return a constant arrayView1d to the initial volume fractions
   */
  arrayView1d< real64 const > const getInitialVolumeFractions() const
  {
    return getPorosityModel().getInitialVolumeFractions();
  }


  /**
   * @brief Create a instantiation of the ReactiveSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of ReactiveSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {

    return ReactiveSolidUpdates< PORO_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >( getSolidModel(),
                                                                            getPorosityModel(),
                                                                            getPermModel(),
                                                                            getSurfaceAreaModel() );
  }
private:
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getSolidModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPermModel;

  SURFACE_AREA_TYPE const & getSurfaceAreaModel() const
  { return this->getParent().template getGroup< SURFACE_AREA_TYPE >( this->m_surfaceAreaModelName ); }

};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_REACTIVESOLID_HPP_ */
