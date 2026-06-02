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

  GEOS_HOST_DEVICE
  virtual void updateSurfaceArea( localIndex const k,
                                  localIndex const q,
                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & initialSurfaceArea,
                                  arraySlice1d< real64, compflow::USD_COMP - 1 > const & surfaceArea ) const override final
  {
    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    real64 const initialPorosity = m_porosityUpdate.getInitialPorosity( k, q );

    for( integer r=0; r < initialSurfaceArea.size(); ++r )
    {
      real64 const volumeFraction_r = m_porosityUpdate.getVolumeFractionForMineral( k, q, r );
      real64 const volumeFraction_r_n = m_porosityUpdate.getVolumeFractionForMineral_n( k, q, r );
      real64 const initialVolumeFraction_r = m_porosityUpdate.getInitialVolumeFractionForMineral( k, q, r );

      if (volumeFraction_r - volumeFraction_r_n < 0) { // dissolution
        surfaceArea[r] = initialSurfaceArea[r] * pow( volumeFraction_r / initialVolumeFraction_r, 2.0/3.0 )
                       * pow( porosity / initialPorosity, 2.0/3.0 );
        if (volumeFraction_r < 1e-6) { // if the mineral is completely dissolved, set surface area to 0 to avoid numerical issues
          surfaceArea[r] = 0.0;
        }
      } else { // precipitation
        surfaceArea[r] = initialSurfaceArea[r] * pow( porosity / initialPorosity, 2.0/3.0 );
      }
    }
  }

private:
  using CoupledSolidUpdates< NullModel, ReactivePorosityBase, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, ReactivePorosityBase, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, ReactivePorosityBase, PERM_TYPE >::m_permUpdate;

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
