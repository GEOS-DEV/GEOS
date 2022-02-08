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
 * @file PorousSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE type of the porosity model
 */
template< typename SOLID_TYPE >
class PorousSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousSolidUpdates( SOLID_TYPE const & solidModel,
                      BiotPorosity const & porosityModel,
                      ConstantPermeability const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >( solidModel, porosityModel, permModel )
  {}

  GEOSX_HOST_DEVICE
  void smallStrainUpdateSinglePhase( localIndex const k,
                                     localIndex const q,
                                     real64 const & initialFluidPressure,
                                     real64 const & fluidPressureOld,
                                     real64 const & deltaFluidPressure,
                                     real64 const ( &strainIncrement )[6],
                                     real64 const & gravityAcceleration,
                                     real64 const ( &gravityVector )[3],
                                     real64 const & solidDensity,
                                     real64 const & initialFluidDensity,
                                     real64 const & fluidDensityOld,
                                     real64 const & fluidDensity,
                                     real64 const & dFluidDensity_dPressure,
                                     real64 ( & totalStress )[6],
                                     real64 ( & dTotalStress_dPressure )[6],
                                     real64 ( & bodyForceIncrement )[3],
                                     real64 ( & dBodyForce_dVolStrainIncrement )[3],
                                     real64 ( & dBodyForce_dPressure )[3],
                                     real64 & fluidMassContentIncrement,
                                     real64 & dFluidMassContent_dPressure,
                                     real64 & dFluidMassContent_dVolStrainIncrement,
                                     DiscretizationOps & stiffness ) const
  {
    // Compute total stress increment and its derivative w.r.t. pressure
    m_solidUpdate.smallStrainUpdate( k,
                                     q,
                                     strainIncrement,
                                     totalStress, // first effective stress increment accumulated
                                     stiffness );

    updateBiotCoefficient( k );

    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
    real64 const initialBiotCoefficient = biotCoefficient; // temporary
    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * ( fluidPressureOld + deltaFluidPressure ) + initialBiotCoefficient * initialFluidPressure );

    dTotalStress_dPressure[0] = -biotCoefficient;
    dTotalStress_dPressure[1] = -biotCoefficient;
    dTotalStress_dPressure[2] = -biotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;

    real64 dPorosity_dPressure;
    real64 dPorosity_dVolStrain;
    m_porosityUpdate.updateFromPressureAndStrain( k,
                                                  q,
                                                  deltaFluidPressure,
                                                  strainIncrement,
                                                  dPorosity_dPressure,
                                                  dPorosity_dVolStrain );

    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    real64 const porosityOld = m_porosityUpdate.getOldPorosity( k, q );
    real64 const porosityInit = m_porosityUpdate.getInitialPorosity( k, q );

    // Compute body force vector and its derivatives w.r.t. to
    // volumetric strain and pressure. The following assumption
    // are made at the moment:
    // 1. dMixtureDens_dVolStrainIncrement is neglected,
    // 2. grains are assumed incompressible
    real64 const mixtureDensity = ( 1.0 - porosity ) * solidDensity + porosity * fluidDensity;
    real64 const initialMixtureDensity = ( 1.0 - porosityInit ) * solidDensity + porosityInit * initialFluidDensity;
    real64 const mixtureDensityIncrement = mixtureDensity - initialMixtureDensity;

    real64 const dMixtureDens_dVolStrainIncrement = 0.0;
    real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -solidDensity + fluidDensity )
                                          + porosity * dFluidDensity_dPressure;

    if( gravityAcceleration > 0.0 )
    {
      LvArray::tensorOps::scaledCopy< 3 >( bodyForceIncrement, gravityVector, mixtureDensityIncrement );
      LvArray::tensorOps::scaledCopy< 3 >( dBodyForce_dVolStrainIncrement, gravityVector, dMixtureDens_dVolStrainIncrement );
      LvArray::tensorOps::scaledCopy< 3 >( dBodyForce_dPressure, gravityVector, dMixtureDens_dPressure );
    }

    // Compute fluid mass contents and derivatives w.r.t. to
    // volumetric strain and pressure
    fluidMassContentIncrement = porosity * fluidDensity - porosityOld * fluidDensityOld;
    dFluidMassContent_dVolStrainIncrement = dPorosity_dVolStrain * fluidDensity;
    dFluidMassContent_dPressure = dPorosity_dPressure * fluidDensity + porosity * dFluidDensity_dPressure;

    //

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
  }

  template< int NUM_MAX_COMPONENTS >
  GEOSX_HOST_DEVICE
  void smallStrainUpdateMultiphase( localIndex const k,
                                    localIndex const q,
                                    localIndex const NP,
                                    localIndex const NC,
                                    real64 const & initialFluidPressure,
                                    real64 const & fluidPressureOld,
                                    real64 const & deltaFluidPressure,
                                    real64 const ( &strainIncrement )[6],
                                    real64 const & gravityAcceleration,
                                    real64 const ( &gravityVector )[3],
                                    real64 const & solidDensity,
                                    real64 const & initialFluidTotalMassDensity,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & fluidPhaseDensity,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & fluidPhaseDensityOld,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & dFluidPhaseDensity_dPressure,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dFluidPhaseDensity_dGlobalCompFraction,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & fluidPhaseCompFrac,
                                    arraySlice2d< real64 const, compflow::USD_PHASE_COMP - 1 > const & fluidPhaseCompFracOld,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & dFluidPhaseCompFrac_dPressure,
                                    arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC -2 > const & dFluidPhaseCompFraction_dGlobalCompFraction,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & fluidPhaseMassDensity,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & fluidPhaseSaturation,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & fluidPhaseSaturationOld,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dFluidPhaseSaturation_dPressure,
                                    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dFluidPhaseSaturation_dGlobalCompDensity,
                                    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dGlobalCompFraction_dGlobalCompDensity,
                                    real64 ( & totalStress )[6],
                                    real64 ( & dTotalStress_dPressure )[6],
                                    real64 ( & bodyForceIncrement )[3],
                                    real64 ( & dBodyForce_dVolStrainIncrement )[3],
                                    real64 ( & dBodyForce_dPressure )[3],
                                    real64 ( & componentMassContentIncrement )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dVolStrainIncrement )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dPressure )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dComponents )[NUM_MAX_COMPONENTS][NUM_MAX_COMPONENTS],
                                    DiscretizationOps & stiffness,
                                    real64 & poreVolumeConstraint,
                                    real64 (&dPoreVolumeConstraint_dPressure ),
                                    real64 (& dPoreVolumeConstraint_dComponents )[NUM_MAX_COMPONENTS] ) const
  {
    // Compute total stress increment and its derivative w.r.t. pressure
    m_solidUpdate.smallStrainUpdate( k,
                                     q,
                                     strainIncrement,
                                     totalStress, // first effective stress increment accumulated
                                     stiffness );

    updateBiotCoefficient( k );

    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
    real64 const initialBiotCoefficient = biotCoefficient; // temporary
    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * ( fluidPressureOld + deltaFluidPressure ) + initialBiotCoefficient * initialFluidPressure );

    dTotalStress_dPressure[0] = -biotCoefficient;
    dTotalStress_dPressure[1] = -biotCoefficient;
    dTotalStress_dPressure[2] = -biotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;

    real64 dPorosity_dPressure;
    real64 dPorosity_dVolStrain;
    m_porosityUpdate.updateFromPressureAndStrain( k,
                                                  q,
                                                  deltaFluidPressure,
                                                  strainIncrement,
                                                  dPorosity_dPressure,
                                                  dPorosity_dVolStrain );

    //real64 const
    real64 porosity = m_porosityUpdate.getPorosity( k, q );
    real64 const porosityOld = m_porosityUpdate.getOldPorosity( k, q );
    real64 const porosityInit = m_porosityUpdate.getInitialPorosity( k, q );

    // Compute body force vector.
    // No derivatives computed. Current assumptions:
    // (  i) dMixtureDens_dVolStrain contribution is neglected
    // ( ii) grains are assumed incompressible
    // (iii) TODO add dMixtureDens_dPressure and dMixtureDens_dGlobalCompDensity
    if( gravityAcceleration > 0.0 )
    {
      // Compute mixture density
      real64 mixtureDensityNew = fluidPhaseSaturation( 0 ) * fluidPhaseMassDensity( 0 );
      for( localIndex i = 1; i < NP; ++i )
      {
        mixtureDensityNew += fluidPhaseSaturation( i ) * fluidPhaseMassDensity( i );
      }
      mixtureDensityNew *= porosity;
      mixtureDensityNew += ( 1.0 - porosity ) * solidDensity;

      real64 mixtureDensityInit = initialFluidTotalMassDensity * porosityInit;
      mixtureDensityInit += ( 1.0 - porosityInit ) * solidDensity;

      real64 const mixtureDensityIncrement = mixtureDensityNew - mixtureDensityInit;
      LvArray::tensorOps::scaledCopy< 3 >( bodyForceIncrement, gravityVector, mixtureDensityIncrement );

      GEOSX_UNUSED_VAR( dBodyForce_dVolStrainIncrement );
      GEOSX_UNUSED_VAR( dBodyForce_dPressure );
    }

    // Compute component mass contents and derivatives w.r.t. to
    // volumetric strain, pressure and components

    // --- temporary work arrays
    real64 dPhaseAmount_dC[NUM_MAX_COMPONENTS];
    real64 dPhaseCompFrac_dC[NUM_MAX_COMPONENTS];

    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( componentMassContentIncrement, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( dComponentMassContent_dVolStrainIncrement, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( dComponentMassContent_dPressure, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS, NUM_MAX_COMPONENTS >( dComponentMassContent_dComponents, 0.0 );

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const phaseAmountNew = porosity * fluidPhaseSaturation( ip ) * fluidPhaseDensity( ip );
      real64 const phaseAmountOld = porosityOld * fluidPhaseSaturationOld( ip ) * fluidPhaseDensityOld( ip );

      real64 const dPhaseAmount_dP = dPorosity_dPressure * fluidPhaseSaturation( ip ) * fluidPhaseDensity( ip )
                                     + porosity * ( dFluidPhaseSaturation_dPressure( ip ) * fluidPhaseDensity( ip )
                                                    + fluidPhaseSaturation( ip ) * dFluidPhaseDensity_dPressure( ip ) );

      // assemble density dependence
      applyChainRule( NC,
                      dGlobalCompFraction_dGlobalCompDensity,
                      dFluidPhaseDensity_dGlobalCompFraction[ip],
                      dPhaseAmount_dC );

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * fluidPhaseSaturation( ip )
                              + fluidPhaseDensity( ip ) * dFluidPhaseSaturation_dGlobalCompDensity( ip, jc );
        dPhaseAmount_dC[jc] *= porosity;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        componentMassContentIncrement[ic] = componentMassContentIncrement[ic]
                                            + phaseAmountNew * fluidPhaseCompFrac( ip, ic )
                                            - phaseAmountOld * fluidPhaseCompFracOld( ip, ic );

        dComponentMassContent_dPressure[ic] = dPhaseAmount_dP * fluidPhaseCompFrac( ip, ic )
                                              + phaseAmountNew * dFluidPhaseCompFrac_dPressure( ip, ic );

        dComponentMassContent_dVolStrainIncrement[ic] = dComponentMassContent_dVolStrainIncrement[ic]
                                                        + fluidPhaseDensity( ip )
                                                        * fluidPhaseSaturation( ip )
                                                        * fluidPhaseCompFrac( ip, ic );

        applyChainRule( NC,
                        dGlobalCompFraction_dGlobalCompDensity,
                        dFluidPhaseCompFraction_dGlobalCompFraction[ip][ic],
                        dPhaseCompFrac_dC );

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dComponentMassContent_dComponents[ic][jc] = dComponentMassContent_dComponents[ic][jc]
                                                      + dPhaseCompFrac_dC[jc] * phaseAmountNew
                                                      + fluidPhaseCompFrac( ip, ic ) * dPhaseAmount_dC[jc];
        }
      }
    }

    // --- Volume balance equation
    poreVolumeConstraint = 1.0;
    dPoreVolumeConstraint_dPressure = 0.0;
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( dPoreVolumeConstraint_dComponents, 0.0 );
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      poreVolumeConstraint = poreVolumeConstraint - fluidPhaseSaturation( ip );
      dPoreVolumeConstraint_dPressure = dPoreVolumeConstraint_dPressure
                                        - dFluidPhaseSaturation_dPressure( ip ) * porosity
                                        - dPorosity_dPressure * fluidPhaseSaturation( ip );

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPoreVolumeConstraint_dComponents[jc] = dPoreVolumeConstraint_dComponents[jc]
                                                - dFluidPhaseSaturation_dGlobalCompDensity( ip, jc )  * porosity;
      }
    }
    poreVolumeConstraint = poreVolumeConstraint * porosity;

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
  }

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficient( k, bulkModulus );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @note If the material model has a strain-dependent material stiffness (e.g.
   * any plasticity, damage, or nonlinear elastic model) then this interface will
   * not work.  Users should instead use one of the interfaces where a strain
   * tensor is provided as input.
   *
   * @param k the element number
   * @param stiffness the stiffness array
   */
  GEOSX_HOST_DEVICE
  void getElasticStiffness( localIndex const k, localIndex const q, real64 ( & stiffness )[6][6] ) const
  {
    m_solidUpdate.getElasticStiffness( k, q, stiffness );
  }

private:

  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_permUpdate;

};

/**
 * @brief PorousSolidBase class used for dispatch of all Porous solids.
 */
class PorousSolidBase
{};

/**
 * @brief Class to represent a porous material for poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE >
class PorousSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousSolidUpdates< SOLID_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( getSolidModel(),
                          getPorosityModel(),
                          getPermModel() );
  }

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getSolidModel().getDensity();
  }

private:
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPermModel;
};



}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_ */
