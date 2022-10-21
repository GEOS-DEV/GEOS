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
                                     real64 const & fluidPressure_n,
                                     real64 const & fluidPressure,
                                     real64 const ( &strainIncrement )[6],
                                     real64 const & gravityAcceleration,
                                     real64 const ( &gravityVector )[3],
                                     real64 const & solidDensity,
                                     real64 const & initialFluidDensity,
                                     real64 const & fluidDensity_n,
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
    // Compute total stress increment and its derivative
    computeTotalStress( k,
                        q,
                        initialFluidPressure,
                        fluidPressure,
                        strainIncrement,
                        totalStress,
                        dTotalStress_dPressure,
                        stiffness );

    // Compute porosity and its derivatives
    real64 const deltaFluidPressure = fluidPressure - fluidPressure_n;
    real64 porosity;
    real64 porosity_n;
    real64 porosityInit;
    real64 dPorosity_dVolStrain;
    real64 dPorosity_dPressure;
    computePorosity( k,
                     q,
                     deltaFluidPressure,
                     strainIncrement,
                     porosity,
                     porosity_n,
                     porosityInit,
                     dPorosity_dVolStrain,
                     dPorosity_dPressure );

    // Compute body force vector and its derivatives
    if( gravityAcceleration > 0.0 )
    {
      computeBodyForce( solidDensity,
                        initialFluidDensity,
                        fluidDensity,
                        dFluidDensity_dPressure,
                        porosityInit,
                        porosity,
                        dPorosity_dVolStrain,
                        dPorosity_dPressure,
                        gravityVector,
                        bodyForceIncrement,
                        dBodyForce_dVolStrainIncrement,
                        dBodyForce_dPressure );
    }

    // Compute fluid mass contents and  its derivatives
    fluidMassContentIncrement = porosity * fluidDensity - porosity_n * fluidDensity_n;
    dFluidMassContent_dVolStrainIncrement = dPorosity_dVolStrain * fluidDensity;
    dFluidMassContent_dPressure = dPorosity_dPressure * fluidDensity + porosity * dFluidDensity_dPressure;

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
  }

  GEOSX_HOST_DEVICE
  void smallStrainUpdateThermalSinglePhase( localIndex const k,
                                            localIndex const q,
                                            real64 const & initialFluidPressure,
                                            real64 const & fluidPressure_n,
                                            real64 const & fluidPressure,
                                            real64 const & initialTemperature,
                                            real64 const & temperature,
                                            real64 const ( &strainIncrement )[6],
                                            real64 ( & totalStress )[6],
                                            DiscretizationOps & stiffness ) const
  {
    // Compute total stress increment and its derivative
    computeTotalStressThermal( k,
                               q,
                               initialFluidPressure,
                               fluidPressure,
                               initialTemperature,
                               temperature,
                               strainIncrement,
                               totalStress,
                               stiffness );

    // Compute porosity
    real64 const deltaFluidPressure = fluidPressure - fluidPressure_n;
    real64 porosity;
    real64 porosity_n;
    real64 porosityInit;
    real64 dPorosity_dVolStrain; // No use. Just to input something
    real64 dPorosity_dPressure; // No use. Just to input something
    computePorosity( k,
                     q,
                     deltaFluidPressure,
                     strainIncrement,
                     porosity,
                     porosity_n,
                     porosityInit,
                     dPorosity_dVolStrain,
                     dPorosity_dPressure );
  }

  template< int NUM_MAX_COMPONENTS >
  GEOSX_HOST_DEVICE
  void smallStrainUpdateMultiphase( localIndex const k,
                                    localIndex const q,
                                    localIndex const NP,
                                    localIndex const NC,
                                    real64 const & initialFluidPressure,
                                    real64 const & fluidPressure_n,
                                    real64 const & fluidPressure,
                                    real64 const ( &strainIncrement )[6],
                                    real64 const & gravityAcceleration,
                                    real64 const ( &gravityVector )[3],
                                    real64 const & solidDensity,
                                    real64 const & initialFluidTotalMassDensity,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & fluidPhaseDensity,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & fluidPhaseDensity_n,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dFluidPhaseDensity,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & fluidPhaseCompFrac,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & fluidPhaseCompFrac_n,
                                    arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC -2 > const & dFluidPhaseCompFrac,
                                    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & fluidPhaseMassDensity,
                                    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dFluidPhaseMassDensity,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & fluidPhaseSaturation,
                                    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & fluidPhaseSaturation_n,
                                    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dFluidPhaseSaturation,
                                    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dGlobalCompFraction_dGlobalCompDensity,
                                    real64 ( & totalStress )[6],
                                    real64 ( & dTotalStress_dPressure )[6],
                                    real64 ( & bodyForceIncrement )[3],
                                    real64 ( & dBodyForce_dVolStrainIncrement )[3],
                                    real64 ( & dBodyForce_dPressure )[3],
                                    real64 ( & dBodyForce_dComponents )[3][NUM_MAX_COMPONENTS],
                                    real64 ( & componentMassContentIncrement )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dVolStrainIncrement )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dPressure )[NUM_MAX_COMPONENTS],
                                    real64 ( & dComponentMassContent_dComponents )[NUM_MAX_COMPONENTS][NUM_MAX_COMPONENTS],
                                    DiscretizationOps & stiffness,
                                    real64 & poreVolumeConstraint,
                                    real64 (&dPoreVolumeConstraint_dPressure ),
                                    real64 (& dPoreVolumeConstraint_dComponents )[1][NUM_MAX_COMPONENTS] ) const
  {
    // Compute total stress increment and its derivatives
    computeTotalStress( k,
                        q,
                        initialFluidPressure,
                        fluidPressure,
                        strainIncrement,
                        totalStress,
                        dTotalStress_dPressure,
                        stiffness );

    // Compute porosity and its derivatives
    real64 const deltaFluidPressure = fluidPressure - fluidPressure_n;
    real64 porosity;
    real64 porosity_n;
    real64 porosityInit;
    real64 dPorosity_dVolStrain;
    real64 dPorosity_dPressure;
    computePorosity( k,
                     q,
                     deltaFluidPressure,
                     strainIncrement,
                     porosity,
                     porosity_n,
                     porosityInit,
                     dPorosity_dVolStrain,
                     dPorosity_dPressure );

    // Compute body force and its derivative
    using Deriv = constitutive::multifluid::DerivativeOffset;

    if( gravityAcceleration > 0.0 )
    {
      // Compute mixture density
      real64 fluidTotalMassDensity = fluidPhaseSaturation( 0 ) * fluidPhaseMassDensity( 0 );
      real64 dFluidTotalMassDensity_dPressure = dFluidPhaseSaturation( 0, Deriv::dP ) * fluidPhaseMassDensity( 0 )
                                                + fluidPhaseSaturation( 0 ) * dFluidPhaseDensity( 0, Deriv::dP );
      for( integer i = 1; i < NP; ++i )
      {
        fluidTotalMassDensity = fluidTotalMassDensity + fluidPhaseSaturation( i ) * fluidPhaseMassDensity( i );
        dFluidTotalMassDensity_dPressure = dFluidTotalMassDensity_dPressure
                                           + dFluidPhaseSaturation( i, Deriv::dP ) * fluidPhaseMassDensity( i )
                                           + fluidPhaseSaturation( i ) * dFluidPhaseDensity( i, Deriv::dP );
      }
      real64 dFluidTotalMassDensity_dComponents[NUM_MAX_COMPONENTS]{};
      real64 dFluidPhaseMassDensity_dC[NUM_MAX_COMPONENTS];
      for( integer ip = 0; ip < NP; ++ip )
      {
        applyChainRule( NC,
                        dGlobalCompFraction_dGlobalCompDensity,
                        dFluidPhaseMassDensity[ip],
                        dFluidPhaseMassDensity_dC,
                        Deriv::dC );
        for( integer jc = 0; jc < NC; ++jc )
        {
          dFluidTotalMassDensity_dComponents[jc] = dFluidTotalMassDensity_dComponents[jc]
                                                   + dFluidPhaseSaturation( ip, Deriv::dC+jc ) * fluidPhaseMassDensity( ip )
                                                   + fluidPhaseSaturation( ip ) * dFluidPhaseMassDensity_dC[jc];
        }
      }
      LvArray::tensorOps::scale< NUM_MAX_COMPONENTS >( dFluidTotalMassDensity_dComponents, porosity );

      computeBodyForce( solidDensity,
                        initialFluidTotalMassDensity,
                        fluidTotalMassDensity,
                        dFluidTotalMassDensity_dPressure,
                        dFluidTotalMassDensity_dComponents,
                        porosityInit,
                        porosity,
                        dPorosity_dVolStrain,
                        dPorosity_dPressure,
                        gravityVector,
                        bodyForceIncrement,
                        dBodyForce_dVolStrainIncrement,
                        dBodyForce_dPressure,
                        dBodyForce_dComponents );
    }

    // Compute component mass contents and their derivatives

    // --- temporary work arrays
    real64 dPhaseAmount_dC[NUM_MAX_COMPONENTS];
    real64 dPhaseCompFrac_dC[NUM_MAX_COMPONENTS];

    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( componentMassContentIncrement, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( dComponentMassContent_dVolStrainIncrement, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS >( dComponentMassContent_dPressure, 0.0 );
    LvArray::tensorOps::fill< NUM_MAX_COMPONENTS, NUM_MAX_COMPONENTS >( dComponentMassContent_dComponents, 0.0 );

    for( integer ip = 0; ip < NP; ++ip )
    {
      real64 const phaseAmount    = porosity    * fluidPhaseSaturation( ip )    * fluidPhaseDensity( ip );
      real64 const phaseAmount_n = porosity_n * fluidPhaseSaturation_n( ip ) * fluidPhaseDensity_n( ip );

      real64 const dPhaseAmount_dP = dPorosity_dPressure * fluidPhaseSaturation( ip ) * fluidPhaseDensity( ip )
                                     + porosity * ( dFluidPhaseSaturation( ip, Deriv::dP ) * fluidPhaseDensity( ip )
                                                    + fluidPhaseSaturation( ip ) * dFluidPhaseDensity( ip, Deriv::dP ) );

      // assemble density dependence
      applyChainRule( NC,
                      dGlobalCompFraction_dGlobalCompDensity,
                      dFluidPhaseDensity[ip],
                      dPhaseAmount_dC,
                      Deriv::dC );

      for( integer jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * fluidPhaseSaturation( ip )
                              + fluidPhaseDensity( ip ) * dFluidPhaseSaturation( ip, Deriv::dC+jc );
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * porosity;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( integer ic = 0; ic < NC; ++ic )
      {
        componentMassContentIncrement[ic] = componentMassContentIncrement[ic]
                                            + phaseAmount * fluidPhaseCompFrac( ip, ic )
                                            - phaseAmount_n * fluidPhaseCompFrac_n( ip, ic );

        dComponentMassContent_dPressure[ic] = dPhaseAmount_dP * fluidPhaseCompFrac( ip, ic )
                                              + phaseAmount * dFluidPhaseCompFrac( ip, ic, Deriv::dP );

        dComponentMassContent_dVolStrainIncrement[ic] = dComponentMassContent_dVolStrainIncrement[ic]
                                                        + dPorosity_dVolStrain
                                                        * fluidPhaseDensity( ip )
                                                        * fluidPhaseSaturation( ip )
                                                        * fluidPhaseCompFrac( ip, ic );

        applyChainRule( NC,
                        dGlobalCompFraction_dGlobalCompDensity,
                        dFluidPhaseCompFrac[ip][ic],
                        dPhaseCompFrac_dC,
                        Deriv::dC );

        for( integer jc = 0; jc < NC; ++jc )
        {
          dComponentMassContent_dComponents[ic][jc] = dComponentMassContent_dComponents[ic][jc]
                                                      + dPhaseCompFrac_dC[jc] * phaseAmount
                                                      + fluidPhaseCompFrac( ip, ic ) * dPhaseAmount_dC[jc];
        }
      }
    }

    // Compute pore volume constraint and its derivatives
    poreVolumeConstraint = 1.0;
    dPoreVolumeConstraint_dPressure = 0.0;
    LvArray::tensorOps::fill< 1, NUM_MAX_COMPONENTS >( dPoreVolumeConstraint_dComponents, 0.0 );
    for( integer ip = 0; ip < NP; ++ip )
    {
      poreVolumeConstraint = poreVolumeConstraint - fluidPhaseSaturation( ip );
      dPoreVolumeConstraint_dPressure = dPoreVolumeConstraint_dPressure
                                        - dFluidPhaseSaturation( ip, Deriv::dP ) * porosity
                                        - dPorosity_dPressure * fluidPhaseSaturation( ip );

      for( integer jc = 0; jc < NC; ++jc )
      {
        dPoreVolumeConstraint_dComponents[0][jc] = dPoreVolumeConstraint_dComponents[0][jc]
                                                   - dFluidPhaseSaturation( ip, Deriv::dC+jc )  * porosity;
      }
    }
    poreVolumeConstraint = poreVolumeConstraint * porosity;

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
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


  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficient( k, bulkModulus );
  }

  GEOSX_HOST_DEVICE
  void computeBodyForce( real64 const & solidDensity,
                         real64 const & initialFluidDensity,
                         real64 const & fluidDensity,
                         real64 const & dFluidDensity_dPressure,
                         real64 const & porosityInit,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const ( &gravityVector )[3],
                         real64 ( & bodyForceIncrement )[3],
                         real64 ( & dBodyForce_dVolStrainIncrement )[3],
                         real64 ( & dBodyForce_dPressure )[3] ) const
  {
    real64 const mixtureDensity = ( 1.0 - porosity ) * solidDensity
                                  + porosity * fluidDensity;
    real64 const initialMixtureDensity = ( 1.0 - porosityInit ) * solidDensity
                                         + porosityInit * initialFluidDensity;
    real64 const mixtureDensityIncrement = mixtureDensity - initialMixtureDensity;

    real64 const dMixtureDens_dVolStrainIncrement = dPorosity_dVolStrain * ( -solidDensity + fluidDensity );
    real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -solidDensity + fluidDensity )
                                          + ( 1.0 - porosity ) * m_porosityUpdate.dGrainDensity_dPressure()
                                          + porosity * dFluidDensity_dPressure;

    LvArray::tensorOps::scaledCopy< 3 >( bodyForceIncrement, gravityVector, mixtureDensityIncrement );
    LvArray::tensorOps::scaledCopy< 3 >( dBodyForce_dVolStrainIncrement, gravityVector, dMixtureDens_dVolStrainIncrement );
    LvArray::tensorOps::scaledCopy< 3 >( dBodyForce_dPressure, gravityVector, dMixtureDens_dPressure );
  }

  template< int NUM_MAX_COMPONENTS >
  GEOSX_HOST_DEVICE
  void computeBodyForce( real64 const & solidDensity,
                         real64 const & initialFluidTotalMassDensity,
                         real64 const & fluidTotalMassDensity,
                         real64 const & dFluidTotalMassDensity_dPressure,
                         real64 const ( &dFluidTotalMassDensity_dComponents)[NUM_MAX_COMPONENTS],
                         real64 const & porosityInit,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const ( &gravityVector )[3],
                         real64 ( & bodyForceIncrement )[3],
                         real64 ( & dBodyForce_dVolStrainIncrement )[3],
                         real64 ( & dBodyForce_dPressure )[3],
                         real64 ( & dBodyForce_dComponents )[3][NUM_MAX_COMPONENTS] ) const
  {
    computeBodyForce( solidDensity,
                      initialFluidTotalMassDensity,
                      fluidTotalMassDensity,
                      dFluidTotalMassDensity_dPressure,
                      porosityInit,
                      porosity,
                      dPorosity_dVolStrain,
                      dPorosity_dPressure,
                      gravityVector,
                      bodyForceIncrement,
                      dBodyForce_dVolStrainIncrement,
                      dBodyForce_dPressure );

    LvArray::tensorOps::Rij_eq_AiBj< 3, NUM_MAX_COMPONENTS >( dBodyForce_dComponents, gravityVector, dFluidTotalMassDensity_dComponents );
  }

  GEOSX_HOST_DEVICE
  void computePorosity( localIndex const k,
                        localIndex const q,
                        real64 const & deltaFluidPressure,
                        real64 const ( &strainIncrement )[6],
                        real64 & porosity,
                        real64 & porosity_n,
                        real64 & porosityInit,
                        real64 & dPorosity_dVolStrain,
                        real64 & dPorosity_dPressure ) const
  {
    m_porosityUpdate.updateFromPressureAndStrain( k,
                                                  q,
                                                  deltaFluidPressure,
                                                  strainIncrement,
                                                  dPorosity_dPressure,
                                                  dPorosity_dVolStrain );

    porosity = m_porosityUpdate.getPorosity( k, q );
    porosity_n = m_porosityUpdate.getPorosity_n( k, q );
    porosityInit = m_porosityUpdate.getInitialPorosity( k, q );
  }

  GEOSX_HOST_DEVICE
  void computeTotalStress( localIndex const k,
                           localIndex const q,
                           real64 const & initialFluidPressure,
                           real64 const & fluidPressure,
                           real64 const ( &strainIncrement )[6],
                           real64 ( & totalStress )[6],
                           real64 ( & dTotalStress_dPressure )[6],
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

    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * fluidPressure + initialBiotCoefficient * initialFluidPressure );

    dTotalStress_dPressure[0] = -biotCoefficient;
    dTotalStress_dPressure[1] = -biotCoefficient;
    dTotalStress_dPressure[2] = -biotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;
  }

  GEOSX_HOST_DEVICE
  void computeTotalStressThermal( localIndex const k,
                                  localIndex const q,
                                  real64 const & initialFluidPressure,
                                  real64 const & fluidPressure,
                                  real64 const & initialTemperature,
                                  real64 const & temperature,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & totalStress )[6],
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

    // Add the contribution of the pore pressure into the total stress
    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * fluidPressure + initialBiotCoefficient * initialFluidPressure );

    // Add the contribution of the thermal expansion into the total stress
    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -3 * thermalExpansionCoefficient * bulkModulus * ( temperature - initialTemperature ) );
  }
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

  /**
   * @brief Const/non-mutable accessor for biot coefficient
   * @return Accessor
   */
  arrayView1d< real64 const > const getBiotCoefficient() const
  {
    return getPorosityModel().getBiotCoefficient();
  }

private:
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPermModel;
};



}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_ */
