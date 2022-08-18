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
 * @file PorousDamageSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSDAMAGESOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSDAMAGESOLID_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/DamagePermeability.hpp"

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
class PorousDamageSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, DamagePermeability >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousDamageSolidUpdates( SOLID_TYPE const & solidModel,
                            BiotPorosity const & porosityModel,
                            DamagePermeability const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, DamagePermeability >( solidModel, porosityModel, permModel )
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
                                     real64 ( & fractureFlowTerm )[3],
                                     real64 ( & dFractureFlowTerm_dPressure )[3],
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

    // Compute fracture flow term and its derivative w.r.t pressure only
    computeFractureFlowTerm( k,
                             q,
                             fluidPressure,
                             fractureFlowTerm,
                             dFractureFlowTerm_dPressure );

    // Compute fluid mass contents and  its derivatives
    fluidMassContentIncrement = porosity * fluidDensity - porosity_n * fluidDensity_n;
    dFluidMassContent_dVolStrainIncrement = dPorosity_dVolStrain * fluidDensity;
    dFluidMassContent_dPressure = dPorosity_dPressure * fluidDensity + porosity * dFluidDensity_dPressure;

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );

    updateMatrixPermeability( k ); 
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

  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, DamagePermeability >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, DamagePermeability >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, DamagePermeability >::m_permUpdate;


  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficient( k, bulkModulus );
    
    // Update the Biot coefficient in the damage model
    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );

    m_solidUpdate.updateBiotCoefficient( k, biotCoefficient );

  }

  GEOSX_HOST_DEVICE
  void updateMatrixPermeability( localIndex const k ) const
  {
    // We tentatively update the permeability using the damage on the first quadrature point
    real64 const damage = fmax( fmin( 1.0, m_solidUpdate.getDamage( k, 0 ) ), 0.0 );

    m_permUpdate.updateDamagePermeability( k, damage );
  }

  // Do we need to consider the damage on the solid density?
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

  // Note: it is a tentative method to compute porosity, and is subject to further changes
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

    real64 const damage = fmax( fmin( 1.0, m_solidUpdate.getDamage( k, q ) ), 0.0 );
    real64 const damage_n = fmax( fmin( 1.0, m_solidUpdate.getOldDamage( k, q ) ), 0.0 );

    porosity = damage + ( 1 - damage ) * m_porosityUpdate.getPorosity( k, q );
    porosity_n = damage_n + ( 1 - damage_n ) * m_porosityUpdate.getPorosity_n( k, q );
    porosityInit = m_porosityUpdate.getInitialPorosity( k, q );

    dPorosity_dVolStrain *= ( 1 - damage );
    dPorosity_dPressure *= ( 1 - damage );
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

    real64 const damagedBiotCoefficient = m_solidUpdate.pressureDamageFunction( k, q ) * biotCoefficient;

    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -damagedBiotCoefficient * fluidPressure + initialBiotCoefficient * initialFluidPressure );

    dTotalStress_dPressure[0] = -damagedBiotCoefficient;
    dTotalStress_dPressure[1] = -damagedBiotCoefficient;
    dTotalStress_dPressure[2] = -damagedBiotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;
  }

  GEOSX_HOST_DEVICE
  void computeFractureFlowTerm( localIndex const k,
                                localIndex const q,
                                real64 const & fluidPressure,
                                real64 ( & fractureFlowTerm )[3],
                                real64 ( & dFractureFlowTerm_dPressure )[3] ) const
  {
    // Compute fracture flow term and its derivative w.r.t pressure
    real64 damageGrad[3]{};
    real64 pressureDamageGrad[3]{};

    m_solidUpdate.getDamageGrad( k, q, damageGrad );

    real64 const damage = m_solidUpdate.getDamage( k, q );

    LvArray::tensorOps::scaledCopy< 3 >( pressureDamageGrad, damageGrad, damage );
    LvArray::tensorOps::scaledCopy< 3 >( fractureFlowTerm, pressureDamageGrad, fluidPressure );

    LvArray::tensorOps::copy< 3 >( dFractureFlowTerm_dPressure, pressureDamageGrad );
  }

};

/**
 * @brief PorousDamageSolidBase class used for dispatch of all Porous damage solids.
 */
class PorousDamageSolidBase
{};

/**
 * @brief Class to represent a porous material for phase-field poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE >
class PorousDamageSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity, DamagePermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousDamageSolidUpdates< SOLID_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousDamageSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousDamageSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + SOLID_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousDamageSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousDamageSolidUpdates.
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
  using CoupledSolid< SOLID_TYPE, BiotPorosity, DamagePermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, DamagePermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, DamagePermeability >::getPermModel;
};



}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSDAMAGESOLID_HPP_ */
