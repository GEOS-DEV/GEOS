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
 * @file PorousReactiveSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_POROUSREACTIVESOLID_HPP_
#define GEOS_CONSTITUTIVE_SOLID_POROUSREACTIVESOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotReactivePorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE type of the porosity model
 */
template< typename SOLID_TYPE,
          typename PERM_TYPE >
class PorousReactiveSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousReactiveSolidUpdates( SOLID_TYPE const & solidModel,
                              BiotReactivePorosity const & porosityModel,
                              PERM_TYPE const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >( solidModel, porosityModel, permModel )
  {}

  GEOS_HOST_DEVICE
  virtual void updateStateReactionsFixedStress( localIndex const k,
                                                localIndex const q,
                                                real64 const & pressure,
                                                real64 const & pressure_k,
                                                real64 const & pressure_n,
                                                real64 const & temperature,
                                                real64 const & temperature_k,
                                                real64 const & temperature_n,
                                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > mineralReactionMolarIncrements ) const override final
  {
    updateBiotCoefficientAndAssignModuli( k );

    m_porosityUpdate.updateFixedStress( k, q,
                                        pressure, pressure_k, pressure_n,
                                        temperature, temperature_k, temperature_n,
                                        mineralReactionMolarIncrements );
  }

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
      real64 const initialVolumeFraction_r = m_porosityUpdate.getInitialVolumeFractionForMineral( k, q, r );
      surfaceArea[r] = initialSurfaceArea[r] * pow( volumeFraction_r / initialVolumeFraction_r, 2.0/3.0 )
                       * pow( porosity / initialPorosity, 2.0/3.0 );
    }
  }

  GEOS_HOST_DEVICE
  void smallStrainUpdateChemoMechanicsFixedStress( localIndex const k,
                                                   localIndex const q,
                                                   real64 const & timeIncrement,
                                                   real64 const & pressure,
                                                   real64 const & pressure_n,
                                                   real64 const & temperature,
                                                   real64 const & temperature_n,
                                                   real64 const & referenceTemperature,
                                                   arraySlice1d< real64 const, compflow::USD_COMP - 1 > mineralReactionMolarIncrements,
                                                   real64 const ( &strainIncrement )[6],
                                                   real64 ( & totalStress )[6],
                                                   DiscretizationOps & stiffness ) const
  {
    GEOS_UNUSED_VAR( referenceTemperature );

    real64 anelasticStrainIncrement = 0.0;

    for( integer r=0; r < mineralReactionMolarIncrements.size(); ++r )
    {
      real64 const molarWeight = m_porosityUpdate.getMolarWeights( r );
      real64 const mineralDensity = m_porosityUpdate.getMineralDensities( r );

      anelasticStrainIncrement -= mineralReactionMolarIncrements[r] * molarWeight/mineralDensity;
    }

    // Compute total stress increment and its derivative
    real64 const deltaPressureFromLastStep = pressure - pressure_n;
    real64 const deltaTemperatureFromLastStep = temperature - temperature_n;
    computeTotalStress( k,
                        q,
                        timeIncrement,
                        pressure,
                        deltaPressureFromLastStep,
                        deltaTemperatureFromLastStep,
                        anelasticStrainIncrement,
                        strainIncrement,
                        totalStress,
                        stiffness );
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
  GEOS_HOST_DEVICE
  inline
  void getElasticStiffness( localIndex const k, localIndex const q, real64 ( & stiffness )[6][6] ) const
  {
    m_solidUpdate.getElasticStiffness( k, q, stiffness );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @param [in] k the element number
   * @param [out] thermalExpansionCoefficient the thermal expansion coefficient
   */
  GEOS_HOST_DEVICE
  inline
  void getThermalExpansionCoefficient( localIndex const k, real64 & thermalExpansionCoefficient ) const
  {
    thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
  }

private:

  using CoupledSolidUpdates< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::m_permUpdate;

  GEOS_HOST_DEVICE
  inline
  void updateBiotCoefficientAndAssignModuli( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficientAndAssignModuli( k, bulkModulus );
  }

  GEOS_HOST_DEVICE
  inline
  void computeTotalStress( localIndex const k,
                           localIndex const q,
                           real64 const & timeIncrement,
                           real64 const & pressure,
                           real64 const & deltaPressureFromLastStep,
                           real64 const & deltaTemperatureFromLastStep,
                           real64 const & anelasticStrainIncrement,
                           real64 const ( &strainIncrement )[6],
                           real64 ( & totalStress )[6],
                           DiscretizationOps & stiffness ) const
  {
    // For now, the model only used for the sequential fixed stress scheme
    // So we ignore the derivatives wrt pressure and temperature
    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );

    real64 mechanicsStrainIncrement[6]{};
    mechanicsStrainIncrement[0] = strainIncrement[0] - thermalExpansionCoefficient * deltaTemperatureFromLastStep;
    mechanicsStrainIncrement[1] = strainIncrement[1] - thermalExpansionCoefficient * deltaTemperatureFromLastStep;
    mechanicsStrainIncrement[2] = strainIncrement[2] - thermalExpansionCoefficient * deltaTemperatureFromLastStep;
    mechanicsStrainIncrement[3] = strainIncrement[3];
    mechanicsStrainIncrement[4] = strainIncrement[4];
    mechanicsStrainIncrement[5] = strainIncrement[5];

    // Add the contributions of pore material stress/pressure
    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );

    // Compute total stress increment and its derivative w.r.t. pressure
    m_solidUpdate.smallStrainUpdate( k,
                                     q,
                                     timeIncrement,
                                     mechanicsStrainIncrement,
                                     totalStress, // first effective stress increment accumulated
                                     stiffness );

    // Compute effective stress increment for the porosity update
    real64 const bulkModulus = stiffness.m_bulkModulus;
    real64 const meanEffectiveStressIncrement = bulkModulus * ( mechanicsStrainIncrement[0] + mechanicsStrainIncrement[1] + mechanicsStrainIncrement[2] );

    m_porosityUpdate.updateMeanEffectiveStressIncrement( k, q, meanEffectiveStressIncrement );

    // Update mineral pressure
    real64 dMineralPres_dMeanEffStressIncre = 0.0;
    m_porosityUpdate.updatePoreMineralPressure( k, q,
                                                deltaPressureFromLastStep,
                                                meanEffectiveStressIncrement,
                                                anelasticStrainIncrement,
                                                dMineralPres_dMeanEffStressIncre );

    real64 const mineralPressure = m_porosityUpdate.getPoreMineralPressure( k );
    real64 const totalPorePressure = pressure + mineralPressure;

    // Add the contributions of pressure to the total stress
    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * totalPorePressure );

    // Add the contributions of mineral pressure to the stiffness
    stiffness.m_bulkModulus = bulkModulus - biotCoefficient * dMineralPres_dMeanEffStressIncre * bulkModulus;
  }

};

/**
 * @brief PorousReactiveSolidBase class used for dispatch of all Porous solids.
 */
class PorousReactiveSolidBase
{};

/**
 * @brief Class to represent a porous material for poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE,
          typename PERM_TYPE >
class PorousReactiveSolid : public CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousReactiveSolidUpdates< SOLID_TYPE, PERM_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousReactiveSolid( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< PERM_TYPE, ConstantPermeability > )   // default case
    {
      return string( "PorousReactive" ) + SOLID_TYPE::catalogName();
    }
    else   // special cases
    {
      return string( "PorousReactive" ) + SOLID_TYPE::catalogName() + PERM_TYPE::catalogName();
    }
  }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousReactiveSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousReactiveSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( getSolidModel(),
                          getPorosityModel(),
                          getPermModel() );
  }

  /**
   * @brief initialize the constitutive models fields.
   */
  virtual void initializeState() const override final;

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getSolidModel().getDensity();
  }

private:
  using CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::getPermModel;
};



}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_POROUSREACTIVESOLID_HPP_ */
