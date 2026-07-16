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
 * @file ReactivePorosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITY_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITY_HPP_

#include "PorosityBase.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{
namespace constitutive
{

class ReactivePorosityUpdates : public PorosityBaseUpdates
{
public:

  ReactivePorosityUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 const > const & porosity_n,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView2d< real64 > const & dPorosity_dTemperature,
                           arrayView2d< real64 const > const & initialPorosity,
                           arrayView1d< real64 const > const & referencePorosity,
                           arrayView3d< real64, reactivefluid::USD_SPECIES > const & volumeFractions,
                           arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialVolumeFractions,
                           arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & volumeFractions_n,
                           integer const numKineticReactions,
                           arrayView1d< real64 const > const & molarWeights,
                           arrayView1d< real64 const > const & mineralDensities ):
    PorosityBaseUpdates( newPorosity,
                         porosity_n,
                         dPorosity_dPressure,
                         dPorosity_dTemperature,
                         initialPorosity,
                         referencePorosity ),
    m_volumeFractions( volumeFractions ),
    m_initialVolumeFractions( initialVolumeFractions ),
    m_volumeFractions_n( volumeFractions_n ),
    m_numKineticReactions( numKineticReactions ),
    m_molarWeights( molarWeights ),
    m_mineralDensities( mineralDensities )
  {}

  GEOS_HOST_DEVICE
  void computePorosity( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements,
                        arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                        arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions_n,
                        real64 & porosity,
                        real64 const & porosity_n,
                        integer const & numKineticReactions,
                        arrayView1d< real64 const > const & molarWeights,
                        arrayView1d< real64 const > const & mineralDensities ) const
  {
    real64 porosityIncrement = 0.0;

    for( integer r=0; r < numKineticReactions; ++r )
    {
      real64 volumeFractionIncrement = -kineticReactionMolarIncrements[r] * molarWeights[r]/mineralDensities[r];
      if (volumeFractionIncrement < 0) { // Dissolution
        real64 const decayRate = -volumeFractionIncrement / volumeFractions_n[r];
        volumeFractions[r] = volumeFractions_n[r] * exp(-decayRate);  
        if (volumeFractions[r] < 1e-6) {
          volumeFractions[r] = 1e-6; // Avoid negative or very small volume fractions
          volumeFractionIncrement = 0.0; // No further change in porosity from this reaction
        }
      } else { // Precipitation
        volumeFractions[r] = volumeFractions_n[r] + volumeFractionIncrement;  
      }
      porosityIncrement -= volumeFractionIncrement;
    }
    
    real64 const phi0 = porosity_n;
    real64 const logPorosity = log(std::max(phi0, 1e-4)) + porosityIncrement / phi0;
    porosity = exp(logPorosity);

    //porosity = porosity_n + porosityIncrement;

    if( porosity < 1e-4 )
    {
      porosity = 1e-4; // Avoid negative or very small porosity values that can cause numerical issues
    }
    else if( porosity > 1.0 )
    {
      porosity = 1.0;
    }

  }

  GEOS_HOST_DEVICE
  void updateFromReactions( localIndex const k,
                            localIndex const q,
                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements ) const
  {
    computePorosity( kineticReactionMolarIncrements,
                     m_volumeFractions[k][q],
                     m_volumeFractions_n[k][q],
                     m_newPorosity[k][q],
                     m_porosity_n[k][q],
                     m_numKineticReactions,
                     m_molarWeights,
                     m_mineralDensities );
  }

  GEOS_HOST_DEVICE
  inline
  real64 getVolumeFractionForMineral( localIndex const k,
                                      localIndex const q,
                                      localIndex const r ) const
  {
    return m_volumeFractions[k][q][r];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getVolumeFractionForMineral_n( localIndex const k,
                                        localIndex const q,
                                        localIndex const r ) const
  {
    return m_volumeFractions_n[k][q][r];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getInitialVolumeFractionForMineral( localIndex const k,
                                             localIndex const q,
                                             localIndex const r ) const
  {
    return m_initialVolumeFractions[k][q][r];
  }

protected:

  arrayView3d< real64, reactivefluid::USD_SPECIES > m_volumeFractions;
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > m_initialVolumeFractions;
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > const m_volumeFractions_n;

  integer const m_numKineticReactions;
  arrayView1d< real64 const > const m_molarWeights;
  arrayView1d< real64 const > const m_mineralDensities;
};


class ReactivePorosity : public PorosityBase
{
public:
  ReactivePorosity( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ReactivePorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void saveConvergedState() const override;

  integer numKineticReactions() const { return m_numKineticReactions; }

  virtual void initializeState() const override;

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * defaultInitialVolumeFractionsString() { return "defaultInitialVolumeFractions"; }
    static constexpr char const * initialVolumeFractionsString() { return "initialVolumeFractions"; }
    static constexpr char const * volumeFractionsString() { return "volumeFractions"; }
    static constexpr char const * volumeFractions_nString() { return "volumeFractions_n"; }
    static constexpr char const * molarWeightsString() { return "molarWeights"; }
    static constexpr char const * mineralDensitiesString() { return "mineralDensities"; }
  } viewKeys;


  using KernelWrapper = ReactivePorosityUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( m_newPorosity,
                          m_porosity_n,
                          m_dPorosity_dPressure,
                          m_dPorosity_dTemperature,
                          m_initialPorosity,
                          m_referencePorosity,
                          m_volumeFractions,
                          m_initialVolumeFractions,
                          m_volumeFractions_n,
                          m_numKineticReactions,
                          m_molarWeights,
                          m_mineralDensities );
  }


private:
  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts );

  array1d< real64 > m_defaultInitialVolumeFractions;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_volumeFractions;
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_initialVolumeFractions;
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_volumeFractions_n;

  integer m_numKineticReactions;
  array1d< real64 > m_molarWeights;
  array1d< real64 > m_mineralDensities;

};


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITY_HPP_
