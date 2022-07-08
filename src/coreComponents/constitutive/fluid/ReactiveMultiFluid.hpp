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
 * @file ReactiveMultiFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP_


#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/chemicalReactions/EquilibriumReactions.hpp"
#include "constitutive/fluid/chemicalReactions/KineticReactions.hpp"

#include <memory>

namespace geosx
{

namespace constitutive
{

class ReactiveMultiFluid : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  ReactiveMultiFluid( string const & name,
                      Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for ReactiveMultiFluid.
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {

private:
  
  KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                 bool const useMass,
                 bool const isThermal,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity,
                 EquilibriumReactions const & equilibriumReactions,
                 KineticReactions const & kineticReactions, 
                 arrayView2d<real64> const &  primarySpeciesConcentration,
                 arrayView2d<real64> const & secondarySpeciesConcentration,
                 arrayView2d<real64> const & primarySpeciesTotalConcentration,
                 arrayView2d<real64> const & kineticReactionRates );        

  EquilibriumReactions::KernelWrapper m_equilibriumReactions;
  
  KineticReactions::KernelWrapper m_kineticReactions;

  arrayView2d<real64>  m_primarySpeciesConcentration;

  arrayView2d<real64>  m_secondarySpeciesConcentration;

  arrayView2d<real64>  m_primarySpeciesTotalConcentration;

  arrayView2d<real64>  m_kineticReactionRates;
  };


  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
  
  };

protected:

  virtual void postProcessInput() override;

private:

  /// Reaction related terms
  integer m_numPrimarySpecies;
  
  integer m_numSecondarySpecies;

  std::unique_ptr< EquilibriumReactions > m_equilibriumReaction;
  
  std::unique_ptr< KineticReactions > m_kineticReactions;

  array2d<real64>  m_primarySpeciesConcentration;

  array2d<real64>  m_secondarySpeciesConcentration;

  array2d<real64>  m_primarySpeciesTotalConcentration;

  array2d<real64>  m_kineticReactionRates;
};

GEOSX_HOST_DEVICE
inline void
ReactiveMultiFluid::KernelWrapper::
  compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
           arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
           arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const
{
  // I am assuming that the primary variable is the concentration of the primary species.
  for(int i=0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesTotalConcentration[i] = composition[i];
  }
  
  m_equilibriumReactions.updateConcentrations( temperature,
                                               primarySpeciesTotalConcentration,
                                               primarySpeciesContentration, 
                                               secondarySpeciesConcentration );

  m_kineticReactions.computeReactionRates( temperature,
                                           primarySpeciesContentration, 
                                           secondarySpeciesConcentration,
                                           kineticReactionRates );                                                                     
}

GEOSX_HOST_DEVICE inline void
ReactiveMultiFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geosx::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  GEOSX_UNUSED_VAR(q);

  compute( pressure,
           temperature,
           composition,
           m_primarySpeciesConcentration[k],
           m_secondarySpeciesConcentration[k],
           m_primarySpeciesTotalConcentration[k],
           m_reactionRates[k] );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP
