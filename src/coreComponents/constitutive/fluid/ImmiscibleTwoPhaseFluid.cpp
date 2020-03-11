/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ImmiscibleTwoPhaseFluid.cpp
 */
#include "ImmiscibleTwoPhaseFluid.hpp"

#include "common/Path.hpp"
#include "managers/ProblemManager.hpp"


using namespace std;

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ImmiscibleTwoPhaseFluid::ImmiscibleTwoPhaseFluid( std::string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::phaseCompressibilityString, &m_phaseCompressibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Phase compressibilities");

  registerWrapper( viewKeyStruct::phaseViscosibilityString, &m_phaseViscosibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Phase viscosity exponential coefficients");

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference pressure");

  registerWrapper( viewKeyStruct::referencePhaseDensityString, &m_referencePhaseDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference phase densities");

  registerWrapper( viewKeyStruct::referencePhaseViscosityString, &m_referencePhaseViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference phase viscosities");
  
}

ImmiscibleTwoPhaseFluid::~ImmiscibleTwoPhaseFluid()
{

}

void ImmiscibleTwoPhaseFluid::DeliverClone( string const & name,
                                                  Group * const parent,
                                                  std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr< ImmiscibleTwoPhaseFluid > newModel = std::make_unique<ImmiscibleTwoPhaseFluid>( name, parent );

  newModel->m_useMass = this->m_useMass;

  newModel->m_componentNames          = this->m_componentNames;
  newModel->m_componentMolarWeight    = this->m_componentMolarWeight;
  newModel->m_phaseNames              = this->m_phaseNames;

  newModel->m_phaseCompressibility    = this->m_phaseCompressibility;
  newModel->m_phaseViscosibility      = this->m_phaseViscosibility;

  newModel->m_referencePhaseDensity   = this->m_referencePhaseDensity;
  newModel->m_referencePhaseViscosity = this->m_referencePhaseViscosity;
  newModel->m_referencePressure       = this->m_referencePressure;
  
  clone = std::move( newModel );

}

 
void ImmiscibleTwoPhaseFluid::PostProcessInput() 
{

#define TWO_PHASE_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOSX_ERROR( "ImmiscibleTwoPhaseFluid: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)"); \
  }

  TWO_PHASE_CHECK_INPUT_LENGTH( m_phaseNames, NUM_PHASES, viewKeyStruct::phaseNamesString )
    
  TWO_PHASE_CHECK_INPUT_LENGTH( m_phaseCompressibility, NUM_PHASES, viewKeyStruct::phaseCompressibilityString )
  TWO_PHASE_CHECK_INPUT_LENGTH( m_phaseViscosibility, NUM_PHASES, viewKeyStruct::phaseViscosibilityString )    
    
  TWO_PHASE_CHECK_INPUT_LENGTH( m_referencePhaseDensity, NUM_PHASES, viewKeyStruct::referencePhaseDensityString )
  TWO_PHASE_CHECK_INPUT_LENGTH( m_referencePhaseViscosity, NUM_PHASES, viewKeyStruct::referencePhaseViscosityString )    
    
#undef TWO_PHASE_CHECK_INPUT_LENGTH

  for (localIndex ip = 0; ip < NUM_PHASES; ++ip)
  {
    GEOSX_ERROR_IF( m_phaseCompressibility[ip] < 0.0,
                   "ImmiscibleTwoPhaseFluid: invalid phase compressibility: " << m_phaseCompressibility[ip] );
    GEOSX_ERROR_IF( m_phaseViscosibility[ip] < 0.0,
                   "ImmiscibleTwoPhaseFluid: invalid phase viscosibility: " << m_phaseViscosibility[ip] );

    GEOSX_ERROR_IF( m_referencePhaseDensity[ip] <= 0.0,
                   "ImmiscibleTwoPhaseFluid: invalid reference phase density: " << m_referencePhaseDensity[ip] );
    GEOSX_ERROR_IF( m_referencePhaseViscosity[ip] <= 0.0,
                   "ImmiscibleTwoPhaseFluid: invalid reference phase viscosity: " << m_referencePhaseViscosity[ip] );
  }
  
}
  
void ImmiscibleTwoPhaseFluid::InitializePostSubGroups( Group * const group )
{
  MultiFluidBase::InitializePostSubGroups(group);
}

void ImmiscibleTwoPhaseFluid::PointUpdate( real64 const & pressure,
                                           real64 const & temperature,
                                           arraySlice1d<real64 const> const & composition,
                                           localIndex const k,
                                           localIndex const q )
{
  m_phaseFraction = 0;
  m_dPhaseFraction_dPressure = 0;
  m_dPhaseFraction_dTemperature = 0;
  m_dPhaseFraction_dGlobalCompFraction = 0;

  m_dPhaseDensity_dTemperature = 0;
  m_dPhaseDensity_dGlobalCompFraction = 0;

  m_dPhaseViscosity_dTemperature = 0;
  m_dPhaseViscosity_dGlobalCompFraction = 0;
  
  m_phaseCompFraction = 0;
  m_dPhaseCompFraction_dPressure = 0;
  m_dPhaseCompFraction_dTemperature = 0;
  m_dPhaseCompFraction_dGlobalCompFraction = 0;

  m_totalDensity = 0;
  m_dTotalDensity_dPressure = 0;
  m_dTotalDensity_dTemperature = 0;
  m_dTotalDensity_dGlobalCompFraction = 0;
  
  Compute( numFluidComponents(),
           numFluidPhases(),
           m_useMass,
           m_phaseNames,
           m_componentMolarWeight,
           pressure,
           temperature,
           composition,
           m_phaseFraction[k][q],
           m_dPhaseFraction_dPressure[k][q],
           m_dPhaseFraction_dTemperature[k][q],
           m_dPhaseFraction_dGlobalCompFraction[k][q],
           m_phaseDensity[k][q],
           m_dPhaseDensity_dPressure[k][q],
           m_dPhaseDensity_dTemperature[k][q],
           m_dPhaseDensity_dGlobalCompFraction[k][q],
           m_phaseViscosity[k][q],
           m_dPhaseViscosity_dPressure[k][q],
           m_dPhaseViscosity_dTemperature[k][q],
           m_dPhaseViscosity_dGlobalCompFraction[k][q],
           m_phaseCompFraction[k][q],
           m_dPhaseCompFraction_dPressure[k][q],
           m_dPhaseCompFraction_dTemperature[k][q],
           m_dPhaseCompFraction_dGlobalCompFraction[k][q],
           m_totalDensity[k][q],
           m_dTotalDensity_dPressure[k][q],
           m_dTotalDensity_dTemperature[k][q],
           m_dTotalDensity_dGlobalCompFraction[k][q],
           m_referencePressure,
           m_referencePhaseDensity,
           m_referencePhaseViscosity,      
           m_phaseCompressibility,
           m_phaseViscosibility );
}

void ImmiscibleTwoPhaseFluid::BatchUpdate( arrayView1d<real64 const> const & pressure,
                                           arrayView1d<real64 const> const & temperature,
                                           arrayView2d<real64 const> const & composition )
{
  m_phaseFraction = 0;
  m_dPhaseFraction_dPressure = 0;
  m_dPhaseFraction_dTemperature = 0;
  m_dPhaseFraction_dGlobalCompFraction = 0;

  m_dPhaseDensity_dTemperature = 0;
  m_dPhaseDensity_dGlobalCompFraction = 0;

  m_dPhaseViscosity_dTemperature = 0;
  m_dPhaseViscosity_dGlobalCompFraction = 0;
  
  m_phaseCompFraction = 0;
  m_dPhaseCompFraction_dPressure = 0;
  m_dPhaseCompFraction_dTemperature = 0;
  m_dPhaseCompFraction_dGlobalCompFraction = 0;

  m_totalDensity = 0;
  m_dTotalDensity_dPressure = 0;
  m_dTotalDensity_dTemperature = 0;
  m_dTotalDensity_dGlobalCompFraction = 0;
  
  MultiFluidBase::BatchUpdateKernel<ImmiscibleTwoPhaseFluid, RAJA::seq_exec>( pressure,
                                                                              temperature,
                                                                              composition,
                                                                              m_referencePressure,
                                                                              m_referencePhaseDensity,
                                                                              m_referencePhaseViscosity,           
                                                                              m_phaseCompressibility,
                                                                              m_phaseViscosibility );
}
 

void ImmiscibleTwoPhaseFluid::Compute( localIndex const GEOSX_UNUSED_PARAM( NC ),
                                       localIndex const GEOSX_UNUSED_PARAM( NP ),
                                       bool const GEOSX_UNUSED_PARAM( useMass ),
                                       arrayView1d<string const> const & GEOSX_UNUSED_PARAM( phaseNames ),
                                       arrayView1d<real64 const> const & GEOSX_UNUSED_PARAM( componentMolarWeight ),
                                       real64 const & pressure,
                                       real64 const & GEOSX_UNUSED_PARAM( temperature ),
                                       arraySlice1d<real64 const> const & GEOSX_UNUSED_PARAM( composition ),
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( phaseFraction ),
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( dPhaseFraction_dPressure ),
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( dPhaseFraction_dTemperature ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( dPhaseFraction_dGlobalCompFraction ),
                                       arraySlice1d<real64> const & phaseDensity,
                                       arraySlice1d<real64> const & dPhaseDensity_dPressure,
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( dPhaseDensity_dTemperature ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( dPhaseDensity_dGlobalCompFraction ),
                                       arraySlice1d<real64> const & phaseViscosity,
                                       arraySlice1d<real64> const & dPhaseViscosity_dPressure,
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( dPhaseViscosity_dTemperature ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( dPhaseViscosity_dGlobalCompFraction ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( phaseCompFraction ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dPressure ),
                                       arraySlice2d<real64> const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dTemperature ),
                                       arraySlice3d<real64> const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dGlobalCompFraction ),
                                       real64 & GEOSX_UNUSED_PARAM( totalDensity ),
                                       real64 & GEOSX_UNUSED_PARAM( dTotalDensity_dPressure ),
                                       real64 & GEOSX_UNUSED_PARAM( dTotalDensity_dTemperature ),
                                       arraySlice1d<real64> const & GEOSX_UNUSED_PARAM( dTotalDensity_dGlobalCompFraction ),
                                       real64 const & referencePressure,
                                       arrayView1d<real64 const> const referencePhaseDensity,
                                       arrayView1d<real64 const> const referencePhaseViscosity,
                                       arrayView1d<real64 const> const phaseCompressibility,
                                       arrayView1d<real64 const> const phaseViscosibility ) 
{
 
  for (localIndex ip = 0; ip < NUM_PHASES; ++ip)
  {
    makeExponentialRelation( ExponentApproximationType::Full,
                             referencePressure,
                             referencePhaseDensity[ip],
                             phaseCompressibility[ip], [&] ( auto relation )
    {
      ComputePhaseProperty( pressure, phaseDensity[ip], dPhaseDensity_dPressure[ip], relation );
    });
    makeExponentialRelation( ExponentApproximationType::Full,
                             referencePressure,
                             referencePhaseViscosity[ip],
                             phaseViscosibility[ip], [&] ( auto relation )
    {
      ComputePhaseProperty( pressure, phaseViscosity[ip], dPhaseViscosity_dPressure[ip], relation );
    });
  }
}

void ImmiscibleTwoPhaseFluid::Compute(real64 const & pressure, real64 const & temperature,
                                      arraySlice1d<double const> const & composition,
                                      arraySlice1d<real64> const & phaseFraction,
                                      arraySlice1d<real64> const & dPhaseFraction_dPressure,
                                      arraySlice1d<real64> const & dPhaseFraction_dTemperature,
                                      arraySlice2d<real64> const & dPhaseFraction_dGlobalCompFraction,
                                      arraySlice1d<real64> const & phaseDensity,
                                      arraySlice1d<real64> const & dPhaseDensity_dPressure,
                                      arraySlice1d<real64> const & dPhaseDensity_dTemperature,
                                      arraySlice2d<real64> const & dPhaseDensity_dGlobalCompFraction,
                                      arraySlice1d<real64> const & phaseViscosity,
                                      arraySlice1d<real64> const & dPhaseViscosity_dPressure,
                                      arraySlice1d<real64> const & dPhaseViscosity_dTemperature,
                                      arraySlice2d<real64> const & dPhaseViscosity_dGlobalCompFraction,
                                      arraySlice2d<real64> const & phaseCompFraction,
                                      arraySlice2d<real64> const & dPhaseCompFraction_dPressure,
                                      arraySlice2d<real64> const & dPhaseCompFraction_dTemperature,
                                      arraySlice3d<real64> const & dPhaseCompFraction_dGlobalCompFraction,
                                      real64 & totalDensity, real64 & dTotalDensity_dPressure,
                                      real64 & dTotalDensity_dTemperature,
                                      arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction) const
{
  Compute( numFluidComponents(), numFluidPhases(), m_useMass, m_phaseNames, m_componentMolarWeight,
           pressure, temperature, composition,
           phaseFraction, dPhaseFraction_dPressure, dPhaseFraction_dTemperature, dPhaseFraction_dGlobalCompFraction,
           phaseDensity, dPhaseDensity_dPressure, dPhaseDensity_dTemperature, dPhaseDensity_dGlobalCompFraction,
           phaseViscosity, dPhaseViscosity_dPressure, dPhaseViscosity_dTemperature, dPhaseViscosity_dGlobalCompFraction,
           phaseCompFraction, dPhaseCompFraction_dPressure, dPhaseCompFraction_dTemperature, dPhaseCompFraction_dGlobalCompFraction,
           totalDensity, dTotalDensity_dPressure, dTotalDensity_dTemperature, dTotalDensity_dGlobalCompFraction,
           m_referencePressure, m_referencePhaseDensity, m_referencePhaseViscosity, m_phaseCompressibility, m_phaseViscosibility );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ImmiscibleTwoPhaseFluid, std::string const &, Group * const )

} //namespace constitutive

} //namespace geosx
