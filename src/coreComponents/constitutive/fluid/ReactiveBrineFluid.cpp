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
 * @file ReactiveBrineFluid.cpp
 */
#include "ReactiveBrineFluid.hpp"

#include "constitutive/fluid/MultiFluidExtrinsicData.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename PHASE > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsThermalFluid"; }
};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiThermalFluid"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename PHASE >
string ReactiveBrineFluid< PHASE > ::catalogName()
{
  return TwoPhaseCatalogNames< PHASE > ::name();
}

template< typename PHASE >
ReactiveBrineFluid< PHASE > ::
ReactiveBrineFluid( string const & name, Group * const parent ):
  ReactiveMultiFluid( name, parent )
{
  registerWrapper( viewKeyStruct::phasePVTParaFilesString(), &m_phasePVTParaFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if( isThermal() )
  {
    getExtrinsicData< extrinsicMeshData::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getExtrinsicData< extrinsicMeshData::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }
}

template< typename PHASE >
bool ReactiveBrineFluid< PHASE > ::isThermal() const
{
  return ( PHASE1::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() &&
           PHASE2::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() );
}


template< typename PHASE >
std::unique_ptr< ConstitutiveBase >
ReactiveBrineFluid< PHASE > ::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ReactiveMultiFluid::deliverClone( name, parent );

  ReactiveBrineFluid & newConstitutiveRelation = dynamicCast< ReactiveBrineFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename PHASE >
integer ReactiveBrineFluid< PHASE > ::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}


template< typename PHASE >
void ReactiveBrineFluid< PHASE > ::postProcessInput()
{
  ReactiveMultiFluid::postProcessInput();

  GEOSX_THROW_IF_NE_MSG( numFluidPhases(), 2,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_NE_MSG( numFluidComponents(), 2,
                         GEOSX_FMT( "{}: invalid number of components", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_NE_MSG( m_phasePVTParaFiles.size(), 2,
                         GEOSX_FMT( "{}: invalid number of values in attribute '{}'", getFullName() ),
                         InputError );

  // NOTE: for now, the names of the phases are still hardcoded here
  // Later, we could read them from the XML file and we would then have a general class here

  string const expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_p1Index = PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );

  string const expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_p2Index = PVTFunctionHelpers::findName( m_phaseNames, expectedGasPhaseNames, viewKeyStruct::phaseNamesString() );

  createPVTModels();
}

template< typename PHASE >
void ReactiveBrineFluid< PHASE > ::createPVTModels()
{

  // TODO: get rid of these external files and move into XML, this is too error prone
  // For now, to support the legacy input, we read all the input parameters at once in the arrays below, and then we create the models
  array1d< array1d< string > > phase1InputParams;
  phase1InputParams.resize( 3 );
  array1d< array1d< string > > phase2InputParams;
  phase2InputParams.resize( 3 );

  // 1) Create the viscosity, density, enthalpy models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );

      if( strs[0] == "DensityFun" )
      {
        if( strs[1] == PHASE1::Density::catalogName() )
        {
          phase1InputParams[PHASE1::InputParamOrder::DENSITY] = strs;
        }
        else if( strs[1] == PHASE2::Density::catalogName() )
        {
          phase2InputParams[PHASE2::InputParamOrder::DENSITY] = strs;
        }
      }
      else if( strs[0] == "ViscosityFun" )
      {
        if( strs[1] == PHASE1::Viscosity::catalogName() )
        {
          phase1InputParams[PHASE1::InputParamOrder::VISCOSITY] = strs;
        }
        else if( strs[1] == PHASE2::Viscosity::catalogName() )
        {
          phase2InputParams[PHASE2::InputParamOrder::VISCOSITY] = strs;
        }
      }
      else if( strs[0] == "EnthalpyFun" )
      {
        if( strs[1] == PHASE1::Enthalpy::catalogName() )
        {
          phase1InputParams[PHASE1::InputParamOrder::ENTHALPY] = strs;
        }
        else if( strs[1] == PHASE2::Enthalpy::catalogName() )
        {
          phase2InputParams[PHASE2::InputParamOrder::ENTHALPY] = strs;
        }
      }
      else
      {
        GEOSX_THROW( GEOSX_FMT( "{}: invalid PVT function type '{}'", getFullName(), strs[0] ), InputError );
      }
    }
    is.close();
  }

  // at this point, we have read the file and we check the consistency of non-thermal models
  GEOSX_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::DENSITY].empty(),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Density::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::DENSITY].empty(),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Density::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::VISCOSITY].empty(),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Viscosity::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::VISCOSITY].empty(),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Viscosity::catalogName() ),
                  InputError );

  // we also detect any inconsistency arising in the enthalpy models
  GEOSX_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::ENTHALPY].empty() &&
                  ( PHASE1::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Enthalpy::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::ENTHALPY].empty() &&
                  ( PHASE2::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Enthalpy::catalogName() ),
                  InputError );

  // then, we are ready to instantiate the phase models
  m_phase1 = std::make_unique< PHASE >( getName() + "_phaseModel1", phase1InputParams, m_componentNames, m_componentMolarWeight );
}

template< typename PHASE >
typename ReactiveBrineFluid< PHASE > ::KernelWrapper
ReactiveBrineFluid< PHASE > ::createKernelWrapper()
{
  return KernelWrapper( *m_phase,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView(),
                        m_numPrimarySpecies,
                        *m_equilibriumReactions,
                        *m_kineticReactions,
                        m_primarySpeciesConcentration.toView(),
                        m_secondarySpeciesConcentration.toView(),
                        m_primarySpeciesTotalConcentration.toView(),
                        m_kineticReactionRates.toView() ) );
}

template< typename PHASE >
ReactiveBrineFluid< PHASE > ::KernelWrapper::
  KernelWrapper( PHASE const & phase,
                 arrayView1d< real64 const > componentMolarWeight,
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
                 integer const numPrimarySpecies,
                 chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                 chemicalReactions::KineticReactions const & kineticReactions,
                 arrayView2d< real64 > const & primarySpeciesConcentration,
                 arrayView2d< real64 > const & secondarySpeciesConcentration,
                 arrayView2d< real64 > const & primarySpeciesTotalConcentration,
                 arrayView2d< real64 > const & kineticReactionRates )
  : ReactiveMultiFluid::KernelWrapper( std::move( componentMolarWeight ),
                                       useMass,
                                       std::move( phaseFraction ),
                                       std::move( phaseDensity ),
                                       std::move( phaseMassDensity ),
                                       std::move( phaseViscosity ),
                                       std::move( phaseEnthalpy ),
                                       std::move( phaseInternalEnergy ),
                                       std::move( phaseCompFraction ),
                                       std::move( totalDensity ),
                                       numPrimarySpecies,
                                       equilibriumReactions,
                                       kineticReactions,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration,
                                       primarySpeciesTotalConcentration,
                                       kineticReactionRates ),
  m_isThermal( isThermal ),
  m_phase( phase.createKernelWrapper() )
  {}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class ReactiveBrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction > >;
template class ReactiveBrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy > >;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveBrinePhillipsFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveBrinePhillipsThermalFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
