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
 * @file ReactiveMultiFluid.cpp
 */
#include "ReactiveMultiFluid.hpp"
#include "ReactiveMultiFluidExtrinsicData.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename PHASE1, typename PHASE2, typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "ReactiveCO2BrinePhillipsFluid"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename PHASE1, typename PHASE2, typename FLASH >
string ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< PHASE1, PHASE2, FLASH >::name();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::
ReactiveMultiFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::phasePVTParaFilesString(), &m_phasePVTParaFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  registerWrapper( viewKeyStruct::flashModelParaFileString(), &m_flashModelParaFile ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Name of the file defining the parameters of the flash model" );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
bool ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::isThermal() const
{
  return ( PHASE1::Enthalpy::catalogName()       != PVTProps::NoOpPVTFunction::catalogName() &&
           PHASE1::InternalEnergy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() &&
           PHASE2::Enthalpy::catalogName()       != PVTProps::NoOpPVTFunction::catalogName() &&
           PHASE2::InternalEnergy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() );
}


template< typename PHASE1, typename PHASE2, typename FLASH >
std::unique_ptr< ConstitutiveBase >
ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  ReactiveMultiFluid & newConstitutiveRelation = dynamicCast< ReactiveMultiFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename PHASE1, typename PHASE2, typename FLASH >
integer ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}


template< typename PHASE1, typename PHASE2, typename FLASH >
void ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  GEOSX_THROW_IF_NE_MSG( numFluidPhases(), 2,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
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

template< typename PHASE1, typename PHASE2, typename FLASH >
void ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::createPVTModels()
{

  // TODO: get rid of these external files and move into XML, this is too error prone
  // For now, to support the legacy input, we read all the input parameters at once in the arrays below, and then we create the models
  array1d< array1d< string > > phase1InputParams;
  phase1InputParams.resize( 4 );
  array1d< array1d< string > > phase2InputParams;
  phase2InputParams.resize( 4 );

  // 1) Create the viscosity, density, enthalpy, and internal energy models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );

      if( strs[0] == toString( SubModelInputNames::DENSITY ) )
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
      else if( strs[0] == toString( SubModelInputNames::VISCOSITY ) )
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

  // then, we are ready to instantiate the phase models
  m_phase1 = std::make_unique< PHASE1 >( getName() + "_phaseModel1", phase1InputParams, m_componentNames, m_componentMolarWeight );
  m_phase2 = std::make_unique< PHASE2 >( getName() + "_phaseModel2", phase2InputParams, m_componentNames, m_componentMolarWeight );

  // 2) Create the flash model
  {
    std::ifstream is( m_flashModelParaFile );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );
      if( strs[0] == "FlashModel" )
      {
        if( strs[1] == FLASH::catalogName() )
        {
          m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                               strs,
                                               m_phaseNames,
                                               m_componentNames,
                                               m_componentMolarWeight );
        }
      }
      else
      {
        GEOSX_THROW( GEOSX_FMT( "{}: invalid flash model type '{}'", getFullName(), strs[0] ), InputError );
      }
    }
    is.close();
  }

  GEOSX_THROW_IF( m_flash == nullptr,
                  GEOSX_FMT( "{}: flash model {} not found in input files", getFullName(), FLASH::catalogName() ),
                  InputError );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
typename ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::KernelWrapper
ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::createKernelWrapper()
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_phase1,
                        *m_phase2,
                        *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
ReactiveMultiFluid< PHASE1, PHASE2, FLASH >::KernelWrapper::
  KernelWrapper( integer const p1Index,
                 integer const p2Index,
                 PHASE1 const & phase1,
                 PHASE2 const & phase2,
                 FLASH const & flash,
                 arrayView1d< geosx::real64 const > componentMolarWeight,
                 bool const useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_p1Index( p1Index ),
  m_p2Index( p2Index ),
  m_phase1( phase1.createKernelWrapper() ),
  m_phase2( phase2.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class ReactiveMultiFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction >,
                                   PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction >,
                                   PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveCO2BrinePhillipsFluid, string const &, Group * const )
} //namespace constitutive

} //namespace geosx
