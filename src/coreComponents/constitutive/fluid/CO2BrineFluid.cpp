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
 * @file CO2BrineFluid.cpp
 */
#include "CO2BrineFluid.hpp"

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
template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                        PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                        PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsThermalFluid"; }
};

template<> class
  TwoPhaseCatalogNames< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                        PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                        PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiThermalFluid"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
string CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                      P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                      FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                               FLASH >::name();
}

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::
CO2BrineFluid( string const & name, Group * const parent ):
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

  m_thermalFlag = ( P1ENTH::catalogName()      != PVTProps::NoOpPVTFunction::catalogName() &&
                    P1INTENERGY::catalogName() != PVTProps::NoOpPVTFunction::catalogName() &&
                    P2ENTH::catalogName()      != PVTProps::NoOpPVTFunction::catalogName() &&
                    P2INTENERGY::catalogName() != PVTProps::NoOpPVTFunction::catalogName() );

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if( m_thermalFlag )
  {
    getExtrinsicData< extrinsicMeshData::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getExtrinsicData< extrinsicMeshData::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }
}

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
std::unique_ptr< ConstitutiveBase >
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  CO2BrineFluid & newConstitutiveRelation = dynamicCast< CO2BrineFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
integer CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                       P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                       FLASH >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}


template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
void CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                    P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                    FLASH >::postProcessInput()
{
  MultiFluidBase::postProcessInput();

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

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
void CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                    P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                    FLASH >::createPVTModels()
{

  // TODO: get rid of these external files and move into XML, this is too error prone

  // 1) Create the viscosity and density models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );

      if( strs[0] == "DensityFun" )
      {
        if( strs[1] == P1DENS::catalogName() )
        {
          m_p1Density = std::make_unique< P1DENS >( getName() + '_' + P1DENS::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2DENS::catalogName() )
        {
          m_p2Density = std::make_unique< P2DENS >( getName() + '_' + P2DENS::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else if( strs[0] == "ViscosityFun" )
      {
        if( strs[1] == P1VISC::catalogName() )
        {
          m_p1Viscosity = std::make_unique< P1VISC >( getName() + '_' + P1VISC::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2VISC::catalogName() )
        {
          m_p2Viscosity = std::make_unique< P2VISC >( getName() + '_' + P2VISC::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else if( strs[0] == "EnthalpyFun" )
      {
        if( strs[1] == P1ENTH::catalogName() )
        {
          m_p1Enthalpy = std::make_unique< P1ENTH >( getName() + '_' + P1ENTH::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2ENTH::catalogName() )
        {
          m_p2Enthalpy = std::make_unique< P2ENTH >( getName() + '_' + P2ENTH::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else if( strs[0] == "InternalEnergyFun" )
      {
        if( strs[1] == P1INTENERGY::catalogName() )
        {
          m_p1IntEnergy = std::make_unique< P1INTENERGY >( getName() + '_' + P1INTENERGY::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2INTENERGY::catalogName() )
        {
          m_p2IntEnergy = std::make_unique< P2INTENERGY >( getName() + '_' + P2INTENERGY::catalogName(), strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else
      {
        GEOSX_THROW( GEOSX_FMT( "{}: invalid PVT function type '{}'", getFullName(), strs[0] ), InputError );
      }
    }
    is.close();
  }

  // at this point, we have read the file, and we can detect any inconsistency arising in the enthalpy/internal energy models
  GEOSX_THROW_IF( m_p1Enthalpy == nullptr && ( P1ENTH::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P1ENTH::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p1IntEnergy == nullptr && ( P1INTENERGY::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P1INTENERGY::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p2Enthalpy == nullptr && ( P2ENTH::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P2ENTH::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p2IntEnergy == nullptr && ( P2INTENERGY::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P2INTENERGY::catalogName() ),
                  InputError );
  // we also check the consistency of non-thermal models
  GEOSX_THROW_IF( m_p1Density == nullptr,
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P1DENS::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p2Density == nullptr,
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P2DENS::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p1Viscosity == nullptr,
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P1VISC::catalogName() ),
                  InputError );
  GEOSX_THROW_IF( m_p2Viscosity == nullptr,
                  GEOSX_FMT( "{}: PVT model {} not found in input files", getFullName(), P2VISC::catalogName() ),
                  InputError );

  // for non-thermal models, the user skips the definition of the enthalpy/internal energy models, so we create NoOp models now
  if( m_p1Enthalpy == nullptr )
  {
    m_p1Enthalpy = std::make_unique< P1ENTH >( getName() + '_' + P1ENTH::catalogName(), string_array(), m_componentNames, m_componentMolarWeight );
  }
  if( m_p2Enthalpy == nullptr )
  {
    m_p2Enthalpy = std::make_unique< P2ENTH >( getName() + '_' + P2ENTH::catalogName(), string_array(), m_componentNames, m_componentMolarWeight );
  }
  if( m_p1IntEnergy == nullptr )
  {
    m_p1IntEnergy = std::make_unique< P1INTENERGY >( getName() + '_' + P1INTENERGY::catalogName(), string_array(), m_componentNames, m_componentMolarWeight );
  }
  if( m_p2IntEnergy == nullptr )
  {
    m_p2IntEnergy = std::make_unique< P2INTENERGY >( getName() + '_' + P2INTENERGY::catalogName(), string_array(), m_componentNames, m_componentMolarWeight );
  }

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
          m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(), strs, m_phaseNames, m_componentNames, m_componentMolarWeight );
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

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
typename CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
                        P2DENS, P2VISC, P2ENTH, P2INTENERGY,
                        FLASH >::KernelWrapper
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::createKernelWrapper()
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_p1Density,
                        *m_p1Viscosity,
                        *m_p1Enthalpy,
                        *m_p1IntEnergy,
                        *m_p2Density,
                        *m_p2Viscosity,
                        *m_p2Enthalpy,
                        *m_p2IntEnergy,
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

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::KernelWrapper::
  KernelWrapper( integer const p1Index,
                 integer const p2Index,
                 P1DENS const & p1Density,
                 P1VISC const & p1Viscosity,
                 P1ENTH const & p1Enthalpy,
                 P1INTENERGY const & p1IntEnergy,
                 P2DENS const & p2Density,
                 P2VISC const & p2Viscosity,
                 P2ENTH const & p2Enthalpy,
                 P2INTENERGY const & p2IntEnergy,
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
  m_p1Density( p1Density.createKernelWrapper() ),
  m_p1Viscosity( p1Viscosity.createKernelWrapper() ),
  m_p1Enthalpy( p1Enthalpy.createKernelWrapper() ),
  m_p1IntEnergy( p1IntEnergy.createKernelWrapper() ),
  m_p2Density( p2Density.createKernelWrapper() ),
  m_p2Viscosity( p2Viscosity.createKernelWrapper() ),
  m_p2Enthalpy( p2Enthalpy.createKernelWrapper() ),
  m_p2IntEnergy( p2IntEnergy.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class CO2BrineFluid< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                              PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                              PVTProps::CO2Solubility >;
template class CO2BrineFluid< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                              PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                              PVTProps::CO2Solubility >;

template class CO2BrineFluid< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                              PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                              PVTProps::CO2Solubility >;
template class CO2BrineFluid< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                              PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                              PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsThermalFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiThermalFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
