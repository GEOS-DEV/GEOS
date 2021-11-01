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
 * @file MultiPhaseMultiComponentFluid.cpp
 */
#include "MultiPhaseMultiComponentFluid.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PVTProps::PhillipsBrineDensity,
                        PVTProps::PhillipsBrineViscosity,
                        PVTProps::SpanWagnerCO2Density,
                        PVTProps::FenghourCO2Viscosity,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PVTProps::EzrokhiBrineDensity,
                        PVTProps::EzrokhiBrineViscosity,
                        PVTProps::SpanWagnerCO2Density,
                        PVTProps::FenghourCO2Viscosity,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "EzrokhiCO2BrineFluid"; }
};
} // end namespace

// provide a definition for catalogName()
template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
string MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::name();
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::
MultiPhaseMultiComponentFluid( string const & name, Group * const parent ):
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

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
std::unique_ptr< ConstitutiveBase >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  MultiPhaseMultiComponentFluid & newConstitutiveRelation = dynamicCast< MultiPhaseMultiComponentFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
integer MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}


template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
void MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::postProcessInput()
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

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
void MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::createPVTModels()
{
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
      else
      {
        GEOSX_THROW( GEOSX_FMT( "{}: invalid PVT function type '{}'", getFullName(), strs[0] ), InputError );
      }
    }
    is.close();
  }

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

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
typename MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::createKernelWrapper()
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_p1Density,
                        *m_p1Viscosity,
                        *m_p2Density,
                        *m_p2Viscosity,
                        *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper::
  KernelWrapper( localIndex const p1Index,
                 localIndex const p2Index,
                 P1DENS const & p1DensityWrapper,
                 P1VISC const & p1ViscosityWrapper,
                 P2DENS const & p2DensityWrapper,
                 P2VISC const & p2ViscosityWrapper,
                 FLASH const & flashWrapper,
                 arrayView1d< geosx::real64 const > componentMolarWeight,
                 bool const useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_p1Index( p1Index ),
  m_p2Index( p2Index ),
  m_p1Density( p1DensityWrapper.createKernelWrapper() ),
  m_p1Viscosity( p1ViscosityWrapper.createKernelWrapper() ),
  m_p2Density( p2DensityWrapper.createKernelWrapper() ),
  m_p2Viscosity( p2ViscosityWrapper.createKernelWrapper() ),
  m_flash( flashWrapper.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use CO2BrinePhillipsFluid alias for this
template class MultiPhaseMultiComponentFluid< PVTProps::PhillipsBrineDensity,
                                              PVTProps::PhillipsBrineViscosity,
                                              PVTProps::SpanWagnerCO2Density,
                                              PVTProps::FenghourCO2Viscosity,
                                              PVTProps::CO2Solubility >;

template class MultiPhaseMultiComponentFluid< PVTProps::EzrokhiBrineDensity,
                                              PVTProps::EzrokhiBrineViscosity,
                                              PVTProps::SpanWagnerCO2Density,
                                              PVTProps::FenghourCO2Viscosity,
                                              PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EzrokhiCO2BrineFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
