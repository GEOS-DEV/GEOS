/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

#include "common/Path.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/FlashModelBase.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionBase.hpp"
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
  TwoPhaseCatalogNames< PVTProps::BrineCO2Density,
                        PVTProps::BrineViscosity,
                        PVTProps::SpanWagnerCO2Density,
                        PVTProps::FenghourCO2Viscosity,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineFluid"; }
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

  MultiPhaseMultiComponentFluid * const newConstitutiveRelation = dynamic_cast< MultiPhaseMultiComponentFluid * >(clone.get());
  newConstitutiveRelation->m_p1Index = this->m_p1Index;
  newConstitutiveRelation->m_p2Index = this->m_p2Index;

  newConstitutiveRelation->createPVTModels();

  return clone;
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
void MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  GEOSX_THROW_IF_NE_MSG( numFluidPhases(), 2,
                         "The number of phases in this model should be equal to 2",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( numFluidComponents(), 2,
                         "The number of components in this model should be equal to 2",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( m_phasePVTParaFiles.size(), 2,
                         "The number of phasePVTParaFiles is not the same as the number of phases!",
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
    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ) )
    {
      string const str( buf );
      string_array const strs = stringutilities::tokenize( str, " " );

      if( strs[0] == "DensityFun" )
      {
        if( strs[1] == P1DENS::catalogName() )
        {
          m_p1Density = std::make_unique< P1DENS >( strs, m_componentNames, m_componentMolarWeight );
          m_p1DensityWrapper.emplace_back( m_p1Density->createKernelWrapper() );
        }
        else if( strs[1] == P2DENS::catalogName() )
        {
          m_p2Density = std::make_unique< P2DENS >( strs, m_componentNames, m_componentMolarWeight );
          m_p2DensityWrapper.emplace_back( m_p2Density->createKernelWrapper() );
        }
      }
      else if( strs[0] == "ViscosityFun" )
      {
        if( strs[1] == P1VISC::catalogName() )
        {
          m_p1Viscosity = std::make_unique< P1VISC >( strs, m_componentNames, m_componentMolarWeight );
          m_p1ViscosityWrapper.emplace_back( m_p1Viscosity->createKernelWrapper() );
        }
        else if( strs[1] == P2VISC::catalogName() )
        {
          m_p2Viscosity = std::make_unique< P2VISC >( strs, m_componentNames, m_componentMolarWeight );
          m_p2ViscosityWrapper.emplace_back( m_p2Viscosity->createKernelWrapper() );
        }
      }
      else
      {
        GEOSX_THROW( "Error: Invalid PVT function: " << strs[0] << ".", InputError );
      }
    }
    is.close();
  }

  GEOSX_THROW_IF( m_p1Density == nullptr, P1DENS::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p2Density == nullptr, P2DENS::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p1Viscosity == nullptr, P1VISC::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p2Viscosity == nullptr, P2VISC::catalogName() << " model not found", InputError );

  // 2) Create the flash model
  {
    std::ifstream is( m_flashModelParaFile );
    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ) )
    {
      string const str( buf );
      string_array const strs = stringutilities::tokenize( str, " " );
      if( strs[0] == "FlashModel" && strs[1] == FLASH::catalogName() )
      {
        m_flash = std::make_unique< FLASH >( strs, m_phaseNames, m_componentNames, m_componentMolarWeight );
        m_flashWrapper.emplace_back( m_flash->createKernelWrapper() );
      }
      else
      {
        GEOSX_THROW( "Error: Invalid flash model: " << strs[0] << ", " << strs[1] << ".", InputError );
      }
    }
    is.close();
  }

  GEOSX_THROW_IF( m_flash == nullptr,
                  FLASH::catalogName() << " not found",
                  InputError );
}

// explicit instantiation of the model template; unfortunately we can't use CO2BrineFluid alias for this
template class MultiPhaseMultiComponentFluid< PVTProps::BrineCO2Density,
                                              PVTProps::BrineViscosity,
                                              PVTProps::SpanWagnerCO2Density,
                                              PVTProps::FenghourCO2Viscosity,
                                              PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
