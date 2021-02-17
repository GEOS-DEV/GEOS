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
 * @file NewMultiPhaseMultiComponentFluid.cpp
 */
#include "NewMultiPhaseMultiComponentFluid.hpp"

#include "common/Path.hpp"
#include "managers/ProblemManager.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace PVTProps;
using namespace stringutilities;

namespace constitutive
{

NewMultiPhaseMultiComponentFluid::NewMultiPhaseMultiComponentFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::phasePVTParaFilesString, &m_phasePVTParaFiles )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  registerWrapper( viewKeyStruct::flashModelParaFileString, &m_flashModelParaFile )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Name of the file defining the parameters of the flash model" );
}

NewMultiPhaseMultiComponentFluid::~NewMultiPhaseMultiComponentFluid()
{}


std::unique_ptr< ConstitutiveBase >
NewMultiPhaseMultiComponentFluid::deliverClone( string const & name,
                                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  NewMultiPhaseMultiComponentFluid * const newConstitutiveRelation = dynamic_cast< NewMultiPhaseMultiComponentFluid * >(clone.get());
  newConstitutiveRelation->m_phaseGasIndex = this->m_phaseGasIndex;
  newConstitutiveRelation->m_phaseLiquidIndex = this->m_phaseLiquidIndex;

  newConstitutiveRelation->createPVTModels();

  return clone;
}

void NewMultiPhaseMultiComponentFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();
  localIndex const numComps = numFluidComponents();
  GEOSX_ERROR_IF( numPhases != 2, "The number of phases in this model should be equal to 2" );
  GEOSX_ERROR_IF( numComps != 2, "The number of components in this model should be equal to 2" );
  GEOSX_ERROR_IF( m_phasePVTParaFiles.size() != 2, "The number of phasePVTParaFiles is not the same as the number of phases!" );

  bool notFound = true;
  for( localIndex i = 0; i < m_phaseNames.size(); ++i )
  {
    if( m_phaseNames[i] == "CO2" || m_phaseNames[i] == "co2" ||
        m_phaseNames[i] == "gas" || m_phaseNames[i] == "Gas" )
    {
      m_phaseGasIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Phase co2/gas is not found!" );

  notFound = true;
  for( localIndex i = 0; i < m_phaseNames.size(); ++i )
  {
    if( m_phaseNames[i] == "Water" || m_phaseNames[i] == "water" ||
        m_phaseNames[i] == "Liquid" || m_phaseNames[i] == "liquid" )
    {
      m_phaseLiquidIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Phase water/liquid is not found!" );

  createPVTModels();
}

void NewMultiPhaseMultiComponentFluid::createPVTModels()
{
  // 1) Create the viscosity and density models
  for( string & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ) )
    {
      string const str( buf );
      string_array const strs = Tokenize( str, " " );
      if( strs[0] == "DensityFun" )
      {
        if( strs[1] == "SpanWagnerCO2Density" )
        {
          m_co2Density = std::make_unique< SpanWagnerCO2Density >( strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == "BrineCO2Density" )
        {
          m_brineDensity = std::make_unique< BrineCO2Density >( strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else if( strs[0] == "ViscosityFun" )
      {
        if( strs[1] == "FenghourCO2Viscosity" )
        {
          m_co2Viscosity = std::make_unique< FenghourCO2Viscosity >( strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == "BrineViscosity" )
        {
          m_brineViscosity = std::make_unique< BrineViscosity >( strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else
      {
        GEOSX_ERROR( "Error: Invalid PVT function: " << strs[0] << "." );
      }
    }
    is.close();
  }

  GEOSX_ERROR_IF( m_co2Density == nullptr, "SpanWagnerCO2Density model not found" );
  GEOSX_ERROR_IF( m_brineDensity == nullptr, "BrineCO2Density model not found" );
  GEOSX_ERROR_IF( m_co2Viscosity == nullptr, "FenghourCO2Viscosity model not found" );
  GEOSX_ERROR_IF( m_brineViscosity == nullptr, "BrineViscosity model not found" );

  // 2) Create the flash model
  {
    std::ifstream is( m_flashModelParaFile );
    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ) )
    {
      string const str( buf );
      string_array const strs = Tokenize( str, " " );
      if( strs[0] == "FlashModel" && strs[1] == "CO2Solubility" )
      {
        m_co2Solubility = std::make_unique< CO2Solubility >( strs, m_phaseNames, m_componentNames, m_componentMolarWeight );
      }
      else
      {
        GEOSX_ERROR( "Error: Invalid flash model: " << strs[0] << ", " << strs[1] << "." );
      }
    }
    is.close();
  }

  GEOSX_ERROR_IF( m_co2Solubility == nullptr, "CO2Solubility model not found" );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, NewMultiPhaseMultiComponentFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
