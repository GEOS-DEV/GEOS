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
 * @file AquiferBoundaryCondition.cpp
 */

#include "AquiferBoundaryCondition.hpp"

#include "mesh/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;

AquiferBoundaryCondition::AquiferBoundaryCondition( string const & name, Group * parent )
  : FieldSpecificationBase( name, parent ),
  m_cumulativeFlux( 0.0 )
{
  registerWrapper( viewKeyStruct::aquiferPorosityString(), &m_porosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer porosity" );

  registerWrapper( viewKeyStruct::aquiferPermeabilityString(), &m_permeability ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer permeability [m^2]" );

  registerWrapper( viewKeyStruct::aquiferInitialPressureString(), &m_initialPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer initial pressure [Pa]" );

  registerWrapper( viewKeyStruct::aquiferWaterViscosityString(), &m_viscosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer water viscosity [Pa.s]" );

  registerWrapper( viewKeyStruct::aquiferWaterDensityString(), &m_density ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer water density [kg.m^-3]" );

  registerWrapper( viewKeyStruct::aquiferTotalCompressibilityString(), &m_totalCompressibility ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer total compressibility (rock and fluid) [Pa^-1]" );

  registerWrapper( viewKeyStruct::aquiferWaterPhaseComponentFractionString(), &m_phaseComponentFraction ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Aquifer water phase component fraction" );

  registerWrapper( viewKeyStruct::aquiferWaterPhaseComponentNamesString(), &m_phaseComponentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Aquifer water phase component names" );

  registerWrapper( viewKeyStruct::aquiferElevationString(), &m_elevation ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer elevation (positive going upward) [m]" );

  registerWrapper( viewKeyStruct::aquiferThicknessString(), &m_thickness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer thickness [m]" );

  registerWrapper( viewKeyStruct::aquiferInnerRadiusString(), &m_innerRadius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer inner radius [m]" );

  registerWrapper( viewKeyStruct::aquiferAngleString(), &m_angle ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Angle subtended by the aquifer boundary from the center of the reservoir [degress]" );

  registerWrapper( viewKeyStruct::pressureInfluenceFunctionNameString(), &m_pressureInfluenceFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the table describing the pressure influence function" );

  registerWrapper( viewKeyStruct::cumulativeFluxString(), &m_cumulativeFlux ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );

}

void AquiferBoundaryCondition::postProcessInput()
{
  GEOSX_THROW_IF( m_permeability <= 0.0,
                  getCatalogName() << " " << getName() << ": the aquifer permeability cannot be equal to zero or negative",
                  InputError );

  GEOSX_THROW_IF( m_pressureInfluenceFunctionName.empty(),
                  getCatalogName() << " " << getName() << ": the pressure influence table name must be specified using the keyword " << viewKeyStruct::pressureInfluenceFunctionNameString(),
                  InputError );

  FunctionManager const & functionManager = FunctionManager::getInstance();
  GEOSX_THROW_IF( !functionManager.hasGroup( m_pressureInfluenceFunctionName ),
                  getCatalogName() << " " << getName() << ": the pressure influence table " << m_pressureInfluenceFunctionName << " could not be found",
                  InputError );

  computeTimeConstant();
  computeInfluxConstant();

  GEOSX_THROW_IF( m_timeConstant <= 0.0,
                  getCatalogName() << " " << getName() << ": the aquifer time constant is equal to zero or negative, the simulation cannot procede",
                  InputError );

  GEOSX_THROW_IF( m_influxConstant <= 0.0,
                  getCatalogName() << " " << getName() << ": the aquifer influx constant is equal to zero or negative, the simulation cannot procede",
                  InputError );

  GEOSX_THROW_IF( m_phaseComponentFraction.size() != m_phaseComponentNames.size(),
                  getCatalogName() << " " << getName() << ": the sizes of "
                                   << viewKeyStruct::aquiferWaterPhaseComponentFractionString() << " and " << viewKeyStruct::aquiferWaterPhaseComponentNamesString()
                                   << " are inconsistent",
                  InputError );

}

void AquiferBoundaryCondition::initializePreSubGroups()
{}

void AquiferBoundaryCondition::computeTimeConstant()
{
  // equation 5.3 of the Eclipse TD
  m_timeConstant = m_viscosity * m_porosity * m_totalCompressibility * m_innerRadius * m_innerRadius;
  m_timeConstant /= m_permeability;
}

void AquiferBoundaryCondition::computeInfluxConstant()
{
  // equation 5.4 of the Eclipse TD
  m_influxConstant = m_thickness * m_angle * m_porosity * m_totalCompressibility * m_innerRadius * m_innerRadius;
}

AquiferBoundaryCondition::KernelWrapper AquiferBoundaryCondition::createKernelWrapper() const
{
  FunctionManager const & functionManager = FunctionManager::getInstance();
  TableFunction const & pressureInfluenceFunction = functionManager.getGroup< TableFunction >( m_pressureInfluenceFunctionName );

  return AquiferBoundaryCondition::KernelWrapper( m_initialPressure,
                                                  m_density,
                                                  m_elevation * m_gravityVector[2],
                                                  m_timeConstant,
                                                  m_influxConstant,
                                                  m_cumulativeFlux,
                                                  pressureInfluenceFunction.createKernelWrapper() );
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, AquiferBoundaryCondition, string const &, Group * const )


} /* namespace geosx */
