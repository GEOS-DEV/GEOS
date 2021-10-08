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
 * @file TableCapillaryPressure.cpp
 */

#include "TableCapillaryPressure.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableCapillaryPressure::TableCapillaryPressure( std::string const & name,
                                                Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::wettingNonWettingCapPresTableNameString(), &m_wettingNonWettingCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table for the pair (wetting phase, non-wetting phase)\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateCapPresTableNameString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateCapPresTableNameString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::wettingIntermediateCapPresTableNameString(), &m_wettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table for the pair (wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingCapPresTableNameString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::nonWettingIntermediateCapPresTableNameString(), &m_nonWettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table for the pair (non-wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingCapPresTableNameString() ) +
                    " to specify the table names" );

  registerWrapper( "capPresWrappers", &m_capPresKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TableCapillaryPressure::postProcessInput()
{
  CapillaryPressureBase::postProcessInput();

  localIndex const numPhases = m_phaseNames.size();
  GEOSX_THROW_IF( numPhases != 2 && numPhases != 3,
                  GEOSX_FMT( "{}: the expected number of fluid phases is either two, or three",
                             getFullName() ),
                  InputError );

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( m_wettingNonWettingCapPresTableName.empty(),
                    GEOSX_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the capillary pressure table for the pair (wetting phase, non-wetting phase)",
                               getFullName(),
                               viewKeyStruct::wettingNonWettingCapPresTableNameString() ),
                    InputError );
  }
  else if( numPhases == 3 )
  {
    GEOSX_THROW_IF( m_wettingIntermediateCapPresTableName.empty() || m_nonWettingIntermediateCapPresTableName.empty(),
                    GEOSX_FMT(
                      "{}: for a three-phase flow simulation, we must use {} to specify the capillary pressure table for the pair (wetting phase, intermediate phase), and {} to specify the capillary pressure table for the pair (non-wetting phase, intermediate phase)",
                      getFullName(),
                      viewKeyStruct::wettingIntermediateCapPresTableNameString(),
                      viewKeyStruct::nonWettingIntermediateCapPresTableNameString()  ),
                    InputError );
  }
}

void TableCapillaryPressure::initializePreSubGroups()
{
  CapillaryPressureBase::initializePreSubGroups();

  localIndex const numPhases = m_phaseNames.size();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( !functionManager.hasGroup( m_wettingNonWettingCapPresTableName ),
                    GEOSX_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_wettingNonWettingCapPresTableName ),
                    InputError );
    TableFunction const & capPresTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingCapPresTableName );
    bool const capPresMustBeIncreasing = ( m_phaseOrder[PhaseType::WATER] < 0 )
      ? true   // pc on the gas phase, function must be increasing
      : false; // pc on the water phase, function must be decreasing
    validateCapillaryPressureTable( capPresTable, capPresMustBeIncreasing );
  }
  else if( numPhases == 3 )
  {
    GEOSX_THROW_IF( !functionManager.hasGroup( m_wettingIntermediateCapPresTableName ),
                    GEOSX_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_wettingIntermediateCapPresTableName ),
                    InputError );
    TableFunction const & capPresTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateCapPresTableName );
    validateCapillaryPressureTable( capPresTableWI, false );

    GEOSX_THROW_IF( !functionManager.hasGroup( m_nonWettingIntermediateCapPresTableName ),
                    GEOSX_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_nonWettingIntermediateCapPresTableName ),
                    InputError );
    TableFunction const & capPresTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateCapPresTableName );
    validateCapillaryPressureTable( capPresTableNWI, true );
  }
}

void TableCapillaryPressure::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  localIndex const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_capPresKernelWrappers.clear();
  if( numPhases == 2 )
  {
    TableFunction const & capPresTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTable.createKernelWrapper() );
  }
  else if( numPhases == 3 )
  {
    TableFunction const & capPresTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTableWI.createKernelWrapper() );
    TableFunction const & capPresTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTableNWI.createKernelWrapper() );
  }
}

void TableCapillaryPressure::validateCapillaryPressureTable( TableFunction const & capPresTable,
                                                             bool const capPresMustBeIncreasing ) const
{
  ArrayOfArraysView< real64 const > coords = capPresTable.getCoordinates();

  GEOSX_THROW_IF_NE_MSG( capPresTable.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                         GEOSX_FMT( "{}: in table '{}' interpolation method must be linear", getFullName(), capPresTable.getName() ),
                         InputError );
  GEOSX_THROW_IF_NE_MSG( capPresTable.numDimensions(), 1,
                         GEOSX_FMT( "{}: table '{}' must have a single independent coordinate", getFullName(), capPresTable.getName() ),
                         InputError );
  GEOSX_THROW_IF_LT_MSG( coords.sizeOfArray( 0 ), 2,
                         GEOSX_FMT( "{}: table `{}` must contain at least two values", getFullName(), capPresTable.getName() ),
                         InputError );

  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const capPres = capPresTable.getValues();

  for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
  {
    // check phase volume fraction
    GEOSX_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                    GEOSX_FMT( "{}: in table '{}' values must be between 0 and 1", getFullName(), capPresTable.getName() ),
                    InputError );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check the monotonicity of the capillary pressure table
    if( capPresMustBeIncreasing )
    {
      GEOSX_THROW_IF( !isZero( capPres[i] ) && (capPres[i] - capPres[i-1]) < -1e-10,
                      GEOSX_FMT( "{}: in table '{}' values must be strictly increasing", getFullName(), capPresTable.getName() ),
                      InputError );
    }
    else
    {
      GEOSX_THROW_IF( !isZero( capPres[i] ) && (capPres[i] - capPres[i-1]) > 1e-10,
                      GEOSX_FMT( "{}: in table '{}' values must be strictly decreasing", getFullName(), capPresTable.getName() ),
                      InputError );
    }

  }
}

TableCapillaryPressure::KernelWrapper::
  KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                 arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac )
  : CapillaryPressureBaseUpdate( phaseTypes,
                                 phaseOrder,
                                 phaseCapPres,
                                 dPhaseCapPres_dPhaseVolFrac ),
  m_capPresKernelWrappers( capPresKernelWrappers )
{}

TableCapillaryPressure::KernelWrapper
TableCapillaryPressure::createKernelWrapper()
{
  createAllTableKernelWrappers();
  return KernelWrapper( m_capPresKernelWrappers,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableCapillaryPressure, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
