/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableCapillaryPressure.cpp
 */

#include "TableCapillaryPressure.hpp"

#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

TableCapillaryPressure::TableCapillaryPressure( std::string const & name,
                                                Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::wettingNonWettingCapPresTableNameString(), &m_wettingNonWettingCapPresTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table [Pa] for the pair (wetting phase, non-wetting phase)\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateCapPresTableNameString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateCapPresTableNameString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::wettingIntermediateCapPresTableNameString(), &m_wettingIntermediateCapPresTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table [Pa] for the pair (wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingCapPresTableNameString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::nonWettingIntermediateCapPresTableNameString(), &m_nonWettingIntermediateCapPresTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Capillary pressure table [Pa] for the pair (non-wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingCapPresTableNameString() ) +
                    " to specify the table names" );

  getWrapperBase( viewKeyStruct::phaseMinVolumeFractionString() )
    .setInputFlag( InputFlags::FALSE );

  registerWrapper( viewKeyStruct::capPresWrappersString(), &m_capPresKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::inverseCapPresWrappersString(), &m_inverseCapPresWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TableCapillaryPressure::postInputInitialization()
{
  CapillaryPressureBase::postInputInitialization();

  integer const numPhases = m_phaseNames.size();
  GEOS_THROW_IF( numPhases != 2 && numPhases != 3,
                 GEOS_FMT( "{}: the expected number of fluid phases is either two, or three",
                           getFullName() ),
                 InputError, getDataContext() );

  // Populate the minimum phase volumes
  m_phaseMinVolumeFraction.resize( numPhases );

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( m_wettingNonWettingCapPresTableName.empty(),
                   GEOS_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the capillary pressure table for the pair (wetting phase, non-wetting phase)",
                             getFullName(),
                             viewKeyStruct::wettingNonWettingCapPresTableNameString() ),
                   InputError, getDataContext() );
  }
  else if( numPhases == 3 )
  {
    GEOS_THROW_IF( m_wettingIntermediateCapPresTableName.empty() || m_nonWettingIntermediateCapPresTableName.empty(),
                   GEOS_FMT( "{}: for a three-phase flow simulation, we must use {} to specify the capillary pressure table "
                             "for the pair (wetting phase, intermediate phase), and {} to specify the capillary pressure table "
                             "for the pair (non-wetting phase, intermediate phase)",
                             getFullName(),
                             viewKeyStruct::wettingIntermediateCapPresTableNameString(),
                             viewKeyStruct::nonWettingIntermediateCapPresTableNameString()  ),
                   InputError, getDataContext() );
  }
}

void TableCapillaryPressure::initializePreSubGroups()
{
  CapillaryPressureBase::initializePreSubGroups();

  integer const numPhases = m_phaseNames.size();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_wettingNonWettingCapPresTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_wettingNonWettingCapPresTableName ),
                   InputError );
    TableFunction const & capPresTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingCapPresTableName );
    bool const capPresMustBeIncreasing = ( m_phaseOrder[PhaseType::WATER] < 0 )
      ? true   // pc on the gas phase, function must be increasing
      : false; // pc on the water phase, function must be decreasing
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(), capPresMustBeIncreasing );
  }
  else if( numPhases == 3 )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_wettingIntermediateCapPresTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_wettingIntermediateCapPresTableName ),
                   InputError, getDataContext() );
    TableFunction const & capPresTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateCapPresTableName );
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableWI, getFullName(), false );

    GEOS_THROW_IF( !functionManager.hasGroup( m_nonWettingIntermediateCapPresTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_nonWettingIntermediateCapPresTableName ),
                   InputError, getDataContext() );
    TableFunction const & capPresTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateCapPresTableName );
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableNWI, getFullName(), true );
  }
}

void TableCapillaryPressure::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_capPresKernelWrappers.clear();
  m_inverseCapPresWrappers.clear();

  if( numPhases == 2 )
  {
    TableFunction const & capPresTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTable.createKernelWrapper() );

    auto const & satArrayView = capPresTable.getCoordinates()[0];
    auto const & capPresArrayView   = capPresTable.getValues();

    std::vector< real64 > satVec( satArrayView.size() );
    std::vector< real64 > pcVec( capPresArrayView.size() );

    std::copy( satArrayView.begin(), satArrayView.end(), satVec.begin() );
    std::copy( capPresArrayView.begin(), capPresArrayView.end(), pcVec.begin() );

    // Reverse both arrays (if original J is decreasing in S)
    std::reverse( pcVec.begin(), pcVec.end() );
    std::reverse( satVec.begin(), satVec.end() );


    auto inverseTable = std::make_shared< TableFunction >( "inverseCapPres", this );

    real64_array invPcVec( pcVec.size() );
    real64_array invSatVec( satVec.size() );
    std::copy( pcVec.begin(), pcVec.end(), invPcVec.data() );
    std::copy( satVec.begin(), satVec.end(), invSatVec.data() );

    array1d< real64_array > coordinates;
    coordinates.emplace_back( std::move( invPcVec ) );


    std::vector< units::Unit > dimUnits = { units::Unknown }; // or actual unit if available

    inverseTable->setTableCoordinates( coordinates, dimUnits );
    inverseTable->setTableValues( std::move( invSatVec ), units::Unknown );
    inverseTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );

    m_inverseCapPresWrappers.emplace_back( inverseTable->createKernelWrapper() );
    m_inverseTables.emplace_back( std::move( inverseTable ) );

    // Populate the end-points from the tables
    TableCapillaryPressureHelpers::populateMinPhaseVolumeFraction( m_phaseOrder.toSliceConst(), capPresTable, m_phaseMinVolumeFraction );
  }
  else if( numPhases == 3 )
  {
    TableFunction const & capPresTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTableWI.createKernelWrapper() );
    m_inverseCapPresWrappers.emplace_back( capPresTableWI.createKernelWrapper() );
    TableFunction const & capPresTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateCapPresTableName );
    m_capPresKernelWrappers.emplace_back( capPresTableNWI.createKernelWrapper() );
    m_inverseCapPresWrappers.emplace_back( capPresTableNWI.createKernelWrapper() );

    // Populate the end-points from the tables
    TableCapillaryPressureHelpers::populateMinPhaseVolumeFraction( m_phaseOrder.toSliceConst(), capPresTableWI, capPresTableNWI, m_phaseMinVolumeFraction );
  }
}


TableCapillaryPressure::KernelWrapper::
  KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
                 arrayView1d< TableFunction::KernelWrapper const > const & inverseCapPresWrappers,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrapped,
                 arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                 arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac )
  : CapillaryPressureBaseUpdate( phaseTypes,
                                 phaseOrder,
                                 phaseTrapped,
                                 phaseCapPres,
                                 dPhaseCapPres_dPhaseVolFrac ),
  m_capPresKernelWrappers( capPresKernelWrappers ),
  m_inverseCapPresWrappers( inverseCapPresWrappers )
{}

TableCapillaryPressure::KernelWrapper
TableCapillaryPressure::createKernelWrapper()
{
  createAllTableKernelWrappers();
  return KernelWrapper( m_capPresKernelWrappers,
                        m_inverseCapPresWrappers,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseTrappedVolFrac,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableCapillaryPressure, std::string const &, Group * const )

} // namespace constitutive

} // namespace geos
