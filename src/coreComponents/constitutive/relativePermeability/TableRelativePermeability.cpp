/*
   1;5202;0c * ------------------------------------------------------------------------------------------------------------
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
 * @file TableRelativePermeability.cpp
 */

#include "TableRelativePermeability.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableRelativePermeability::TableRelativePermeability( std::string const & name,
                                                      Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::waterOilRelPermTableNamesString(), &m_waterOilRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (water phase, oil phase)\n"
                    "The expected format is \"{ waterPermTableName, oilPermTableName }\", in that order" );

  registerWrapper( viewKeyStruct::gasOilRelPermTableNamesString(), &m_gasOilRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (gas phase, oil phase)\n"
                    "The expected format is \"{ gasPermTableName, oilPermTableName }\", in that order" );

  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setSizedFromParent( 0 );

  registerWrapper( "waterOilRelPermTableWrappers", &m_waterOilRelPermTableKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( "gasOilRelPermTableWrappers", &m_gasOilRelPermTableKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TableRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::OIL] < 0,
                  getFullName() << ": reference oil phase has not been defined and must be included in model",
                  InputError );

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::WATER] >= 0 && !(m_waterOilRelPermTableNames.size() == 2),
                  getFullName() << ": since water is present you must define two tables for the oil-water relperms: "
                                   "one for the oil phase and one for the water phase",
                  InputError );

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::GAS] >= 0 && !(m_gasOilRelPermTableNames.size() == 2),
                  getFullName() << ": since gas is present you must define two tables for the oil-gas relperms: "
                                   "one for the oil phase and one for the gas phase",
                  InputError );

}

void TableRelativePermeability::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();
  createAllTableKernelWrappers();
}

void TableRelativePermeability::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  m_phaseMinVolumeFraction.resize( MAX_NUM_PHASES );

  if( m_waterOilRelPermTableKernelWrappers.empty() && m_gasOilRelPermTableKernelWrappers.empty() )
  {

    // check water-oil relperms
    for( localIndex ip = 0; ip < m_waterOilRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_waterOilRelPermTableNames[ip] );
      real64 const minVolPhaseFrac = validateRelativePermeabilityTable( relPermTable );
      if( ip == 0 ) // water
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::WATER]] = minVolPhaseFrac;
      }
      else if( ip == 1 ) // oil
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::OIL]] = minVolPhaseFrac;
      }
      else
      {
        GEOSX_THROW( getFullName() << ": there should be only two table names for the water-oil pair", InputError );
      }
      m_waterOilRelPermTableKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }

    // check gas-oil relperms
    for( localIndex ip = 0; ip < m_gasOilRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_gasOilRelPermTableNames[ip] );
      real64 const minVolPhaseFrac = validateRelativePermeabilityTable( relPermTable );
      if( ip == 0 ) // gas
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::GAS]] = minVolPhaseFrac;
      }
      else if( ip == 1 ) // oil
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::OIL]] = minVolPhaseFrac;
      }
      else
      {
        GEOSX_THROW( getFullName() << ": there should be only two table names for the gas-oil pair", InputError );
      }

      m_gasOilRelPermTableKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }
  }
}

real64 TableRelativePermeability::validateRelativePermeabilityTable( TableFunction const & relPermTable ) const
{
  ArrayOfArraysView< real64 const > coords = relPermTable.getCoordinates();
  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const relPerm = relPermTable.getValues();
  real64 minVolFraction = phaseVolFrac[0];

  GEOSX_THROW_IF_NE_MSG( relPermTable.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                         getFullName() << ": the interpolation method for the tables must be linear",
                         InputError );

  GEOSX_THROW_IF_NE_MSG( relPermTable.numDimensions(), 1,
                         getFullName() << ": the table must contain one vector of phase volume fraction, and one vector of relative permeabilities",
                         InputError );
  GEOSX_THROW_IF( coords.sizeOfArray( 0 ) < 2,
                  getFullName() << ": the table must contain at least two values",
                  InputError );

  // note that the TableFunction class has already checked that coords.sizeOfArray( 0 ) == relPerm.size()
  for( localIndex i = 0; i < coords.sizeOfArray( 0 ); ++i )
  {

    // check phase volume fraction
    GEOSX_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                    getFullName() << ": the phase volume fraction (i.e. saturation) must be between 0 and 1",
                    InputError );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check phase relative permeability
    if( i == 0 )
    {
      GEOSX_THROW_IF( !isZero( relPerm[i] ),
                      getFullName() << ": the first relative permeability value must be equal to zero",
                      InputError );
    }
    else
    {
      GEOSX_THROW_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-10,
                      getFullName() << ": the relative permeability must be strictly increasing",
                      InputError );

      if( isZero( relPerm[i-1] ) && !isZero( relPerm[i] ) )
      {
        minVolFraction = phaseVolFrac[i-1];
      }
    }

  }
  return minVolFraction;
}

TableRelativePermeability::KernelWrapper::
  KernelWrapper( arrayView1d< geosx::TableFunction::KernelWrapper const > const & waterOilRelPermTableKernelWrappers,
                 arrayView1d< geosx::TableFunction::KernelWrapper const > const & gasOilRelPermTableKernelWrappers,
                 arrayView1d< geosx::real64 const > const & phaseMinVolumeFraction,
                 arrayView1d< geosx::integer const > const & phaseTypes,
                 arrayView1d< geosx::integer const > const & phaseOrder,
                 arrayView3d< geosx::real64, relperm::USD_RELPERM > const & phaseRelPerm,
                 arrayView4d< geosx::real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
  : RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac ),
  m_waterOilRelPermTableKernelWrappers( waterOilRelPermTableKernelWrappers ),
  m_gasOilRelPermTableKernelWrappers( gasOilRelPermTableKernelWrappers ),
  m_phaseMinVolumeFraction( phaseMinVolumeFraction )
{}

TableRelativePermeability::KernelWrapper TableRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_waterOilRelPermTableKernelWrappers,
                        m_gasOilRelPermTableKernelWrappers,
                        m_phaseMinVolumeFraction,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeability, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
