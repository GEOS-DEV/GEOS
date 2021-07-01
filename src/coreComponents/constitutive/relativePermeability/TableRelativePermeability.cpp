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
}

TableRelativePermeability::~TableRelativePermeability()
{}

void TableRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::OIL] < 0,
                  "TableRelativePermeability: reference oil phase has not been defined and must be included in model",
                  InputError );

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::WATER] >= 0 && !(m_waterOilRelPermTableNames.size() == 2),
                  "TableRelativePermeability: since water is present you must define two tables for the oil-water relperms: "
                  "one for the oil phase and one for the water phase",
                  InputError );

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::GAS] >= 0 && !(m_gasOilRelPermTableNames.size() == 2),
                  "TableRelativePermeability: since gas is present you must define two tables for the oil-gas relperms: "
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

  m_phaseMinVolumeFraction.resize( PhaseType::MAX_NUM_PHASES );

  if( m_waterOilRelPermTableKernelWrappers.empty() && m_gasOilRelPermTableKernelWrappers.empty() )
  {

    // check water-oil relperms
    for( localIndex ip = 0; ip < m_waterOilRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction const >( m_waterOilRelPermTableNames[ip] );
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
        GEOSX_THROW( "There should be only two table names for the water-oil pair", InputError );
      }
      m_waterOilRelPermTableKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }

    // check gas-oil relperms
    for( localIndex ip = 0; ip < m_gasOilRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction const >( m_gasOilRelPermTableNames[ip] );
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
        GEOSX_THROW( "There should be only two table names for the gas-oil pair", InputError );
      }

      m_gasOilRelPermTableKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }
  }
}

real64 TableRelativePermeability::validateRelativePermeabilityTable( TableFunction const & relPermTable ) const
{
  ArrayOfArraysView< real64 const > coords = relPermTable.getCoordinates();
  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  array1d< real64 > const & relPerm = relPermTable.getValues();
  real64 minVolFraction = phaseVolFrac[0];

  GEOSX_THROW_IF( relPermTable.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                  "The interpolation method for the relative permeability tables must be linear",
                  InputError );

  GEOSX_THROW_IF( coords.size() != 1,
                  "The relative permeability table must contain one vector of phase volume fraction, and one vector of relative permeabilities",
                  InputError );
  GEOSX_THROW_IF( coords.sizeOfArray( 0 ) < 2,
                  "The relative permeability table must contain at least two values",
                  InputError );

  // note that the TableFunction class has already checked that coords.sizeOfArray( 0 ) == relPerm.size()
  for( localIndex i = 0; i < coords.sizeOfArray( 0 ); ++i )
  {

    // check phase volume fraction
    GEOSX_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                    "In the relative permeability table, the phase volume fraction (i.e., saturation) must be between 0 and 1",
                    InputError );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check phase relative permeability
    if( i == 0 )
    {
      GEOSX_THROW_IF( !isZero( relPerm[i] ),
                      "In the relative permeability table, the first relative permeability value must be equal to zero",
                      InputError );
    }
    else
    {
      GEOSX_THROW_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-10,
                      "In the relative permeability table, the relative permeability must be strictly increasing",
                      InputError );

      if( isZero( relPerm[i-1] ) && !isZero( relPerm[i] ) )
      {
        minVolFraction = phaseVolFrac[i-1];
      }
    }

  }
  return minVolFraction;
}

std::unique_ptr< ConstitutiveBase >
TableRelativePermeability::deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase >
  clone = RelativePermeabilityBase::deliverClone( name, parent );
  TableRelativePermeability & tableRelPerm = dynamicCast< TableRelativePermeability & >( *clone );

  // TODO: see if it is simpler to just add a copy constructor for the wrappers, and let GEOSX copy everything automatically.
  tableRelPerm.createAllTableKernelWrappers();
  return clone;
}

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
