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
 * @file TableRelativePermeability.cpp
 */

#include "TableRelativePermeability.hpp"

#include "managers/Functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableRelativePermeability::TableRelativePermeability( std::string const & name,
                                                      Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::relPermTableNamesString, &m_relPermTableNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "List of relative permeability tables" );
}

TableRelativePermeability::~TableRelativePermeability()
{}


void TableRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();
}

void TableRelativePermeability::InitializePreSubGroups( Group * const )
{
  FunctionManager const & functionManager = FunctionManager::Instance();
  for( localIndex ip = 0; ip < m_relPermTableNames.size(); ++ip )
  {
    TableFunction const & relPermTable = *functionManager.GetGroup< TableFunction const >( m_relPermTableNames[ip] );
    ValidateRelativePermeabilityTable( relPermTable );
    m_relPermTableKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
  }
}

void TableRelativePermeability::ValidateRelativePermeabilityTable( TableFunction const & relPermTable ) const
{
  ArrayOfArrays< real64 > const & coords = relPermTable.getCoordinates();
  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  array1d< real64 > const & relPerm = relPermTable.getValues();

  GEOSX_ERROR_IF( coords.size() != 1,
                  "The relative permeability table must contain one vector of phase volume fraction, and one vector of relative permeabilities" );
  GEOSX_ERROR_IF( coords.sizeOfArray( 0 ) < 2,
                  "The relative permeability table must contain at least two values" );

  // note that the TableFunction class has already checked that coords.sizeOfArray( 0 ) == relPerm.size()
  for( localIndex i = 0; i < coords.sizeOfArray( 0 ); ++i )
  {

    // check phase volume fraction
    GEOSX_ERROR_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                    "In the relative permeability table, the phase volume fraction (i.e., saturation) must be between 0 and 1" );

    if( i >= 1 )
    {
      GEOSX_ERROR_IF( phaseVolFrac[i] - phaseVolFrac[i-1] < 1e-10,
                      "In the relative permeability table, the phase volume fraction (i.e., saturation) must be strictly increasing" );
    }

    // check phase relative permeability
    if( i == 0 )
    {
      GEOSX_ERROR_IF( !isZero( relPerm[i] ),
                      "In the relative permeability table, the first relative permeability value must be equal to zero" );
    }
    else
    {
      GEOSX_ERROR_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-10,
                      "In the relative permeability table, the relative permeability must be strictly increasing" );
    }
  }
}

TableRelativePermeability::KernelWrapper TableRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_relPermTableKernelWrappers,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeability, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
