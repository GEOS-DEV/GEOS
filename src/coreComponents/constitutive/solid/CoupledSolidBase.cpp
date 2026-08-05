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
 * @file CoupledSolidBase.cpp
 */

#include "CoupledSolidBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CoupledSolidBase::CoupledSolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the porosity model." );

  registerWrapper( viewKeyStruct::permeabilityModelNameString(), &m_permeabilityModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the permeability model." );

  registerWrapper( viewKeyStruct::solidInternalEnergyModelNameString(), &m_solidInternalEnergyModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the solid internal energy model." );

  registerWrapper( viewKeyStruct::surfaceAreaModelNameString(), &m_surfaceAreaModelName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the reactive surface area model." );
}

} /* namespace constitutive */

} /* namespace geos */
