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
 * @file CoupledSolidBase.cpp
 */

#include "CoupledSolidBase.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "ElasticOrthotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "Damage.hpp"
#include "DamageSpectral.hpp"
#include "DamageVolDev.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

CoupledSolidBase::CoupledSolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_solidModelName(),
  m_porosityModelName(),
  m_permeabilityModelName(),
  m_solidInternalEnergyModelName()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the porosity model." );

  registerWrapper( viewKeyStruct::permeabilityModelNameString(), &m_permeabilityModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the permeability model." );

  registerWrapper( viewKeyStruct::solidInternalEnergyModelNameString(), &m_solidInternalEnergyModelName ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the solid internal energy model." );
}

CoupledSolidBase::~CoupledSolidBase() = default;

}
} /* namespace geosx */
