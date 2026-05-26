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
 * @file CrystalHealMPMEvent.cpp
 */

#include "CrystalHealMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

CrystalHealMPMEvent::CrystalHealMPMEvent( const string & name,
                                          Group * const parent ):
  MPMEventBase( name, parent ),
  m_targetRegion( "mat1" ),
  m_healType( 0 ),
  m_markedParticlesToHeal( 0 ),
  m_strengthKnockdown( 0.0 )
{
  registerWrapper( viewKeyStruct::targetRegionString(), &m_targetRegion ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Particle region to perform heal on" );

  registerWrapper( viewKeyStruct::healTypeString(), &m_healType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of heal algorithm to perform" );

  registerWrapper( viewKeyStruct::markedParticlesToHealString(), &m_markedParticlesToHeal ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_markedParticlesToHeal ).
    setDescription( "Flag whether identification of particles to heal has been performed" );

  registerWrapper( viewKeyStruct::strengthKnockdownString(), &m_strengthKnockdown ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_strengthKnockdown ).
    setDescription( "Knockdown factor to apply to fully damaged material after healing" );
}

CrystalHealMPMEvent::~CrystalHealMPMEvent()
{}

void CrystalHealMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( m_strengthKnockdown < 0.0 || m_strengthKnockdown > 1.0,  "Crystal heal MPM Event " << getName() << " strength knockdown should be between 0.0 and 1.0");
}

REGISTER_CATALOG_ENTRY( MPMEventBase, CrystalHealMPMEvent, string const &, Group * const )

} /* namespace geos */