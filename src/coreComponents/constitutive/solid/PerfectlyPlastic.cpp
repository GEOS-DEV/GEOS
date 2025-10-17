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
 *  @file PerfectlyPlastic.cpp
 */

#include "PerfectlyPlastic.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

PerfectlyPlastic::PerfectlyPlastic( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent )
{
  // register default values
  registerWrapper( viewKeyStruct::defaultYieldStressString(), &m_defaultYieldStress ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default yield stress" );

  // register fields
  registerField< fields::solid::yieldStress >( &m_yieldStress );
}


void PerfectlyPlastic::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultYieldStress < 0.0, "Negative yield stress detected", InputError );

  getField< fields::solid::yieldStress >().setApplyDefaultValue( m_defaultYieldStress );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PerfectlyPlastic, std::string const &, Group * const )
}
} /* namespace geos */
