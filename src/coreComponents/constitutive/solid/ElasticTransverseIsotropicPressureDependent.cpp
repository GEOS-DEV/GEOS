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
 *  @file ElasticTransverseIsotropicPressureDependent.cpp
 */

#include "ElasticTransverseIsotropicPressureDependent.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticTransverseIsotropicPressureDependent::ElasticTransverseIsotropicPressureDependent( string const & name, Group * const parent ):
  ElasticTransverseIsotropic( name, parent ),
  m_dc11dp(),
  m_dc13dp(),
  m_dc33dp(),
  m_dc44dp(),
  m_dc66dp()
{
  registerWrapper( viewKeyStruct::dc11dpString(), &m_dc11dp ).
    setApplyDefaultValue( 0 ).
    setDescription( "Elastic Stiffness Field C11 pressure derivative" );

  registerWrapper( viewKeyStruct::dc13dpString(), &m_dc13dp ).
    setApplyDefaultValue( 0 ).
    setDescription( "Elastic Stiffness Field C13 pressure derivative" );

  registerWrapper( viewKeyStruct::dc33dpString(), &m_dc33dp ).
    setApplyDefaultValue( 0 ).
    setDescription( "Elastic Stiffness Field C33 pressure derivative" );

  registerWrapper( viewKeyStruct::dc44dpString(), &m_dc44dp ).
    setApplyDefaultValue( 0 ).
    setDescription( "Elastic Stiffness Field C44 pressure derivative" );

  registerWrapper( viewKeyStruct::dc66dpString(), &m_dc66dp ).
    setApplyDefaultValue( 0 ).
    setDescription( "Elastic Stiffness Field C66 pressure derivative" );
}

ElasticTransverseIsotropicPressureDependent::~ElasticTransverseIsotropicPressureDependent()
{}

void ElasticTransverseIsotropicPressureDependent::postProcessInput()
{
  ElasticTransverseIsotropic::postProcessInput();
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticTransverseIsotropicPressureDependent, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
