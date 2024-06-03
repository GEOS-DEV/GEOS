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
 * @file ConstantViscosity.cpp
 */

#include "ConstantViscosity.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ConstantViscosity::ConstantViscosity( string const & name,
                                      ComponentProperties const & componentProperties,
                                      integer const phaseIndex,
                                      ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  Parameters const * parameters = dynamic_cast< Parameters const * >(&modelParameters);
  m_constantPhaseViscosity = parameters->m_constantPhaseViscosity[phaseIndex];
}

ConstantViscosityUpdate::ConstantViscosityUpdate( real64 const constantPhaseViscosity ):
  m_constantPhaseViscosity( constantPhaseViscosity )
{}

ConstantViscosity::KernelWrapper
ConstantViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_constantPhaseViscosity );
}

void ConstantViscosity::Parameters::registerParameters( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::constantPhaseViscosityString(), &m_constantPhaseViscosity ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Constant phase viscosity" );
}

void ConstantViscosity::Parameters::postProcessInput( MultiFluidBase const * fluid )
{
  integer const numPhase = fluid->numFluidPhases();

  if( m_constantPhaseViscosity.empty() )
  {
    m_constantPhaseViscosity.resize( numPhase );
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      m_constantPhaseViscosity[ip] = ConstantViscosity::defaultViscosity;
    }
  }
  GEOS_THROW_IF_NE_MSG( m_constantPhaseViscosity.size(), numPhase,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", fluid->getFullName(),
                                  viewKeyStruct::constantPhaseViscosityString() ),
                        InputError );
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
