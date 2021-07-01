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

#include "MultiFluidPVTPackageWrapper.hpp"

#include <map>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

MultiFluidPVTPackageWrapper::MultiFluidPVTPackageWrapper( string const & name,
                                                          Group * const parent )
  :
  MultiFluidBase( name, parent ),
  m_fluid( nullptr )
{ }

MultiFluidPVTPackageWrapper::~MultiFluidPVTPackageWrapper() = default;

namespace
{

pvt::PHASE_TYPE getPVTPackagePhaseType( string const & name )
{
  static std::map< string, pvt::PHASE_TYPE > const phaseTypes{
    { "gas", pvt::PHASE_TYPE::GAS },
    { "oil", pvt::PHASE_TYPE::OIL },
    { "water", pvt::PHASE_TYPE::LIQUID_WATER_RICH }
  };
  auto const it = phaseTypes.find( name );
  GEOSX_ERROR_IF( it == phaseTypes.end(), "Fluid phase not supported by PVTPackage: " << name );
  return it->second;
}

} // namespace

void MultiFluidPVTPackageWrapper::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  m_phaseTypes.resize( numFluidPhases() );
  std::transform( m_phaseNames.begin(), m_phaseNames.end(), m_phaseTypes.begin(), getPVTPackagePhaseType );
}

void MultiFluidPVTPackageWrapper::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

std::unique_ptr< ConstitutiveBase >
MultiFluidPVTPackageWrapper::deliverClone( string const & name,
                                           Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  MultiFluidPVTPackageWrapper & model = dynamicCast< MultiFluidPVTPackageWrapper & >( *clone );
  model.m_phaseTypes = m_phaseTypes;

  return clone;
}

} //namespace constitutive

} //namespace geosx
