/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
  * @file CompressibleSinglePhaseFluid.cpp
  */

#include "CompressibleSinglePhaseFluid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

static ExponentApproximationType stringToExponentType( string const & model )
{
  if (model == "linear")
  {
    return ExponentApproximationType::Linear;
  }
  else if (model == "quadratic")
  {
    return ExponentApproximationType::Quadratic;
  }
  else if (model == "exponential")
  {
    return ExponentApproximationType::Full;
  }
  GEOSX_ERROR("Model type not supported: " << model);

  // otherwise compilers complain about reaching the end of non-void function
  return ExponentApproximationType::Full;
}

CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( std::string const & name, Group * const parent ):
  SingleFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid compressibility");

  registerWrapper( viewKeyStruct::viscosibilityString, &m_viscosibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid viscosity exponential coefficient");

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference pressure");

  registerWrapper( viewKeyStruct::referenceDensityString, &m_referenceDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid density");

  registerWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid viscosity");

  registerWrapper( viewKeyStruct::densityModelString, &m_densityModelString, false )->
    setApplyDefaultValue("linear")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Type of density model (linear, quadratic, exponential)");

  registerWrapper( viewKeyStruct::viscosityModelString, &m_viscosityModelString, false )->
    setApplyDefaultValue("linear")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Type of viscosity model (linear, quadratic, exponential)");
}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

void CompressibleSinglePhaseFluid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::AllocateConstitutiveData(parent, numConstitutivePointsPerParentIndex);

  m_density = m_referenceDensity;
  m_viscosity = m_referenceViscosity;
}

void
CompressibleSinglePhaseFluid::DeliverClone( string const & name,
                                            Group * const parent,
                                            std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<CompressibleSinglePhaseFluid>( name, parent );
  }
  SingleFluidBase::DeliverClone( name, parent, clone );
  CompressibleSinglePhaseFluid * const newConstitutiveRelation = dynamic_cast<CompressibleSinglePhaseFluid *>(clone.get());


  newConstitutiveRelation->m_compressibility      = this->m_compressibility;
  newConstitutiveRelation->m_viscosibility        = this->m_viscosibility;
  newConstitutiveRelation->m_referencePressure    = this->m_referencePressure;
  newConstitutiveRelation->m_referenceDensity     = this->m_referenceDensity;
  newConstitutiveRelation->m_referenceViscosity   = this->m_referenceViscosity;
  newConstitutiveRelation->m_densityModelString   = this->m_densityModelString;
  newConstitutiveRelation->m_viscosityModelString = this->m_viscosityModelString;
  newConstitutiveRelation->m_densityModelType     = this->m_densityModelType;
  newConstitutiveRelation->m_viscosityModelType   = this->m_viscosityModelType;

}

void CompressibleSinglePhaseFluid::PostProcessInput()
{
  SingleFluidBase::PostProcessInput();

  GEOSX_ERROR_IF( m_compressibility < 0.0, "An invalid value of fluid compressibility ("
                                          << m_compressibility << ") is specified" );

  GEOSX_ERROR_IF( m_viscosibility < 0.0, "An invalid value of fluid viscosibility ("
                                        << m_viscosibility << ") is specified" );

  GEOSX_ERROR_IF( m_referenceDensity <= 0.0, "An invalid value of fluid reference density ("
                                            << m_referenceDensity << ") is specified" );

  GEOSX_ERROR_IF( m_referenceViscosity <= 0.0, "An invalid value of fluid reference viscosity ("
                                              << m_compressibility << ") is specified" );

  m_densityModelType   = stringToExponentType( m_densityModelString );
  m_viscosityModelType = stringToExponentType( m_viscosityModelString );

  real64 dRho_dP;
  real64 dVisc_dP;
  Compute( m_referencePressure, m_referenceDensity, dRho_dP, m_referenceViscosity, dVisc_dP );
  this->getWrapper< array2d<real64> >(viewKeyStruct::dDens_dPresString)->setDefaultValue( dRho_dP );
  this->getWrapper< array2d<real64> >(viewKeyStruct::dVisc_dPresString)->setDefaultValue( dVisc_dP );
}

void CompressibleSinglePhaseFluid::PointUpdateViscosity( real64 const & pressure, localIndex const k, localIndex const q )
{
  if (pressure > m_referencePressure)
    makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto relation )
    {
      Compute( pressure, m_viscosity[k][q], m_dViscosity_dPressure[k][q], relation );
    } );
}

void CompressibleSinglePhaseFluid::PointUpdateDensity( real64 const & pressure, localIndex const k, localIndex const q )
{
  if (pressure > m_referencePressure)
    makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto relation )
    {
      Compute( pressure, m_density[k][q], m_dDensity_dPressure[k][q], relation );
    } );
}

void CompressibleSinglePhaseFluid::PointUpdatePressureExplicit( real64 & pressure, localIndex const k, localIndex const q )
{
  real64 pressureCap = 1e8;
  pressure = m_density[k][q] < m_referenceDensity ? 0 : (pressure <= 0.5 * pressureCap
                                                        ? (1 - m_referenceDensity / m_density[k][q]) / m_compressibility + m_referencePressure
                                                        : 0.5 * pressureCap + ( pressure - 0.5 * pressureCap) / (1.0 / m_compressibility - 0.5 * pressureCap) * 0.5 * pressureCap);

  // In explicit solver, density is calculated from mass, not pressure; below may only be used for pressure-initialization
  //  real64 pressureCap = 1e8;
  //  m_density[k][q] = pressure <= 0.5 * pressureCap ? m_referenceDensity / (1 - pressure * m_compressibility)
  //                                                  : m_referenceDensity / (1 - 0.5 * pressureCap * m_compressibility - (2.0 * pressure / pressureCap - 1.0) * (1.0 - 0.5 * pressureCap * m_compressibility));

  //  if (m_density[k][q] > m_referenceDensity)
  //    Compute( pressure, m_density[k][q], m_viscosity[k][q], m_dViscosity_dPressure[k][q] );
  //  else
  //    pressure = 0;
  //  real64 pressureCap = 1e8;
  //  pressure = min(pressure, pressureCap);
}

void CompressibleSinglePhaseFluid::PointUpdatePressure( real64 & pressure, real64 const & mass, real64 const & volume, real64 const & poroRef, real64 const & totalCompressibility)
{
  real64 value = mass / ( m_referenceDensity * volume * poroRef );
  if (value > 1)
    switch( m_densityModelType )
    {
      case ExponentApproximationType::Full:
      {
        pressure = std::log( value ) / totalCompressibility + m_referencePressure;
        break;
      }
      case ExponentApproximationType::Linear:
      {
        pressure = ( value - 1) / totalCompressibility + m_referencePressure;
        break;
      }
      default:
      {
        GEOSX_ERROR( "PointUpdatePressure() ExponentApproximationType is invalid: Only Full or Linear is available!" );
      }
    }
  else
    pressure = 0;
}

void CompressibleSinglePhaseFluid::PointUpdate( real64 const & pressure, localIndex const k, localIndex const q )
{
  Compute( pressure, m_density[k][q], m_dDensity_dPressure[k][q], m_viscosity[k][q], m_dViscosity_dPressure[k][q] );
}

void CompressibleSinglePhaseFluid::BatchUpdate( arrayView1d<double const> const & pressure )
{
  makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto relation )
  {
    SingleFluidBase::BatchDensityUpdateKernel<CompressibleSinglePhaseFluid>( pressure, relation );
  } );
  makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto relation )
  {
    SingleFluidBase::BatchDensityUpdateKernel<CompressibleSinglePhaseFluid>( pressure, relation );
  } );
}

void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                            real64 & density, real64 & dDensity_dPressure,
                                            real64 & viscosity, real64 & dViscosity_dPressure ) const
{
  makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto relation )
  {
    Compute( pressure, density, dDensity_dPressure, relation );
  } );
  makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto relation )
  {
    Compute( pressure, viscosity, dViscosity_dPressure, relation );
  } );
}

//void CompressibleSinglePhaseFluid::Compute( real64 & pressure,
//                                            real64 const & density,
//                                            real64 & viscosity, real64 & dViscosity_dPressure ) const
//{
//  makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto relation )
//  {
//    Inverse( density, pressure, relation );
//  } );
//  makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto relation )
//  {
//    Compute( pressure, viscosity, dViscosity_dPressure, relation );
//  } );
//}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
