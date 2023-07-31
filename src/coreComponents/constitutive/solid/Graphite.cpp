/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file Graphite.cpp
 */

#include "Graphite.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

    /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// State variable: The material direction for each element/particle
  array2d< real64 > m_materialDirection;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Material parameter: The value of the failure strength
  real64 m_failureStrength;

  /// Material parameter: The value of crack speed
  real64 m_crackSpeed;

Graphite::Graphite( string const & name, Group * const parent ):
  ElasticTransverseIsotropicPressureDependent( name, parent ),
  m_plasticStrain(),
  m_relaxation(),
  m_damage(),
  m_jacobian(),
  m_materialDirection(),
  m_lengthScale(),
  m_failureStrength(),
  m_crackSpeed(),
  m_damagedMaterialFrictionalSlope(),
  m_distortionShearResponseX2(),
  m_distortionShearResponseY1(),
  m_distortionShearResponseY2(),
  m_distortionShearResponseM1(),
  m_inPlaneShearResponseX2(),
  m_inPlaneShearResponseY1(),
  m_inPlaneShearResponseY2(),
  m_inPlaneShearResponseM1(),
  m_coupledShearResponseX2(),
  m_coupledShearResponseY1(),
  m_coupledShearResponseY2(),
  m_coupledShearResponseM1(),
  m_maximumPlasticStrain()
{
  // register default values
  registerWrapper( viewKeyStruct::failureStrengthString(), &m_failureStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crack speed" );

  registerWrapper( viewKeyStruct::damagedMaterialFrictionalSlopeString(), &m_damagedMaterialFrictionalSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Damaged material frictional slope" );

  // register fields
  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::LEVEL_0).
    setDescription( "Strength scale" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point damage values" );

  //Need to register inputs for failure envelope
}


Graphite::~Graphite()
{}


void Graphite::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticTransverseIsotropicPressureDependent::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
}


void Graphite::postProcessInput()
{
  ElasticTransverseIsotropicPressureDependent::postProcessInput();

  GEOS_THROW_IF( m_failureStrength <= 0.0, "Maximum theoretical strength must be greater than 0", InputError );
  GEOS_THROW_IF( m_crackSpeed <= 0.0, "Crack speed must be a positive number.", InputError );
}


void Graphite::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Graphite, std::string const &, Group * const )
}
} /* namespace geos */
