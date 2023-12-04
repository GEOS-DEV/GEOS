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
 * @file ContactBase.cpp
 */

#include "ContactBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ContactBase::ContactBase( string const & name,
                          Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_apertureTable( nullptr )
{
  registerWrapper( viewKeyStruct::penaltyStiffnessString(), &m_penaltyStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the penetration penalty stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::shearStiffnessString(), &m_shearStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the shear elastic stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::apertureToleranceString(), &m_apertureTolerance ).
    setApplyDefaultValue( 1.0e-9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to be used to avoid floating point errors in expressions involving aperture. "
                    "For example in the case of dividing by the actual aperture (not the effective aperture "
                    "that results from the aperture function) this value may be used to avoid the 1/0 error. "
                    "Note that this value may have some physical significance in its usage, as it may be used "
                    "to smooth out highly nonlinear behavior associated with 1/0 in addition to avoiding the "
                    "1/0 error." );

  registerWrapper( viewKeyStruct::displacementJumpThresholdString(), &m_displacementJumpThreshold ).
    setApplyDefaultValue( std::numeric_limits< real64 >::epsilon() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "A threshold valued to determine whether a fracture is open or not." );


  registerWrapper( viewKeyStruct::apertureTableNameString(), &m_apertureTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the aperture table" );

  registerWrapper( viewKeyStruct::useApertureModelString(), &m_useApertureModel ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag for using the nonlinear aperture model instead of the aperture model."
                    "For now, 1 is for the Barton--Bandis model, and 2 is for the exponential model." );

  registerWrapper( viewKeyStruct::refNormalStressString(), &m_refNormalStress ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the reference normal stress in the aperture model." );
}

ContactBase::~ContactBase()
{}


void ContactBase::postProcessInput()
{

  GEOS_THROW_IF( m_apertureTableName.empty() && m_useApertureModel == 0,
                 getFullName() << ": the aperture table name " << m_apertureTableName << " is empty", InputError );

  GEOS_THROW_IF( m_useApertureModel >= 1 && m_refNormalStress <= 0.0,
                 getFullName() << ": the reference normal stress is not given or has a wrong value", InputError );

}

void ContactBase::initializePreSubGroups()
{}


void ContactBase::allocateConstitutiveData( Group & parent,
                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  if( m_useApertureModel == 0 )
  {
    FunctionManager & functionManager = FunctionManager::getInstance();

    GEOS_THROW_IF( !functionManager.hasGroup( m_apertureTableName ),
                  getFullName() << ": the aperture table named " << m_apertureTableName << " could not be found",
                  InputError );

    TableFunction & apertureTable = functionManager.getGroup< TableFunction >( m_apertureTableName );
    validateApertureTable( apertureTable );

    ArrayOfArraysView< real64 > coords = apertureTable.getCoordinates();
    arraySlice1d< real64 const > apertureValues = coords[0];
    array1d< real64 > & hydraulicApertureValues = apertureTable.getValues();

    localIndex const n = apertureValues.size()-1;
    real64 const slope = ( hydraulicApertureValues[n] - hydraulicApertureValues[n-1] ) / ( apertureValues[n] - apertureValues[n-1] );
    real64 const apertureTransition = ( hydraulicApertureValues[n] - slope * apertureValues[n] ) / ( 1.0 - slope );

    // if the aperture transition is larger than the last coordinates, we enlarge the table
    // this check is necessary to ensure that the coordinates are strictly increasing
    if( apertureTransition > apertureValues[apertureValues.size()-1] )
    {
      GEOS_LOG_RANK_0( GEOS_FMT ( "Adding aperture transition for table {}:", m_apertureTableName ) );
      std::ostringstream s_orig;
      for( localIndex i = 0; i < apertureValues.size(); i++ )
        s_orig << "[ " << apertureValues[i] << ", " << hydraulicApertureValues[i] << " ] ";
      GEOS_LOG_RANK_0( GEOS_FMT ( "    Original table = {}", s_orig.str()));

      coords.emplaceBack( 0, apertureTransition );
      hydraulicApertureValues.emplace_back( apertureTransition );
      // if the aperture transition is larger than 0, we keep enlarging the table
      // this check is necessary to ensure that the coordinates are strictly increasing
      if( apertureTransition > 0 )
      {
        coords.emplaceBack( 0, apertureTransition*10e9 );
        hydraulicApertureValues.emplace_back( apertureTransition*10e9 );
        apertureTable.reInitializeFunction();
      }

      std::ostringstream s_mod;
      for( localIndex i = 0; i < apertureValues.size(); i++ )
        s_mod << "[ " << apertureValues[i] << ", " << hydraulicApertureValues[i] << " ] ";
      GEOS_LOG_RANK_0( GEOS_FMT ( "    Modified table = {}", s_mod.str()));
    }

    m_apertureTable = &apertureTable;
  }
}


void ContactBase::validateApertureTable( TableFunction const & apertureTable ) const
{
  ArrayOfArraysView< real64 const > const coords = apertureTable.getCoordinates();
  arrayView1d< real64 const > const & hydraulicApertureValues = apertureTable.getValues();

  GEOS_THROW_IF( coords.size() > 1,
                 getFullName() << ": Aperture limiter table cannot be greater than a 1D table.",
                 InputError );

  arraySlice1d< real64 const > apertureValues = coords[0];
  localIndex const size = apertureValues.size();

  GEOS_THROW_IF( coords( 0, size-1 ) > 0.0 || coords( 0, size-1 ) < 0.0,
                 getFullName() << ": Invalid aperture limiter table. Last coordinate must be zero!",
                 InputError );

  GEOS_THROW_IF( apertureValues.size() < 2,
                 getFullName() << ": Invalid aperture limiter table. Must have more than two points specified",
                 InputError );

  localIndex const n = apertureValues.size()-1;
  real64 const slope = ( hydraulicApertureValues[n] - hydraulicApertureValues[n-1] ) / ( apertureValues[n] - apertureValues[n-1] );

  GEOS_THROW_IF( slope >= 1.0,
                 getFullName() << ": Invalid aperture table. The slope of the last two points >= 1 is invalid.",
                 InputError );
}



ContactBaseUpdates ContactBase::createKernelWrapper() const
{
  return ContactBaseUpdates( m_penaltyStiffness,
                             m_shearStiffness,
                             m_displacementJumpThreshold,
                             *m_apertureTable,
                             m_useApertureModel,
                             m_refNormalStress );
}

} /* end namespace constitutive */

} /* end namespace geos */
