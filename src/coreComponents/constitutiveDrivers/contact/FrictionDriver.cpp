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

#include "FrictionDriver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutiveDrivers/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionBase.hpp"
#include "constitutive/contact/FrictionSelector.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

FrictionDriver::FrictionDriver( const string & name, Group * const parent )
  : ConstitutiveDriver( name, parent )
{
  registerWrapper( viewKeyStruct::frictionNameString(), &m_frictionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Friction model to test" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of sample step to take in both jumps and traction increments" );

  registerWrapper( viewKeyStruct::jumpFunctionString(), &m_jumpFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing jump function along world x-axis" );

  registerWrapper( viewKeyStruct::dJumpFunctionString(), &m_dJumpFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing deltaDisplacementJump function along world x-axis" );

  registerWrapper( viewKeyStruct::tractionFunctionString(), &m_tractionFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing traction function along world x-axis" );

  registerWrapper( viewKeyStruct::thetaString(), &m_theta ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "x-Tilt angle in degree" );

  registerWrapper( viewKeyStruct::phiString(), &m_phi ).
    setInputFlag( InputFlags::INVALID ).
    setDescription( "y-Tilt angle in degree" );

  //first batch of parameters
  registerWrapper( viewKeyStruct::normalDispTolFac(), &m_normalDispTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Displacement Tolerance (scale as inverse of average Young modulus)." );

  registerWrapper( viewKeyStruct::normalTractionTolFac(), &m_normalTracTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Traction Tolerance" );

  registerWrapper( viewKeyStruct::slidingTolFac(), &m_slidingTolFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "tangential Displacement Tolerance" );

  registerWrapper( viewKeyStruct::iterPenNFac(), &m_iterPenNFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "normal Penalty Factor" );

  registerWrapper( viewKeyStruct::iterPenTFac(), &m_iterPenTFac ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "tangential Penatly Factor" );


  //geometry
  registerWrapper( viewKeyStruct::area(), &m_area ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Face area" );

  registerWrapper( viewKeyStruct::volume(), &m_volume ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' volume" );

  registerWrapper( viewKeyStruct::bulk(), &m_bulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' bulk Modulus" );

  registerWrapper( viewKeyStruct::shear(), &m_shearModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Neighboring cells' shear Modulus" );

  //algo tune
  registerWrapper( viewKeyStruct::simultaneous(), &m_isSimultaneous ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "isSimultaneous" );



  addLogLevel< logInfo::LogOutput >();
}

void FrictionDriver::postInputInitialization()
{
  ConstitutiveDriver::postInputInitialization();

  // Check that the functions exist
  FunctionManager & functionManager = FunctionManager::getInstance();
  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_jumpFunctionName ),
                 GEOS_FMT( "Jump function with name '{}' not found", m_jumpFunctionName ),
                 getWrapperDataContext( viewKeyStruct::jumpFunctionString() ) );

  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_dJumpFunctionName ),
                 GEOS_FMT( "dJump function with name '{}' not found", m_dJumpFunctionName ),
                 getWrapperDataContext( viewKeyStruct::dJumpFunctionString() ) );

  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_tractionFunctionName ),
                 GEOS_FMT( "Traction function with name '{}' not found", m_tractionFunctionName ),
                 getWrapperDataContext( viewKeyStruct::tractionFunctionString() ) );

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  // initialize functions
  TableFunction & jumpFunction = functionManager.getGroup< TableFunction >( m_jumpFunctionName );
  TableFunction & dJumpFunction = functionManager.getGroup< TableFunction >( m_dJumpFunctionName );
  TableFunction & tractionFunction = functionManager.getGroup< TableFunction >( m_tractionFunctionName );

  jumpFunction.initializeFunction();
  dJumpFunction.initializeFunction();
  tractionFunction.initializeFunction();

  // TODO: Maybe we should take the maximum extent of jumpFunction and tractionFunction
  ArrayOfArraysView< real64 > coordinates = jumpFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];

  // Allocate the data
  allocateTable( numCols, minTime, maxTime );

  // set input columns
  initializeTable();
}

bool FrictionDriver::execute()
{
  FrictionBase & baseFriction = getFriction();

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Friction Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Friction ............... " << m_frictionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type ................... " << baseFriction.getCatalogName() );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps .................. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << m_outputFile );

  // create a dummy discretization with one quadrature point for
  // storing constitutive data
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  integer const numRows = m_table.size( 0 );
  discretization.resize( numRows ); // numRows elements
  baseFriction.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  constitutiveUpdatePassThru( baseFriction, [&]( auto & selectedFrictionModel )
  {
    using FRICTION_TYPE = TYPEOFREF( selectedFrictionModel );
    runTest< FRICTION_TYPE >( selectedFrictionModel, m_table );
  } );

  return false;
}

void FrictionDriver::getColumnNames( string_array & columnNames ) const
{
  columnNames.emplace_back( "time" );
  columnNames.emplace_back( "traction,normal" );
  columnNames.emplace_back( "traction,tangent1" );
  columnNames.emplace_back( "traction,tangent2" );
  columnNames.emplace_back( "displacement jump,normal" );
  columnNames.emplace_back( "displacement jump,tangent1" );
  columnNames.emplace_back( "displacement jump,tangent2" );
  columnNames.emplace_back( "delta displacement jump,normal" );
  columnNames.emplace_back( "delta displacement jump,tangent1" );
  columnNames.emplace_back( "delta displacement jump,tangent2" );
  columnNames.emplace_back( "encoded constaint (0:converged, 1:stick & gn>0 (opening), 2: interpenetration, 3: stick & gt>lim (disp-sliding), 4: tau>taulim (trac-sliding) )" );
  columnNames.emplace_back( "fracture state (0:stick, 1:slip , 2: new slip, 3: open)" );
  columnNames.emplace_back( "newtraction,normal" );
  columnNames.emplace_back( "newtraction,tangent1" );
  columnNames.emplace_back( "newtraction,tangent2" );
  columnNames.emplace_back( "iterative penalty, normal" );
  columnNames.emplace_back( "iterative penalty, tangent" );

  if( dynamic_cast< CoulombFriction const * >(&getFriction()) != nullptr )
  {
    columnNames.emplace_back( "tau limit" );
  }
}

void FrictionDriver::initializeTable()
{
  integer const numRows = m_table.size( 0 );

  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction const & jumpFunction = functionManager.getGroup< TableFunction >( m_jumpFunctionName );
  TableFunction const & dJumpFunction = functionManager.getGroup< TableFunction >( m_dJumpFunctionName );
  TableFunction const & tractionFunction = functionManager.getGroup< TableFunction >( m_tractionFunctionName );

  real64 const cos_theta = cos( m_theta * M_PI/180.0 );
  real64 const sin_theta = sin( m_theta * M_PI/180.0 );

  real64 const cos_phi = cos( m_phi * M_PI/180.0 );
  real64 const sin_phi = sin( m_phi * M_PI/180.0 );

  for( integer index = 0; index < numRows; ++index )
  {
    real64 const time = m_table( index, TIME );

    real64 const traction = tractionFunction.evaluate( &time );
    m_table( index, NTRAC ) = traction*sin_theta;
    m_table( index, STRAC0 ) = traction*cos_theta*cos_phi;
    m_table( index, STRAC1 ) = traction*cos_theta*sin_phi;

    real64 const jump = jumpFunction.evaluate( &time );
    m_table( index, NJUMP ) = jump*sin_theta;
    m_table( index, SLIP0 ) = jump*cos_theta*cos_phi;
    m_table( index, SLIP1 ) = jump*cos_theta*sin_phi;

    real64 const dJump = dJumpFunction.evaluate( &time );
    m_table( index, NDJUMP ) = dJump*sin_theta;
    m_table( index, DSLIP0 ) = dJump*cos_theta*cos_phi;
    m_table( index, DSLIP1 ) = dJump*cos_theta*sin_phi;

    m_table( index, CC )    = 0;
    m_table( index, FS )    = fields::contact::FractureState::Stick;
  }

  if( CoulombFriction const * coulombFriction = dynamic_cast< CoulombFriction const * >(&getFriction()) )
  {
    real64 const cohesion = coulombFriction->getCohesion();
    real64 const frictionCoeff = coulombFriction->getFrictionCoeff();
    for( integer index = 0; index < numRows; ++index )
    {
      real64 const normal_traction = m_table( index, NTRAC );
      m_table( index, TLIM ) = cohesion - normal_traction * frictionCoeff;
    }
  }
}

FrictionBase & FrictionDriver::getFriction()
{
  return getConstitutiveManager().getGroup< FrictionBase >( m_frictionName );
}

FrictionBase const & FrictionDriver::getFriction() const
{
  return getConstitutiveManager().getGroup< FrictionBase >( m_frictionName );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        FrictionDriver,
                        string const &, dataRepository::Group * const )

}
