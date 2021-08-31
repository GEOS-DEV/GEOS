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

#include "BlackOilFluid.hpp"

#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "math/extrapolation/Extrapolation.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

BlackOilFluid::BlackOilFluid( string const & name,
                              Group * const parent )
  :
  BlackOilFluidBase( name, parent )
{}

void BlackOilFluid::postProcessInput()
{
  BlackOilFluidBase::postProcessInput();

  GEOSX_THROW_IF( numFluidPhases() != 3,
                  "BlackOilFluid model named " << getName() << ": this model only supports three-phase flow ",
                  InputError );
}

std::unique_ptr< ConstitutiveBase >
BlackOilFluid::deliverClone( string const & name,
                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  BlackOilFluid & model = dynamicCast< BlackOilFluid & >( *clone );
  model.m_phaseTypes = m_phaseTypes;
  model.m_phaseOrder = m_phaseOrder;
  model.m_hydrocarbonPhaseOrder = m_hydrocarbonPhaseOrder;

  model.m_PVTO.Rs = m_PVTO.Rs;
  model.m_PVTO.bubblePressure = m_PVTO.bubblePressure;
  model.m_PVTO.saturatedBo = m_PVTO.saturatedBo;
  model.m_PVTO.saturatedViscosity = m_PVTO.saturatedViscosity;
  model.m_PVTO.undersaturatedPressure2d = m_PVTO.undersaturatedPressure2d;
  model.m_PVTO.undersaturatedBo2d = m_PVTO.undersaturatedBo2d;
  model.m_PVTO.undersaturatedViscosity2d = m_PVTO.undersaturatedViscosity2d;
  model.m_PVTO.surfaceMassDensity = m_PVTO.surfaceMassDensity;
  model.m_PVTO.surfaceMoleDensity = m_PVTO.surfaceMoleDensity;

  model.createAllKernelWrappers();
  return clone;
}

void BlackOilFluid::readInputDataFromTableFunctions()
{
  GEOSX_ERROR( "BlackOilFluid: this option is not implemented yet, please provide PVT files in standard Eclipse format" );
}

void BlackOilFluid::readInputDataFromPVTFiles()
{
  GEOSX_THROW_IF( m_formationVolFactorTableNames.size() > 0.0 || m_viscosityTableNames.size() > 0.0,
                  "BlackOilFluid model named " << getName() << ": input is redundant (user provided both TableFunction names and pvt files)",
                  InputError );

  using PT = BlackOilFluid::PhaseType;

  PVTProps::BlackOilTables boTables;
  boTables.buildAllTables( PT::OIL,
                           PT::WATER,
                           PT::GAS,
                           m_phaseTypes,
                           m_tableFiles );

  // water data
  array1d< array1d< real64 > > data( 1 );
  data[0] = boTables.getWaterTable();
  fillWaterData( data );

  // gas data
  fillHydrocarbonData( PT::GAS, boTables.getGasTable() );

  // for the Black-Oil model, the oil PVT is treated differently from gas
  fillPVTOData( boTables.getOilTable(),
                m_surfacePhaseMassDensity[m_phaseOrder[PT::OIL]],
                m_componentMolarWeight[m_phaseOrder[PT::OIL]],
                m_surfacePhaseMassDensity[m_phaseOrder[PT::GAS]],
                m_componentMolarWeight[m_phaseOrder[PT::GAS]] );

  // check consistency
  checkTableConsistency();

  // create missing undersaturated branches
  createUndersaturatedProperties();

  // extend existing branches
  extendUndersaturatedProperties();

  // refine table branches
  refineUndersaturatedTables( 100 );

  // check consistency
  checkTableConsistency();
}


void BlackOilFluid::fillPVTOData( array1d< array1d< real64 > > const & oilTable,
                                  real64 const oilSurfaceMassDensity,
                                  real64 const oilSurfaceMolecularWeight,
                                  real64 const gasSurfaceMassDensity,
                                  real64 const gasSurfaceMolecularWeight )
{

  // (standard) format of oilTable:
  // if oilTable.size() == 4 (saturated case):   Rs, bubble point pressure, Bo, viscosity
  // if ollTable.size() == 3 (unsaturated case):     unsaturated pressure,  Bo, viscosity


  // Step 1: count the number of saturated points by looping through oilTable, and resize tables accordingly

  localIndex numSaturatedPoints = 0;
  for( localIndex i = 0; i < oilTable.size(); ++i )
  {
    if( oilTable[i].size() == 4 )
    {
      numSaturatedPoints++;
    }
  }
  m_PVTO.numSaturatedPoints = numSaturatedPoints;

  // store data in a new structure
  m_PVTO.Rs.resize( numSaturatedPoints );
  m_PVTO.bubblePressure.resize( numSaturatedPoints );
  m_PVTO.saturatedBo.resize( numSaturatedPoints );
  m_PVTO.saturatedViscosity.resize( numSaturatedPoints );
  m_PVTO.undersaturatedPressure.resize( numSaturatedPoints );
  m_PVTO.undersaturatedBo.resize( numSaturatedPoints );
  m_PVTO.undersaturatedViscosity.resize( numSaturatedPoints );


  // Step 2: copy values from the tables

  localIndex iSat = 0;
  for( localIndex i = 0; i < oilTable.size(); ++i )
  {

    // Step 2.1: get the saturated values

    m_PVTO.Rs[iSat] = oilTable[i][0];
    m_PVTO.bubblePressure[iSat] = oilTable[i][1];
    m_PVTO.saturatedBo[iSat] = oilTable[i][2];
    m_PVTO.saturatedViscosity[iSat] = oilTable[i][3];

    // Step 2.2: get the undersaturated values

    // Step 2.2.1: start from pressure = P - Pbub = 0

    m_PVTO.undersaturatedPressure[iSat].emplace_back( 0 );
    m_PVTO.undersaturatedBo[iSat].emplace_back( m_PVTO.saturatedBo[iSat] );
    m_PVTO.undersaturatedViscosity[iSat].emplace_back( m_PVTO.saturatedViscosity[iSat] );

    // Step 2.2.2: loop over the undersaturated values

    // Note: these undersaturated values may be absent, in which case an undersaturated branch is created in createUndersaturatedProperties

    localIndex branchSize = 0;
    localIndex j = i + 1;
    while( oilTable[j].size() == 3 )
    {
      m_PVTO.undersaturatedPressure[iSat].emplace_back( oilTable[j][0] - m_PVTO.bubblePressure[iSat] );
      m_PVTO.undersaturatedBo[iSat].emplace_back( oilTable[j][1] );
      m_PVTO.undersaturatedViscosity[iSat].emplace_back( oilTable[j][2] );
      branchSize++;
      j++;
      if( j == oilTable.size())
      {
        break;
      }
    }

    // Step 2.2.3: go to the next saturated part
    i = i + branchSize;
    iSat++;
  }


  // Step 3: add 1 atm value if does not exist yet

  if( !isZero( m_PVTO.Rs[0] ) )
  {
    real64 const refPressure = 101325.0;
    real64 const viscosity = extrapolation::logExtrapolation( m_PVTO.bubblePressure[1],
                                                              m_PVTO.bubblePressure[0],
                                                              m_PVTO.saturatedViscosity[1],
                                                              m_PVTO.saturatedViscosity[0],
                                                              refPressure );


    m_PVTO.Rs.emplace( 0, 0. );
    m_PVTO.bubblePressure.emplace( 0, refPressure );
    m_PVTO.saturatedBo.emplace( 0, 1.0 );
    m_PVTO.saturatedViscosity.emplace( 0, viscosity );
    m_PVTO.numSaturatedPoints++;

    // Note: the additional undersaturated values of this branch will be created in createUndersaturatedProperties

    array1d< real64 > tmp( 1 );
    tmp[0] = 0.;
    m_PVTO.undersaturatedPressure.emplace( 0, tmp );
    tmp[0] = 1.;
    m_PVTO.undersaturatedBo.emplace( 0, tmp );
    tmp[0] = viscosity;
    m_PVTO.undersaturatedViscosity.emplace( 0, tmp );
  }


  // Step 4: find the max relative pressure (later used to extend the tables)

  m_PVTO.maxRelativePressure = 0;
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    real64 const maxP = *( std::max_element( m_PVTO.undersaturatedPressure[i].begin(), m_PVTO.undersaturatedPressure[i].end() ) );
    m_PVTO.maxRelativePressure = std::max( maxP, m_PVTO.maxRelativePressure );
  }


  // Step 5: save surface densities (for now, this info is passed twice to the kernels)

  m_PVTO.surfaceMassDensity.resize( BlackOilFluidUpdate::NC_BO-1 );
  m_PVTO.surfaceMoleDensity.resize( BlackOilFluidUpdate::NC_BO-1 );

  using PT = BlackOilFluid::PhaseType;

  m_PVTO.surfaceMassDensity[PT::OIL] = oilSurfaceMassDensity;
  m_PVTO.surfaceMoleDensity[PT::OIL] = oilSurfaceMassDensity / oilSurfaceMolecularWeight;

  m_PVTO.surfaceMassDensity[PT::GAS] = gasSurfaceMassDensity;
  m_PVTO.surfaceMoleDensity[PT::GAS] = gasSurfaceMassDensity / gasSurfaceMolecularWeight;

}

void BlackOilFluid::createUndersaturatedProperties()
{
  localIndex iBranchUp = m_PVTO.numSaturatedPoints - 1;

  // This function creates undersaturated properties ONLY IF it is not provided by the user
  // Specifically, if undersaturated data is missing, we construct undersaturated branches by interpolation, starting from second highest Rs

  for( localIndex iCurrent = m_PVTO.numSaturatedPoints - 2; iCurrent >= 0; --iCurrent )
  {
    localIndex iBranchLow = iCurrent;

    // if only the saturated part is present, then we populate the undersaturated part
    // otherwise, we do not do anything for this saturated value
    if( m_PVTO.undersaturatedPressure[iCurrent].size() == 1 )
    {

      // search for a lower undersaturated branch already created
      for( localIndex iLow = iCurrent; iLow >= 0; --iLow )
      {
        if( m_PVTO.undersaturatedPressure[iLow].size() > 1 )
        {
          iBranchLow = iLow;
          break;
        }
      }

      // if there is no undersaturated branch created yet, we start from the top branch
      if( iBranchLow == iCurrent )
      {
        iBranchLow = iBranchUp;
      }


      // Step 1: collect all the pressure values in the undersaturated pressure tables on these branches

      std::vector< real64 > allPressures;
      for( localIndex i = 0; i < m_PVTO.undersaturatedPressure[iBranchUp].size(); ++i )
      {
        allPressures.emplace_back( m_PVTO.undersaturatedPressure[iBranchUp][i] );
      }
      for( localIndex i = 0; i < m_PVTO.undersaturatedPressure[iBranchLow].size(); ++i )
      {
        allPressures.emplace_back( m_PVTO.undersaturatedPressure[iBranchLow][i] );
      }
      std::sort( allPressures.begin(), allPressures.end() );
      allPressures.erase( std::unique( allPressures.begin(), allPressures.end() ), allPressures.end() );

      // Step 2: put these pressures into an array1d (for views later), we should ultimately remove this step

      array1d< real64 > pres( allPressures.size() );
      for( localIndex i = 0; i < pres.size(); ++i )
      {
        pres[i] = allPressures[i];
      }

      // Step 3: then populate the undersaturated branch with two-step interpolation

      array1d< real64 > BoInterpUp( pres.size() );
      array1d< real64 > BoInterpLow( pres.size() );
      array1d< real64 > viscInterpUp( pres.size() );
      array1d< real64 > viscInterpLow( pres.size() );

      real64 const dRsUp = LvArray::math::abs( m_PVTO.Rs[iBranchUp] - m_PVTO.Rs[iCurrent] );
      real64 const dRsLow = LvArray::math::abs( m_PVTO.Rs[iCurrent] - m_PVTO.Rs[iBranchLow] );

      for( localIndex i = 1; i < pres.size(); ++i )
      {

        // pressure
        m_PVTO.undersaturatedPressure[iCurrent].emplace_back( pres[i] );

        arrayView1d< real64 const > presUp = m_PVTO.undersaturatedPressure[iBranchUp].toViewConst();
        arrayView1d< real64 const > presLow = m_PVTO.undersaturatedPressure[iBranchLow].toViewConst();

        // Bo
        arrayView1d< real64 const > const & BoUp = m_PVTO.undersaturatedBo[iBranchUp];
        arrayView1d< real64 const > const & BoLow = m_PVTO.undersaturatedBo[iBranchLow];
        interpolation::linearInterpolation( presUp, BoUp,
                                            pres.toViewConst(), BoInterpUp.toView() );
        interpolation::linearInterpolation( presLow, BoLow,
                                            pres.toViewConst(), BoInterpLow.toView() );
        real64 const BoSlopeUp  = ( BoInterpUp[i] - BoInterpUp[i - 1] ) / ( pres[i] - pres[i - 1] );
        real64 const BoSlopeLow = ( BoInterpLow[i] - BoInterpLow[i - 1] ) / ( pres[i] - pres[i - 1] );
        real64 const BoSlope    = ( 1 / dRsUp * BoSlopeUp + 1 / dRsLow * BoSlopeLow ) / ( 1 / dRsUp + 1 / dRsLow );
        m_PVTO.undersaturatedBo[iCurrent].emplace_back( m_PVTO.undersaturatedBo[iCurrent][i - 1] + BoSlope * ( pres[i] - pres[i - 1] ) );

        // viscosity
        arrayView1d< real64 const > const & viscUp = m_PVTO.undersaturatedViscosity[iBranchUp];
        arrayView1d< real64 const > const & viscLow = m_PVTO.undersaturatedViscosity[iBranchLow];
        interpolation::linearInterpolation( presUp, viscUp,
                                            pres.toViewConst(), viscInterpUp.toView() );
        interpolation::linearInterpolation( presLow, viscLow,
                                            pres.toViewConst(), viscInterpLow.toView() );
        real64 const viscSlopeUp  = ( viscInterpUp[i] - viscInterpUp[i - 1] ) / ( pres[i] - pres[i - 1] );
        real64 const viscSlopeLow = ( viscInterpLow[i] - viscInterpLow[i - 1] ) / ( pres[i] - pres[i - 1] );
        real64 const viscSlope    = ( 1 / dRsUp * viscSlopeUp + 1 / dRsLow * viscSlopeLow ) / ( 1 / dRsUp + 1 / dRsLow );
        m_PVTO.undersaturatedViscosity[iCurrent].emplace_back( m_PVTO.undersaturatedViscosity[iCurrent][i - 1] + viscSlope * ( pres[i] - pres[i - 1] ) );
      }

    }
    iBranchUp = iCurrent;
  }
}

void BlackOilFluid::extendUndersaturatedProperties()
{
  // for all the branches in the undersaturated tables
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    real64 const deltaPresExtended = m_PVTO.maxRelativePressure - m_PVTO.undersaturatedPressure[i].back();

    array1d< real64 > const & presUndersat = m_PVTO.undersaturatedPressure[i];
    array1d< real64 > const & BoUndersat = m_PVTO.undersaturatedBo[i];
    array1d< real64 > const & viscUndersat = m_PVTO.undersaturatedViscosity[i];

    // if the (desired) max pressure is above the largest undersaturated pressure currently in the table
    if( LvArray::math::abs( deltaPresExtended ) > 0 )
    {
      localIndex const branchSize = m_PVTO.undersaturatedPressure[i].size();

      // then we compute an extrapolated value of Bo and viscosity
      real64 const Bo = extrapolation::linearExtrapolation( presUndersat[branchSize - 2], presUndersat[branchSize - 1],
                                                            BoUndersat[branchSize - 2], BoUndersat[branchSize - 1],
                                                            m_PVTO.maxRelativePressure );
      real64 const visc = extrapolation::linearExtrapolation( presUndersat[branchSize - 2], presUndersat[branchSize - 1],
                                                              viscUndersat[branchSize - 2], viscUndersat[branchSize - 1],
                                                              m_PVTO.maxRelativePressure );

      // then add thes values to the undersaturated tables
      m_PVTO.undersaturatedPressure[i].emplace_back( m_PVTO.maxRelativePressure );
      m_PVTO.undersaturatedBo[i].emplace_back( Bo );
      m_PVTO.undersaturatedViscosity[i].emplace_back( visc );

    }
  }
}



void BlackOilFluid::refineUndersaturatedTables( localIndex const numRefinedPressurePoints )
{
  // Step 1: find the maximum pressure found in the table

  real64 maxPresInTable = 0;
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < m_PVTO.undersaturatedPressure[i].size(); ++j )
    {
      if( maxPresInTable < m_PVTO.undersaturatedPressure[i][j] )
      {
        maxPresInTable = m_PVTO.undersaturatedPressure[i][j];
      }
    }
  }

  // Step 2: populate the refinedPres table from 0 to the max pressure found above

  array1d< real64 > refinedPres( numRefinedPressurePoints );
  for( localIndex i = 0; i < numRefinedPressurePoints; ++i )
  {
    refinedPres[i] = i * ( maxPresInTable / ( numRefinedPressurePoints - 1 ) );
  }


  // Step 3: interpolate in the Bo and viscosity tables to get the refined undersaturated values

  array1d< real64 > refinedProp( refinedPres.size() );
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    // Bo
    interpolation::linearInterpolation( m_PVTO.undersaturatedPressure[i].toViewConst(), m_PVTO.undersaturatedBo[i].toViewConst(),
                                        refinedPres.toViewConst(), refinedProp.toView() );
    m_PVTO.undersaturatedBo[i] = refinedProp;

    // viscosity
    interpolation::linearInterpolation( m_PVTO.undersaturatedPressure[i].toViewConst(), m_PVTO.undersaturatedViscosity[i].toViewConst(),
                                        refinedPres.toViewConst(), refinedProp.toView() );
    m_PVTO.undersaturatedViscosity[i] = refinedProp;

    // copy pressure
    m_PVTO.undersaturatedPressure[i] = refinedPres;
  }

  // Step 4: copy undersaturated data (array1d< array1d< real64 > >) into the final 2D arrays that will be used in the kernels

  m_PVTO.undersaturatedPressure2d.resize( m_PVTO.numSaturatedPoints, refinedPres.size() );
  m_PVTO.undersaturatedBo2d.resize( m_PVTO.numSaturatedPoints, refinedPres.size() );
  m_PVTO.undersaturatedViscosity2d.resize( m_PVTO.numSaturatedPoints, refinedPres.size() );
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < refinedPres.size(); ++j )
    {
      m_PVTO.undersaturatedPressure2d[i][j] = m_PVTO.undersaturatedPressure[i][j];
      m_PVTO.undersaturatedBo2d[i][j] = m_PVTO.undersaturatedBo[i][j];
      m_PVTO.undersaturatedViscosity2d[i][j] = m_PVTO.undersaturatedViscosity[i][j];
    }
  }
}


void BlackOilFluid::checkTableConsistency() const
{
  using PT = BlackOilFluid::PhaseType;

  // check for the presence of one bubble point
  GEOSX_THROW_IF( m_PVTO.undersaturatedPressure[m_PVTO.numSaturatedPoints - 1].size() <= 1,
                  "BlackOilFluid model named " << getName() << ": At least one bubble pressure is required in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );

  // check for saturated region
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints - 1; ++i )
  {
    // Rs must increase with Pb
    GEOSX_THROW_IF( ( m_PVTO.Rs[i + 1] - m_PVTO.Rs[i] ) <= 0,
                    "BlackOilFluid model named " << getName() << ": Rs must increase with Pb in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
    // Bo must increase with Pb
    GEOSX_THROW_IF( ( m_PVTO.saturatedBo[i + 1] - m_PVTO.saturatedBo[i] ) <= 0,
                    "BlackOilFluid model named " << getName() << ": Bo must increase with Pb in saturated region in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
    // Viscosity must decrease with Pb
    GEOSX_THROW_IF( ( m_PVTO.saturatedViscosity[i + 1] - m_PVTO.saturatedViscosity[i] ) >= 0,
                    "BlackOilFluid model named " << getName() << ": Viscosity must decrease with Pb in saturated region in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
  }

  // check for under-saturated branches
  for( localIndex i = 0; i < m_PVTO.numSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < m_PVTO.undersaturatedPressure[i].size() - 1; ++j )
    {
      // Pressure
      GEOSX_THROW_IF( ( m_PVTO.undersaturatedPressure[i][j + 1] - m_PVTO.undersaturatedPressure[i][j] ) <= 0,
                      "BlackOilFluid model named " << getName() << ": P must decrease in undersaturated region in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
      // Bo must decrease with P
      GEOSX_THROW_IF( ( m_PVTO.undersaturatedBo[i][j + 1] - m_PVTO.undersaturatedBo[i][j] ) >= 0,
                      "BlackOilFluid model named " << getName() << ": Bo must decrease with P in undersaturated region in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
      // Viscosity must increase with Pb
      GEOSX_THROW_IF( ( m_PVTO.undersaturatedViscosity[i][j + 1] - m_PVTO.undersaturatedViscosity[i][j] ) < -1e-10,
                      "BlackOilFluid model named " << getName() << ": Viscosity must increase with P in undersaturated region in " << m_tableFiles[m_phaseOrder[PT::OIL]], InputError );
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
