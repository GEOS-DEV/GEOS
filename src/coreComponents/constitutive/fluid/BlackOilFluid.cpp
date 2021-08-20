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

void BlackOilFluid::createAllKernelWrappers()
{
  BlackOilFluidBase::createAllKernelWrappers();

  // create the kernel wrapper
  m_PVTO.m_kernelWrapper.create( m_PVTO.m_Rs.toViewConst(),
                                 m_PVTO.m_bubblePressure.toViewConst(),
                                 m_PVTO.m_saturatedBo.toViewConst(),
                                 m_PVTO.m_saturatedViscosity.toViewConst(),
                                 m_PVTO.m_undersaturatedPressure2d.toViewConst(),
                                 m_PVTO.m_undersaturatedBo2d.toViewConst(),
                                 m_PVTO.m_undersaturatedViscosity2d.toViewConst(),
                                 m_PVTO.m_surfaceMassDensity.toViewConst(),
                                 m_PVTO.m_surfaceMoleDensity.toViewConst() );
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

  model.m_PVTO.m_Rs = m_PVTO.m_Rs;
  model.m_PVTO.m_bubblePressure = m_PVTO.m_bubblePressure;
  model.m_PVTO.m_saturatedBo = m_PVTO.m_saturatedBo;
  model.m_PVTO.m_saturatedViscosity = m_PVTO.m_saturatedViscosity;
  model.m_PVTO.m_undersaturatedPressure2d = m_PVTO.m_undersaturatedPressure2d;
  model.m_PVTO.m_undersaturatedBo2d = m_PVTO.m_undersaturatedBo2d;
  model.m_PVTO.m_undersaturatedViscosity2d = m_PVTO.m_undersaturatedViscosity2d;
  model.m_PVTO.m_surfaceMassDensity = m_PVTO.m_surfaceMassDensity;
  model.m_PVTO.m_surfaceMoleDensity = m_PVTO.m_surfaceMoleDensity;

  model.createAllKernelWrappers();
  return clone;
}

void BlackOilFluid::useProvidedTableFunctions()
{
  GEOSX_ERROR( "BlackOilFluid: this option is not implemented yet" );
}

void BlackOilFluid::readInputDataFromPVTFiles()
{
  using PT = BlackOilFluid::PhaseType;

  PVTProps::BlackOilTables boTables;
  boTables.buildAllTables( PT::OIL,
                           PT::WATER,
                           PT::GAS,
                           m_phaseTypes,
                           m_tableFiles );

  // water data
  array1d< array1d< real64 > > data( 1 );
  data[0] = boTables.m_waterTable;
  fillWaterData( data );

  // gas data
  fillHydrocarbonData( PT::GAS, boTables.m_gasTable );

  // for the Black-Oil model, the oil PVT is treated differently from gas
  fillPVTOData( boTables.m_oilTable,
                m_surfacePhaseMassDensity[m_phaseOrder[PT::OIL]],
                m_componentMolarWeight[m_phaseOrder[PT::OIL]],
                m_surfacePhaseMassDensity[m_phaseOrder[PT::GAS]],
                m_componentMolarWeight[m_phaseOrder[PT::GAS]] );

  // check consistency
  checkTableConsistency();

  // create missing undersaturated branches
  createUnderSaturatedProperties();

  // extend existing branches
  extendUnderSaturatedProperties();

  // refine table branches
  refineTableAndCopy( 100 );

  // check consistency
  checkTableConsistency();
}

void BlackOilFluid::fillPVTOData( array1d< array1d< real64 > > const & PVT,
                                  real64 oilSurfaceMassDensity,
                                  real64 oilSurfaceMolecularWeight,
                                  real64 gasSurfaceMassDensity,
                                  real64 gasSurfaceMolecularWeight )
{
  // count number of saturated points
  localIndex nSaturatedPoints = 0;
  for( localIndex i = 0; i < PVT.size(); ++i )
  {
    if( PVT[i].size() == 4 )
    {
      nSaturatedPoints++;
    }
  }
  m_PVTO.m_nSaturatedPoints = nSaturatedPoints;

  // store data in a new structure
  m_PVTO.m_Rs.resize( nSaturatedPoints );
  m_PVTO.m_bubblePressure.resize( nSaturatedPoints );
  m_PVTO.m_saturatedBo.resize( nSaturatedPoints );
  m_PVTO.m_saturatedViscosity.resize( nSaturatedPoints );
  m_PVTO.m_undersaturatedPressure.resize( nSaturatedPoints );
  m_PVTO.m_undersaturatedBo.resize( nSaturatedPoints );
  m_PVTO.m_undersaturatedViscosity.resize( nSaturatedPoints );

  localIndex iSat = 0;
  for( localIndex i = 0; i < PVT.size(); ++i )
  {
    // saturated part
    m_PVTO.m_Rs[iSat] = PVT[i][0];
    m_PVTO.m_bubblePressure[iSat] = PVT[i][1];
    m_PVTO.m_saturatedBo[iSat] = PVT[i][2];
    m_PVTO.m_saturatedViscosity[iSat] = PVT[i][3];

    // add unsaturated properties: pressure = P - Pbub
    m_PVTO.m_undersaturatedPressure[iSat].emplace_back( 0 );
    m_PVTO.m_undersaturatedBo[iSat].emplace_back( m_PVTO.m_saturatedBo[iSat] );
    m_PVTO.m_undersaturatedViscosity[iSat].emplace_back( m_PVTO.m_saturatedViscosity[iSat] );

    localIndex branchSize = 0;
    localIndex j = i + 1;
    while( PVT[j].size() == 3 )
    {
      m_PVTO.m_undersaturatedPressure[iSat].emplace_back( PVT[j][0] - m_PVTO.m_bubblePressure[iSat] );
      m_PVTO.m_undersaturatedBo[iSat].emplace_back( PVT[j][1] );
      m_PVTO.m_undersaturatedViscosity[iSat].emplace_back( PVT[j][2] );
      branchSize++;
      j++;
      if( j == PVT.size())
      {
        break;
      }
    }
    i = i + branchSize;
    iSat++;
  }

  // add 1 atm value if does not exist yet
  if( !isZero( m_PVTO.m_Rs[0] ) )
  {
    real64 const refPressure = 101325.0;
    real64 const viscosity = logExtrapolation( m_PVTO.m_bubblePressure[1],
                                               m_PVTO.m_saturatedViscosity[1],
                                               m_PVTO.m_bubblePressure[0],
                                               m_PVTO.m_saturatedViscosity[0],
                                               refPressure );

    m_PVTO.m_Rs.emplace( 0, 0. );
    m_PVTO.m_bubblePressure.emplace( 0, refPressure );
    m_PVTO.m_saturatedBo.emplace( 0, 1.0 );
    m_PVTO.m_saturatedViscosity.emplace( 0, viscosity );
    m_PVTO.m_nSaturatedPoints++;

    array1d< real64 > tmp( 1 );
    tmp[0] = 0.;
    m_PVTO.m_undersaturatedPressure.emplace( 0, tmp );
    tmp[0] = 1.;
    m_PVTO.m_undersaturatedBo.emplace( 0, tmp );
    tmp[0] = viscosity;
    m_PVTO.m_undersaturatedViscosity.emplace( 0, tmp );
  }

  m_PVTO.m_maxRelativePressure = 0;
  m_PVTO.m_minRelativePressure = 1e8;

  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {
    real64 const maxP = *( std::max_element( m_PVTO.m_undersaturatedPressure[i].begin(), m_PVTO.m_undersaturatedPressure[i].end() ) );
    real64 const minP = *( std::min_element( m_PVTO.m_undersaturatedPressure[i].begin(), m_PVTO.m_undersaturatedPressure[i].end() ) );
    m_PVTO.m_maxRelativePressure = std::max( maxP, m_PVTO.m_maxRelativePressure );
    m_PVTO.m_minRelativePressure = std::min( minP, m_PVTO.m_minRelativePressure );
  }

  // densities
  m_PVTO.m_surfaceMassDensity.resize( NC_BO-1 );
  m_PVTO.m_surfaceMoleDensity.resize( NC_BO-1 );

  using PT = BlackOilFluid::PhaseType;

  m_PVTO.m_surfaceMassDensity[PT::OIL] = oilSurfaceMassDensity;
  m_PVTO.m_surfaceMoleDensity[PT::OIL] = oilSurfaceMassDensity / oilSurfaceMolecularWeight;

  m_PVTO.m_surfaceMassDensity[PT::GAS] = gasSurfaceMassDensity;
  m_PVTO.m_surfaceMoleDensity[PT::GAS] = gasSurfaceMassDensity / gasSurfaceMolecularWeight;

}

void BlackOilFluid::extendUnderSaturatedProperties()
{
  // extrapolate m_undersaturated properties up to max pressure
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {
    real64 const dPext = m_PVTO.m_maxRelativePressure - m_PVTO.m_undersaturatedPressure[i].back();

    array1d< real64 > const & Pusat = m_PVTO.m_undersaturatedPressure[i];
    array1d< real64 > const & Bousat = m_PVTO.m_undersaturatedBo[i];
    array1d< real64 > const & viscusat = m_PVTO.m_undersaturatedViscosity[i];

    if( std::fabs( dPext ) > 0 )
    {
      localIndex const branchSize = m_PVTO.m_undersaturatedPressure[i].size();
      real64 const Bo = linearExtrapolation( Pusat[branchSize - 2],
                                             Bousat[branchSize - 2],
                                             Pusat[branchSize - 1],
                                             Bousat[branchSize - 1],
                                             m_PVTO.m_maxRelativePressure );
      real64 const visc = linearExtrapolation( Pusat[branchSize - 2],
                                               viscusat[branchSize - 2],
                                               Pusat[branchSize - 1],
                                               viscusat[branchSize - 1],
                                               m_PVTO.m_maxRelativePressure );
      m_PVTO.m_undersaturatedBo[i].emplace_back( Bo );
      m_PVTO.m_undersaturatedViscosity[i].emplace_back( visc );
      m_PVTO.m_undersaturatedPressure[i].emplace_back( m_PVTO.m_maxRelativePressure );
    }
  }
}

void BlackOilFluid::createUnderSaturatedProperties()
{
  localIndex upperBranchIndex = m_PVTO.m_nSaturatedPoints - 1;

  // construct undersaturated branches by interpolation, start from second highest Rs

  for( localIndex iCurrent = m_PVTO.m_nSaturatedPoints - 2; iCurrent-- > 0; )
  {
    localIndex lowerBranchIndex = iCurrent;
    if( m_PVTO.m_undersaturatedPressure[iCurrent].size() == 1 )  // only saturated part is present then create missing part
    {
      for( localIndex iLower = iCurrent; iLower-- > 0; )
      {
        if( m_PVTO.m_undersaturatedPressure[iLower].size() > 1 ) // search for lower under-saturated branch
        {
          lowerBranchIndex = iLower;
          break;
        }
      }

      if( lowerBranchIndex == iCurrent )
      {
        lowerBranchIndex = upperBranchIndex;
      }

      real64 dRs_up = std::fabs( m_PVTO.m_Rs[upperBranchIndex] - m_PVTO.m_Rs[iCurrent] );
      real64 dRs_dn = std::fabs( m_PVTO.m_Rs[iCurrent] - m_PVTO.m_Rs[lowerBranchIndex] );

      // generate merge of pressures
      std::vector< real64 > pTarget; // TODO: change to array1d
      for( localIndex i = 0; i < m_PVTO.m_undersaturatedPressure[upperBranchIndex].size(); ++i )
      {
        pTarget.emplace_back( m_PVTO.m_undersaturatedPressure[upperBranchIndex][i] );
      }
      for( localIndex i = 0; i < m_PVTO.m_undersaturatedPressure[lowerBranchIndex].size(); ++i )
      {
        pTarget.emplace_back( m_PVTO.m_undersaturatedPressure[lowerBranchIndex][i] );
      }
      std::sort( pTarget.begin(), pTarget.end() );
      pTarget.erase( std::unique( pTarget.begin(), pTarget.end() ), pTarget.end() );

      // shift pressures up and down
      array1d< real64 > x_up( m_PVTO.m_undersaturatedPressure[upperBranchIndex].size() );
      array1d< real64 > x_dn( m_PVTO.m_undersaturatedPressure[lowerBranchIndex].size() );
      for( localIndex i = 0; i < m_PVTO.m_undersaturatedPressure[upperBranchIndex].size(); ++i )
      {
        x_up[i] = m_PVTO.m_undersaturatedPressure[upperBranchIndex][i];
      }

      for( localIndex i = 0; i < m_PVTO.m_undersaturatedPressure[lowerBranchIndex].size(); ++i )
      {
        x_dn[i] = m_PVTO.m_undersaturatedPressure[lowerBranchIndex][i];
      }

      // create branch
      array1d< real64 > BoInterp_up, BoInterp_dn, viscInterp_up, viscInterp_dn;
      array1d< real64 > ppTarget( pTarget.size() );
      for( localIndex i = 0; i < localIndex( pTarget.size() ); ++i ) // TODO: remove cast
      {
        ppTarget[i] = pTarget[i];
      }
      for( localIndex i = 1; i < ppTarget.size(); ++i )
      {

        m_PVTO.m_undersaturatedPressure[iCurrent].emplace_back( ppTarget[i] );

        // Bo
        array1d< real64 > const & Bo_up = m_PVTO.m_undersaturatedBo[upperBranchIndex];
        array1d< real64 > const & Bo_dn = m_PVTO.m_undersaturatedBo[lowerBranchIndex];
        BoInterp_up.resize( ppTarget.size() );
        BoInterp_dn.resize( ppTarget.size() );
        interpolation1( x_up, Bo_up, ppTarget, BoInterp_up );
        interpolation1( x_dn, Bo_dn, ppTarget, BoInterp_dn );
        real64 const BoSlope_up = ( BoInterp_up[i] - BoInterp_up[i - 1] ) / ( ppTarget[i] - ppTarget[i - 1] );
        real64 const BoSlope_dn = ( BoInterp_dn[i] - BoInterp_dn[i - 1] ) / ( ppTarget[i] - ppTarget[i - 1] );
        real64 const BoSlope = ( 1 / dRs_up * BoSlope_up + 1 / dRs_dn * BoSlope_dn ) / ( 1 / dRs_up + 1 / dRs_dn );
        m_PVTO.m_undersaturatedBo[iCurrent].emplace_back( m_PVTO.m_undersaturatedBo[iCurrent][i - 1] + BoSlope * ( ppTarget[i] - ppTarget[i - 1] ) );

        // visc
        array1d< real64 > const & visc_up = m_PVTO.m_undersaturatedViscosity[upperBranchIndex];
        array1d< real64 > const & visc_dn = m_PVTO.m_undersaturatedViscosity[lowerBranchIndex];
        viscInterp_up.resize( ppTarget.size() );
        viscInterp_dn.resize( ppTarget.size() );
        interpolation1( x_up, visc_up, ppTarget, viscInterp_up );
        interpolation1( x_dn, visc_dn, ppTarget, viscInterp_dn );
        real64 const viscSlope_up = ( viscInterp_up[i] - viscInterp_up[i - 1] ) / ( ppTarget[i] - ppTarget[i - 1] );
        real64 const viscSlope_dn = ( viscInterp_dn[i] - viscInterp_dn[i - 1] ) / ( ppTarget[i] - ppTarget[i - 1] );
        real64 const viscSlope = ( 1 / dRs_up * viscSlope_up + 1 / dRs_dn * viscSlope_dn ) / ( 1 / dRs_up + 1 / dRs_dn );
        m_PVTO.m_undersaturatedViscosity[iCurrent].emplace_back( m_PVTO.m_undersaturatedViscosity[iCurrent][i - 1] + viscSlope * ( ppTarget[i] - ppTarget[i - 1] ) );
      }

    }
    upperBranchIndex = iCurrent;
  }
}

void BlackOilFluid::refineTableAndCopy( localIndex const nLevels )
{
  std::vector< real64 > refineP; // TODO:change to array1d
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < m_PVTO.m_undersaturatedPressure[i].size(); ++j )
    {
      refineP.emplace_back( m_PVTO.m_undersaturatedPressure[i][j] );
    }
  }
  std::sort( refineP.begin(), refineP.end() );
  refineP.erase( std::unique( refineP.begin(), refineP.end() ), refineP.end() );

  // linspace
  real64 const maxVal = *std::max_element( refineP.begin(), refineP.end() );
  array1d< real64 > refinePP( refineP.size() );
  for( localIndex i = 0; i < localIndex( refineP.size() ); ++i )
  {
    refinePP[i] = refineP[i];
  }

  refinePP.resize( nLevels );
  for( localIndex i = 0; i != nLevels; ++i )
  {
    refinePP[i] = i * ( maxVal / ( nLevels - 1 ) );
  }

  localIndex const ppSize = refinePP.size();
  array1d< real64 > propRefined( ppSize );
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {

    // Bo
    interpolation1( m_PVTO.m_undersaturatedPressure[i], m_PVTO.m_undersaturatedBo[i], refinePP, propRefined );
    m_PVTO.m_undersaturatedBo[i].resize( ppSize );
    for( localIndex j = 0; j < ppSize; ++j )
    {
      m_PVTO.m_undersaturatedBo[i][j] = propRefined[j];
    }

    // viscosity
    interpolation1( m_PVTO.m_undersaturatedPressure[i], m_PVTO.m_undersaturatedViscosity[i], refinePP, propRefined );
    m_PVTO.m_undersaturatedViscosity[i].resize( ppSize );
    for( localIndex j = 0; j < ppSize; ++j )
    {
      m_PVTO.m_undersaturatedViscosity[i][j] = propRefined[j];
    }

    for( localIndex j = 0; j < ppSize; ++j )
    {
      m_PVTO.m_undersaturatedPressure[i].resize( ppSize );
      m_PVTO.m_undersaturatedPressure[i][j] = refinePP[j];
    }
  }

  // copy undersaturated into array2d
  m_PVTO.m_undersaturatedPressure2d.resize( m_PVTO.m_nSaturatedPoints, ppSize );
  m_PVTO.m_undersaturatedBo2d.resize( m_PVTO.m_nSaturatedPoints, ppSize );
  m_PVTO.m_undersaturatedViscosity2d.resize( m_PVTO.m_nSaturatedPoints, ppSize );
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < ppSize; ++j )
    {
      m_PVTO.m_undersaturatedPressure2d[i][j] = m_PVTO.m_undersaturatedPressure[i][j];
      m_PVTO.m_undersaturatedBo2d[i][j] = m_PVTO.m_undersaturatedBo[i][j];
      m_PVTO.m_undersaturatedViscosity2d[i][j] = m_PVTO.m_undersaturatedViscosity[i][j];
    }
  }
}

void BlackOilFluid::checkTableConsistency() const
{
  // check for the presence of one bubble point
  GEOSX_THROW_IF( m_PVTO.m_undersaturatedPressure[m_PVTO.m_nSaturatedPoints - 1].size() <= 1,
                  "At least one bubble pressure is required", InputError );

  // check for saturated region
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints - 1; ++i )
  {
    // Rs must increase with Pb
    GEOSX_THROW_IF( ( m_PVTO.m_Rs[i + 1] - m_PVTO.m_Rs[i] ) <= 0,
                    "Rs must increase with Pb", InputError );
    // Bo must increase with Pb
    GEOSX_THROW_IF( ( m_PVTO.m_saturatedBo[i + 1] - m_PVTO.m_saturatedBo[i] ) <= 0,
                    "Bo must increase with Pb in saturated region", InputError );
    // Viscosity must decrease with Pb
    GEOSX_THROW_IF( ( m_PVTO.m_saturatedViscosity[i + 1] - m_PVTO.m_saturatedViscosity[i] ) >= 0,
                    "Viscosity must decrease with Pb in saturated region", InputError );
  }

  // check for under-saturated branches
  for( localIndex i = 0; i < m_PVTO.m_nSaturatedPoints; ++i )
  {
    for( localIndex j = 0; j < m_PVTO.m_undersaturatedPressure[i].size() - 1; ++j )
    {
      // Pressure
      GEOSX_THROW_IF( ( m_PVTO.m_undersaturatedPressure[i][j + 1] - m_PVTO.m_undersaturatedPressure[i][j] ) <= 0,
                      "P must decrease in undersaturated region", InputError );
      // Bo must decrease with P
      GEOSX_THROW_IF( ( m_PVTO.m_undersaturatedBo[i][j + 1] - m_PVTO.m_undersaturatedBo[i][j] ) >= 0,
                      "Bo must decrease with P in undersaturated region", InputError );
      // Viscosity must increase with Pb
      GEOSX_THROW_IF( ( m_PVTO.m_undersaturatedViscosity[i][j + 1] - m_PVTO.m_undersaturatedViscosity[i][j] ) < -1e-10,
                      "Viscosity must increase with P in undersaturated region", InputError );
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
