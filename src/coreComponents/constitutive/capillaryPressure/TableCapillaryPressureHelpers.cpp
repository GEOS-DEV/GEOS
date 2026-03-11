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
 * @file TableCapillaryPressureHelpers.cpp
 */

#include "TableCapillaryPressureHelpers.hpp"
#include "CapillaryPressureBase.hpp"

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

void
TableCapillaryPressureHelpers::validateCapillaryPressureTable( TableFunction const & capPresTable,
                                                               string const &,
                                                               bool const capPresMustBeIncreasing )
{
  ArrayOfArraysView< real64 const > coords = capPresTable.getCoordinates();

  GEOS_THROW_IF_NE_MSG( capPresTable.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                        GEOS_FMT( "in table '{}' interpolation method must be linear", capPresTable.getName() ),
                        InputError, capPresTable.getDataContext() );
  GEOS_THROW_IF_NE_MSG( capPresTable.numDimensions(), 1,
                        GEOS_FMT( "table '{}' must have a single independent coordinate", capPresTable.getName() ),
                        InputError, capPresTable.getDataContext() );
  GEOS_THROW_IF_LT_MSG( coords.sizeOfArray( 0 ), 2,
                        GEOS_FMT( "table `{}` must contain at least two values", capPresTable.getName() ),
                        InputError, capPresTable.getDataContext() );

  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const capPres = capPresTable.getValues();

  for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
  {
    // check phase volume fraction
    GEOS_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                   GEOS_FMT( "in table '{}' values must be between 0 and 1", capPresTable.getName() ),
                   InputError, capPresTable.getDataContext() );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check the monotonicity of the capillary pressure table
    if( capPresMustBeIncreasing )
    {
      GEOS_THROW_IF( !isZero( capPres[i] ) && (capPres[i] - capPres[i-1]) < -1e-15,
                     "values must be strictly increasing (i.e. |Delta Pc| > 1e-15 between two non-zero values)",
                     InputError, capPresTable.getDataContext() );
    }
    else
    {
      GEOS_THROW_IF( !isZero( capPres[i] ) && (capPres[i] - capPres[i-1]) > 1e-15,
                     "values must be strictly decreasing  (i.e. |Delta Pc| > 1e-15 between two non-zero values)",
                     InputError, capPresTable.getDataContext() );
    }
  }
}

void TableCapillaryPressureHelpers::populateMinPhaseVolumeFraction(
  arraySlice1d< integer const > const phaseOrder,
  TableFunction const & capPresTable,
  arraySlice1d< real64 > minPhaseVolumeFraction )
{
  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = phaseOrder[PT::WATER];
  integer const ipOil   = phaseOrder[PT::OIL];
  integer const ipGas   = phaseOrder[PT::GAS];

  ArrayOfArraysView< real64 const > coords = capPresTable.getCoordinates();
  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  real64 const minSaturation = phaseVolFrac[0];
  real64 const maxSaturation = phaseVolFrac[phaseVolFrac.size()-1];
  if( ipWater < 0 )
  {
    minPhaseVolumeFraction[ipGas] = minSaturation;
    minPhaseVolumeFraction[ipOil] = 1.0 - maxSaturation;
  }
  else
  {
    minPhaseVolumeFraction[ipWater] = minSaturation;
    if( 0 <= ipGas )
    {
      minPhaseVolumeFraction[ipGas] = 1.0 - maxSaturation;
    }
    else
    {
      minPhaseVolumeFraction[ipOil] = 1.0 - maxSaturation;
    }
  }
}

void TableCapillaryPressureHelpers::populateMinPhaseVolumeFraction(
  arraySlice1d< integer const > const phaseOrder,
  TableFunction const & capPresTableWettingIntermediate,
  TableFunction const & capPresTableNonWettingIntermediate,
  arraySlice1d< real64 > minPhaseVolumeFraction )
{
  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = phaseOrder[PT::WATER];
  integer const ipOil   = phaseOrder[PT::OIL];
  integer const ipGas   = phaseOrder[PT::GAS];

  ArrayOfArraysView< real64 const > coords = capPresTableWettingIntermediate.getCoordinates();
  arraySlice1d< real64 const > wettingSaturation = coords[0];
  real64 const minWettingSaturation = wettingSaturation[0];
  coords = capPresTableNonWettingIntermediate.getCoordinates();
  arraySlice1d< real64 const > nonWettingSaturation = coords[0];
  real64 const minNonWettingSaturation = nonWettingSaturation[0];

  minPhaseVolumeFraction[ipWater] = minWettingSaturation;
  minPhaseVolumeFraction[ipGas] = minNonWettingSaturation;
  minPhaseVolumeFraction[ipOil] = 0.0;
}

void
TableCapillaryPressureHelpers::validateCapillaryPressureTable( const geos::TableFunction & capPresTable,
                                                               const geos::string & fullConstitutiveName,
                                                               const bool capPresMustBeIncreasing,
                                                               geos::real64 & phaseMax,
                                                               geos::real64 & phaseMin,
                                                               geos::real64 & phaseCapPresMinEndPoint,
                                                               geos::real64 & phaseCapPresMaxEndPoint )
{

  TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, fullConstitutiveName, capPresMustBeIncreasing );
  ArrayOfArraysView< real64 const > coords = capPresTable.getCoordinates();
  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const capPres = capPresTable.getValues();

  phaseMin = phaseVolFrac[0];
  phaseCapPresMinEndPoint = capPres[0];
  phaseMax = phaseVolFrac[phaseVolFrac.size()-1];
  phaseCapPresMaxEndPoint = capPres[phaseVolFrac.size()-1];


  if( capPresMustBeIncreasing )
  {

    for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
    {
      if( isZero( capPres[i - 1] ) && !isZero( capPres[i] ))
      {
        phaseMin = phaseVolFrac[i - 1];
        phaseCapPresMinEndPoint = capPres[i - 1];
      }
    }
  }
  else
  {
    for( localIndex i = coords.sizeOfArray( 0 )-2; i>0; --i )
    {
      if( isZero( capPres[i + 1] ) && !isZero( capPres[i] ))
      {
        phaseMax = phaseVolFrac[i + 1];
        phaseCapPresMaxEndPoint = capPres[i + 1];
      }
    }

  }
}


} // namespace constitutive

} // namespace geos
