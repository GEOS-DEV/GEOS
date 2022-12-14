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

#include "KilloughHysteresis.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

KilloughHysteresis::KernelKilloughHysteresisBase::KernelKilloughHysteresisBase( const arrayView1d< const geosx::real64 > & landParam,
                                                                                real64 const & jerauldParamA,
                                                                                real64 const & jerauldParamB,
                                                                                real64 const & killoughCurvatureParam,
                                                                                const arrayView3d< geosx::real64, cappres::USD_CAPPRES > & phaseTrappedVolFrac ):
  m_jerauldParam_a( jerauldParamA ),
  m_jerauldParam_b( jerauldParamB ),
  m_killoughCurvatureParamRelPerm( killoughCurvatureParam ),
  m_landParam( landParam )
{}



real64 KilloughHysteresis::KernelKilloughHysteresisBase::getJerauldParamA() const
{
  return m_jerauldParam_a;
}

real64 KilloughHysteresis::KernelKilloughHysteresisBase::getJerauldParamB() const
{
  return m_jerauldParam_b;
}

real64 KilloughHysteresis::KernelKilloughHysteresisBase::getCurvatureParam() const
{
  return m_killoughCurvatureParamRelPerm;
}


void KilloughHysteresis::setRelPermParameters( const geosx::real64 & jerauldA,
                                               const geosx::real64 & jerauldB,
                                               const geosx::real64 & relpermCurv )
{
  m_jerauldParam_a = jerauldA;
  m_jerauldParam_b = jerauldB;
  m_killoughCurvatureParamRelPerm = relpermCurv;
}


void KilloughHysteresis::postProcessInput()
{
  GEOSX_THROW_IF( m_jerauldParam_a < 0,
                  GEOSX_FMT( "{}: the parameter {} must be positive",
                             getCatalogName(),
                             viewKeyStruct::jerauldParameterAString() ),
                  InputError );

  GEOSX_THROW_IF( m_jerauldParam_b < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getCatalogName(),
                             viewKeyStruct::jerauldParameterBString() ),
                  InputError );

  GEOSX_THROW_IF( m_killoughCurvatureParamRelPerm < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getCatalogName(),
                             viewKeyStruct::killoughCurvatureParameterString() ),
                  InputError );

}

KilloughHysteresis::KernelKilloughHysteresisBase
KilloughHysteresis::createKernelWrapper( const arrayView1d< const geosx::real64 > & landParam,
                                         const arrayView3d< geosx::real64, relperm::USD_RELPERM > & phaseTrappedVolFrac ) const
{
  return KilloughHysteresis::KernelKilloughHysteresisBase( landParam,
                                                           m_jerauldParam_a,
                                                           m_jerauldParam_b,
                                                           m_killoughCurvatureParamRelPerm,
                                                           phaseTrappedVolFrac );
}

//TODO
void KilloughHysteresis::computeLandCoefficient( KilloughHysteresis::HysteresisCurve_t const & hcurve,
                                                 real64 & landParam )
{

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase
  if( hcurve.isWetting() )
  {
    real64 const Scrd = hcurve.oppositeBoundSat;
    real64 const Smxd = hcurve.drainageExtremaSat;
    real64 const Smxi = hcurve.imbibitionExtremaSat;
    real64 const Swc = Scrd;
    GEOSX_THROW_IF(  (Smxi - Smxd) > 0,
                     GEOSX_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                                "Crossing relative permeability curves.\n",
                                getCatalogName(),
                                Smxi,
                                Smxd ),
                     InputError );

    landParam = ( Smxd - Swc ) / LvArray::math::max( KilloughHysteresis::KernelKilloughHysteresisBase::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
  }
  else
  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Smx =  hcurve.oppositeBoundSat;
    real64 const Scrd = hcurve.drainageExtremaSat;
    real64 const Scri = hcurve.imbibitionExtremaSat;
    GEOSX_THROW_IF( (Scrd - Scri) > 0,
                    GEOSX_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                               "Crossing relative permeability curves.\n",
                               getCatalogName(),
                               Scrd,
                               Scri ),
                    InputError );

    landParam = ( Smx - Scrd ) / LvArray::math::max( KilloughHysteresis::KernelKilloughHysteresisBase::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
  }
}



}//end namespace
}//end namespace
