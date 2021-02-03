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

/**
 * @file MultiPhaseMultiComponentFluid.cpp
 */
#include "MultiPhaseMultiComponentFluid.hpp"

#include "common/Path.hpp"
#include "managers/ProblemManager.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "PVTFunctions/FlashModelBase.hpp"
#include "PVTFunctions/PVTFunctionBase.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace PVTProps;
using namespace stringutilities;

namespace constitutive
{

MultiPhaseMultiComponentFluid::MultiPhaseMultiComponentFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::phasePVTParaFilesString, &m_phasePVTParaFiles )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "List of the names of the files including PVT function parameters" );

  registerWrapper( viewKeyStruct::flashModelParaFileString, &m_flashModelParaFile )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "name of the filen including flash calculation function parameters" );

}

MultiPhaseMultiComponentFluid::~MultiPhaseMultiComponentFluid()
{}


std::unique_ptr< ConstitutiveBase >
MultiPhaseMultiComponentFluid::deliverClone( string const & name,
                                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  MultiPhaseMultiComponentFluid * const newConstitutiveRelation = dynamic_cast< MultiPhaseMultiComponentFluid * >(clone.get());


  newConstitutiveRelation->m_useMass = this->m_useMass;

  newConstitutiveRelation->m_componentNames       = this->m_componentNames;
  newConstitutiveRelation->m_componentMolarWeight = this->m_componentMolarWeight;

  newConstitutiveRelation->m_phaseNames           = this->m_phaseNames;

  newConstitutiveRelation->m_phasePVTParaFiles       = this->m_phasePVTParaFiles;

  newConstitutiveRelation->m_flashModelParaFile = this->m_flashModelParaFile;

  //  newConstitutiveRelation->CreatePVTModels();

  newConstitutiveRelation->m_phaseDensityFuns = this->m_phaseDensityFuns;
  newConstitutiveRelation->m_phaseViscosityFuns = this->m_phaseViscosityFuns;

  newConstitutiveRelation->m_flashModel = this->m_flashModel;

  return clone;
}

void MultiPhaseMultiComponentFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF( m_phasePVTParaFiles.size() != NP, "The number of phasePVTParaFiles is not the same as the number of phases!" );

  createPVTModels();

}

void MultiPhaseMultiComponentFluid::initializePostSubGroups( Group * const group )
{
  MultiFluidBase::initializePostSubGroups( group );

  //  CreatePVTModels();

}


void MultiPhaseMultiComponentFluid::createPVTModels()
{
  for( string & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );

    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ))
    {
      string const str( buf );
      string_array const strs = Tokenize( str, " " );

      if( strs[0] == "DensityFun" )
      {
        m_phaseDensityFuns.emplace_back( PVTFunction::CatalogInterface::factory( strs[ 1 ],
                                                                                 strs,
                                                                                 m_componentNames,
                                                                                 m_componentMolarWeight ) );
      }
      else if( strs[0] == "ViscosityFun" )
      {
        m_phaseViscosityFuns.emplace_back( PVTFunction::CatalogInterface::factory( strs[ 1 ],
                                                                                   strs,
                                                                                   m_componentNames,
                                                                                   m_componentMolarWeight ) );
      }
      else
      {
        GEOSX_ERROR( "Error: Invalid PVT function: " << strs[0] << "." );
      }
    }

    is.close();
  }

  {
    std::ifstream is( m_flashModelParaFile );

    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ))
    {
      string const str( buf );
      string_array const strs = Tokenize( str, " " );

      if( strs[0] == "FlashModel" )
      {
        m_flashModel = FlashModel::CatalogInterface::factory( strs[1],
                                                              strs,
                                                              m_phaseNames,
                                                              m_componentNames,
                                                              m_componentMolarWeight );
      }
      else
      {
        GEOSX_ERROR( "Error: Not flash model: " << strs[0] << "." );
      }
    }

    is.close();
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MultiPhaseMultiComponentFluid, string const &, Group * const )

void MultiPhaseMultiComponentFluidUpdate::compute( real64 pressure,
                                                   real64 temperature,
                                                   arraySlice1d< real64 const > const & composition,
                                                   arraySlice1d< real64 > const & phaseFraction,
                                                   arraySlice1d< real64 > const & phaseDensity,
                                                   arraySlice1d< real64 > const & phaseMassDensity,
                                                   arraySlice1d< real64 > const & phaseViscosity,
                                                   arraySlice2d< real64 > const & phaseCompFraction,
                                                   real64 & totalDensity ) const
{
  GEOSX_UNUSED_VAR( pressure )
  GEOSX_UNUSED_VAR( temperature )
  GEOSX_UNUSED_VAR( composition )
  GEOSX_UNUSED_VAR( phaseFraction )
  GEOSX_UNUSED_VAR( phaseDensity )
  GEOSX_UNUSED_VAR( phaseMassDensity )
  GEOSX_UNUSED_VAR( phaseViscosity )
  GEOSX_UNUSED_VAR( phaseCompFraction )
  GEOSX_UNUSED_VAR( totalDensity )
  GEOSX_ERROR( "Not implemented" );
}

void MultiPhaseMultiComponentFluidUpdate::compute( real64 pressure,
                                                   real64 temperature,
                                                   arraySlice1d< real64 const > const & composition,
                                                   arraySlice1d< real64 > const & phaseFraction,
                                                   arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                                                   arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                                                   arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                                                   arraySlice1d< real64 > const & phaseDensity,
                                                   arraySlice1d< real64 > const & dPhaseDensity_dPressure,
                                                   arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
                                                   arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                                                   arraySlice1d< real64 > const & phaseMassDensity,
                                                   arraySlice1d< real64 > const & dPhaseMassDensity_dPressure,
                                                   arraySlice1d< real64 > const & dPhaseMassDensity_dTemperature,
                                                   arraySlice2d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
                                                   arraySlice1d< real64 > const & phaseViscosity,
                                                   arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
                                                   arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
                                                   arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                                                   arraySlice2d< real64 > const & phaseCompFraction,
                                                   arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                                                   arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                                                   arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                                                   real64 & totalDensity,
                                                   real64 & dTotalDensity_dPressure,
                                                   real64 & dTotalDensity_dTemperature,
                                                   arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const
{
  CompositionalVarContainer< 1 > phaseFrac {
    phaseFraction,
    dPhaseFraction_dPressure,
    dPhaseFraction_dTemperature,
    dPhaseFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseDens {
    phaseDensity,
    dPhaseDensity_dPressure,
    dPhaseDensity_dTemperature,
    dPhaseDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseMassDens {
    phaseMassDensity,
    dPhaseMassDensity_dPressure,
    dPhaseMassDensity_dTemperature,
    dPhaseMassDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseVisc {
    phaseViscosity,
    dPhaseViscosity_dPressure,
    dPhaseViscosity_dTemperature,
    dPhaseViscosity_dGlobalCompFraction
  };

  CompositionalVarContainer< 2 > phaseCompFrac {
    phaseCompFraction,
    dPhaseCompFraction_dPressure,
    dPhaseCompFraction_dTemperature,
    dPhaseCompFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 0 > totalDens {
    totalDensity,
    dTotalDensity_dPressure,
    dTotalDensity_dTemperature,
    dTotalDensity_dGlobalCompFraction
  };

#if defined(__CUDACC__)
  // For some reason nvcc thinks these aren't used.
  GEOSX_UNUSED_VAR( phaseFrac, phaseDens, phaseMassDens, phaseVisc, phaseCompFrac, totalDens );
#endif

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  localIndex const NC = numComponents();
  localIndex const NP = numPhases();

  stackArray1d< EvalVarArgs, maxNumComp > C( NC );

  if( m_useMass )
  {
    stackArray1d< EvalVarArgs, maxNumComp > X( NC );
    EvalVarArgs totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      X[ic].m_var = composition[ic];
      X[ic].m_der[ic+1] = 1.0;

      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      C[ic] = X[ic] * mwInv; // this is molality (units of mole/mass)
      totalMolality += C[ic];
    }

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      C[ic] /= totalMolality;
    }
  }
  else
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      C[ic].m_var = composition[ic];
      C[ic].m_der[ic+1] = 1.0;
    }
  }

  EvalVarArgs P =  pressure;
  P.m_der[0] = 1.0;

  static real64 TK = 273.15;
  EvalVarArgs T =  temperature - TK;

  stackArray1d< EvalVarArgs, maxNumPhase > phaseFractionTemp( NP );
  stackArray2d< EvalVarArgs, maxNumPhase * maxNumComp > phaseCompFractionTemp( NP, NC );

  //phaseFractionTemp and phaseCompFractionTemp all are mole fraction,
  //w.r.t mole fraction or mass fraction (useMass)
  m_flashModel->partition( P, T, C, phaseFractionTemp, phaseCompFractionTemp );

  stackArray1d< EvalVarArgs, maxNumPhase > phaseDensityTemp( NP );
  stackArray1d< EvalVarArgs, maxNumPhase > phaseMassDensityTemp( NP );
  stackArray1d< EvalVarArgs, maxNumPhase > phaseViscosityTemp( NP );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // molarDensity or massDensity (useMass)
    m_phaseDensityFuns[ip]->evaluation( P, T, phaseCompFractionTemp[ip], phaseDensityTemp[ip], m_useMass );
    m_phaseViscosityFuns[ip]->evaluation( P, T, phaseCompFractionTemp[ip], phaseViscosityTemp[ip] );
  }

  if( m_useMass )
  {
    stackArray1d< EvalVarArgs, maxNumPhase > phaseMW( NP );
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      // copy phaseDens into phaseMassDens
      phaseMassDensityTemp[ip] = phaseDensityTemp[ip];

      // compute the molecular weight to get the mass phase (component) fractions
      EvalVarArgs molarDens;
      m_phaseDensityFuns[ip]->evaluation( P, T, phaseCompFractionTemp[ip], molarDens, 0 );
      phaseMW[ip] =  phaseDensityTemp[ip] /  molarDens;
    }

    EvalVarArgs totalMass = 0.0;
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFractionTemp[ip] *= phaseMW[ip];
      totalMass += phaseFractionTemp[ip];
    }

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFractionTemp[ip] /= totalMass;
    }

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {

        real64 compMW = m_componentMolarWeight[ic];

        phaseCompFractionTemp[ip][ic] = phaseCompFractionTemp[ip][ic] * compMW /  phaseMW[ip];

      }
    }
  }
  else
  {
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      // recompute the mass density
      EvalVarArgs massDens;
      m_phaseDensityFuns[ip]->evaluation( P, T, phaseCompFractionTemp[ip], massDens, 1 );

      // copy phaseDens into phaseMassDens
      phaseMassDensityTemp[ip] = massDens;
    }
  }

  EvalVarArgs totalDensityTemp = 0.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    totalDensityTemp += phaseFractionTemp[ip] / phaseDensityTemp[ip];
  }
  totalDensityTemp  = 1.0 / totalDensityTemp;

  //transfer data
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    phaseFrac.value[ip] = phaseFractionTemp[ip].m_var;
    phaseFrac.dPres[ip] = phaseFractionTemp[ip].m_der[0];
    phaseFrac.dTemp[ip] = 0.0;

    phaseDens.value[ip] = phaseDensityTemp[ip].m_var;
    phaseDens.dPres[ip] = phaseDensityTemp[ip].m_der[0];
    phaseDens.dTemp[ip] = 0.0;

    phaseMassDens.value[ip] = phaseMassDensityTemp[ip].m_var;
    phaseMassDens.dPres[ip] = phaseMassDensityTemp[ip].m_der[0];
    phaseMassDens.dTemp[ip] = 0.0;

    phaseVisc.value[ip] = phaseViscosityTemp[ip].m_var;
    phaseVisc.dPres[ip] = phaseViscosityTemp[ip].m_der[0];
    phaseVisc.dTemp[ip] = 0.0;

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      phaseFrac.dComp[ip][ic]     = phaseFractionTemp[ip].m_der[ic+1];
      phaseDens.dComp[ip][ic]     = phaseDensityTemp[ip].m_der[ic+1];
      phaseMassDens.dComp[ip][ic] = phaseMassDensityTemp[ip].m_der[ic+1];
      phaseVisc.dComp[ip][ic]     = phaseViscosityTemp[ip].m_der[ic+1];

      phaseCompFrac.value[ip][ic] = phaseCompFractionTemp[ip][ic].m_var;
      phaseCompFrac.dPres[ip][ic] = phaseCompFractionTemp[ip][ic].m_der[0];
      phaseCompFrac.dTemp[ip][ic] = 0.0;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        phaseCompFrac.dComp[ip][ic][jc] = phaseCompFractionTemp[ip][ic].m_der[jc+1];
      }
    }
  }

  totalDens.value = totalDensityTemp.m_var;
  totalDens.dPres = totalDensityTemp.m_der[0];
  totalDens.dTemp = 0.0;

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    totalDens.dComp[ic] = totalDensityTemp.m_der[ic+1];
  }
}

} //namespace constitutive

} //namespace geosx
