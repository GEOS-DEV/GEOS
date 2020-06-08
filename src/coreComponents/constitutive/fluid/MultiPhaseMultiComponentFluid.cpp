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
 * @file MultiPhaseMultiComponentFluid.cpp
 */
#include "MultiPhaseMultiComponentFluid.hpp"

#include "common/Path.hpp"
#include "managers/ProblemManager.hpp"
#include "constitutive/fluid/PVTFunctions/CO2SolubilityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/BrineCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/BrineViscosityFunction.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace PVTProps;
using namespace stringutilities;

namespace constitutive
{

MultiPhaseMultiComponentFluid::MultiPhaseMultiComponentFluid( std::string const & name, Group * const parent ):
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


void MultiPhaseMultiComponentFluid::DeliverClone( string const & name,
                                                  Group * const parent,
                                                  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< MultiPhaseMultiComponentFluid >( name, parent );
  }
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

}

void MultiPhaseMultiComponentFluid::PostProcessInput()
{
  MultiFluidBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF( m_phasePVTParaFiles.size() != NP, "The number of phasePVTParaFiles is not the same as the number of phases!" );

  CreatePVTModels();

}

void MultiPhaseMultiComponentFluid::InitializePostSubGroups( Group * const group )
{
  MultiFluidBase::InitializePostSubGroups( group );

  //  CreatePVTModels();

}


void MultiPhaseMultiComponentFluid::CreatePVTModels()
{
  string flashModelParaFile;

  ProblemManager const * const problemManager = this->GetGroupByPath< ProblemManager >( "/" );
  if( problemManager != nullptr )
  {
    string inputFileName = problemManager->getInputFileName();
    if( inputFileName.empty())
    {
      inputFileName = problemManager->getRestartFileName();
    }
    string inputFileDir;
    splitPath( inputFileName, inputFileDir, inputFileName );

    for( std::string & filename : m_phasePVTParaFiles )
    {
      if( !isAbsolutePath( filename ))
      {
        getAbsolutePath( inputFileDir + '/' + filename, filename );
      }
    }

    flashModelParaFile = inputFileDir + '/' + m_flashModelParaFile;

  }

  for( std::string & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );

    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ))
    {
      std::string str( buf );
      string_array strs = Tokenize( str, " " );

      if( streq( strs[0], "DensityFun" ))
      {
        m_phaseDensityFuns.push_back( PVTFunction::CatalogInterface::Factory( strs[1], strs, m_componentNames, m_componentMolarWeight ));
      }
      else if( streq( strs[0], "ViscosityFun" ))
      {
        m_phaseViscosityFuns.push_back( PVTFunction::CatalogInterface::Factory( strs[1], strs, m_componentNames, m_componentMolarWeight ));
      }
      else
        GEOSX_ERROR( "Error: Invalid PVT function: " << strs[0] << "." );
    }

    is.close();

  }

  {

    std::ifstream is( flashModelParaFile );

    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    while( is.getline( buf, buf_size ))
    {
      std::string str( buf );
      string_array strs = Tokenize( str, " " );

      if( streq( strs[0], "FlashModel" ))
      {

        m_flashModel = ( FlashModel::CatalogInterface::Factory( strs[1],
                                                                strs,
                                                                m_phaseNames,
                                                                m_componentNames,
                                                                m_componentMolarWeight ) );
      }
      else
        GEOSX_ERROR( "Error: Not flash model: " << strs[0] << "." );
    }

    is.close();

  }
}

void MultiPhaseMultiComponentFluid::PointUpdate( real64 const & pressure,
                                                 real64 const & temperature,
                                                 arraySlice1d< real64 const > const & composition,
                                                 localIndex const k,
                                                 localIndex const q )
{
  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  Compute( NC, NP, m_useMass,
           m_phaseNames,
           m_componentMolarWeight,
           pressure,
           temperature,
           composition,
           m_phaseFraction[k][q],
           m_dPhaseFraction_dPressure[k][q],
           m_dPhaseFraction_dTemperature[k][q],
           m_dPhaseFraction_dGlobalCompFraction[k][q],
           m_phaseDensity[k][q],
           m_dPhaseDensity_dPressure[k][q],
           m_dPhaseDensity_dTemperature[k][q],
           m_dPhaseDensity_dGlobalCompFraction[k][q],
           m_phaseViscosity[k][q],
           m_dPhaseViscosity_dPressure[k][q],
           m_dPhaseViscosity_dTemperature[k][q],
           m_dPhaseViscosity_dGlobalCompFraction[k][q],
           m_phaseCompFraction[k][q],
           m_dPhaseCompFraction_dPressure[k][q],
           m_dPhaseCompFraction_dTemperature[k][q],
           m_dPhaseCompFraction_dGlobalCompFraction[k][q],
           m_totalDensity[k][q],
           m_dTotalDensity_dPressure[k][q],
           m_dTotalDensity_dTemperature[k][q],
           m_dTotalDensity_dGlobalCompFraction[k][q],
           m_phaseDensityFuns,
           m_phaseViscosityFuns,
           m_flashModel );
}

void MultiPhaseMultiComponentFluid::BatchUpdate( arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( pressure ),
                                                 arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( temperature ),
                                                 arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( composition ) )
{}

void MultiPhaseMultiComponentFluid::Compute( localIndex const NC, localIndex const NP, bool const useMass,
                                             arrayView1d< string const > const & GEOSX_UNUSED_PARAM( phaseNames ),
                                             arrayView1d< real64 const > const & componentMolarWeight,
                                             real64 const & pressure,
                                             real64 const & temperature,
                                             arraySlice1d< real64 const > const & composition,
                                             arraySlice1d< real64 > const & phaseFraction,
                                             arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                                             arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                                             arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                                             arraySlice1d< real64 > const & phaseDensity,
                                             arraySlice1d< real64 > const & dPhaseDensity_dPressure,
                                             arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
                                             arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                                             arraySlice1d< real64 > const & phaseViscosity,
                                             arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
                                             arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
                                             arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                                             arraySlice2d< real64 > const & phaseCompFraction,
                                             arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                                             arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                                             arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                                             real64 & totalDensity, real64 & dTotalDensity_dPressure,
                                             real64 & dTotalDensity_dTemperature,
                                             arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction,
                                             array1d< std::shared_ptr< PVTProps::PVTFunction > > const & phaseDensityFuns,
                                             array1d< std::shared_ptr< PVTProps::PVTFunction > > const & phaseViscosityFuns,
                                             std::shared_ptr< PVTProps::FlashModel > const & flashModel )
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


  localIndex constexpr maxNumComp = MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumPhase = MAX_NUM_PHASES;

  stackArray1d< EvalVarArgs, maxNumComp > C( NC );

  if( useMass )
  {

    stackArray1d< EvalVarArgs, maxNumComp > X( NC );

    EvalVarArgs totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      X[ic].m_var = composition[ic];
      X[ic].m_der[ic+1] = 1.0;

      realT const mwInv = 1.0 / componentMolarWeight[ic];
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
  flashModel->Partition( P, T, C, phaseFractionTemp, phaseCompFractionTemp );

  stackArray1d< EvalVarArgs, maxNumPhase > phaseDensityTemp( NP );
  stackArray1d< EvalVarArgs, maxNumPhase > phaseViscosityTemp( NP );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // molarDensity or massDensity (useMass)
    phaseDensityFuns[ip]->Evaluation( P, T, phaseCompFractionTemp[ip], phaseDensityTemp[ip], useMass );

    phaseViscosityFuns[ip]->Evaluation( P, T, phaseCompFractionTemp[ip], phaseViscosityTemp[ip] );

  }

  if( useMass )
  {

    stackArray1d< EvalVarArgs, maxNumPhase > phaseMW( NP );

    for( localIndex ip = 0; ip < NP; ++ip )
    {

      EvalVarArgs molarPhaseDensity;

      phaseDensityFuns[ip]->Evaluation( P, T, phaseCompFractionTemp[ip], molarPhaseDensity, 0 );

      phaseMW[ip] =  phaseDensityTemp[ip] /  molarPhaseDensity;

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

        realT compMW = componentMolarWeight[ic];

        phaseCompFractionTemp[ip][ic] = phaseCompFractionTemp[ip][ic] * compMW /  phaseMW[ip];

      }
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

    phaseVisc.value[ip] = phaseViscosityTemp[ip].m_var;
    phaseVisc.dPres[ip] = phaseViscosityTemp[ip].m_der[0];
    phaseVisc.dTemp[ip] = 0.0;

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      phaseFrac.dComp[ip][ic] = phaseFractionTemp[ip].m_der[ic+1];
      phaseDens.dComp[ip][ic] = phaseDensityTemp[ip].m_der[ic+1];
      phaseVisc.dComp[ip][ic] = phaseViscosityTemp[ip].m_der[ic+1];

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

void MultiPhaseMultiComponentFluid::Compute( real64 const & GEOSX_UNUSED_PARAM( pressure ), real64 const & GEOSX_UNUSED_PARAM( temperature ),
                                             arraySlice1d< double const > const & GEOSX_UNUSED_PARAM( composition ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( phaseFraction ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseFraction_dPressure ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseFraction_dTemperature ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseFraction_dGlobalCompFraction ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( phaseDensity ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseDensity_dPressure ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseDensity_dTemperature ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseDensity_dGlobalCompFraction ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( phaseViscosity ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseViscosity_dPressure ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseViscosity_dTemperature ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseViscosity_dGlobalCompFraction ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( phaseCompFraction ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dPressure ),
                                             arraySlice2d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dTemperature ),
                                             arraySlice3d< real64 > const & GEOSX_UNUSED_PARAM( dPhaseCompFraction_dGlobalCompFraction ),
                                             real64 & GEOSX_UNUSED_PARAM( totalDensity ),
                                             real64 & GEOSX_UNUSED_PARAM( dTotalDensity_dPressure ),
                                             real64 & GEOSX_UNUSED_PARAM( dTotalDensity_dTemperature ),
                                             arraySlice1d< real64 > const & GEOSX_UNUSED_PARAM( dTotalDensity_dGlobalCompFraction ) ) const
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MultiPhaseMultiComponentFluid, std::string const &, Group * const )

} //namespace constitutive

} //namespace geosx
