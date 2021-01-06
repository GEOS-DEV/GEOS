/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ReactiveFluidBase.cpp
 */

#include "ReactiveFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ReactiveFluidBase::ReactiveFluidBase( std::string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{

  registerWrapper( viewKeyStruct::basisSpeciesNamesString, &m_basisSpeciesNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Basis speciese name list" );


  registerWrapper( viewKeyStruct::inputDependentSpeciesNamesString, &m_inputDependentSpeciesNames )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Input dependent speciese name list" );


  registerWrapper( viewKeyStruct::logActH2OString, &m_logActH2O )->
    setApplyDefaultValue( -1.65676E-02 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "logActH2O" );

  registerWrapper( viewKeyStruct::logFO2gString, &m_logFO2g )->
    setApplyDefaultValue( -0.7 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "logFO2g" );

  registerWrapper( viewKeyStruct::dependentConcString, &m_dependentConc );

  registerWrapper( viewKeyStruct::dDependentConc_dConcString, &m_dDependentConc_dConc );

  registerWrapper( viewKeyStruct::kineticReactionRateString, &m_kineticReactionRate );

  registerWrapper( viewKeyStruct::dKineticReactionRate_dConcString, &m_dKineticReactionRate_dConc );

  registerWrapper( viewKeyStruct::kineticSpeciesReactionRateString, &m_kineticSpeciesReactionRate );

  registerWrapper( viewKeyStruct::dKineticSpeciesReactionRate_dConcString, &m_dKineticSpeciesReactionRate_dConc );


  registerWrapper( viewKeyStruct::totalConcString, &m_totalConc );

  registerWrapper( viewKeyStruct::dTotalConc_dConcString, &m_dTotalConc_dConc );

  registerWrapper( viewKeyStruct::concentrationActString, &m_concentrationAct );

}

ReactiveFluidBase::~ReactiveFluidBase() = default;

void ReactiveFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

}

void ReactiveFluidBase::allocateConstitutiveData( Group * const parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


} //namespace constitutive

} //namespace geosx
