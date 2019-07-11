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

ReactiveFluidBase::ReactiveFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent )
{

  RegisterViewWrapper( viewKeyStruct::basisSpeciesNamesString, &m_basisSpeciesNames, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Basis speciese name list");

  RegisterViewWrapper( viewKeyStruct::logActH2OString, &m_logActH2O, false )->
    setApplyDefaultValue(-1.65676E-02)->
    setInputFlag(InputFlags::OPTIONAL)->    
    setDescription("logActH2O");

  RegisterViewWrapper( viewKeyStruct::logFO2gString, &m_logFO2g, false )->
    setApplyDefaultValue(-0.7)->
    setInputFlag(InputFlags::OPTIONAL)->    
    setDescription("logFO2g");    

  RegisterViewWrapper( viewKeyStruct::dependentConcString, &m_dependentConc, false );

  RegisterViewWrapper( viewKeyStruct::dDependentConc_dConcString, &m_dDependentConc_dConc, false );  

  
}

ReactiveFluidBase::~ReactiveFluidBase() = default;

void ReactiveFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  bool HplusNotFound = 1;
  bool H2OFound = 0;  

  localIndex const NBasis = numBasisSpecies();

  m_isHplus.resize(NBasis);
  m_isHplus = 0;
  
  for(localIndex id = 0; id < NBasis; ++id)
    {

      if(m_basisSpeciesNames[id] == "H+")
	{
	  HplusNotFound = 0;
	  m_isHplus[id] = 1;
	}
      
      if(m_basisSpeciesNames[id] == "H2O")
	H2OFound = 1;      

    }

  GEOS_ERROR_IF( HplusNotFound, "ReactiveFluidBase: H+ is not specified in basisSpeciesNames" );

  GEOS_ERROR_IF( H2OFound, "ReactiveFluidBase: H2O cannot be specified in basisSpeciesNames" );  
  
}

void ReactiveFluidBase::AllocateConstitutiveData( ManagedGroup * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

}
  
void
ReactiveFluidBase::DeliverClone( string const & name,
                               ManagedGroup * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOS_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  ReactiveFluidBase * const newConstitutiveRelation = dynamic_cast<ReactiveFluidBase *>(clone.get());

  newConstitutiveRelation->m_basisSpeciesNames = m_basisSpeciesNames;
  newConstitutiveRelation->m_dependentSpeciesNames = m_dependentSpeciesNames;
  newConstitutiveRelation->m_stochMatrix = m_stochMatrix;
  newConstitutiveRelation->m_dependentConc = m_dependentConc;
  newConstitutiveRelation->m_dDependentConc_dConc = m_dDependentConc_dConc;
  newConstitutiveRelation->m_logFO2g = m_logFO2g;
  newConstitutiveRelation->m_logActH2O = m_logActH2O;
}

} //namespace constitutive

} //namespace geosx
