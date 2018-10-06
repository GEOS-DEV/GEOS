/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
  * @file CompositionalMultiphaseFluid.cpp
  */

#include "MultiFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


MultiFluidBase::MultiFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent ),
    m_useMassFractions( false )
{
  RegisterViewWrapper( viewKeysMultiFluidBase.phases.Key(), &m_phases, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.componentNames.Key(), &m_componentNames, false );

  RegisterViewWrapper( viewKeysMultiFluidBase.phaseFraction.Key(), &m_phaseFraction, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseFraction_dPressure.Key(), &m_dPhaseFraction_dPressure, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseFraction_dTemperature.Key(), &m_dPhaseFraction_dTemperature, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseFraction_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeysMultiFluidBase.phaseVolumeFraction.Key(), &m_phaseVolumeFraction, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseVolumeFraction_dPressure.Key(), &m_dPhaseVolumeFraction_dPressure, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseVolumeFraction_dTemperature.Key(), &m_dPhaseVolumeFraction_dTemperature, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseVolumeFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseVolumeFraction_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeysMultiFluidBase.phaseDensity.Key(), &m_phaseDensity, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseDensity_dPressure.Key(), &m_dPhaseDensity_dPressure, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseDensity_dTemperature.Key(), &m_dPhaseDensity_dTemperature, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseDensity_dGlobalCompMoleFraction.Key(), &m_dPhaseDensity_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeysMultiFluidBase.phaseCompFraction.Key(), &m_phaseCompFraction, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseCompFraction_dPressure.Key(), &m_dPhaseCompFraction_dPressure, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseCompFraction_dTemperature.Key(), &m_dPhaseCompFraction_dTemperature, false );
  RegisterViewWrapper( viewKeysMultiFluidBase.dPhaseCompFraction_dGlobalCompFraction.Key(), &m_dPhaseCompFraction_dGlobalCompFraction, false );
}

void MultiFluidBase::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                               localIndex const numPts )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numPts );

  localIndex const size = parent->size();
  localIndex const numPhase = numFluidPhases();
  localIndex const numComp = numFluidComponents();

  this->resize( size );

  m_phaseFraction.resize( size, numPts, numPhase );
  m_dPhaseFraction_dPressure.resize( size, numPts, numPhase );
  m_dPhaseFraction_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseFraction_dGlobalCompFraction.resize( size, numPts, numPhase, numComp );

  m_phaseVolumeFraction.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dPressure.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dGlobalCompFraction.resize( size, numPts, numPhase, numComp );

  m_phaseDensity.resize( size, numPts, numPhase );
  m_dPhaseDensity_dPressure.resize( size, numPts, numPhase );
  m_dPhaseDensity_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseDensity_dGlobalCompFraction.resize( size, numPts, numPhase, numComp );

  m_phaseCompFraction.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompFraction_dPressure.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompFraction_dTemperature.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompFraction_dGlobalCompFraction.resize( size, numPts, numPhase, numComp, numComp );
}

MultiFluidBase::~MultiFluidBase()
{

}

void MultiFluidBase::FillDocumentationNode()
{
  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->GetCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Multi-component multiphase fluid model" );

  docNode->AllocateChildNode( viewKeysMultiFluidBase.phases.Key(),
                              viewKeysMultiFluidBase.phases.Key(),
                              -1,
                              "string_array",
                              "string_array",
                              "List of fluid phases",
                              "List of fluid phases",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeysMultiFluidBase.componentNames.Key(),
                              viewKeysMultiFluidBase.componentNames.Key(),
                              -1,
                              "string_array",
                              "string_array",
                              "List of component names",
                              "List of component names",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );
}

void MultiFluidBase::ReadXML_PostProcess()
{
    GEOS_ERROR_IF(numFluidComponents() > MAX_NUM_COMPONENTS,
      "MultiFluidBase: Number of fluid components exceeds the maximum of " << MAX_NUM_COMPONENTS);
}

localIndex MultiFluidBase::numFluidComponents() const
{
  return integer_conversion<localIndex>(m_componentNames.size());
}

localIndex MultiFluidBase::numFluidPhases() const
{
  return integer_conversion<localIndex>(m_phases.size());
}

bool MultiFluidBase::getMassFractionFlag() const
{
  return m_useMassFractions;
}

void MultiFluidBase::setMassFractionFlag(bool flag)
{
  m_useMassFractions = flag;
}

} //namespace constitutive

} //namespace geosx
