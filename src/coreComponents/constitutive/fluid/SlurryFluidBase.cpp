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
  * @file SlurryFluidBase.cpp
  */

#include "SlurryFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SlurryFluidBase::SlurryFluidBase( std::string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ), m_isNewtonianFluid(1)
{

  registerWrapper( viewKeyStruct::componentNamesString, &m_componentNames, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of fluid component names");

  registerWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Default value for density.");  

  registerWrapper( viewKeyStruct::defaultCompressibilityString, &m_defaultCompressibility, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Default value for compressibility.");  

  registerWrapper( viewKeyStruct::defaultViscosityString, &m_defaultViscosity, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Default value for viscosity.");  
  
  
  registerWrapper( viewKeyStruct::flowBehaviorIndexString, &m_nIndices, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flow behavior index");

  registerWrapper( viewKeyStruct::flowConsistencyIndexString, &m_Ks, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flow consistency index");      

  registerWrapper( viewKeyStruct::densityString, &m_density, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dDens_dPresString, &m_dDens_dPres, false );
  registerWrapper( viewKeyStruct::dDens_dProppantConcString, &m_dDens_dProppantConc, false );
  registerWrapper( viewKeyStruct::dDens_dCompConcString, &m_dDens_dCompConc, false );  

  registerWrapper( viewKeyStruct::fluidDensityString, &m_fluidDensity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dFluidDens_dPresString, &m_dFluidDens_dPres, false );
  registerWrapper( viewKeyStruct::dFluidDens_dCompConcString, &m_dFluidDens_dCompConc, false );    

    registerWrapper( viewKeyStruct::fluidViscosityString, &m_fluidViscosity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dFluidVisc_dPresString, &m_dFluidVisc_dPres, false );
  registerWrapper( viewKeyStruct::dFluidVisc_dCompConcString, &m_dFluidVisc_dCompConc, false );    

  registerWrapper( viewKeyStruct::componentDensityString, &m_componentDensity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dCompDens_dPresString, &m_dCompDens_dPres, false );
  registerWrapper( viewKeyStruct::dCompDens_dCompConcString, &m_dCompDens_dCompConc, false );    
  
  
  registerWrapper( viewKeyStruct::viscosityString, &m_viscosity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dVisc_dPresString, &m_dVisc_dPres, false );
  registerWrapper( viewKeyStruct::dVisc_dProppantConcString, &m_dVisc_dProppantConc, false );
  registerWrapper( viewKeyStruct::dVisc_dCompConcString, &m_dVisc_dCompConc, false );  
  
}

SlurryFluidBase::~SlurryFluidBase() = default;

void SlurryFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  localIndex const NC = numFluidComponents();

  GEOSX_ERROR_IF( m_defaultDensity.size() != NC,
                  "The number of flow behavior indices is not the same as the component number" );

  GEOSX_ERROR_IF( m_defaultCompressibility.size() != NC,
                  "The number of flow behavior indices is not the same as the component number" );

  GEOSX_ERROR_IF( m_defaultViscosity.size() != NC,
                  "The number of flow behavior indices is not the same as the component number" );

  /*
  GEOSX_ERROR_IF( m_nIndices.size() != NC, "The number of flow behavior indices is not the same as the component number" );

  GEOSX_ERROR_IF( m_Ks.size() != NC, "The number of flow consistency indices is not the same as the component number" );  
  */

  
}

localIndex SlurryFluidBase::numFluidComponents() const
{
  return integer_conversion<localIndex>(m_componentNames.size());
}

string const & SlurryFluidBase::componentName(localIndex ic) const
{
  GEOSX_ERROR_IF( ic >= numFluidComponents(), "Index " << ic << " exceeds number of fluid components" );
  return m_componentNames[ic];
}

  
void SlurryFluidBase::AllocateConstitutiveData( Group * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  localIndex const NC = numFluidComponents();
  
  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDens_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDens_dProppantConc.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDens_dCompConc.resize( parent->size(), numConstitutivePointsPerParentIndex, NC );  

  m_componentDensity.resize( parent->size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dCompConc.resize( parent->size(), numConstitutivePointsPerParentIndex, NC, NC );  

  
  m_fluidDensity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dCompConc.resize( parent->size(), numConstitutivePointsPerParentIndex , NC );  

  m_fluidViscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dCompConc.resize( parent->size(), numConstitutivePointsPerParentIndex , NC );  
  
  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dProppantConc.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dCompConc.resize( parent->size(), numConstitutivePointsPerParentIndex, NC );  
  
}


void
SlurryFluidBase::DeliverClone( string const & name,
                               Group * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOSX_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  SlurryFluidBase * const newConstitutiveRelation = dynamic_cast<SlurryFluidBase *>(clone.get());

  newConstitutiveRelation->m_componentNames  = this->m_componentNames;

  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_defaultCompressibility = m_defaultCompressibility;
  newConstitutiveRelation->m_defaultViscosity = m_defaultViscosity;    

  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_dDens_dPres = m_dDens_dPres;
  newConstitutiveRelation->m_dDens_dProppantConc = m_dDens_dProppantConc;
  newConstitutiveRelation->m_dDens_dCompConc = m_dDens_dCompConc;  

  newConstitutiveRelation->m_componentDensity = m_componentDensity;
  newConstitutiveRelation->m_dCompDens_dPres = m_dCompDens_dPres;
  newConstitutiveRelation->m_dCompDens_dCompConc = m_dCompDens_dCompConc;  
  
  newConstitutiveRelation->m_fluidDensity = m_fluidDensity;  
  newConstitutiveRelation->m_dFluidDens_dPres = m_dFluidDens_dPres;
  newConstitutiveRelation->m_dFluidDens_dCompConc = m_dFluidDens_dCompConc;

  newConstitutiveRelation->m_fluidViscosity = m_fluidViscosity;  
  newConstitutiveRelation->m_dFluidVisc_dPres = m_dFluidVisc_dPres;
  newConstitutiveRelation->m_dFluidVisc_dCompConc = m_dFluidVisc_dCompConc;
  
  newConstitutiveRelation->m_viscosity = m_viscosity;
  newConstitutiveRelation->m_dVisc_dPres = m_dVisc_dPres;
  newConstitutiveRelation->m_dVisc_dProppantConc = m_dVisc_dProppantConc;
  newConstitutiveRelation->m_dVisc_dCompConc = m_dVisc_dCompConc;    

  newConstitutiveRelation->m_nIndices   = this->m_nIndices;
  newConstitutiveRelation->m_Ks   = this->m_Ks;

  newConstitutiveRelation->m_isNewtonianFluid   = this->m_isNewtonianFluid;    
  
}

  
} //namespace constitutive

} //namespace geosx
