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

/**
 * @file BiotPorosity.cpp
 */

#include "BiotPorosity.hpp"
#include "PorosityFields.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


BiotPorosity::BiotPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultGrainBulkModulusString(), &m_defaultGrainBulkModulus ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Grain bulk modulus" );

  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default thermal expansion coefficient" );

  registerWrapper( viewKeyStruct::useUniaxialFixedStressString(), &m_useUniaxialFixedStress ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag enabling uniaxial approximation in fixed stress update" );

  registerWrapper( viewKeyStruct::defaultBiotCoefficientString(), &m_defaultBiotCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Default Biot coefficient. If not specified it will be homogeneous for the material." );

  registerField( fields::porosity::biotCoefficient{}, &m_biotCoefficient ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Biot coefficient" );

  registerField( fields::porosity::grainBulkModulus{}, &m_grainBulkModulus ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Grain Bulk modulus." );

  registerField( fields::porosity::thermalExpansionCoefficient{}, &m_thermalExpansionCoefficient );

  registerField( fields::porosity::meanTotalStressIncrement_k{}, &m_meanTotalStressIncrement_k );

  registerField( fields::porosity::averageMeanTotalStressIncrement_k{}, &m_averageMeanTotalStressIncrement_k );

  registerWrapper( viewKeyStruct::solidBulkModulusString(), &m_bulkModulus ).
    setApplyDefaultValue( 1e-6 ).
    setDescription( "Solid bulk modulus" );

  registerWrapper( viewKeyStruct::solidShearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( 1e-6 ).
    setDescription( "Solid shear modulus" );
}

void BiotPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_meanTotalStressIncrement_k.resize( 0, numConstitutivePointsPerParentIndex );
}

void BiotPorosity::postProcessInput()
{
  PorosityBase::postProcessInput();

  GEOS_ERROR_IF( (m_defaultBiotCoefficient > -1.0) & (m_defaultGrainBulkModulus > -1.0), 
                 "Both the defaultBiotCoefficient and the defaultGrainBulkModulus were specified. Only one of them can be specified" ); 

  GEOS_ERROR_IF( (m_defaultBiotCoefficient < 0.0) & (m_defaultGrainBulkModulus < 0.0), 
                 "It seems that neither the defaultBiotCoefficient and the defaultGrainBulkModulus were not specified. One of them mus be specified" );                  

  getWrapper< array1d< real64 > >( fields::porosity::thermalExpansionCoefficient::key() ).
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );

  // set results as array default values
  getWrapper< array1d< real64 > >( fields::porosity::biotCoefficient::key() ).
    setApplyDefaultValue( m_defaultBiotCoefficient );

  // set results as array default values
  getWrapper< array1d< real64 > >( fields::porosity::grainBulkModulus::key() ).
    setApplyDefaultValue( m_defaultGrainBulkModulus );
}

void BiotPorosity::initializeState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView1d< real64 const > referencePorosity = m_referencePorosity;
  arrayView2d< real64 >             newPorosity = m_newPorosity;
  arrayView2d< real64 >             porosity_n  = m_porosity_n;
  arrayView2d< real64 >         initialPorosity = m_initialPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      newPorosity[k][q]     = referencePorosity[k];
      porosity_n[k][q]      = referencePorosity[k];
      initialPorosity[k][q] = referencePorosity[k];
    }
  } );
}

void BiotPorosity::initializeBiotCoefficient( arrayView1d< real64 const> const bulkModulus ) const
{
  localIndex const numE = numElem();

  arrayView1d< real64 >  biotCoefficient = m_biotCoefficient.toView();
  arrayView1d< real64 >  grainBulkModulus = m_grainBulkModulus.toView();
  /// Note: this only works for linearelasticity but since this assumption is made in many other places
  /// we can do this for now. I will have to be removed / modified so that the biotcoefficient is always computed
  /// based on the solid model.
  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    if ( grainBulkModulus[k] > -1.0 )
    {
      biotCoefficient[k] = 1 - bulkModulus[k] / grainBulkModulus[k];
    }
    else if ( biotCoefficient[k] > -1.0)
    {
      grainBulkModulus[k] = bulkModulus[k] / (1 - grainBulkModulus[k]);
    }
  } );
}

void BiotPorosity::saveConvergedState() const
{
  PorosityBase::saveConvergedState();
  m_meanTotalStressIncrement_k.zero();
  m_averageMeanTotalStressIncrement_k.zero();
}

void BiotPorosity::ignoreConvergedState() const
{
  PorosityBase::ignoreConvergedState();
  m_meanTotalStressIncrement_k.zero();
  m_averageMeanTotalStressIncrement_k.zero();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, BiotPorosity, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
