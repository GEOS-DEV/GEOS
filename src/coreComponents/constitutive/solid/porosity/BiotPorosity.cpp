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

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


BiotPorosity::BiotPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::grainBulkModulusString(), &m_grainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Grain bulk modulus" );

  registerWrapper( viewKeyStruct::biotCoefficientString(), &m_biotCoefficient ).
    setDescription( "Biot coefficient." );
}

void BiotPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void BiotPorosity::postProcessInput()
{
  PorosityBase::postProcessInput();
  // TODO valdate input
}

void BiotPorosity::initializeState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView1d< real64 const > referencePorosity = m_referencePorosity;
  arrayView2d< real64 >             newPorosity = m_newPorosity;
  arrayView2d< real64 >             oldPorosity = m_oldPorosity;
  arrayView2d< real64 >         initialPorosity = m_initialPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      newPorosity[k][q]     = referencePorosity[k];
      oldPorosity[k][q]     = referencePorosity[k];
      initialPorosity[k][q] = referencePorosity[k];
    }
  } );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, BiotPorosity, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geosx */
