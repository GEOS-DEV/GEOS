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
 * @file CohesiveZoneBase.cpp
 */

#include "CohesiveZoneBase.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

CohesiveZoneBase::CohesiveZoneBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newNormalStress( 0, 3 ),
  m_oldNormalStress( 0, 3 ),
  m_newShearStress( 0, 3 ),
  m_oldShearStress( 0, 3 )
{
  string const labels[3] = { "X", "Y", "Z" };

  registerWrapper( viewKeyStruct::newNormalStressString(), &m_newNormalStress ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setApplyDefaultValue( 0.0 ). // default to zero initial stress
    setDescription( "Current normal stress" ).
    setDimLabels( 1, labels );

  registerWrapper( viewKeyStruct::newShearStressString(), &m_newShearStress ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setApplyDefaultValue( 0.0 ). // default to zero initial stress
    setDescription( "Current shear stress" ).
    setDimLabels( 1, labels );

  registerWrapper( viewKeyStruct::oldNormalStressString(), &m_oldNormalStress ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setApplyDefaultValue( 0.0 ). // default to zero initial stress
    setDescription( "Previous normal stress" ).
    setDimLabels( 1, labels );

  registerWrapper( viewKeyStruct::oldShearStressString(), &m_oldShearStress ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setApplyDefaultValue( 0.0 ). // default to zero initial stress
    setDescription( "Previous shear stress" ).
    setDimLabels( 1, labels );
}


void CohesiveZoneBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();
}


void CohesiveZoneBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_newNormalStress.resize( 0, 3 );
  m_newShearStress.resize( 0, 3 );
  m_oldNormalStress.resize( 0, 3 );
  m_oldShearStress.resize( 0, 3 );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CohesiveZoneBase, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
