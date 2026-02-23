/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BartonBandis.cpp
 */


#include "BartonBandisStressPathDriven.hpp"


namespace geos
{

namespace constitutive
{

using namespace dataRepository;

BartonBandisStressPathDriven::BartonBandisStressPathDriven( string const & name, Group * const parent ):
  HydraulicApertureBase( name, parent )
{
  
  registerWrapper( viewKeyStruct::biotString(), &m_biot ).
    setApplyDefaultValue( 1.0 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Biot coefficient." );
  registerWrapper( viewKeyStruct::poissonString(), &m_poisson ).
    setApplyDefaultValue( 0.3 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Poisson ratio." );
  registerWrapper( viewKeyStruct::normalStiffnessString(), &m_normalStiffness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal stiffness: Kni." );
  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 1.0e5 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure: p_0." );
  registerWrapper( viewKeyStruct::referenceTotalStressString(), &m_referenceTotalStress ).
    setApplyDefaultValue( { 85.0e6, 85.0e6, 105.0e6 } ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Total stress at reference state: sigmaT_0." );
    
}


BartonBandisStressPathDrivenUpdates BartonBandisStressPathDriven::createKernelWrapper() const
{
  return KernelWrapper( m_aperture0, 
                        m_biot, 
                        m_poisson, 
                        m_normalStiffness, 
                        m_referencePressure, 
                        m_referenceTotalStress);
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BartonBandisStressPathDriven, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
