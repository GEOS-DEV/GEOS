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

#include "OedometricStressPath.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
// GEOS_PHYSICSSOLVERS_OEDOMETRICSTRESSPATH_HPP_
namespace geos
{

using namespace dataRepository;

OedometricStressPath::OedometricStressPath( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );

  registerWrapper( viewKeysStruct::biotString(), &m_biot ).
    setApplyDefaultValue( 1.0 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Biot coefficient." );
  registerWrapper( viewKeysStruct::poissonString(), &m_poisson ).
    setApplyDefaultValue( 0.3 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( " Poisson ratio." );
  registerWrapper( viewKeysStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 1.0e5 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( " Reference pressure: p_0." );
  registerWrapper( viewKeysStruct::referenceTotalStressString(), &m_referenceTotalStress ).
    setApplyDefaultValue( { 85.0e6, 85.0e6, 105.0e6 } ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( " Total stress at reference state: sigmaT_0." );

  m_referenceEffectiveStress[0] = m_referenceTotalStress[0] - m_biot*m_referencePressure; 
  m_referenceEffectiveStress[1] = m_referenceTotalStress[1] - m_biot*m_referencePressure; 
  m_referenceEffectiveStress[2] = m_referenceTotalStress[2] - m_biot*m_referencePressure; 


  addLogLevel< logInfo::NonlinearSolver >();
}

void OedometricStressPath::postInputInitialization()
{
  if( getLogLevel() > 0 )
  {
    print(); 
  }
}

void OedometricStressPath::print() const
{
  TableData tableData;
  tableData.addRow( "Log level", getLogLevel());
  tableData.addRow( "Biot", m_biot );
  tableData.addRow( "Poisson", m_poisson );
  tableData.addRow( "Reference Pressure", m_referencePressure );
  tableData.addRow( "Reference total stress", m_referenceTotalStress[0] );
  TableLayout const tableLayout = TableLayout( GEOS_FMT( "{}: oedometric", getParent().getName() ),
                                               { "Parameter", "Value" } );
  TableTextFormatter const tableFormatter( tableLayout );
  GEOS_LOG_RANK_0( tableFormatter.toString( tableData ));
}


real64 OedometricStressPath::computeFractureStress( real64 const pressure,
                                                    R1Tensor const & normal ) const
{  
  auto dot_product = [](real64 const (&u)[3], R1Tensor const &v ) -> real64
  {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  };

  auto matmul = [](real64 const (&u)[3], R1Tensor const &v, real64 (&r)[3]) -> void
  {
    r[0] = u[0]*v[0];
    r[1] = u[1]*v[1];
    r[2] = u[2]*v[2];
  }; 

  real64 const deltaSigmaZ = m_biot * (pressure - m_referencePressure);
  real64 const poisson_deltaSigma = deltaSigmaZ * m_poisson/(1.0 - m_poisson);
  // matrix diagonal
  real64 effectiveStress[3] = { m_referenceEffectiveStress[0] - poisson_deltaSigma,
                                m_referenceEffectiveStress[1] - poisson_deltaSigma,
                                m_referenceEffectiveStress[2] - deltaSigmaZ };

  real64 effectiveStressOnFracture[3] = {0.0, 0.0, 0.0}; // sigma_c
  matmul(effectiveStress, normal, effectiveStressOnFracture);
  real64 normalComponentOfStressOnFracture = dot_product(effectiveStressOnFracture, normal); // sigmaN_N

  return normalComponentOfStressOnFracture;
}



REGISTER_CATALOG_ENTRY( Group, OedometricStressPath, string const &, Group * const )

} /* namespace geos */
