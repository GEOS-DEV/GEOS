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

#include "constitutive/prescribedStressPathGeomechanics/PrescribedStressPathGeomechanicsBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

PrescribedStressPathGeomechanicsBase::PrescribedStressPathGeomechanicsBase( string const & name,
                                                                            Group * const parent ):
  ConstitutiveBase( name, parent )
{
}


real64 PrescribedStressPathGeomechanicsBase::computeFractureStress( real64 const pressure, 
                                                                    R1Tensor const & normal ) const
{
  GEOS_ERROR( "PrescribedStressPathGeomechanicsBase::computeFractureStress called!. Should be overridden." );
  return 0.0;
}

}/* namespace constitutive */

} /* namespace geos */
