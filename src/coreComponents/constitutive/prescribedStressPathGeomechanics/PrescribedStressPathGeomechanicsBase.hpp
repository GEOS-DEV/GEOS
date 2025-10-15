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

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PRESCRIBEDSTRESSPATHGEOMECHANICSBASE_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PRESCRIBEDSTRESSPATHGEOMECHANICSBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{
/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class PrescribedStressPathGeomechanicsBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  PrescribedStressPathGeomechanicsBase( string const & name,
                                        Group * const parent );


  virtual real64 computeFractureStress( real64 const pressure, R1Tensor const & normal) const;
  
};

}/* namespace constitutive */

} /* namespace geos */



#endif /* GEOS_CONSTITUTIVE_PERMEABILITY_PRESCRIBEDSTRESSPATHGEOMECHANICSBASE_HPP_ */


