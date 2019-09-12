/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file solidSelector.hpp
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_

#include "LinearElasticIsotropic.hpp"

namespace geosx
{
namespace constitutive
{


template< typename LAMBDA >
bool constitutiveUpdatePassThru( constitutive::ConstitutiveBase * const constitutiveRelation,
                                 LAMBDA && lambda )
{
  bool rval = true;
  if( dynamic_cast<LinearElasticIsotropic * >( constitutiveRelation ) )
  {
    lambda( static_cast<LinearElasticIsotropic & >( *constitutiveRelation) );
  }
  else
  {
    rval = false;
  }

  return rval;
}

}
}

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_ */
