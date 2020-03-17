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

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_

#include "LinearElasticIsotropic.hpp"
#include "LinearViscoElasticIsotropic.hpp"
#include "LinearElasticAnisotropic.hpp"
#include "BilinearElasticIsotropic.hpp"
#include "HardeningElasticIsotropic.hpp"
#include "CycLiqCPSP.hpp"
#include "NonlinearElasticDuncanChangEB.hpp"

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
//  else if( dynamic_cast<LinearViscoElasticIsotropic * >( constitutiveRelation ) )
//  {
//    lambda( static_cast<LinearViscoElasticIsotropic & >( *constitutiveRelation) );
//  }
  else if( dynamic_cast<BilinearElasticIsotropic * >( constitutiveRelation ) )
  {
    lambda( static_cast<BilinearElasticIsotropic & >( *constitutiveRelation) );
  }
  else if( dynamic_cast<HardeningElasticIsotropic * >( constitutiveRelation ) )
  {
    lambda( static_cast<HardeningElasticIsotropic & >( *constitutiveRelation) );
  }
  else if( dynamic_cast<CycLiqCPSP * >( constitutiveRelation ) )
  {
    lambda( static_cast<CycLiqCPSP & >( *constitutiveRelation) );
  }
  else if( dynamic_cast<NonlinearElasticDuncanChangEB * >( constitutiveRelation ) )
  {
    lambda( static_cast<NonlinearElasticDuncanChangEB & >( *constitutiveRelation) );
  }
  else if( dynamic_cast<LinearElasticAnisotropic * >( constitutiveRelation ) )
  {
#if !defined(__CUDA_ARCH__)
    lambda( static_cast<LinearElasticAnisotropic & >( *constitutiveRelation) );
#else
    GEOSX_ERROR( "Cannot call kernel using constitutiveUpdatePassThru. "
                 "Too many parameters in LinearElasticAnisotropic::KernelWrapper");
#endif
  }
  else
  {
    rval = false;
  }

  return rval;
}

}
}

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_ */
