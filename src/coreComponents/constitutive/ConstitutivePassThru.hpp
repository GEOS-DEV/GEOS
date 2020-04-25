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

#include "Dummy.hpp"
#include "solid/LinearElasticIsotropic.hpp"
#include "solid/LinearElasticAnisotropic.hpp"
#include "solid/LinearElasticTransverseIsotropic.hpp"

namespace geosx
{
namespace constitutive
{

template< typename BASETYPE >
struct ConstitutivePassThru;


template<>
struct ConstitutivePassThru<SolidBase>
{
  template< typename LAMBDA >
  static
  inline
  bool Execute( ConstitutiveBase * const constitutiveRelation,
                LAMBDA && lambda )
  {
    bool rval = true;
    if( dynamic_cast< LinearElasticIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticIsotropic * >( constitutiveRelation) );
    }
    else if( dynamic_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation) );
    }
    else if( dynamic_cast< LinearElasticAnisotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticAnisotropic * >( constitutiveRelation) );
    }
    else
    {
      rval = false;
    }
    return rval;
  }
};


template<>
struct ConstitutivePassThru<Dummy>
{
  template< typename LAMBDA >
  static
  inline
  bool Execute( ConstitutiveBase * const constitutiveRelation,
                LAMBDA && lambda )
  {
    bool rval = true;
    if( dynamic_cast< Dummy * >( constitutiveRelation ) )
    {
      lambda( static_cast< Dummy * >( constitutiveRelation ) );
    }
    else
    {
      lambda( static_cast< ConstitutiveBase * >( constitutiveRelation ) );
      rval = false;
    }
    return rval;
  }
};



template< typename LAMBDA >
static
inline
bool ConstitutiveBasePassThru( ConstitutiveBase * const constitutiveRelation,
                               LAMBDA && lambda )
{
  bool rval = true;
  if( dynamic_cast< SolidBase * >( constitutiveRelation ) )
  {
    ConstitutivePassThru<SolidBase>::Execute( static_cast< SolidBase *>(constitutiveRelation),
                                              std::forward<LAMBDA&&>(lambda) );
  }
  else if( dynamic_cast< Dummy * >( constitutiveRelation ) )
  {
    ConstitutivePassThru<Dummy>::Execute( static_cast< Dummy *>(constitutiveRelation),
                                          std::forward<LAMBDA&&>(lambda) );
  }
  return rval;
}

}
}

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_ */
