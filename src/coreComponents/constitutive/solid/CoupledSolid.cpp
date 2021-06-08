/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CoupledSolid.cpp
 */

#include "CoupledSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "NullModel.hpp"
#include "porosity/BiotPorosity.hpp"
#include "porosity/PressurePorosity.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::CoupledSolid( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_solidModel( nullptr ),
  m_porosityModel( nullptr ),
  m_solidModelName(),
  m_porosityModelName()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Name of the solid model." );

    registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Name of the porosity model." );
}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::~CoupledSolid()
{}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::postProcessInput()
{
  m_solidModel = &this->getParent().getGroup< SOLID_TYPE >( m_solidModelName );
  m_porosityModel = &this->getParent().getGroup< PORO_TYPE >( m_porosityModelName );
}

// Register all CoupleSolid model types.
typedef CoupledSolid< NullModel, PressurePorosity > CompressibleRock;
typedef CoupledSolid< ElasticIsotropic, BiotPorosity > PoroElasticIsotropic;
typedef CoupledSolid< ElasticTransverseIsotropic, BiotPorosity > PoroElasticTransverseIsotropic;
typedef CoupledSolid< DruckerPrager, BiotPorosity > PoroDruckerPrager;
typedef CoupledSolid< DruckerPragerExtended, BiotPorosity > PoroDruckerPragerExtended;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRock, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPragerExtended, string const &, Group * const )

}
} /* namespace geosx */
