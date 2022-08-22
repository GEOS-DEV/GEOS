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
 * @file DuvautLionsSolid.cpp
 */

#include "DuvautLionsSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "ElasticOrthotropic.hpp"
#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
DuvautLionsSolid< SOLID_TYPE >::DuvautLionsSolid( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_solidModelName(),
  m_solidUpdate(),
  m_relaxationTime()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::relaxationTimeString(), &m_relaxationTime ).
    setApplyDefaultValue( -1 ).
    setDescription( "Relaxation time" );
}

template< typename SOLID_TYPE >
DuvautLionsSolid< SOLID_TYPE >::~DuvautLionsSolid() = default;

// Register all DuvautLionsSolid model types. Uncomment the ones listed as needed. 
typedef DuvautLionsSolid< ElasticIsotropic > ViscoElasticIsotropic;
//typedef DuvautLionsSolid< ElasticTransverseIsotropic > ViscoElasticTransverseIsotropic;
//typedef DuvautLionsSolid< ElasticOrthotropic > ViscoElasticOrthotropic;
//typedef DuvautLionsSolid< DelftEgg > ViscoDelftEgg;
//typedef DuvautLionsSolid< DruckerPrager > ViscoDruckerPrager;
//typedef DuvautLionsSolid< DruckerPragerExtended > ViscoDruckerPragerExtended;
//typedef DuvautLionsSolid< ModifiedCamClay > ViscoModifiedCamClay;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoElasticIsotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoElasticTransverseIsotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoElasticOrthotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDelftEgg, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDruckerPrager, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDruckerPragerExtended, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDamageElasticIsotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDamageSpectralElasticIsotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDamageVolDevElasticIsotropic, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoModifiedCamClay, string const &, Group * const )

}
} /* namespace geosx */
