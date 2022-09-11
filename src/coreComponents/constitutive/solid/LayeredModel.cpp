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
 * @file LayeredModel.cpp
 */
#include "LayeredModel.hpp"
#include "ElasticIsotropic.hpp"
//#include "ElasticTransverseIsotropic.hpp"
//#include "ElasticOrthotropic.hpp"
//#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
//#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE1 , 
          typename SOLID_TYPE2 >
LayeredModel< SOLID_TYPE1 , SOLID_TYPE2 >::LayeredModel( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_solidModelNameLayer1(),
  m_solidModelNameLayer2()
  //m_solidUpdate(),
  //m_relaxationTime()
{
  registerWrapper( viewKeyStruct::solidModelNameLayer1String(), &m_solidModelNameLayer1 ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model of layer 1." );

  registerWrapper( viewKeyStruct::solidModelNameLayer2String(), &m_solidModelNameLayer2 ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model of layer 2." );

  // registerWrapper( viewKeyStruct::relaxationTimeString(), &m_relaxationTime ).
  //   setInputFlag( dataRepository::InputFlags::REQUIRED ).
  //   setDescription( "Relaxation time" );
}

template< typename SOLID_TYPE1 , 
          typename SOLID_TYPE2 >
LayeredModel< SOLID_TYPE1 , SOLID_TYPE2 >::~LayeredModel() = default;


template< typename SOLID_TYPE1 , 
          typename SOLID_TYPE2 >
void LayeredModel< SOLID_TYPE1 , SOLID_TYPE2 >::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{

  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

// Register all LayeredModel model types. Uncomment the ones listed as needed. 
typedef LayeredModel< ElasticIsotropic , DruckerPrager > LayeredElasticIsotropicDruckerPrager;
//typedef LayeredModel< ElasticTransverseIsotropic > ViscoElasticTransverseIsotropic;
//typedef LayeredModel< ElasticOrthotropic > ViscoElasticOrthotropic;
//typedef LayeredModel< DelftEgg > ViscoDelftEgg;
//typedef LayeredModel< DruckerPrager > ViscoDruckerPrager;
//typedef LayeredModel< DruckerPragerExtended > ViscoDruckerPragerExtended;
//typedef LayeredModel< ModifiedCamClay > ViscoModifiedCamClay;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LayeredElasticIsotropicDruckerPrager, string const &, Group * const )
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
