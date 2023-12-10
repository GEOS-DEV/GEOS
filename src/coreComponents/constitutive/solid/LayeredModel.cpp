
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

#include "LayeredModel.hpp"
#include "DuvautLionsSolid.hpp"
#include "ElasticIsotropic.hpp"
//#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"

namespace geos
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
{

  registerWrapper( viewKeyStruct::solidModelNameLayer1String(), &m_solidModelNameLayer1 ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model of layer 1."  );

  registerWrapper( viewKeyStruct::solidModelNameLayer2String(), &m_solidModelNameLayer2 ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model of layer 2."  );
}

template< typename SOLID_TYPE1 , 
          typename SOLID_TYPE2 >
LayeredModel< SOLID_TYPE1 , SOLID_TYPE2 >::~LayeredModel() = default;

// template< typename BASE >
// void DuvautLionsSolid< BASE >::postProcessInput()
// {
//   BASE::postProcessInput();
// }

template< typename SOLID_TYPE1 , 
           typename SOLID_TYPE2 >
 void LayeredModel< SOLID_TYPE1 , SOLID_TYPE2 >::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex ) 
 {
   std::cout<<"I'm in LayeredModel::allocateConstitutive data!"<<std::endl;
   SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
   getSolidModelLayer2().allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
   getSolidModelLayer1().allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

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


}
} /* namespace geos */
