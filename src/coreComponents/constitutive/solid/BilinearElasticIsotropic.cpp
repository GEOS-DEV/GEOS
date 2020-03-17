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

/*
 *  @file BilinearElasticIsotropic.cpp
 *
 *  Created on: Feb 13, 2020
 *      Author: ron
 */

#include "BilinearElasticIsotropic.hpp"

namespace geosx {
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive {

BilinearElasticIsotropic::BilinearElasticIsotropic( std::string const & name, Group * const parent ):
		  SolidBase( name, parent ),
		  m_defaultYoungsModulus_1(),
		  m_defaultYoungsModulus_2(),
		  m_defaultPoissonRatio_1(),
		  m_defaultPoissonRatio_2(),
		  m_defaultCriticalStress()

{
	registerWrapper( viewKeyStruct::defaultYoungsModulus_1String, &m_defaultYoungsModulus_1, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Elastic Youngs Modulus Parameter 1");

	registerWrapper( viewKeyStruct::defaultYoungsModulus_2String, &m_defaultYoungsModulus_2, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Elastic Youngs Modulus Parameter 2");

	registerWrapper( viewKeyStruct::defaultPoissonRatio_1String, &m_defaultPoissonRatio_1, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Elastic Poisson Ratio Parameter 1");

	registerWrapper( viewKeyStruct::defaultPoissonRatio_2String, &m_defaultPoissonRatio_2, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Elastic Poisson Ratio Parameter 2");

	registerWrapper( viewKeyStruct::defaultCriticalStressString, &m_defaultCriticalStress, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Elastic Critical Stress Parameter");

	registerWrapper( viewKeyStruct::youngsModulus_1String, &m_youngsModulus_1, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("Elastic Youngs Modulus 1 Field");
	registerWrapper( viewKeyStruct::youngsModulus_2String, &m_youngsModulus_2, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("Elastic Youngs Modulus 2 Field");
	registerWrapper( viewKeyStruct::poissonRatio_1String, &m_poissonRatio_1, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("Elastic Poisson Ratio 1 Field");
	registerWrapper( viewKeyStruct::poissonRatio_2String, &m_poissonRatio_2, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("Elastic Poisson Ratio 2 Field");
	registerWrapper( viewKeyStruct::criticalStressString, &m_criticalStress, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("Elastic Critical Stress Field");

}

BilinearElasticIsotropic::~BilinearElasticIsotropic()
{}

void
BilinearElasticIsotropic::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<BilinearElasticIsotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  BilinearElasticIsotropic * const newConstitutiveRelation = dynamic_cast<BilinearElasticIsotropic *>(clone.get());

  newConstitutiveRelation->m_defaultYoungsModulus_1 = m_defaultYoungsModulus_1;
  newConstitutiveRelation->m_youngsModulus_1 = m_youngsModulus_1;
  newConstitutiveRelation->m_defaultYoungsModulus_2 = m_defaultYoungsModulus_2;
  newConstitutiveRelation->m_youngsModulus_2 = m_youngsModulus_2;
  newConstitutiveRelation->m_defaultPoissonRatio_1 = m_defaultPoissonRatio_1;
  newConstitutiveRelation->m_poissonRatio_1 = m_poissonRatio_1;
  newConstitutiveRelation->m_defaultPoissonRatio_2 = m_defaultPoissonRatio_2;
  newConstitutiveRelation->m_poissonRatio_2 = m_poissonRatio_2;
  newConstitutiveRelation->m_defaultCriticalStress = m_defaultCriticalStress;
  newConstitutiveRelation->m_criticalStress = m_criticalStress;
  newConstitutiveRelation->m_stress = m_stress;
  //newConstitutiveRelation->m_stiffness = m_stiffness;
}

void BilinearElasticIsotropic::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_youngsModulus_1.resize( parent->size() );
  m_youngsModulus_2.resize( parent->size() );
  m_poissonRatio_1.resize( parent->size() );
  m_poissonRatio_2.resize( parent->size() );
  m_criticalStress.resize( parent->size() );

  m_youngsModulus_1 = m_defaultYoungsModulus_1;
  m_youngsModulus_2 = m_defaultYoungsModulus_2;
  m_poissonRatio_1 = m_defaultPoissonRatio_1;
  m_poissonRatio_2 = m_defaultPoissonRatio_2;
  m_criticalStress = m_defaultCriticalStress;

}

void BilinearElasticIsotropic::PostProcessInput()
{
  m_youngsModulus_1 = m_defaultYoungsModulus_1;
  m_youngsModulus_2 = m_defaultYoungsModulus_2;
  m_poissonRatio_1 = m_defaultPoissonRatio_1;
  m_poissonRatio_2 = m_defaultPoissonRatio_2;
  m_criticalStress = m_defaultCriticalStress;
}

void BilinearElasticIsotropic::StateUpdatePoint( localIndex const k,
                                               localIndex const q,
                                               R2SymTensor const & D,
                                               R2Tensor const & Rot,
                                               real64 const GEOSX_UNUSED_ARG( dt ),
                                               integer const GEOSX_UNUSED_ARG( updateStiffnessFlag ) )
{
	real64 const * const restrict Mdata = m_stress[k][q].Data();
	real64 shearModulus, bulkModulus;
	if((Mdata[0]-Mdata[2]) * (Mdata[0]-Mdata[2]) + (Mdata[2]-Mdata[5]) * (Mdata[2]-Mdata[5]) + (Mdata[5]-Mdata[0]) * (Mdata[5]-Mdata[0])
			+ 6 * (Mdata[1]*Mdata[1] + Mdata[3]*Mdata[3] + Mdata[4]*Mdata[4]) < 2 * m_criticalStress[k] *  m_criticalStress[k])
	{
		shearModulus = m_youngsModulus_1[k] / (2 * (1 + m_poissonRatio_1[k]));
		bulkModulus = m_youngsModulus_1[k] / (3 * (1 - 2 * m_poissonRatio_1[k]));

	}
	else
	{
		shearModulus = m_youngsModulus_2[k] / (2 * (1 + m_poissonRatio_2[k]));
		bulkModulus = m_youngsModulus_2[k] / (3 * (1 - 2 * m_poissonRatio_2[k]));
	}
  real64 meanStresIncrement = D.Trace();
  R2SymTensor temp = D;
  temp.PlusIdentity( -meanStresIncrement / 3.0 );
  temp *= 2.0 * shearModulus;
  meanStresIncrement *= bulkModulus;
  temp.PlusIdentity( meanStresIncrement );

  m_stress[k][q] += temp;

  temp.QijAjkQlk( m_stress[k][q], Rot );
  m_stress[k][q] = temp;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BilinearElasticIsotropic, std::string const &, Group * const )
} /* namespace constitutive */
} /* namespace geosx */
