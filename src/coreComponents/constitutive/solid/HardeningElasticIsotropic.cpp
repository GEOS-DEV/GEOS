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
 *  @file HardeningElasticIsotropic.cpp
 *
 *  Created on: Feb 16, 2020
 *      Author: ron
 */

#include "HardeningElasticIsotropic.hpp"
#include <iostream>

namespace geosx {
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive {

HardeningElasticIsotropic::HardeningElasticIsotropic( std::string const & name, Group * const parent ):
				  SolidBase( name, parent ),
			      m_defaultYoungsModulus(),
				  m_defaultHardeningFactor(),
				  m_defaultReferenceVoidRatio(),
				  m_defaultPoissonRatio(),
				  m_voidRatio{}

		{
			registerWrapper( viewKeyStruct::defaultYoungsModulusString, &m_defaultYoungsModulus, 0 )->
			    setApplyDefaultValue(-1)->
			    setInputFlag(InputFlags::OPTIONAL)->
			    setDescription("Elastic Youngs Modulus Parameter");

			registerWrapper( viewKeyStruct::defaultHardeningFactorString, &m_defaultHardeningFactor, 0 )->
			    setApplyDefaultValue(-1)->
			    setInputFlag(InputFlags::OPTIONAL)->
			    setDescription("Hardening Factor Parameter");

			registerWrapper( viewKeyStruct::defaultReferenceVoidRatioString, &m_defaultReferenceVoidRatio, 0 )->
			    setApplyDefaultValue(-1)->
			    setInputFlag(InputFlags::OPTIONAL)->
			    setDescription("Reference Void Ratio Parameter 1");

			registerWrapper( viewKeyStruct::defaultPoissonRatioString, &m_defaultPoissonRatio, 0 )->
			    setApplyDefaultValue(-1)->
			    setInputFlag(InputFlags::OPTIONAL)->
			    setDescription("Elastic Poisson Ratio Parameter");

		    registerWrapper( viewKeyStruct::voidRatioString, &m_voidRatio, 0 )->
			    setPlotLevel(PlotLevel::LEVEL_0)->
			    setDescription("Void Ratio");

		    //registerWrapper( viewKeyStruct::stressString, &m_strain, 0 )->
		    //  setPlotLevel(PlotLevel::LEVEL_0)->
		    //  setDescription("Strain Deviator");

		}

HardeningElasticIsotropic::~HardeningElasticIsotropic()
{}

void
HardeningElasticIsotropic::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<HardeningElasticIsotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  HardeningElasticIsotropic * const newConstitutiveRelation = dynamic_cast<HardeningElasticIsotropic *>(clone.get());

  newConstitutiveRelation->m_defaultYoungsModulus = m_defaultYoungsModulus;
  newConstitutiveRelation->m_youngsModulus = m_youngsModulus;
  newConstitutiveRelation->m_defaultHardeningFactor = m_defaultHardeningFactor;
  newConstitutiveRelation->m_hardeningFactor = m_hardeningFactor;
  newConstitutiveRelation->m_defaultReferenceVoidRatio = m_defaultReferenceVoidRatio;
  newConstitutiveRelation->m_referenceVoidRatio = m_referenceVoidRatio;
  newConstitutiveRelation->m_defaultPoissonRatio = m_defaultPoissonRatio;
  newConstitutiveRelation->m_poissonRatio = m_poissonRatio;

  newConstitutiveRelation->m_voidRatio = m_voidRatio;
  //newConstitutiveRelation->m_strain = m_strain;
}

void HardeningElasticIsotropic::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_youngsModulus.resize( parent->size() );
  m_hardeningFactor.resize( parent->size() );
  m_referenceVoidRatio.resize( parent->size() );
  m_poissonRatio.resize( parent->size() );
  m_voidRatio.resize( parent->size(), numConstitutivePointsPerParentIndex );
  //m_strain.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_youngsModulus = m_defaultYoungsModulus;
  m_hardeningFactor = m_defaultHardeningFactor;
  m_referenceVoidRatio = m_defaultReferenceVoidRatio;
  m_poissonRatio = m_defaultPoissonRatio;
  m_voidRatio = m_defaultReferenceVoidRatio;
}

void HardeningElasticIsotropic::PostProcessInput()
{
	m_youngsModulus = m_defaultYoungsModulus;
	m_hardeningFactor = m_defaultHardeningFactor;
	m_referenceVoidRatio = m_defaultReferenceVoidRatio;
	m_poissonRatio = m_defaultPoissonRatio;

	m_voidRatio = m_defaultReferenceVoidRatio;
}

void HardeningElasticIsotropic::StateUpdatePoint( localIndex const k,
                                               localIndex const q,
                                               R2SymTensor const & D,
                                               R2Tensor const & Rot,
                                               real64 const GEOSX_UNUSED_ARG( dt ),
                                               integer const GEOSX_UNUSED_ARG( updateStiffnessFlag ) )
{
	  real64 meanStresIncrement = D.Trace();
	if(m_voidRatio[k][q] > 0)
	{
		m_voidRatio[k][q] += meanStresIncrement * (1 + m_voidRatio[k][q]);
	}
	real64 youngsModulus, shearModulus, bulkModulus;
	youngsModulus = m_youngsModulus[k] + m_hardeningFactor[k] * (m_voidRatio[k][q] - m_referenceVoidRatio[k]);
	shearModulus = youngsModulus / (2 * (1 + m_poissonRatio[k]));
	bulkModulus = youngsModulus / (3 * (1 - 2 * m_poissonRatio[k]));
	R2SymTensor temp = D;
	temp.PlusIdentity( -meanStresIncrement / 3.0 );
	temp *= 2.0 * shearModulus;
	meanStresIncrement *= bulkModulus;
	temp.PlusIdentity( meanStresIncrement );

	m_stress[k][q] += temp;

	temp.QijAjkQlk( m_stress[k][q], Rot );
	m_stress[k][q] = temp;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HardeningElasticIsotropic, std::string const &, Group * const )
} /* namespace constitutive */
} /* namespace geosx */
