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
 *  @file HardeningElasticIsotropic.hpp
 *
 *  Created on: Feb 16, 2020
 *      Author: ron
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_HARDENINGELASTICISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_HARDENINGELASTICISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

namespace constitutive
{
/**
 * @class HardeningElasticIsotropic
 *
 * Class to provide a hardening (linear with ) elastic isotropic material response.
 */

class HardeningElasticIsotropic : public SolidBase
{
public:
	HardeningElasticIsotropic( string const & name, Group * const parent );

	virtual ~HardeningElasticIsotropic() override;


	virtual void
	DeliverClone( string const & name,
	              Group * const parent,
	              std::unique_ptr<ConstitutiveBase> & clone ) const override;


	virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
	                                       localIndex const numConstitutivePointsPerParentIndex ) override;

	static constexpr auto m_catalogNameString = "HardeningElasticIsotropic";
	static std::string CatalogName() { return m_catalogNameString; }
	virtual string GetCatalogName() override { return CatalogName(); }

    virtual void StateUpdatePoint( localIndex const k,
	                               localIndex const q,
	                               R2SymTensor const & D,
	                               R2Tensor const & Rot,
	                               real64 const dt,
	                               integer const updateStiffnessFlag ) override;

	struct viewKeyStruct : public SolidBase::viewKeyStruct
	  {
	    static constexpr auto defaultYoungsModulusString  = "defaultYoungsModulus";
	    static constexpr auto defaultHardeningFactorString = "defaultHardeningFactor";
	    static constexpr auto defaultReferenceVoidRatioString =  "defaultReferenceVoidRatio";
	    static constexpr auto defaultPoissonRatioString =  "defaultPoissonRatio" ;

	    static constexpr auto voidRatioString  = "voidRatio";
	  };

	real64 constrainedModulus(localIndex k) const
	{
		return ( (m_youngsModulus[k] - m_hardeningFactor[k] * m_referenceVoidRatio[k]) * ( 1 - m_poissonRatio[k] )
				/ (1 + m_poissonRatio[k]) / (1 - 2 * m_poissonRatio[k]));
	}

	  class KernelWrapper
	  {
	  public:
	    KernelWrapper( arrayView1d<real64 const> const & youngsModulus,
	                   arrayView1d<real64 const> const & hardeningFactor,
	                   arrayView1d<real64 const> const & referenceVoidRatio,
	                   arrayView1d<real64 const> const & poissonRatio) :
	      m_youngsModulus(youngsModulus),
		  m_hardeningFactor(hardeningFactor),
		  m_referenceVoidRatio(referenceVoidRatio),
		  m_poissonRatio(poissonRatio)
	    {}

	    /**
	     * accessor to return the stiffness at a given element
	     * @param k the element number
	     * @param c the stiffness array
	     */
	    GEOSX_HOST_DEVICE inline
	    void GetStiffness( localIndex const k, real64 (&c)[6][6] ) const
	    {
	    	real64 const G = (m_youngsModulus[k] - m_hardeningFactor[k] * m_referenceVoidRatio[k]) / (2 * (1 + m_poissonRatio[k]));
			real64 const Lame = m_poissonRatio[k] *(m_youngsModulus[k] - m_hardeningFactor[k] * m_referenceVoidRatio[k])
					/ (1 + m_poissonRatio[k]) / (1 - 2 * m_poissonRatio[k]);

	        memset( c, 0, sizeof( c ) );

	        c[0][0] = Lame + 2 * G;
	        c[0][1] = Lame;
	        c[0][2] = Lame;

	        c[1][0] = Lame;
	        c[1][1] = Lame + 2 * G;
	        c[1][2] = Lame;

	        c[2][0] = Lame;
	        c[2][1] = Lame;
	        c[2][2] = Lame + 2 * G;

	        c[3][3] = G;

	        c[4][4] = G;

	        c[5][5] = G;
	    }

	  private:
	    arrayView1d<real64 const> const m_youngsModulus;
	    arrayView1d<real64 const> const m_hardeningFactor;
	    arrayView1d<real64 const> const m_referenceVoidRatio;
	    arrayView1d<real64 const> const m_poissonRatio;

	  };

	  KernelWrapper createKernelWrapper() const
	  { return KernelWrapper( m_youngsModulus,
			  m_hardeningFactor,
			  m_referenceVoidRatio,
			  m_poissonRatio); }

protected:
  virtual void PostProcessInput() override;

	private:
	  /// reference pressure parameter
	    //real64 m_referencePressure;

	  /// scalar compressibility parameter
	    //real64 m_compressibility;

	    real64 m_defaultYoungsModulus;
	    real64 m_defaultHardeningFactor;
	    real64 m_defaultReferenceVoidRatio;
	    real64 m_defaultPoissonRatio;

	    array1d<real64> m_youngsModulus;
	    array1d<real64> m_hardeningFactor;
	    array1d<real64> m_referenceVoidRatio;
	    array1d<real64> m_poissonRatio;

	    array2d<real64> m_voidRatio;
	    //array2d<R2SymTensor> m_strain;
};

} /* namespace constitutive */
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_HARDENINGELASTICISOTROPIC_HPP_ */
