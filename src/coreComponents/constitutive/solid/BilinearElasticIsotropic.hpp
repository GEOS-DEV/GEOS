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
 *  @file BilinearElasticIsotropic.hpp
 *
 *  Created on: Feb 13, 2020
 *      Author: ron
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_BILINEARELASTICISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_BILINEARELASTICISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class BilinearElasticIsotropic
 *
 * Class to provide a bilinear elastic isotropic material response.
 */

class BilinearElasticIsotropic : public SolidBase
{
public:
	BilinearElasticIsotropic( string const & name, Group * const parent );

	virtual ~BilinearElasticIsotropic() override;

	virtual void
	DeliverClone( string const & name,
	              Group * const parent,
	              std::unique_ptr<ConstitutiveBase> & clone ) const override;


	virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
	                                       localIndex const numConstitutivePointsPerParentIndex ) override;

	static constexpr auto m_catalogNameString = "BilinearElasticIsotropic";
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
	    static constexpr auto defaultYoungsModulus_1String  = "defaultYoungsModulus_1";
	    static constexpr auto defaultYoungsModulus_2String =  "defaultYoungsModulus_2" ;
	    static constexpr auto defaultPoissonRatio_1String = "defaultPoissonRatio_1";
	    static constexpr auto defaultPoissonRatio_2String =  "defaultPoissonRatio_2" ;
	    static constexpr auto defaultCriticalStressString =  "defaultCriticalStress" ;

	    static constexpr auto youngsModulus_1String  = "youngsModulus_1";
	    static constexpr auto youngsModulus_2String =  "youngsModulus_2" ;
	    static constexpr auto poissonRatio_1String = "poissonRatio_1";
	    static constexpr auto poissonRatio_2String =  "poissonRatio_2" ;
	    static constexpr auto criticalStressString =  "criticalStress" ;

	  };

	real64 constrainedModulus(localIndex k) const { return ( m_youngsModulus_1[k] * ( 1 - m_poissonRatio_1[k] )/ (1 + m_poissonRatio_1[k]) / (1 - 2 * m_poissonRatio_1[k])); }
	real64 compressibility()
	{
		real64 bulkModulus;
		bulkModulus = m_defaultYoungsModulus_1 / (3 * (1 - 2 * m_defaultPoissonRatio_1));
		return 1 / bulkModulus;
	}
	real64 compressibility() const
	{
		real64 bulkModulus;
		bulkModulus = m_defaultYoungsModulus_1 / (3 * (1 - 2 * m_defaultPoissonRatio_1));
		return 1 / bulkModulus;
	}

	  class KernelWrapper
	  {
	  public:
	    KernelWrapper( arrayView1d<real64 const> const & youngsModulus_1,
	                   arrayView1d<real64 const> const & youngsModulus_2,
	                   arrayView1d<real64 const> const & poissonRatio_1,
	                   arrayView1d<real64 const> const & poissonRatio_2,
	                   arrayView1d<real64 const> const & criticalStress) :
	      m_youngsModulus_1(youngsModulus_1),
		  m_youngsModulus_2(youngsModulus_2),
		  m_poissonRatio_1(poissonRatio_1),
		  m_poissonRatio_2(poissonRatio_2),
		  m_criticalStress(criticalStress)
	    {}

	    /**
	     * accessor to return the stiffness at a given element
	     * @param k the element number
	     * @param c the stiffness array
	     */
	    GEOSX_HOST_DEVICE inline
	    void GetStiffness( localIndex const k, real64 (&c)[6][6] ) const
	    {
	    	real64 const G = m_youngsModulus_1[k] / (2 * (1 + m_poissonRatio_1[k]));
			real64 const Lame = m_poissonRatio_1[k] * m_youngsModulus_1[k] / (1 + m_poissonRatio_1[k]) / (1 - 2 * m_poissonRatio_1[k]);

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
	    arrayView1d<real64 const> const m_youngsModulus_1;
	    arrayView1d<real64 const> const m_youngsModulus_2;
	    arrayView1d<real64 const> const m_poissonRatio_1;
	    arrayView1d<real64 const> const m_poissonRatio_2;
	    arrayView1d<real64 const> const m_criticalStress;

	  };

	  KernelWrapper createKernelWrapper() const
	  { return KernelWrapper( m_youngsModulus_1,
			  m_youngsModulus_2,
			  m_poissonRatio_1,
			  m_poissonRatio_2,
			  m_criticalStress); }

protected:
    virtual void PostProcessInput() override;

private:
  /// reference pressure parameter
    //real64 m_referencePressure;

  /// scalar compressibility parameter
    //real64 m_compressibility;

    real64 m_defaultYoungsModulus_1;
    real64 m_defaultYoungsModulus_2;
    real64 m_defaultPoissonRatio_1;
    real64 m_defaultPoissonRatio_2;
    real64 m_defaultCriticalStress;

    array1d<real64> m_youngsModulus_1;
    array1d<real64> m_youngsModulus_2;
    array1d<real64> m_poissonRatio_1;
    array1d<real64> m_poissonRatio_2;
    array1d<real64> m_criticalStress;
};

} /* namespace constitutive */
} /* namespace geosx */

#endif /*GEOSX_CONSTITUTIVE_SOLID_BILINEARELASTICISOTROPIC_HPP_ */
