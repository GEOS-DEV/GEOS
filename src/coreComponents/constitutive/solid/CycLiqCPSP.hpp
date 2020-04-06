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
 *  @file CycLiqCPSP.hpp
 *
 *  Created on: Feb 27, 2020
 *      Author: ron
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

namespace constitutive
{
/**
 * @class CycLiqCPSP
 *
 * Class to provide a CycLiqCPSP material response.
 */
class CycLiqCPSP : public SolidBase
{
public:
	CycLiqCPSP( string const & name, Group * const parent );

	virtual ~CycLiqCPSP() override;

	virtual void
	DeliverClone( string const & name,
			Group * const parent,
			std::unique_ptr<ConstitutiveBase> & clone ) const override;

	virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
			localIndex const numConstitutivePointsPerParentIndex ) override;

	static constexpr auto m_catalogNameString = "CycLiqCPSP";
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
		static constexpr auto defaultG0String = "defaultG0";
		static constexpr auto defaultKappaString = "defaultKappa";
		static constexpr auto defaultHString = "defaultH";
		static constexpr auto defaultMString = "defaultM";
		static constexpr auto defaultDre1String = "defaultDre1";
		static constexpr auto defaultDre2String = "defaultDre2";
		static constexpr auto defaultDirString = "defaultDir";
		static constexpr auto defaultEtaString = "defaultEta";
		static constexpr auto defaultRdrString = "defaultRdr";
		static constexpr auto defaultNpString = "defaultNp";
		static constexpr auto defaultNdString = "defaultNd";
		static constexpr auto defaultLamdacString = "defaultLamdac";
		static constexpr auto defaultE0String = "defaultE0";
		static constexpr auto defaultKsiString = "defaultKsi";
		static constexpr auto defaultEinString = "defaultEin";
        static constexpr auto defaultInitialTimeString = "defaultInitialTime";

		static constexpr auto G0String = "G0";
		static constexpr auto kappaString = "kappa";
		static constexpr auto hString = "h";
		static constexpr auto MString = "M";
		static constexpr auto dre1String = "dre1";
		static constexpr auto dre2String = "dre2";
		static constexpr auto dirString = "dir";
		static constexpr auto etaString = "eta";
		static constexpr auto rdrString = "rdr";
		static constexpr auto npString = "np";
		static constexpr auto ndString = "nd";
		static constexpr auto lamdacString = "lamdac";
		static constexpr auto e0String = "e0";
		static constexpr auto ksiString = "ksi";
		static constexpr auto einString = "ein";

		static constexpr auto strainString = "strain";
		static constexpr auto epsvirString = "epsvir";
		static constexpr auto epsvreString = "epsvre";
		static constexpr auto gammamonoString = "gammamono";
		static constexpr auto epsvcString = "epsvc";
		static constexpr auto etamString = "etam";
		static constexpr auto alphaString = "alpha";
		static constexpr auto initialTimeString = "initialTime";
	};

	// 如何更新？？
	real64 constrainedModulus(localIndex k) const
	{
		/*
		localIndex numConstitutivePointsPerParentIndex = m_stress[k].size();
		R2SymTensor stress;
		for(localIndex i = 0; i < numConstitutivePointsPerParentIndex; ++i)
		{
			stress += m_stress[k][i];
		}
		stress /= numConstitutivePointsPerParentIndex;
		real64 p = stress.Trace() / 3.0;
		*/
		real64 p = 1e10;

		real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat );
		real64 K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat );
		G = 3.0 / 8.0 * K;
		return ( K + 4 / 3.0 * G );
	}

	real64 compressibility()
	{
		real64 bulkModulus;
		real64 p = 1e10;
		bulkModulus = (1 + m_defaultEin) / m_defaultKappa * pat * sqrt( p / pat );
		return 1 / bulkModulus;
	}
	real64 compressibility() const
	{
		real64 bulkModulus;
		real64 p = 1e10;
		bulkModulus = (1 + m_defaultEin) / m_defaultKappa * pat * sqrt( p / pat );
		return 1 / bulkModulus;
	}

	class KernelWrapper
	{
		public:
		KernelWrapper( arrayView1d<real64 const> const & G0,
		               arrayView1d<real64 const> const & ein,
		               arrayView1d<real64 const> const & kappa) :
		            	   m_G0(G0),
						   m_ein(ein),
						   m_kappa(kappa)
		{}

	    /**
	     * accessor to return the stiffness at a given element
	     * @param k the element number
	     * @param c the stiffness array
	     */
		 GEOSX_HOST_DEVICE inline
		 void GetStiffness( localIndex const k, real64 (&c)[6][6] ) const // 如何更新？？
		 {
			 real64 const G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( 1000 / pat );
		 	 real64 const Lame = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( 1000 / pat ) - 2.0/3.0 * G;

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
		 arrayView1d<real64 const> const m_G0;
		 arrayView1d<real64 const> const m_ein;
		 arrayView1d<real64 const> const m_kappa;
		 real64 const pat = 101.0;

	};

	KernelWrapper createKernelWrapper() const
	{ return KernelWrapper( m_G0,
			m_ein,
			m_kappa); }

	protected:
	virtual void PostProcessInput() override;

	private:
		real64 m_defaultG0;
		real64 m_defaultKappa;
		real64 m_defaultH;
		real64 m_defaultM;
		real64 m_defaultDre1;
		real64 m_defaultDre2;
		real64 m_defaultDir;
		real64 m_defaultEta;
		real64 m_defaultRdr;
		real64 m_defaultNp;
		real64 m_defaultNd;
		real64 m_defaultLamdac;
		real64 m_defaultE0;
		real64 m_defaultKsi;
		real64 m_defaultEin;
        real64 m_defaultInitialTime;

		array1d<real64> m_G0;
		array1d<real64> m_kappa;
		array1d<real64> m_h;
		array1d<real64> m_M;
		array1d<real64> m_dre1;
		array1d<real64> m_dre2;
		array1d<real64> m_dir;
		array1d<real64> m_eta;
		array1d<real64> m_rdr;
		array1d<real64> m_np;
		array1d<real64> m_nd;
		array1d<real64> m_lamdac;
		array1d<real64> m_e0;
		array1d<real64> m_ksi;
		array1d<real64> m_ein;

		array2d<R2SymTensor> m_strain;
		array2d<real64> m_epsvir;
		array2d<real64> m_epsvre;
		array2d<real64> m_gammamono;
		array2d<real64> m_epsvc;
		array2d<real64> m_etam;
	    array2d<R2SymTensor> m_alpha;
	    array2d<real64> m_initialTime;

	    const real64 root23 = sqrt( 2.0 / 3.0 );
	    const real64 root32 = sqrt( 3.0 / 2.0 );
	    const real64 one3 = 1.0 / 3.0;
	   	const real64 two3 = 2.0 / 3.0;
	   	const real64 pat = 101.0;
	   	const real64 pmin = 0.5;
	   	const real64 pcut = 0.5;
	   	const real64 tolerance = 1.0e-8*pcut;

};

} /* namespace constitutive */
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_ */
