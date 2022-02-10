/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file CycLiqCPSP.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_

#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class CycLiqCPSPUpdates
 *
 * Class to provide CycLiqCPSP material updates that may be
 * called from a kernel function.
 */
class CycLiqCPSPUpdates : public SolidBaseUpdates
{

public:
	CycLiqCPSPUpdates(arrayView1d< real64 const > const & G0,
            arrayView1d< real64 const > const & kappa,
            arrayView1d< real64 const > const & h,
            arrayView1d< real64 const > const & dre1,
            arrayView1d< real64 const > const & dre2,
            arrayView1d< real64 const > const & dir,
            arrayView1d< real64 const > const & eta,
            arrayView1d< real64 const > const & rdr,
            arrayView1d< real64 const > const & np,
			arrayView1d< real64 const > const & nd,
            arrayView1d< real64 const > const & M,
            arrayView1d< real64 const > const & lamdac,
            arrayView1d< real64 const > const & e0,
            arrayView1d< real64 const > const & ksi,
            arrayView1d< real64 const > const & ein,
			arrayView2d< integer > const & initialCycles,
			arrayView2d< real64 > const & epsvir,
			arrayView2d< real64 > const & epsvre,
			arrayView2d< real64 > const & gammamono,
			arrayView2d< real64 > const & epsvc,
			arrayView2d< real64 > const & etam,
			arrayView3d< real64, solid::STRESS_USD > const & alpha,
			arrayView3d< real64, solid::STRESS_USD > const & strain,
            arrayView3d< real64, solid::STRESS_USD > const & newStress,
            arrayView3d< real64, solid::STRESS_USD > const & oldStress):
        SolidBaseUpdates( newStress, oldStress ),
	    m_G0( G0 ),
	    m_kappa( kappa ),
	    m_h( h ),
	    m_dre1( dre1 ),
	    m_dre2( dre2 ),
	    m_dir( dir ),
	    m_eta( eta ),
	    m_rdr( rdr ),
	    m_np( np ),
	    m_nd( nd ),
	    m_M( M ),
	    m_lamdac( lamdac ),
	    m_e0( e0 ),
	    m_ksi( ksi ),
	    m_ein( ein ),
	    m_initialCycles( initialCycles ),
	    m_epsvir( epsvir ),
	    m_epsvre( epsvre ),
	    m_gammamono( gammamono ),
	    m_epsvc( epsvc ),
	    m_etam( etam ),
	    m_alpha( alpha ),
	    m_strain( strain )
	{}

	/// Deleted default constructor
	CycLiqCPSPUpdates() = delete;

	/// Default copy constructor
	CycLiqCPSPUpdates( CycLiqCPSPUpdates const & ) = default;

	/// Default move constructor
	CycLiqCPSPUpdates( CycLiqCPSPUpdates && ) = default;

	/// Deleted copy assignment operator
	CycLiqCPSPUpdates & operator=( CycLiqCPSPUpdates const & ) = delete;

	/// Deleted move assignment operator
	CycLiqCPSPUpdates & operator=( CycLiqCPSPUpdates && ) =  delete;

	/// Use the "isotropic" form of inner product compression
	using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

	/// Use base version of saveConvergedState
	using SolidBaseUpdates::saveConvergedState;


	GEOSX_HOST_DEVICE
	virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
	                                                    localIndex const q,
	                                                    real64 const ( &totalStrain )[6],
	                                                    real64 ( &stress )[6] ) const override final;

	GEOSX_HOST_DEVICE
	virtual void smallStrainNoStateUpdate( localIndex const k,
	                                         localIndex const q,
	                                         real64 const ( &totalStrain )[6],
	                                         real64 ( &stress )[6],
	                                         real64 ( &stiffness )[6][6] ) const override final;

	GEOSX_HOST_DEVICE
	virtual void smallStrainNoStateUpdate( localIndex const k,
	                                         localIndex const q,
	                                         real64 const ( &totalStrain )[6],
	                                         real64 ( &stress )[6],
	                                         DiscretizationOps & stiffness ) const final;

	GEOSX_HOST_DEVICE
	virtual void smallStrainUpdate_StressOnly( localIndex const k,
	                                             localIndex const q,
	                                             real64 const ( &strainIncrement )[6],
	                                             real64 ( &stress )[6] ) const override;

	GEOSX_HOST_DEVICE
	virtual void smallStrainUpdate( localIndex const k,
	                                  localIndex const q,
	                                  real64 const ( &strainIncrement )[6],
	                                  real64 ( &stress )[6],
	                                  real64 ( &stiffness )[6][6] ) const override;

	GEOSX_HOST_DEVICE
	virtual void smallStrainUpdate( localIndex const k,
	                                  localIndex const q,
	                                  real64 const ( &strainIncrement )[6],
	                                  real64 ( &stress )[6],
	                                  DiscretizationOps & stiffness ) const;

	GEOSX_HOST_DEVICE
	virtual void getElasticStiffness( localIndex const k,
	                                    localIndex const q,
	                                    real64 ( &stiffness )[6][6] ) const override;

protected:

  /// A reference to the ArrayView holding the bulk modulus for each element.
  // arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the shear modulus for each element.
  // arrayView1d< real64 const > const m_shearModulus;

  arrayView1d< real64 const > const m_G0;
  arrayView1d< real64 const > const m_kappa;
  arrayView1d< real64 const > const m_h;
  arrayView1d< real64 const > const m_dre1;
  arrayView1d< real64 const > const m_dre2;
  arrayView1d< real64 const > const m_dir;
  arrayView1d< real64 const > const m_eta;
  arrayView1d< real64 const > const m_rdr;
  arrayView1d< real64 const > const m_np;
  arrayView1d< real64 const > const m_nd;
  arrayView1d< real64 const > const m_M;
  arrayView1d< real64 const > const m_lamdac;
  arrayView1d< real64 const > const m_e0;
  arrayView1d< real64 const > const m_ksi;
  arrayView1d< real64 const > const m_ein;
  arrayView2d< integer > const m_initialCycles;
  arrayView2d< real64 > const m_epsvir;
  arrayView2d< real64 > const m_epsvre;
  arrayView2d< real64 > const m_gammamono;
  arrayView2d< real64 > const m_epsvc;
  arrayView2d< real64 > const m_etam;
  arrayView3d< real64, solid::STRESS_USD > const m_alpha;
  arrayView3d< real64, solid::STRESS_USD > const m_strain;

  /// some constant numbers
   const real64 root23 = sqrt( 2.0 / 3.0 );
   const real64 root32 = sqrt( 3.0 / 2.0 );
   const real64 one3 = 1.0 / 3.0;
   const real64 two3 = 2.0 / 3.0;
   /// The atmospheric pressure
   const real64 pat = 101.0;
   /// The minimum effective spherical stress
   const real64 pmin = 0.5;
   const real64 pcut = 0.5;
   const real64 tolerance = 1.0e-8*pcut;

};


/// To be corrected by ron
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::getElasticStiffness( localIndex const k,
                                                   localIndex const q,
                                                   real64 ( & stiffness )[6][6] ) const
{
  GEOSX_UNUSED_VAR( q );
  real64 const p = 1e3;
  real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat )*1000.0;
  real64 const K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat )*1000.0;
  G = 0.3835 * K;
  real64 const lambda = K - 2.0/3.0 * G;

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  stiffness[0][0] = lambda + 2*G;
  stiffness[0][1] = lambda;
  stiffness[0][2] = lambda;

  stiffness[1][0] = lambda;
  stiffness[1][1] = lambda + 2*G;
  stiffness[1][2] = lambda;

  stiffness[2][0] = lambda;
  stiffness[2][1] = lambda;
  stiffness[2][2] = lambda + 2*G;

  stiffness[3][3] = G;
  stiffness[4][4] = G;
  stiffness[5][5] = G;
}


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void CycLiqCPSPUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                                     localIndex const q,
                                                                     real64 const ( &totalStrain )[6],
                                                                     real64 ( & stress )[6] ) const
{
	++ m_initialCycles( k, q );
	if(m_initialCycles( k, q ) < 0)
	{
	     real64 p = 1e3;
	     real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat )*1000.0;
	     real64 K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat )*1000.0;
	     G = 0.3835 * K;
	     real64 const lambda = K - 2.0 / 3.0 * G;
	     real64 const volStrain = totalStrain[0] + totalStrain[1] + totalStrain[2] ;
	     real64 const twoG = 2.0 * G;

	     stress[0] = lambda * volStrain + twoG * totalStrain[0];
	     stress[1] = lambda * volStrain + twoG * totalStrain[1];
	     stress[2] = lambda * volStrain + twoG * totalStrain[2];
	     stress[3] = G * totalStrain[3];
	     stress[4] = G * totalStrain[4];
	     stress[5] = G * totalStrain[5];
	}
	else
	{
		/*
	     real64 p = 1e3;
	     real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat )*1000.0;
	     real64 K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat )*1000.0;
	     G = 0.3835 * K;
	     real64 const lambda = K - 2.0 / 3.0 * G;
	     real64 const volStrain = totalStrain[0] + totalStrain[1] + totalStrain[2] ;
	     real64 const twoG = 2.0 * G;

	     stress[0] = lambda * volStrain + twoG * totalStrain[0];
	     stress[1] = lambda * volStrain + twoG * totalStrain[1];
	     stress[2] = lambda * volStrain + twoG * totalStrain[2];
	     stress[3] = G * totalStrain[3];
	     stress[4] = G * totalStrain[4];
	     stress[5] = G * totalStrain[5];
	     */
		real64 Mfc = m_M[k];
		real64 Mdc = m_M[k];
		real64 sinphi = 3.0 * Mfc / (Mfc + 6.0);
		real64 tanphi = sinphi / sqrt(1.0 - sinphi * sinphi);
		real64 Mfo = 2 * sqrt(3.0) * tanphi / sqrt(3.0 + 4.0 * tanphi * tanphi);
		real64 epsvir_n = m_epsvir[k][q];
		real64 epsvre_n = m_epsvre[k][q];
		real64 epsvc_n = m_epsvc[k][q];

		real64 dev_strain[6] = {0};
		real64 dev_strain_n[6] = {0};
	    real64 ddev_strain_p[6] = {0};
	    real64 dev_stress[6] = {0};
	    real64 dev_stress_n[6] = {0};
	    real64 normal[6] = {0};
	    real64 pass[6] = {0};
	    real64 alpha_ns[6] = {0};
	    real64 rbar[6] = {0};
	    real64 rbar0[6] = {0};
	    real64 rbar1[6] = {0};
	    real64 r1[6] = {0};
	    real64 rd[6] = {0};

	    real64 stress_n[6] = {0};
	    LvArray::tensorOps::copy< 6 >( stress_n, m_oldStress[k][q] );
	    LvArray::tensorOps::scale< 6 >(stress_n, 1.0 / 1000);
	    real64 strain_n[6] = {0};
        LvArray::tensorOps::copy< 6 >( strain_n, m_strain[k][q] );
        real64 strain_nplus1[6] = {0};
        LvArray::tensorOps::copy< 6 >( strain_nplus1, strain_n );
        LvArray::tensorOps::add< 6 >( strain_nplus1, totalStrain );
        real64 alpha_n[6] = {0};
        LvArray::tensorOps::copy< 6 >( alpha_n, m_alpha[k][q] );
        real64 stress_pass[6] = {0};
        real64 alpha_nplus1[6] = {0};
        real64 r[6] = {0};
        real64 r_nplus1[6] = {0};

    	real64 phi = 0.0;
    	real64 phi_n = 0.0;
    	real64 trace = 0.0;
    	real64 trace_n = 0.0;
    	real64 dtracep = 0.0;
    	real64 iconv = 0.0;
    	real64 lambdamax = 0.0;
    	real64 lambdamin = 0.0;
    	real64 dlambda;
    	real64 chi = 0.0;
    	real64 epsvre_ns, epsvir_ns, epsvc_ns;
    	real64 sin3theta, beta0, beta1, Fb0, Fb1, Fb = 0.0, gtheta;
    	real64 intm;
    	real64 beta;
    	real64 ec, psi;
    	int isub, sub1, sub2, sub;
    	real64 N, H;
    	real64 loadindex;
    	real64 rou;
    	real64 roubar;
    	real64 eta_nplus1, Dre_n, Dir_n, D_;

    	//compression as positive
    	LvArray::tensorOps::scale< 6 >(stress_n, -1);
    	LvArray::tensorOps::scale< 6 >(strain_n, -1);
    	LvArray::tensorOps::scale< 6 >(strain_nplus1, -1);

    	//compute the deviatoric stress of last step
    	real64 p_n = LvArray::tensorOps::symTrace< 3 >( stress_n ) * one3;
    	LvArray::tensorOps::copy< 6 >( dev_stress_n, stress_n );
    	LvArray::tensorOps::symAddIdentity< 3 >( dev_stress_n, -p_n );
    	if( p_n < pmin )
    	{
    		p_n = pmin;
    	}

    	//compute the deviatoric strains
    	trace = LvArray::tensorOps::symTrace< 3 >( strain_nplus1 );
    	trace_n = LvArray::tensorOps::symTrace< 3 >( strain_n );
    	LvArray::tensorOps::copy< 6 >( dev_strain_n, strain_n );
    	LvArray::tensorOps::symAddIdentity< 3 >( dev_strain_n, -trace_n * one3 );
    	LvArray::tensorOps::copy< 6 >( dev_strain, strain_nplus1 );
    	LvArray::tensorOps::symAddIdentity< 3 >( dev_strain, -trace * one3 );

    	real64 en = (1 + m_ein[k]) * exp(-trace_n) - 1.0;//void ratio
    	real64 temp[6] = {0};
    	real64 p0;
    	real64 epsvc0;

    	// --------------(I)Initialize-------------------------------------
    	LvArray::tensorOps::copy< 6 >( alpha_nplus1, alpha_n );
    	real64 epsvir_nplus1 = epsvir_n;
    	real64 epsvre_nplus1 = epsvre_n;
    	real64 epsvc_nplus1 = epsvc_n + trace - trace_n;
    	real64 etamplus1 = m_etam[k][q];
    	real64 lambda = 0.0;
        p0 = p_n;
        epsvc0 = - 2 * m_kappa[k] / (1 + m_ein[k]) * ( sqrt( p0 / pat) - sqrt( pmin / pat ));
        LvArray::tensorOps::copy< 6 >( r, dev_stress_n);
        LvArray::tensorOps::scale< 6 >(r, 1.0 / p_n);
    	ec = m_e0[k] - m_lamdac[k] * pow( p_n / pat , m_ksi[k]);
    	psi = en - ec;

    	// --------------Trial Before Substep-------------------------------------
    	real64 p_nplus1;
    	if (epsvc_nplus1 > epsvc0)
    	{
    		p_nplus1 = pat * pow( sqrt( p0 / pat ) + ( 1 + m_ein[k] ) / 2.0 / m_kappa[k] * epsvc_nplus1, 2);
    	}
    	else
    	{
    	  	p_nplus1 = pmin;
    	}
    	real64 K,G;
    	K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( ( p_n + p_nplus1 ) / 2.0 / pat);
    	G = m_G0[k] * pat * ( pow( 2.97 - m_ein[k] , 2 ) / ( 1 + m_ein[k]))
    		* sqrt( ( p_n + p_nplus1 ) / 2.0 / pat);
    	LvArray::tensorOps::copy< 6 >( dev_stress, dev_stress_n);
    	LvArray::tensorOps::copy< 6 >( temp, dev_strain);
    	LvArray::tensorOps::scale< 6 >(temp, 2.0 * G);
    	LvArray::tensorOps::add< 6 >( dev_stress, temp );
    	LvArray::tensorOps::copy< 6 >( temp, dev_strain_n);
    	LvArray::tensorOps::scale< 6 >(temp, -2.0 * G);
    	LvArray::tensorOps::add< 6 >( dev_stress, temp );
    	LvArray::tensorOps::copy< 6 >( r_nplus1, dev_stress);
    	LvArray::tensorOps::scale< 6 >(r_nplus1, 1.0 / p_nplus1);
    	LvArray::tensorOps::copy< 6 >(temp, r_nplus1);
    	LvArray::tensorOps::subtract< 6 >( temp, r );

    	sub1=(int)( root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(temp, temp) ) / 0.05 ) + 1;
    	LvArray::tensorOps::copy< 6 >(temp, dev_strain );
    	LvArray::tensorOps::subtract< 6 >( temp, dev_strain_n );
    	sub2 =(int)( sqrt( two3 *  LvArray::tensorOps::symDoubleContraction< 3 >(temp, temp) ) / 0.001 ) + 1;
    	sub = sub1;
    	if (sub2 > sub1)
    	{
    		sub = sub2;
    	}
    	if (sub > 100)
    	{
    	  	sub = 100;
    	}
    	LvArray::tensorOps::copy< 6 >( alpha_ns, alpha_n);
    	epsvir_ns = epsvir_n;
    	epsvre_ns = epsvre_n;
    	epsvc_ns = epsvc_n;

    	// --------------Trial Substep End-------------------------------------

    	// --------------(II)Elastic Predict-------------------------------------
    	real64 gammamonos = m_gammamono[k][q];
    	real64 eta_n;
    	for (isub = 0; isub < sub; isub++)
    	{
    		// --------------(I)Initialize-------------------------------------
    		LvArray::tensorOps::copy< 6 >( alpha_nplus1, alpha_ns);
    		epsvir_nplus1 = epsvir_ns;
    		epsvre_nplus1 = epsvre_ns;
    		epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub;
    		LvArray::tensorOps::copy< 6 >( r, dev_stress_n);
    		LvArray::tensorOps::scale< 6 >(r, 1.0 / p_n);

    		for (int i = 0; i < 6; ++i)
    		{
    			ddev_strain_p[i] = 0.0;
    		}
    		dtracep = 0.0;
    		lambda = 0.0;
    		p0 = p_n;
    		epsvc0 = - 2 * m_kappa[k] / ( 1 + m_ein[k] ) * ( sqrt( p0 / pat ) - sqrt( pmin / pat ));
    		if (epsvc_nplus1 < epsvc0)
    		{
    			p_nplus1 = pmin;
    		}
    		else
    		{
    			p_nplus1 = pat * pow( sqrt( p0 / pat ) + ( 1 + m_ein[k] ) / 2.0 / m_kappa[k] * epsvc_nplus1, 2);
    		}
    		K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( ( p_n + p_nplus1 ) / 2.0 / pat);
    		G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k]))
    				* sqrt( ( p_n + p_nplus1 ) / 2.0 / pat);
    		LvArray::tensorOps::copy< 6 >( dev_stress, dev_stress_n);
    		LvArray::tensorOps::copy< 6 >(temp, dev_strain);
    		LvArray::tensorOps::scale< 6 >( temp, 2.0 * G / sub );
    		LvArray::tensorOps::add< 6 >( dev_stress, temp );
    		LvArray::tensorOps::copy< 6 >(temp, dev_strain_n);
    		LvArray::tensorOps::scale< 6 >( temp, -2.0 * G / sub );
    		LvArray::tensorOps::add< 6 >( dev_stress, temp );
    		LvArray::tensorOps::copy< 6 >(r_nplus1, dev_stress);
    		LvArray::tensorOps::scale< 6 >( r_nplus1, 1.0 / p_nplus1 );

    		eta_n = root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r, r) );
    		if ( LvArray::tensorOps::symDoubleContraction< 3 >(r, r) < tolerance )
    		{
    		  	r1[0] = -1.0 / sqrt(6.0);
    		  	r1[1] = -1.0 / sqrt(6.0);
    		  	r1[2] = 2.0 / sqrt(6.0);
    		  	r1[3] = 0;
    		  	r1[4] = 0;
    		  	r1[5] = 0;
    		}
    		else
    		{
    			LvArray::tensorOps::copy< 6 >(r1, r);
    			LvArray::tensorOps::scale< 6 >( r1, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r, r) ) );
    		}
    		LvArray::tensorOps::symAijBjk< 3 >( temp, r1, r1 );
    		LvArray::tensorOps::symAijBjk< 3 >( pass, temp, r1 );
    		sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );
    		if ( sin3theta > 1.0 )
    		{
    			sin3theta = 1.0;
    		}
    	    else if ( sin3theta < -1.0 )
    	    {
    	    	sin3theta = -1.0;
    	    }
    		gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    				/ Mfo * ( 1 - sin3theta * sin3theta ) );
    		if (eta_n / gtheta > etamplus1)
    		{
    			etamplus1 = eta_n / gtheta;
    		}
    		if (eta_n / gtheta > Mfc * exp( - m_np[k] * psi ) - tolerance)
    		{
    			etamplus1 = eta_n / gtheta;
    		}
    		beta0 = 0.0;
    		beta1 = 1.0;
    		LvArray::tensorOps::copy< 6 >( rbar0, alpha_ns );
    		LvArray::tensorOps::copy< 6 >( temp, r );
    		LvArray::tensorOps::scale< 6 >( temp, beta0 );
    		LvArray::tensorOps::add< 6 >( rbar0, temp );
    		LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    		LvArray::tensorOps::scale< 6 >( temp, -beta0 );
    		LvArray::tensorOps::add< 6 >( rbar0, temp );

    		LvArray::tensorOps::copy< 6 >( rbar1, alpha_ns );
    		LvArray::tensorOps::copy< 6 >( temp, r );
    		LvArray::tensorOps::scale< 6 >( temp, beta1 );
    		LvArray::tensorOps::add< 6 >( rbar1, temp );
    		LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    		LvArray::tensorOps::scale< 6 >( temp, -beta1 );
    		LvArray::tensorOps::add< 6 >( rbar1, temp );

    		LvArray::tensorOps::copy< 6 >( temp, r );
    		LvArray::tensorOps::subtract< 6 >( temp, alpha_ns );

    		if (LvArray::tensorOps::symDoubleContraction< 3 >(r, r) < tolerance &&
    				LvArray::tensorOps::symDoubleContraction< 3 >(alpha_ns, alpha_ns) < tolerance)
    		{
    			normal[0] = 2.0 / sqrt(6.0);
    			normal[1] = -1.0 / sqrt(6.0);
    			normal[2] = -1.0 / sqrt(6.0);
    			normal[3] = 0;
    			normal[4] = 0;
    			normal[5] = 0;
    			LvArray::tensorOps::copy< 6 >( rbar, normal );
    			LvArray::tensorOps::scale< 6 >( rbar, root23 * Mfc * exp( -m_np[k] * psi ) );
    			beta = 1.0e20;
    		}
    		else if (sqrt(LvArray::tensorOps::symDoubleContraction< 3 >(temp, temp)) < tolerance)
    		{
    			LvArray::tensorOps::copy< 6 >( normal, r1 );
    			LvArray::tensorOps::copy< 6 >( rbar, normal );
    			LvArray::tensorOps::scale< 6 >( rbar, root23 * etamplus1 * sin3theta );
    			beta = 1.0e20;
    		}
    		else
    		{
    			if ( LvArray::tensorOps::symDoubleContraction< 3 >(rbar0, rbar0) < tolerance )
    			{
    				beta0 = 0.01;
    				LvArray::tensorOps::copy< 6 >( rbar0, alpha_ns );
    				LvArray::tensorOps::copy< 6 >( temp, r );
    				LvArray::tensorOps::scale< 6 >( temp, beta0 );
    				LvArray::tensorOps::add< 6 >( rbar0, temp );
    				LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    				LvArray::tensorOps::scale< 6 >( temp, -beta0 );
    				LvArray::tensorOps::add< 6 >( rbar0, temp );
    			}
    			LvArray::tensorOps::copy< 6 >( normal, rbar0 );
    			LvArray::tensorOps::scale< 6 >( normal, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar0, rbar0) ) );
    			LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    			LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    			sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );
    			if (sin3theta > 1.0)
    				sin3theta = 1.0;
    			else if (sin3theta < -1.0)
    				sin3theta = -1.0;
    			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    					/ Mfo * ( 1 - sin3theta * sin3theta ) );
    			LvArray::tensorOps::copy< 6 >( temp, normal );
    			LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    			LvArray::tensorOps::subtract< 6 >( temp, rbar0 );
    			Fb0 = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    			LvArray::tensorOps::copy< 6 >( normal, rbar1);
    			LvArray::tensorOps::scale< 6 >( normal, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar1, rbar1) ) );
    			LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    			LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    			sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );
    			if (sin3theta > 1.0)
    				sin3theta = 1.0;
    			else if (sin3theta < -1.0)
    				sin3theta = -1.0;
    			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    					/ Mfo * ( 1 - sin3theta * sin3theta ) );

    			LvArray::tensorOps::copy< 6 >( temp, normal );
    			LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    			LvArray::tensorOps::subtract< 6 >( temp, rbar1 );
    			Fb1 = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    			if (fabs(Fb0) <= 1.0e-5)
    			{
    				LvArray::tensorOps::copy< 6 >( rbar, rbar0 );
    				beta = beta0;
    			}
    			else if (fabs(Fb1) <= 1.0e-5)
    			{
    				LvArray::tensorOps::copy< 6 >( rbar, rbar1 );
    				beta = beta1;
    			}
    			else
    			{
    				while (Fb0 * Fb1 > 0)
    				{
    					beta0 = beta1;
    					beta1 = 2 * beta1;
    					LvArray::tensorOps::copy< 6 >( rbar0, alpha_ns );
    					LvArray::tensorOps::copy< 6 >( temp, r );
    					LvArray::tensorOps::scale< 6 >( temp, beta0 );
    					LvArray::tensorOps::add< 6 >( rbar0, temp );
    					LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    					LvArray::tensorOps::scale< 6 >( temp, -beta0 );
    					LvArray::tensorOps::add< 6 >( rbar0, temp );

    					LvArray::tensorOps::copy< 6 >( rbar1, alpha_ns );
    					LvArray::tensorOps::copy< 6 >( temp, r );
    					LvArray::tensorOps::scale< 6 >( temp, beta1 );
    					LvArray::tensorOps::add< 6 >( rbar1, temp );
    					LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    					LvArray::tensorOps::scale< 6 >( temp, -beta1 );
    					LvArray::tensorOps::add< 6 >( rbar1, temp );

    					LvArray::tensorOps::copy< 6 >( normal, rbar0);
    					LvArray::tensorOps::scale< 6 >( normal, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar0, rbar0) ) );

    					LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    					LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    					sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );
    					if (sin3theta > 1.0)
    						sin3theta = 1.0;
    					else if (sin3theta < -1.0)
    						sin3theta = -1.0;
    			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );
    					LvArray::tensorOps::copy< 6 >( temp, normal );
    					LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    					LvArray::tensorOps::subtract< 6 >( temp, rbar0 );
    					Fb0 = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    					if(LvArray::tensorOps::symDoubleContraction< 3 >(rbar1, rbar1) > 1.0e20)
    					{
    						std::cout<<beta0;
    						std::cout<<"\n";
    					}
    					LvArray::tensorOps::copy< 6 >( normal, rbar1);
    					LvArray::tensorOps::scale< 6 >( normal, 1.0 /sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar1, rbar1)));

    					LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    					LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    					sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );
    					if (sin3theta > 1.0)
    						sin3theta = 1.0;
    					else if (sin3theta < -1.0)
    						sin3theta = -1.0;
    			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );

    			        LvArray::tensorOps::copy< 6 >( temp, normal );
    			        LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    			        LvArray::tensorOps::subtract< 6 >( temp, rbar1 );
    			        Fb1 = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );

    			   	}
    				if (fabs(Fb0) <= 1.0e-5)
    			    {
    					LvArray::tensorOps::copy< 6 >( rbar, rbar0 );
    			    	beta = beta0;
    			    }
    			    else if (fabs(Fb1) <= 1.0e-5)
    			    {
    			    	LvArray::tensorOps::copy< 6 >( rbar, rbar1 );
    			    	beta = beta1;
    			    }
    				else
    				{
    				    beta = beta1 - Fb1 * ( beta1 - beta0 ) / ( Fb1 - Fb0 );
    				    LvArray::tensorOps::copy< 6 >( rbar, alpha_ns );
    					LvArray::tensorOps::copy< 6 >( temp, r );
    					LvArray::tensorOps::scale< 6 >( temp, beta );
    					LvArray::tensorOps::add< 6 >( rbar, temp );
    					LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    					LvArray::tensorOps::scale< 6 >( temp, -beta );
    					LvArray::tensorOps::add< 6 >( rbar, temp );
    					LvArray::tensorOps::copy< 6 >( normal, rbar);
    					LvArray::tensorOps::scale< 6 >( normal, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar, rbar) ) );

    					LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    					LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    					sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );

    					if (sin3theta > 1.0)
    						sin3theta = 1.0;
    					else if (sin3theta < -1.0)
    						sin3theta = -1.0;
    			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );

    			        LvArray::tensorOps::copy< 6 >( temp, normal );
    			        LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    			        LvArray::tensorOps::subtract< 6 >( temp, rbar );
    			        Fb = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    					intm = 1;
    				    while (fabs(Fb) > 1.0e-6)
    				    {
    						if (Fb * Fb1 < 0)
    						{
    							beta0 = beta1;
    							Fb0 = Fb1;
    							beta1 = beta;
    							Fb1 = Fb;
    						}
    						else
    						{
    							Fb0 = Fb1 * Fb0 / ( Fb1 + Fb );
    							beta1 = beta;
    							Fb1 = Fb;
    						}
    						beta = beta1 - Fb1 * ( beta1 - beta0 ) / ( Fb1 - Fb0 );
    						LvArray::tensorOps::copy< 6 >( rbar, alpha_ns );
    						LvArray::tensorOps::copy< 6 >( temp, r );
    						LvArray::tensorOps::scale< 6 >( temp, beta );
    						LvArray::tensorOps::add< 6 >( rbar, temp );
    						LvArray::tensorOps::copy< 6 >( temp, alpha_ns );
    						LvArray::tensorOps::scale< 6 >( temp, -beta );
    						LvArray::tensorOps::add< 6 >( rbar, temp );
    						LvArray::tensorOps::copy< 6 >( normal, rbar);
    						LvArray::tensorOps::scale< 6 >( normal, 1.0 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(rbar, rbar) ) );

    						LvArray::tensorOps::symAijBjk< 3 >( temp, normal, normal );
    						LvArray::tensorOps::symAijBjk< 3 >( pass, temp, normal );
    						sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );

    						if (sin3theta > 1.0)
    							sin3theta = 1.0;
    						else if (sin3theta < -1.0)
    							sin3theta = -1.0;
    						gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    								/ Mfo * ( 1 - sin3theta * sin3theta ) );
    				        LvArray::tensorOps::copy< 6 >( temp, normal );
    				        LvArray::tensorOps::scale< 6 >( temp, root23 * etamplus1 * gtheta );
    				        LvArray::tensorOps::subtract< 6 >( temp, rbar );
    				        Fb = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    						intm = intm + 1;
    				    }
    				}
    			}
    		}
    		LvArray::tensorOps::scale< 6 >( normal, root32 );
    		N = LvArray::tensorOps::symDoubleContraction< 3 >( r,normal );
    		// --------------(III)Loading/Unloading-------------------------------------
    		LvArray::tensorOps::copy< 6 >( temp, dev_stress );
    		LvArray::tensorOps::subtract< 6 >( temp, dev_stress_n );
    		phi = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal ) - ( p_nplus1 - p_n ) * N;

    		LvArray::tensorOps::copy< 6 >( temp, r_nplus1 );
    		LvArray::tensorOps::subtract< 6 >( temp, r );
    		phi_n = LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );


    		// --------------(IV)Unloading-------------------------------------
    		if (phi < tolerance || phi_n < tolerance)
    		{
    			gammamonos = 0.0;
    			LvArray::tensorOps::copy< 6 >( alpha_nplus1, r );
    			//epsvirpr=epsvir_n;
    		}
    		// --------------(V)Loading-------------------------------------
    		else if (phi > tolerance && phi_n > tolerance)
    		{
    			epsvc_nplus1 = epsvc_ns + ( trace - trace_n ) / sub - dtracep;
    			loadindex = 0.0;
    			lambda = 0.0;
    			lambdamin = 0.0;
    			lambdamax = 0.0;
    			LvArray::tensorOps::copy< 6 >( temp, r );
    			LvArray::tensorOps::subtract< 6 >( temp, alpha_ns );
    			rou = root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(temp, temp) );

    			LvArray::tensorOps::copy< 6 >( temp, rbar );
    			LvArray::tensorOps::subtract< 6 >( temp, alpha_ns );
    			roubar = root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(temp, temp) );

    			if (roubar > tolerance)
    			{
    				H = two3 * m_h[k] * G * gtheta * exp( -m_np[k] * psi ) * ( Mfc * exp( -m_np[k] * psi )
    				/ ( etamplus1 + tolerance ) * roubar / ( rou + tolerance ) - 1.0 );
    				if (H < tolerance && H >= 0)
    				{
    					H = tolerance;
    				}
    				if (H > -tolerance && H < 0)
    				{
    					H = -tolerance;
    				}
    				eta_nplus1 = root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1) );
    				LvArray::tensorOps::copy< 6 >( rd, rbar );
    				LvArray::tensorOps::scale< 6 >( rd,  Mdc * exp( m_nd[k] * psi ) / (etamplus1 + tolerance ) );
    				LvArray::tensorOps::copy< 6 >( temp, rd );
    				LvArray::tensorOps::subtract< 6 >( temp, r );
    				Dre_n = m_dre1[k] * root23 * LvArray::tensorOps::symDoubleContraction< 3 >( temp, normal );
    				if (epsvir_ns > tolerance)
    					chi = - m_dir[k] * epsvre_ns / epsvir_ns;
    				else
    					chi = 0.0;
    				if (chi > 1.0)
    					chi = 1.0;
    				if (Dre_n > 0.0)
    				{
    					Dre_n = pow( -m_dre2[k] * chi, 2) / p_n;
    					if (-epsvre_ns < tolerance)
    						Dre_n = 0.0;
    				}
    				if (Dre_n > 0.0)
    				{
    					if (psi >= 0)
    					{
    						LvArray::tensorOps::copy< 6 >( temp, rd );
    						LvArray::tensorOps::subtract< 6 >( temp, r );
    						Dir_n = m_dir[k] * exp( m_nd[k] * psi- m_eta[k] * epsvir_ns )
    								* ( root23 * LvArray::tensorOps::symDoubleContraction< 3 >( temp , normal ) ) * exp(chi);
    					}
    					else
    					{
    						LvArray::tensorOps::copy< 6 >( temp, rd );
    						LvArray::tensorOps::subtract< 6 >( temp, r );
    						Dir_n = m_dir[k] * exp ( m_nd[k] * psi- m_eta[k] * epsvir_ns ) * ( root23 * LvArray::tensorOps::symDoubleContraction< 3 >( temp , normal )* exp(chi)
    								+ pow( m_rdr[k] * ( 1 - exp( m_nd[k] * psi ) ) / ( m_rdr[k] * ( 1 - exp( m_nd[k] * psi ) ) + gammamonos ) , 2 ) );
    					}
    				}
    				else
    				{
    					if (psi >= 0)
    					{
    						Dir_n = 0.0;
    					}
    					else
    					{
    						Dir_n = m_dir[k] * exp ( m_nd[k] * psi- m_eta[k] * epsvir_ns ) * ( pow( m_rdr[k] * ( 1
    								- exp( m_nd[k] * psi ) ) / ( m_rdr[k] * ( 1 - exp( m_nd[k] * psi ) ) + gammamonos ) , 2 ) );
    					}
    				}
    				D_ = Dir_n + Dre_n;
                    //prem=p_n;
    				int wr = 1;
    				iconv = 0.0;
    				do
    				{
    					if ( fabs(iconv) < tolerance )//original: iconv==0.0
    					{
    						dlambda = phi / fabs( H + 2 * G - K * D_ * N);
    						lambda += dlambda;
    					}
    					else
    					{
    						lambda = 0.5 * ( lambdamax + lambdamin );
    					}
    					loadindex = H * lambda;
    					dtracep = lambda * D_;
    					LvArray::tensorOps::copy< 6 >( ddev_strain_p, normal );
    					LvArray::tensorOps::scale< 6 >( ddev_strain_p, lambda );
    					epsvc_nplus1 = epsvc_ns + ( trace - trace_n ) / sub - dtracep;
    					if (epsvc_nplus1 < epsvc0)
    					{
    						p_nplus1 = pmin;
    						epsvc_nplus1 = epsvc0;
    					}
    					else
    					{
    						p_nplus1 = pat * pow( sqrt( p0 / pat )+( 1 + m_ein[k] ) / 2.0
    								/ m_kappa[k] * epsvc_nplus1 , 2);
    					}
    					G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ), 2 ) / ( 1 + m_ein[k] ) )
    							* sqrt( ( p_n + p_nplus1 ) / 2.0 / pat );

    					LvArray::tensorOps::copy< 6 >( dev_stress, dev_stress_n );
    					LvArray::tensorOps::copy< 6 >( temp, dev_strain );
    					LvArray::tensorOps::scale< 6 >( temp, 2 * G / sub );
    					LvArray::tensorOps::add< 6 >( dev_stress, temp );
    					LvArray::tensorOps::copy< 6 >( temp, dev_strain_n );
    					LvArray::tensorOps::scale< 6 >( temp, -2 * G / sub );
    					LvArray::tensorOps::add< 6 >( dev_stress, temp );
    					LvArray::tensorOps::copy< 6 >( temp, ddev_strain_p );
    					LvArray::tensorOps::scale< 6 >( temp, -2 * G );
    					LvArray::tensorOps::add< 6 >( dev_stress, temp );

    					LvArray::tensorOps::copy< 6 >( temp, dev_stress );
    					LvArray::tensorOps::subtract< 6 >( temp, dev_stress_n );

    					phi = LvArray::tensorOps::symDoubleContraction< 3 >( temp , normal ) - ( p_nplus1 - p_n ) * N - loadindex;
    					//phi=doublecontraction(dev_stress-dev_stress_n,normal) / root32 -(p_nplus1-p_n)*N / root32-loadindex;
    					if (phi < -tolerance)
    					{
    						iconv = 1.0;
    						lambdamax = lambda;
    					}
    					if (phi > tolerance && fabs( iconv - 1.0 ) < tolerance)//original: iconv==1.0
    						lambdamin = lambda;
    					wr = wr + 1;
    					epsvir_nplus1 = lambda * Dir_n + epsvir_ns;
    					epsvre_nplus1 = lambda* Dre_n + epsvre_ns;
    				}while (fabs(phi) > tolerance);
    				gammamonos = gammamonos + lambda;
    			}
    		}
    		LvArray::tensorOps::copy< 6 >( r_nplus1, dev_stress );
    		LvArray::tensorOps::scale< 6 >( r_nplus1, 1.0 / p_nplus1 );

    		eta_nplus1 = root32 * sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1) );

    		if (eta_nplus1 >= Mfc * exp( -m_np[k] * psi ) / ( 1.0 + Mfc / 3.0 ) - tolerance)
    		{
    			LvArray::tensorOps::copy< 6 >( r1, r_nplus1 );
    			LvArray::tensorOps::scale< 6 >( r1, 1 / sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1) ) );
    			LvArray::tensorOps::symAijBjk< 3 >( temp, r1, r1 );
    			LvArray::tensorOps::symAijBjk< 3 >( pass, temp, r1 );
    			sin3theta = -sqrt(6.0) * LvArray::tensorOps::symTrace< 3 >( pass );

    			if (sin3theta > 1.0)
    				sin3theta = 1.0;
    			else if (sin3theta < -1.0)
    				sin3theta = -1.0;
    			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
    					/ Mfo * ( 1 - sin3theta * sin3theta ) );
    			LvArray::tensorOps::scale< 6 >( r1, root23 * Mfc * exp( -m_np[k] * psi ) * gtheta );

    			if (LvArray::tensorOps::symDoubleContraction< 3 >(r1, r1) -LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1) < tolerance)
    			{
    				intm = sqrt( LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1) / LvArray::tensorOps::symDoubleContraction< 3 >(r1, r1) ) + tolerance;
    				LvArray::tensorOps::scale< 6 >( dev_stress, 1.0 / intm );
    			}
    			LvArray::tensorOps::copy< 6 >( r_nplus1, dev_stress );
    			LvArray::tensorOps::scale< 6 >( r_nplus1, 1.0 / p_nplus1 );
    			eta_nplus1 = root32 * sqrt(LvArray::tensorOps::symDoubleContraction< 3 >(r_nplus1, r_nplus1));
    		}
    		LvArray::tensorOps::copy< 6 >( alpha_ns, alpha_nplus1 );
    		epsvir_ns = epsvir_nplus1;
    		epsvre_ns = epsvre_nplus1;
            epsvc_ns = epsvc_ns - epsvc_nplus1 + ( trace - trace_n ) / sub - dtracep;
    		epsvc_nplus1 = epsvc_ns;

    		p_n = p_nplus1;
    		LvArray::tensorOps::copy< 6 >( dev_stress_n, dev_stress );

    	}
    	LvArray::tensorOps::copy< 6 >( stress_pass, dev_stress );
    	LvArray::tensorOps::symAddIdentity< 3 >( stress_pass, p_n );
    	LvArray::tensorOps::copy< 6 >( temp, stress_pass );
    	LvArray::tensorOps::scale< 6 >( temp, -1000 );
    	LvArray::tensorOps::subtract< 6 >( temp, m_oldStress[k][q] );// from positive to negative
    	LvArray::tensorOps::copy< 6 >( stress, temp );

    	LvArray::tensorOps::copy< 6 >( m_strain[ k ][ q ], strain_nplus1 );
    	LvArray::tensorOps::scale< 6 >( m_strain[ k ][ q ], -1 );
    	m_epsvir[k][q] = epsvir_ns;
    	m_epsvre[k][q] = epsvre_ns;
    	m_gammamono[k][q] = gammamonos;
    	m_epsvc[k][q] = epsvc_ns;
    	m_etam[k][q] = etamplus1 ;
    	LvArray::tensorOps::copy< 6 >( m_alpha[k][q], alpha_ns );
	}
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                        localIndex const q,
                                                        real64 const ( &totalStrain )[6],
                                                        real64 ( & stress )[6],
                                                        real64 ( & stiffness )[6][6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  getElasticStiffness( k, q, stiffness );
}


/// To be corrected by ron
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                        localIndex const q,
                                                        real64 const ( &totalStrain )[6],
                                                        real64 ( & stress )[6],
                                                        DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  real64 const p = 1e3;
  real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat )*1000.0;
  real64 const K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat )*1000.0;
  G = 0.3835 * K;
  stiffness.m_bulkModulus = K;
  stiffness.m_shearModulus = G;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                            localIndex const q,
                                                            real64 const ( &strainIncrement )[6],
                                                            real64 ( & stress )[6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress  = incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 real64 ( & stiffness )[6][6] ) const
{
  smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  getElasticStiffness( k, q, stiffness );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CycLiqCPSPUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  real64 const p = 1e3;
  real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat )*1000.0;
  real64 const K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat )*1000.0;
  G = 0.3835 * K;
  stiffness.m_bulkModulus = K;
  stiffness.m_shearModulus = G;
}

class CycLiqCPSP : public SolidBase
{
public:

  /// Alias for CycLiqCPSPUpdates
  using KernelWrapper = CycLiqCPSPUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CycLiqCPSP( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CycLiqCPSP() override;


  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CycLiqCPSP";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  /// Keys for data specified in this class.
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default parameters
    static constexpr char const * defaultG0String() { return "defaultG0"; }
    static constexpr char const * defaultKappaString() { return "defaultKappa"; }
    static constexpr char const * defaultHString() { return "defaultH"; }
    static constexpr char const * defaultMString() { return "defaultM"; }
    static constexpr char const * defaultDre1String() { return "defaultDre1"; }
    static constexpr char const * defaultDre2String() { return "defaultDre2"; }
    static constexpr char const * defaultDirString() { return "defaultDir"; }
    static constexpr char const * defaultEtaString() { return "defaultEta"; }
    static constexpr char const * defaultRdrString() { return "defaultRdr"; }
    static constexpr char const * defaultNpString() { return "defaultNp"; }
    static constexpr char const * defaultNdString() { return "defaultNd"; }
    static constexpr char const * defaultLamdacString() { return "defaultLamdac"; }
    static constexpr char const * defaultE0String() { return "defaultE0"; }
    static constexpr char const * defaultKsiString() { return "defaultKsi"; }
    static constexpr char const * defaultEinString() { return "defaultEin"; }
    static constexpr char const * defaultInitialCyclesString() { return "defaultInitialCycles"; }

    /// string/key for parameters
    static constexpr char const * G0String() { return "G0"; }
    static constexpr char const * kappaString() { return "kappa"; }
    static constexpr char const * hString() { return "h"; }
    static constexpr char const * MString() { return "M"; }
    static constexpr char const * dre1String() { return "dre1"; }
    static constexpr char const * dre2String() { return "dre2"; }
    static constexpr char const * dirString() { return "dir"; }
    static constexpr char const * etaString() { return "eta"; }
    static constexpr char const * rdrString() { return "rdr"; }
    static constexpr char const * npString() { return "np"; }
    static constexpr char const * ndString() { return "nd"; }
    static constexpr char const * lamdacString() { return "lamdac"; }
    static constexpr char const * e0String() { return "e0"; }
    static constexpr char const * ksiString() { return "ksi"; }
    static constexpr char const * einString() { return "ein"; }

    /// string/key for state parameters
    static constexpr char const * strainString() { return "strain"; }
    static constexpr char const * epsvirString() { return "epsvir"; }
    static constexpr char const * epsvreString() { return "epsvre"; }
    static constexpr char const * gammamonoString() { return "gammamono"; }
    static constexpr char const * epsvcString() { return "epsvc"; }
    static constexpr char const * etamString() { return "etam"; }
    static constexpr char const * alphaString() { return "alpha"; }
    static constexpr char const * initialCyclesString() { return "initialCycles"; }
  };

  /**
   * @brief Create a instantiation of the
   *        CycLiqCPSPUpdates class that refers to the
   *        data in this.
   * @return An instantiation of CycLiqCPSPUpdates.
   */
  CycLiqCPSPUpdates createKernelUpdates() const
  {
    return CycLiqCPSPUpdates( m_G0,
    		m_kappa,
			m_h,
    		m_dre1,
    		m_dre2,
    		m_dir,
    		m_eta,
			m_rdr,
    		m_np,
    		m_nd,
    		m_M,
    		m_lamdac,
    		m_e0,
    		m_ksi,
    		m_ein,
    		m_initialCycles,
    		m_epsvir,
    		m_epsvre,
    	    m_gammamono,
    		m_epsvc,
    		m_etam,
    		m_alpha,
    		m_strain,
			m_newStress,
			m_oldStress );
  }


  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
    					m_G0,
    		    		m_kappa,
    					m_h,
    		    		m_dre1,
    		    		m_dre2,
    		    		m_dir,
    		    		m_eta,
    					m_rdr,
    		    		m_np,
    		    		m_nd,
    		    		m_M,
    		    		m_lamdac,
    		    		m_e0,
    		    		m_ksi,
    		    		m_ein,
    		    		m_initialCycles,
    		    		m_epsvir,
    		    		m_epsvre,
    		    	    m_gammamono,
    		    		m_epsvc,
    		    		m_etam,
    		    		m_alpha,
    		    		m_strain,
                          m_newStress,
                          m_oldStress );
  }

  protected:

    /// Post-process XML data
    virtual void postProcessInput() override;

    /// The default value of the [CycLiqCPSP parameters] elastic shear modulus for any new allocations.
    real64 m_defaultG0;

    /// The default value of the [CycLiqCPSP parameters] elastic bulk modulus for any new allocations.
    real64 m_defaultKappa;

    /// The default value of the [CycLiqCPSP parameters] plastic modulus for any new allocations.
    real64 m_defaultH;

    /// The default value of the [CycLiqCPSP parameters] reversible dilatancy generation for any new allocations.
    real64 m_defaultDre1;

    /// The default value of the [CycLiqCPSP parameters] reversible dilatancy release for any new allocations.
    real64 m_defaultDre2;

    /// The default value of the [CycLiqCPSP parameters] irreversible dilatancy for any new allocations.
    real64 m_defaultDir;

    /// The default value of the [CycLiqCPSP parameters] irreversible dilatancy rate decrease for any new allocations.
    real64 m_defaultEta;

    /// The default value of the [CycLiqCPSP parameters] reference shear strain for any new allocations.
    real64 m_defaultRdr;

    /// The default value of the [CycLiqCPSP parameters] plasticity state dependency for any new allocations.
    real64 m_defaultNp;

    /// The default value of the [CycLiqCPSP parameters] dilatancy state dependency for any new allocations.
    real64 m_defaultNd;

    /// The default value of the [CycLiqCPSP parameters] critical state (CS) stress ratio for any new allocations.
    real64 m_defaultM;

    /// The default value of the [CycLiqCPSP parameters] CS in e-p space for any new allocations.
    real64 m_defaultLamdac;

    /// The default value of the [CycLiqCPSP parameters] CS in e-p space for any new allocations.
    real64 m_defaultE0;

    /// The default value of the [CycLiqCPSP parameters] CS in e-p space for any new allocations.
    real64 m_defaultKsi;

    /// The default value of the initial void ratio for any new allocations.
    real64 m_defaultEin;

    /// The default value of the Cycles of elastic stage for any new allocations.
    integer m_defaultInitialCycles;

    /// The [CycLiqCPSP parameters] elastic shear modulus for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_G0;

    /// The [CycLiqCPSP parameters] elastic bulk modulus for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_kappa;

    /// The [CycLiqCPSP parameters] plastic modulus for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_h;

    /// The [CycLiqCPSP parameters] reversible dilatancy generation for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_dre1;

    /// The [CycLiqCPSP parameters] reversible dilatancy release for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_dre2;

    /// The [CycLiqCPSP parameters] irreversible dilatancy for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_dir;

    /// The [CycLiqCPSP parameters] irreversible dilatancy rate decrease for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_eta;

    /// The [CycLiqCPSP parameters] reference shear strain for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_rdr;

    /// The [CycLiqCPSP parameters] plasticity state dependency for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_np;

    /// The [CycLiqCPSP parameters] dilatancy state dependency for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_nd;

    /// The [CycLiqCPSP parameters] critical state (CS) stress ratio for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_M;

    /// The [CycLiqCPSP parameters] CS in e-p space for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_lamdac;

    /// The [CycLiqCPSP parameters] CS in e-p space for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_e0;

    /// The [CycLiqCPSP parameters] CS in e-p space for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_ksi;

    /// The initial void ratio for each upper level dimension (i.e. cell) of *this
    array1d< real64 > m_ein;

    /// The the Cycles of elastic stage at a quadrature point.
    array2d< integer > m_initialCycles;

    /// The material stress at a quadrature point.
    array3d< real64, solid::STRESS_PERMUTATION > m_strain;

    array2d< real64 > m_epsvir;

    array2d< real64 > m_epsvre;

    /// The shear strain generated after the previous stress reversal at a quadrature point.
    array2d< real64 > m_gammamono;

    array2d< real64 > m_epsvc;

    array2d< real64 > m_etam;

    array3d< real64, solid::STRESS_PERMUTATION > m_alpha;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_CYCLIQCPSP_HPP_ */
