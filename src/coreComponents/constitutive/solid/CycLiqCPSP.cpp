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
 *  @file CycLiqCPSP.cpp
 *
 *  Created on: Feb 27, 2020
 *      Author: ron
 */

#include "CycLiqCPSP.hpp"
#include <math.h>

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{

CycLiqCPSP::CycLiqCPSP( std::string const & name, Group * const parent ):
		  SolidBase( name, parent ),
		  m_defaultG0(),
		  m_defaultKappa(),
		  m_defaultH(),
		  m_defaultM(),
		  m_defaultDre1(),
		  m_defaultDre2(),
		  m_defaultDir(),
		  m_defaultEta(),
		  m_defaultRdr(),
		  m_defaultNp(),
		  m_defaultNd(),
		  m_defaultLamdac(),
		  m_defaultE0(),
		  m_defaultKsi(),
		  m_defaultEin(),
          m_defaultInitialTime(),
		  m_strain{},
		  m_epsvir{},
		  m_epsvre{},
		  m_gammamono{},
		  m_epsvc{},
		  m_etam{},
		  m_alpha{},
		  m_initialTime{}
{
	registerWrapper( viewKeyStruct::defaultG0String, &m_defaultG0, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : G0");
	registerWrapper( viewKeyStruct::defaultKappaString, &m_defaultKappa, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : kappa");
	registerWrapper( viewKeyStruct::defaultHString, &m_defaultH, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : h");
	registerWrapper( viewKeyStruct::defaultMString, &m_defaultM, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : M");
	registerWrapper( viewKeyStruct::defaultDre1String, &m_defaultDre1, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : dre1");
	registerWrapper( viewKeyStruct::defaultDre2String, &m_defaultDre2, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : dre2");
	registerWrapper( viewKeyStruct::defaultDirString, &m_defaultDir, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : dir");
	registerWrapper( viewKeyStruct::defaultEtaString, &m_defaultEta, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : eta");
	registerWrapper( viewKeyStruct::defaultRdrString, &m_defaultRdr, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : rdr");
	registerWrapper( viewKeyStruct::defaultNpString, &m_defaultNp, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : np");
	registerWrapper( viewKeyStruct::defaultNdString, &m_defaultNd, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : nd");
	registerWrapper( viewKeyStruct::defaultLamdacString, &m_defaultLamdac, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : lamdac");
	registerWrapper( viewKeyStruct::defaultE0String, &m_defaultE0, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : e0");
	registerWrapper( viewKeyStruct::defaultKsiString, &m_defaultKsi, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : ksi");
	registerWrapper( viewKeyStruct::defaultEinString, &m_defaultEin, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : ein");
	registerWrapper( viewKeyStruct::defaultInitialTimeString, &m_defaultInitialTime, 0 )->
	    setApplyDefaultValue(-1)->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("Default for CycLiqCPSP Model Parameter : initial time");

	registerWrapper( viewKeyStruct::G0String, &m_G0, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : G0");
	registerWrapper( viewKeyStruct::kappaString, &m_kappa, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : kappa");
	registerWrapper( viewKeyStruct::hString, &m_h, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : h");
	registerWrapper( viewKeyStruct::MString, &m_M, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : M");
	registerWrapper( viewKeyStruct::dre1String, &m_dre1, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : dre1");
	registerWrapper( viewKeyStruct::dre2String, &m_dre2, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : dre2");
	registerWrapper( viewKeyStruct::dirString, &m_dir, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : dir");
	registerWrapper( viewKeyStruct::etaString, &m_eta, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : eta");
	registerWrapper( viewKeyStruct::rdrString, &m_rdr, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : rdr");
	registerWrapper( viewKeyStruct::npString, &m_np, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : np");
	registerWrapper( viewKeyStruct::ndString, &m_nd, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : nd");
	registerWrapper( viewKeyStruct::lamdacString, &m_lamdac, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : lamdac");
	registerWrapper( viewKeyStruct::e0String, &m_e0, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : e0");
	registerWrapper( viewKeyStruct::ksiString, &m_ksi, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : ksi");
	registerWrapper( viewKeyStruct::einString, &m_ein, 0 )->
	    setApplyDefaultValue(-1)->
	    setDescription("CycLiqCPSP Model Parameter : ein");

    registerWrapper( viewKeyStruct::strainString, &m_strain, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("strain");
    registerWrapper( viewKeyStruct::epsvirString, &m_epsvir, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("epsvir");
    registerWrapper( viewKeyStruct::epsvreString, &m_epsvre, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("epsvre");
    registerWrapper( viewKeyStruct::gammamonoString, &m_gammamono, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("gammamono");
    registerWrapper( viewKeyStruct::epsvcString, &m_epsvc, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("epsvc");
    registerWrapper( viewKeyStruct::etamString, &m_etam, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("etam");
    registerWrapper( viewKeyStruct::alphaString, &m_alpha, 0 )->
	    setPlotLevel(PlotLevel::LEVEL_0)->
	    setDescription("alpha");
    registerWrapper( viewKeyStruct::initialTimeString, &m_initialTime, 0 )->
  	    setPlotLevel(PlotLevel::LEVEL_0)->
  	    setDescription("initialTime");
}

CycLiqCPSP::~CycLiqCPSP()
{}

void
CycLiqCPSP::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<CycLiqCPSP>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  CycLiqCPSP * const newConstitutiveRelation = dynamic_cast<CycLiqCPSP *>(clone.get());

  newConstitutiveRelation->m_defaultG0 = m_defaultG0;
  newConstitutiveRelation->m_G0 = m_G0;
  newConstitutiveRelation->m_defaultKappa = m_defaultKappa;
  newConstitutiveRelation->m_kappa = m_kappa;
  newConstitutiveRelation->m_defaultH = m_defaultH;
  newConstitutiveRelation->m_h = m_h;
  newConstitutiveRelation->m_defaultM = m_defaultM;
  newConstitutiveRelation->m_M = m_M;
  newConstitutiveRelation->m_defaultDre1 = m_defaultDre1;
  newConstitutiveRelation->m_dre1 = m_dre1;
  newConstitutiveRelation->m_defaultDre2 = m_defaultDre2;
  newConstitutiveRelation->m_dre2 = m_dre2;
  newConstitutiveRelation->m_defaultDir = m_defaultDir;
  newConstitutiveRelation->m_dir = m_dir;
  newConstitutiveRelation->m_defaultEta = m_defaultEta;
  newConstitutiveRelation->m_eta = m_eta;
  newConstitutiveRelation->m_defaultRdr = m_defaultRdr;
  newConstitutiveRelation->m_rdr = m_rdr;
  newConstitutiveRelation->m_defaultNp = m_defaultNp;
  newConstitutiveRelation->m_np = m_np;
  newConstitutiveRelation->m_defaultNd = m_defaultNd;
  newConstitutiveRelation->m_nd = m_nd;
  newConstitutiveRelation->m_defaultLamdac = m_defaultLamdac;
  newConstitutiveRelation->m_lamdac = m_lamdac;
  newConstitutiveRelation->m_defaultE0 = m_defaultE0;
  newConstitutiveRelation->m_e0 = m_e0;
  newConstitutiveRelation->m_defaultKsi = m_defaultKsi;
  newConstitutiveRelation->m_ksi = m_ksi;
  newConstitutiveRelation->m_defaultEin = m_defaultEin;
  newConstitutiveRelation->m_ein = m_ein;
  newConstitutiveRelation->m_defaultInitialTime = m_defaultInitialTime;

  newConstitutiveRelation->m_strain = m_strain;
  newConstitutiveRelation->m_epsvir = m_epsvir;
  newConstitutiveRelation->m_epsvre = m_epsvre;
  newConstitutiveRelation->m_gammamono = m_gammamono;
  newConstitutiveRelation->m_epsvc = m_epsvc;
  newConstitutiveRelation->m_etam = m_etam;
  newConstitutiveRelation->m_alpha = m_alpha;
  newConstitutiveRelation->m_stress = m_stress;
  newConstitutiveRelation->m_initialTime = m_initialTime;
}

void CycLiqCPSP::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_G0.resize( parent->size() );
  m_kappa.resize( parent->size() );
  m_h.resize( parent->size() );
  m_M.resize( parent->size() );
  m_dre1.resize( parent->size() );
  m_dre2.resize( parent->size() );
  m_dir.resize( parent->size() );
  m_eta.resize( parent->size() );
  m_rdr.resize( parent->size() );
  m_np.resize( parent->size() );
  m_nd.resize( parent->size() );
  m_lamdac.resize( parent->size() );
  m_e0.resize( parent->size() );
  m_ksi.resize( parent->size() );
  m_ein.resize( parent->size() );

  m_strain.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_epsvir.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_epsvre.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_gammamono.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_epsvc.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_etam.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_alpha.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_initialTime.resize( parent->size(), numConstitutivePointsPerParentIndex );
  //m_strain.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_G0 = m_defaultG0;
  m_kappa = m_defaultKappa;
  m_h = m_defaultH;
  m_M = m_defaultM;
  m_dre1 = m_defaultDre1;
  m_dre2 = m_defaultDre2;
  m_dir = m_defaultDir;
  m_eta = m_defaultEta;
  m_rdr = m_defaultRdr;
  m_np = m_defaultNp;
  m_nd = m_defaultNd;
  m_lamdac = m_defaultLamdac;
  m_e0 = m_defaultE0;
  m_ksi = m_defaultKsi;
  m_ein = m_defaultEin;
  m_initialTime = -m_defaultInitialTime;
}

void CycLiqCPSP::PostProcessInput()
{
	  m_G0 = m_defaultG0;
	  m_kappa = m_defaultKappa;
	  m_h = m_defaultH;
	  m_M = m_defaultM;
	  m_dre1 = m_defaultDre1;
	  m_dre2 = m_defaultDre2;
	  m_dir = m_defaultDir;
	  m_eta = m_defaultEta;
	  m_rdr = m_defaultRdr;
	  m_np = m_defaultNp;
	  m_nd = m_defaultNd;
	  m_lamdac = m_defaultLamdac;
	  m_e0 = m_defaultE0;
	  m_ksi = m_defaultKsi;
	  m_ein = m_defaultEin;
	  m_initialTime = -m_defaultInitialTime;
}

void CycLiqCPSP::StateUpdatePoint( localIndex const k,
                                               localIndex const q,
                                               R2SymTensor const & D,
                                               R2Tensor const & Rot,
                                               real64 const dt,
                                               integer const GEOSX_UNUSED_ARG( updateStiffnessFlag ) )
{
	if(dt < std::numeric_limits<real64>::max())
	{
		m_initialTime[k][q] += dt;
	}
if(m_initialTime[k][q] < 0)
{
	m_clearDisplacement = 1;

	 //m_strain[k][q] += D;
     real64 p = 1e6;
     real64 G = m_G0[k] * pat * ( pow( ( 2.97 - m_ein[k] ) , 2 ) / ( 1 + m_ein[k])) * sqrt( p / pat );
     real64 K = (1 + m_ein[k]) / m_kappa[k] * pat * sqrt( p / pat );
     G = 3.0 / 8.0 * K;
     real64 meanStresIncrement = D.Trace();
     R2SymTensor temp = D;
     temp.PlusIdentity( -meanStresIncrement / 3.0 );
     temp *= 2.0 * G;
     meanStresIncrement *= K;
     temp.PlusIdentity( meanStresIncrement );
     m_stress[k][q] += temp;
     temp.QijAjkQlk( m_stress[k][q], Rot );
     m_stress[k][q] = temp;
}

else
{
	m_clearDisplacement = -1;

	real64 Mfc = m_M[k];
	real64 Mdc = m_M[k];
	real64 sinphi = 3.0 * Mfc / (Mfc + 6.0);
	real64 tanphi = sinphi / sqrt(1.0 - sinphi * sinphi);
	real64 Mfo = 2 * sqrt(3.0) * tanphi / sqrt(3.0 + 4.0 * tanphi * tanphi);
	real64 epsvir_n = m_epsvir[k][q];
    real64 epsvre_n = m_epsvre[k][q];
    real64 epsvc_n = m_epsvc[k][q];

    R2SymTensor dev_strain;
    R2SymTensor dev_strain_n;
    R2SymTensor ddev_strain_p;
    R2SymTensor dev_stress;
    R2SymTensor dev_stress_n;
	R2SymTensor normal;
	R2SymTensor pass;
	R2SymTensor alpha_ns;
	R2SymTensor rbar;
	R2SymTensor rbar0;
	R2SymTensor rbar1;
	R2SymTensor r1;
	R2SymTensor rd;
	R2SymTensor stress_n = m_stress[k][q];
	R2SymTensor strain_n = m_strain[k][q];
	R2SymTensor strain_nplus1 = strain_n;
	strain_nplus1 += D;
	R2SymTensor alpha_n = m_alpha[k][q];
	R2SymTensor stress_pass;//, ZeroTensor;
	R2SymTensor alpha_nplus1;
	R2SymTensor r;
	R2SymTensor r_nplus1;

	for(int i = 0; i < 6; ++i)
	{
		dev_strain.Data()[i] = 0;
	}

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
	stress_n *= (-1);
	strain_n *= (-1);
	strain_nplus1 *= (-1);

	//compute the deviatoric stress of last step
	real64 p_n = stress_n.Trace() * one3;
	dev_stress_n = stress_n;
	dev_stress_n.PlusIdentity( -p_n );
	if( p_n < pmin )
	{
		p_n = pmin;
	}

	//compute the deviatoric strains
	trace = strain_nplus1.Trace();
	trace_n = strain_n.Trace();
	dev_strain_n = strain_n;
	dev_strain_n.PlusIdentity( -trace_n * one3 );
	dev_strain = strain_nplus1;
	dev_strain.PlusIdentity( -trace * one3 );

	real64 en = (1 + m_ein[k]) * exp(-trace_n) - 1.0;//void ratio
	R2SymTensor temp;
	real64 p0;
	real64 epsvc0;

	/*
	if(!isInitialize[k][q])
	{
		p0 = p_n;
		epsvc0 = -2 * m_kappa[k] / (1 + m_ein[k]) * (sqrt( p0 / pat ) - sqrt( pmin / pat ));
		r = dev_stress_n* ( 1.0 / p_n );
		r1 = r *(1.0 * sqrt( r.doubleContraction(r, r) ));
	    temp.AijBjk( r1, r1 );
        pass.AijBjk( temp, r1 );
	    sin3theta = -sqrt(6.0) * pass.Trace();
		if (sin3theta>1.0)
			sin3theta=1.0;
		else if (sin3theta<-1.0)
			sin3theta=-1.0;
		gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
				/ Mfo * ( 1 - sin3theta * sin3theta ) );
        m_etam[k][q] = root32 * sqrt( r.doubleContraction(r, r) ) / gtheta;
    	isInitialize[k][q] = true;
		//trace_inp = trace_n;
	}
*/
	// --------------(I)Initialize-------------------------------------
	alpha_nplus1 = alpha_n;
	real64 epsvir_nplus1 = epsvir_n;
	real64 epsvre_nplus1 = epsvre_n;
	real64 epsvc_nplus1 = epsvc_n + trace - trace_n;
	real64 etamplus1 = m_etam[k][q];
	real64 lambda = 0.0;
    p0 = p_n;
    epsvc0 = - 2 * m_kappa[k] / (1 + m_ein[k]) * ( sqrt( p0 / pat) - sqrt( pmin / pat ));
	r = dev_stress_n* ( 1.0 / p_n );
	ec = m_e0[k] - m_lamdac[k] * pow( p_n / pat , m_ksi[k]);
	psi = en - ec;
	for (int i = 0; i < 6; ++i)
	{
		ddev_strain_p.Data()[i] = 0.0;
	}

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
	dev_stress = dev_stress_n;
	dev_stress += 2.0 * G * dev_strain;
	dev_stress -= 2.0 * G * dev_strain_n;

	r_nplus1 = dev_stress * ( 1.0 / p_nplus1);
	temp = r_nplus1;
	temp -= r;
	sub1=(int)( root32 * sqrt( temp.doubleContraction(temp, temp) ) / 0.05 ) + 1;
	temp = dev_strain;
	temp -= dev_strain_n;
	sub2 =(int)( sqrt( two3 * temp.doubleContraction(temp, temp) ) / 0.001 ) + 1;
	sub = sub1;
	if (sub2 > sub1)
	{
		sub = sub2;
	}
	if (sub > 100)
	{
	  	sub = 100;
	}
	alpha_ns = alpha_n;
	epsvir_ns=epsvir_n;
	epsvre_ns=epsvre_n;
	epsvc_ns=epsvc_n;
	// --------------Trial Substep End-------------------------------------

	// --------------(II)Elastic Predict-------------------------------------
	real64 gammamonos = m_gammamono[k][q];
	real64 eta_n;
	for (isub = 0; isub < sub; isub++)
	{
		// --------------(I)Initialize-------------------------------------
		alpha_nplus1 = alpha_ns;
		epsvir_nplus1 = epsvir_ns;
		epsvre_nplus1 = epsvre_ns;
		epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub;
		r = dev_stress_n * (1.0 / p_n);
		for (int i = 0; i < 6; ++i)
		{
			ddev_strain_p.Data()[i] = 0.0;
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
		dev_stress = dev_stress_n;
		dev_stress += 2.0 * G / sub * dev_strain;
		dev_stress -= 2.0 * G / sub * dev_strain_n;
		r_nplus1 = dev_stress * ( 1.0 / p_nplus1 );
		eta_n = root32 * sqrt( r.doubleContraction(r, r) );

		if ( r.doubleContraction(r, r) < tolerance )
		{
		  	r1.Data()[0] = -1.0 / sqrt(6.0);
		  	r1.Data()[2] = -1.0 / sqrt(6.0);
		  	r1.Data()[5] = 2.0 / sqrt(6.0);
		  	r1.Data()[1] = 0;
		  	r1.Data()[3] = 0;
		  	r1.Data()[4] = 0;
		}
		else
		{
			r1 = r *(1.0 * sqrt( r.doubleContraction(r, r) ));
		}
        temp.AijBjk( r1, r1 );
        pass.AijBjk( temp, r1 );
		sin3theta = -sqrt(6.0) * pass.Trace();
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
		rbar0 = alpha_ns;
		rbar0 += beta0 * r;
		rbar0 -= beta0 * alpha_ns;
		rbar1 = alpha_ns;
		rbar1 += beta1 * r;
		rbar1 -= beta1 * alpha_ns;
		temp = r;
		temp -= alpha_ns;
		//std::cout<<temp.doubleContraction(temp, temp);
		//std::cout<<"r-alpha\n";
		if (r.doubleContraction(r, r) < tolerance && alpha_ns.doubleContraction(alpha_ns, alpha_ns) < tolerance)
		{
			normal.Data()[0] = 2.0 / sqrt(6.0);
			normal.Data()[2] = -1.0 / sqrt(6.0);
			normal.Data()[5] = -1.0 / sqrt(6.0);
			normal.Data()[1] = 0;
			normal.Data()[2] = 0;
			normal.Data()[3] = 0;
			rbar = root23 * Mfc * exp( -m_np[k] * psi ) * normal;
			beta = 1.0e20;
		}
		else if (sqrt(temp.doubleContraction(temp, temp)) < tolerance)
		{
			normal = r1;
			rbar = root23 * etamplus1 * sin3theta * normal;
			beta = 1.0e20;
		}
		else
		{
			if ( rbar0.doubleContraction(rbar0, rbar0) < tolerance )
			{
				beta0 = 0.01;
				rbar0 = alpha_ns;
				rbar0 += beta0 * r;
				rbar0 -= beta0 * alpha_ns;
			}
			normal = rbar0 * ( 1.0 / sqrt( rbar0.doubleContraction(rbar0, rbar0) ));
	        temp.AijBjk( normal, normal );
	        pass.AijBjk( temp, normal );
			sin3theta = -sqrt(6.0) * pass.Trace();
			if (sin3theta > 1.0)
				sin3theta = 1.0;
			else if (sin3theta < -1.0)
				sin3theta = -1.0;
			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
					/ Mfo * ( 1 - sin3theta * sin3theta ) );
			temp = normal * root23 * etamplus1 * gtheta;
			temp -= rbar0;
			Fb0 = temp.doubleContraction(temp,normal);
			normal = rbar1 * ( 1.0 / sqrt( rbar1.doubleContraction(rbar1, rbar1) ));
	        temp.AijBjk( normal, normal );
	        pass.AijBjk( temp, normal );
			sin3theta = -sqrt(6.0) * pass.Trace();
			if (sin3theta > 1.0)
				sin3theta = 1.0;
			else if (sin3theta < -1.0)
				sin3theta = -1.0;
			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
					/ Mfo * ( 1 - sin3theta * sin3theta ) );
			temp = normal * root23 * etamplus1 * gtheta;
			temp -= rbar1;
			Fb1 = temp.doubleContraction(temp,normal);
			if (abs(Fb0) <= 1.0e-5)
			{
				rbar = rbar0;
				beta = beta0;
			}
			else if (abs(Fb1) <= 1.0e-5)
			{
				rbar = rbar1;
				beta = beta1;
			}
			else
			{
				while (Fb0 * Fb1 > 0)
				{
					//std::cout<<"here1\n";
					beta0 = beta1;
					beta1 = 2 * beta1;
					rbar0 = alpha_ns;
					rbar0 += beta0 * r;
					rbar0 -= beta0 * alpha_ns;
					rbar1 = alpha_ns;
					rbar1 += beta1 * r;
					rbar1 -= beta1 * alpha_ns;
			        normal = rbar0 * ( 1.0 / sqrt( rbar0.doubleContraction(rbar0, rbar0) ));
			        temp.AijBjk( normal, normal );
			        pass.AijBjk( temp, normal );
			        sin3theta = -sqrt(6.0) * pass.Trace();
					if (sin3theta > 1.0)
						sin3theta = 1.0;
					else if (sin3theta < -1.0)
						sin3theta = -1.0;
			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );
					temp = normal * root23 * etamplus1 * gtheta;
					temp -= rbar0;
					Fb0 = temp.doubleContraction(temp,normal);
					if(rbar1.doubleContraction(rbar1, rbar1) > 1.0e20)
					{
						std::cout<<beta0;
						std::cout<<"\n";
						std::cout<<beta1;
						std::cout<<"\n";
						std::cout<<alpha_ns.doubleContraction(alpha_ns,alpha_ns);
						std::cout<<"\n";
						std::cout<<r.doubleContraction(r,r);
						std::cout<<"\n";
					}
			        normal = rbar1 * ( 1.0 /sqrt( rbar1.doubleContraction(rbar1, rbar1) ));
			        temp.AijBjk( normal, normal );
			        pass.AijBjk( temp, normal );
			        sin3theta = -sqrt(6.0) * pass.Trace();
					if (sin3theta > 1.0)
						sin3theta = 1.0;
					else if (sin3theta < -1.0)
						sin3theta = -1.0;
			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );
					temp = normal * root23 * etamplus1 * gtheta;
					temp -= rbar1;
					Fb1 = temp.doubleContraction(temp,normal);
			   	}
				if (abs(Fb0) <= 1.0e-5)
			    {
			    	rbar = rbar0;
			    	beta = beta0;
			    }
			    else if (abs(Fb1) <= 1.0e-5)
			    {
			    	rbar = rbar1;
			    	beta = beta1;
			    }
				else
				{
				    beta = beta1 - Fb1 * ( beta1 - beta0 ) / ( Fb1 - Fb0 );
				    rbar = alpha_ns;
				    rbar += beta * r;
				    rbar -= beta * alpha_ns;
				    normal = rbar * ( 1.0 / sqrt( rbar.doubleContraction(rbar, rbar) ));
			        temp.AijBjk( normal, normal );
			        pass.AijBjk( temp, normal );
			        sin3theta = -sqrt(6.0) * pass.Trace();
					if (sin3theta > 1.0)
						sin3theta = 1.0;
					else if (sin3theta < -1.0)
						sin3theta = -1.0;
			        gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
			        		/ Mfo * ( 1 - sin3theta * sin3theta ) );
					temp = normal * root23 * etamplus1 * gtheta;
					temp -= rbar;
					Fb = temp.doubleContraction(temp,normal);
					intm = 1;
				    while (abs(Fb) > 1.0e-6)
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
					    rbar = alpha_ns;
					    rbar += beta * r;
					    rbar -= beta * alpha_ns;
				        normal = rbar * ( 1.0 / sqrt(rbar.doubleContraction(rbar, rbar)));
				        temp.AijBjk( normal, normal );
				        pass.AijBjk( temp, normal );
						sin3theta = -sqrt(6.0) * pass.Trace();
						if (sin3theta > 1.0)
							sin3theta = 1.0;
						else if (sin3theta < -1.0)
							sin3theta = -1.0;
						gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
								/ Mfo * ( 1 - sin3theta * sin3theta ) );
						temp = normal * root23 * etamplus1 * gtheta;
						temp -= rbar;
						Fb = temp.doubleContraction(temp,normal);
						intm = intm + 1;
				    }
				}
			}
		}
		normal = root32 * normal;
		N = r.doubleContraction( r,normal );

		// --------------(III)Loading/Unloading-------------------------------------
		temp = dev_stress;
		temp -= dev_stress_n;
		phi = temp.doubleContraction( temp, normal ) - ( p_nplus1 - p_n ) * N;
		temp = r_nplus1;
		temp -= r;
		phi_n = temp.doubleContraction( temp, normal );

		// --------------(IV)Unloading-------------------------------------
		if (phi < tolerance || phi_n < tolerance)
		{
			gammamonos = 0.0;
			alpha_nplus1 = r;
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
			temp = r;
			temp -= alpha_ns;
			rou = root32 * sqrt( temp.doubleContraction(temp, temp) );
			temp = rbar;
			temp -= alpha_ns;
			roubar = root32 * sqrt( temp.doubleContraction(temp, temp) );

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
				eta_nplus1 = root32 * sqrt( r_nplus1.doubleContraction(r_nplus1, r_nplus1) );
				rd = Mdc * exp( m_nd[k] * psi ) / (etamplus1 + tolerance ) * rbar;
				temp = rd;
				temp -= r;
				Dre_n = m_dre1[k] * root23 * temp.doubleContraction( temp, normal );
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
						temp = rd;
						temp -= r;
						Dir_n = m_dir[k] * exp( m_nd[k] * psi- m_eta[k] * epsvir_ns )
								* ( root23 * temp.doubleContraction( temp , normal ) ) * exp(chi);
					}
					else
					{
						temp = rd;
						temp -= r;
						Dir_n = m_dir[k] * exp ( m_nd[k] * psi- m_eta[k] * epsvir_ns ) * ( root23 * temp.doubleContraction( temp , normal )* exp(chi)
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
					if ( abs(iconv) < tolerance )//original: iconv==0.0
					{
						//dlambda=phi/(H+2*G-K*D_*N);
						dlambda = phi / abs( H + 2 * G - K * D_ * N);
						lambda += dlambda;
					}
					else
					{
						lambda = 0.5 * ( lambdamax + lambdamin );
					}
					loadindex = H * lambda;
					dtracep = lambda * D_;
					ddev_strain_p = lambda * normal;
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
					dev_stress = dev_stress_n;
					dev_stress += 2 * G / sub * dev_strain;
					dev_stress -= 2 * G / sub * dev_strain_n;
					dev_stress -= 2 * G * ddev_strain_p;

					temp = dev_stress;
					temp -= dev_stress_n;
					phi = temp.doubleContraction( temp , normal ) - ( p_nplus1 - p_n ) * N - loadindex;
					//phi=doublecontraction(dev_stress-dev_stress_n,normal) / root32 -(p_nplus1-p_n)*N / root32-loadindex;
					if (phi < -tolerance)
					{
						iconv = 1.0;
						lambdamax = lambda;
					}
					if (phi > tolerance && abs( iconv - 1.0 ) < tolerance)//original: iconv==1.0
						lambdamin = lambda;
					wr = wr + 1;
					epsvir_nplus1 = lambda * Dir_n + epsvir_ns;
					epsvre_nplus1 = lambda* Dre_n + epsvre_ns;
				}while (abs(phi) > tolerance);
				gammamonos = gammamonos + lambda;
			}
		}
		r_nplus1 = dev_stress * ( 1.0 / p_nplus1 );
		eta_nplus1 = root32 * sqrt( r_nplus1.doubleContraction(r_nplus1, r_nplus1) );

		if (eta_nplus1 >= Mfc * exp( -m_np[k] * psi ) / ( 1.0 + Mfc / 3.0 ) - tolerance)
		{
			r1 = r_nplus1 * ( 1 / sqrt( r_nplus1.doubleContraction(r_nplus1, r_nplus1) ));
	        temp.AijBjk( r1, r1 );
	        pass.AijBjk( temp, r1 );
			sin3theta = -sqrt(6.0) * pass.Trace();
			if (sin3theta > 1.0)
				sin3theta = 1.0;
			else if (sin3theta < -1.0)
				sin3theta = -1.0;
			gtheta = 1 / ( 1 + Mfc / 6.0 * ( sin3theta + sin3theta * sin3theta ) + ( Mfc - Mfo )
					/ Mfo * ( 1 - sin3theta * sin3theta ) );
			r1 = root23 * Mfc * exp( -m_np[k] * psi ) * gtheta * r1;
			if (r1.doubleContraction(r1, r1) -r_nplus1.doubleContraction(r_nplus1, r_nplus1) < tolerance)
			{
				intm = sqrt( r_nplus1.doubleContraction(r_nplus1, r_nplus1) / r1.doubleContraction(r1, r1) ) + tolerance;
				dev_stress = dev_stress * ( 1.0 / intm );
			}
			r_nplus1 = dev_stress * (1.0 / p_nplus1 );
			eta_nplus1 = root32 * sqrt(r_nplus1.doubleContraction(r_nplus1, r_nplus1));
		}

		alpha_ns = alpha_nplus1;
		epsvir_ns = epsvir_nplus1;
		epsvre_ns = epsvre_nplus1;
        epsvc_ns = epsvc_ns - epsvc_nplus1 + ( trace - trace_n ) / sub - dtracep;
		epsvc_nplus1 = epsvc_ns;

		p_n = p_nplus1;
		dev_stress_n = dev_stress;
	}
	stress_pass = dev_stress;
	stress_pass.PlusIdentity(p_n);
	temp = stress_pass * (-1);
	temp -= m_stress[k][q]; // from positive to negative

	m_stress[k][q] += temp;
	temp.QijAjkQlk( m_stress[k][q], Rot );
	m_stress[k][q] = temp;
	m_strain[k][q] = strain_nplus1 * (-1);

	m_epsvir[k][q] = epsvir_ns;
    m_epsvre[k][q] = epsvre_ns;
    m_gammamono[k][q] = gammamonos;
    m_epsvc[k][q] = epsvc_ns;
    m_etam[k][q] = etamplus1;
	m_alpha[k][q] = alpha_ns;
}
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CycLiqCPSP, std::string const &, Group * const )
} /* namespace constitutive */
} /* namespace geosx */
