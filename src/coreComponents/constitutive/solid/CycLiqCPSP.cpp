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
 *  @file CycLiqCPSP.cpp
 */

#include "CycLiqCPSP.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

CycLiqCPSP::CycLiqCPSP( string const & name, Group * const parent ):
		SolidBase( name, parent ),
		  m_defaultG0(),
		  m_defaultKappa(),
		  m_defaultH(),
		  m_defaultDre1(),
		  m_defaultDre2(),
		  m_defaultDir(),
		  m_defaultEta(),
		  m_defaultRdr(),
		  m_defaultNp(),
		  m_defaultNd(),
		  m_defaultM(),
		  m_defaultLamdac(),
		  m_defaultE0(),
		  m_defaultKsi(),
		  m_defaultEin(),
		  m_defaultInitialCycles(),
		  m_G0(),
		  m_kappa(),
		  m_h(),
		  m_dre1(),
		  m_dre2(),
		  m_dir(),
		  m_eta(),
		  m_rdr(),
		  m_np(),
		  m_nd(),
		  m_M(),
		  m_lamdac(),
		  m_e0(),
		  m_ksi(),
		  m_ein(),
		  m_initialCycles(),
		  m_strain(),
		  m_epsvir(),
		  m_epsvre(),
		  m_gammamono(),
		  m_epsvc(),
		  m_etam(),
		  m_alpha()
{
	  registerWrapper( viewKeyStruct::defaultG0String(), &m_defaultG0 ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : G0" );
	  registerWrapper( viewKeyStruct::defaultKappaString(), &m_defaultKappa ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : kappa" );
	  registerWrapper( viewKeyStruct::defaultHString(), &m_defaultH ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : h" );
	  registerWrapper( viewKeyStruct::defaultDre1String(), &m_defaultDre1 ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : dre1" );
	  registerWrapper( viewKeyStruct::defaultDre2String(), &m_defaultDre2 ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : dre2" );
	  registerWrapper( viewKeyStruct::defaultDirString(), &m_defaultDir ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : dir" );
	  registerWrapper( viewKeyStruct::defaultEtaString(), &m_defaultEta ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : eta" );
	  registerWrapper( viewKeyStruct::defaultRdrString(), &m_defaultRdr ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : rdr" );
	  registerWrapper( viewKeyStruct::defaultNpString(), &m_defaultNp ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : np" );
	  registerWrapper( viewKeyStruct::defaultNdString(), &m_defaultNd ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : nd" );
	  registerWrapper( viewKeyStruct::defaultMString(), &m_defaultM ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : M" );
	  registerWrapper( viewKeyStruct::defaultLamdacString(), &m_defaultLamdac ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : lamdac" );
	  registerWrapper( viewKeyStruct::defaultE0String(), &m_defaultE0 ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : e0" );
	  registerWrapper( viewKeyStruct::defaultKsiString(), &m_defaultKsi ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : ksi" );
	  registerWrapper( viewKeyStruct::defaultEinString(), &m_defaultEin ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : ein" );
	  registerWrapper( viewKeyStruct::defaultInitialCyclesString(), &m_defaultInitialCycles ).
	    setApplyDefaultValue( -1 ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Default for CycLiqCPSP Model Parameter : initial Cycles" );

	  registerWrapper( viewKeyStruct::G0String(), &m_G0 ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : G0" );
	  registerWrapper( viewKeyStruct::kappaString(), &m_kappa ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : kappa" );
	  registerWrapper( viewKeyStruct::hString(), &m_h ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : h" );
	  registerWrapper( viewKeyStruct::dre1String(), &m_dre1 ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : dre1" );
	  registerWrapper( viewKeyStruct::dre2String(), &m_dre2 ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : dre2" );
	  registerWrapper( viewKeyStruct::dirString(), &m_dir ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : dir" );
	  registerWrapper( viewKeyStruct::etaString(), &m_eta ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : eta" );
	  registerWrapper( viewKeyStruct::rdrString(), &m_rdr ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : rdr" );
	  registerWrapper( viewKeyStruct::npString(), &m_np ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : np" );
	  registerWrapper( viewKeyStruct::ndString(), &m_nd ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : nd" );
	  registerWrapper( viewKeyStruct::MString(), &m_M ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : M" );
	  registerWrapper( viewKeyStruct::lamdacString(), &m_lamdac ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : lamdac" );
	  registerWrapper( viewKeyStruct::e0String(), &m_e0 ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : e0" );
	  registerWrapper( viewKeyStruct::ksiString(), &m_ksi ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : ksi" );
	  registerWrapper( viewKeyStruct::einString(), &m_ein ).
	    setApplyDefaultValue( -1 ).
	    setDescription( "CycLiqCPSP Model Parameter : ein" );
	  registerWrapper( viewKeyStruct::initialCyclesString(), &m_initialCycles ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "initialCycles" );

	  registerWrapper( viewKeyStruct::strainString(), &m_strain ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "strain" );
	  registerWrapper( viewKeyStruct::epsvirString(), &m_epsvir ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "epsvir" );
	  registerWrapper( viewKeyStruct::epsvreString(), &m_epsvre ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "epsvre" );
	  registerWrapper( viewKeyStruct::gammamonoString(), &m_gammamono ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "gammamono" );
	  registerWrapper( viewKeyStruct::epsvcString(), &m_epsvc ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "epsvc" );
	  registerWrapper( viewKeyStruct::etamString(), &m_etam ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "etam" );
	  registerWrapper( viewKeyStruct::alphaString(), &m_alpha ).
	    setPlotLevel( PlotLevel::LEVEL_0 ).
	    setDescription( "alpha" );
}

CycLiqCPSP::~CycLiqCPSP()
{}

void CycLiqCPSP::postProcessInput()
{
  SolidBase::postProcessInput();

	this->getWrapper< array1d< real64 > >( viewKeyStruct::G0String() ).
	    setApplyDefaultValue( m_defaultG0 );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::kappaString() ).
	    setApplyDefaultValue( m_defaultKappa );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::hString() ).
	    setApplyDefaultValue( m_defaultH );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::dre1String() ).
	    setApplyDefaultValue( m_defaultDre1 );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::dre2String() ).
	    setApplyDefaultValue( m_defaultDre2 );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::dirString() ).
	    setApplyDefaultValue( m_defaultDir );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::etaString() ).
	    setApplyDefaultValue( m_defaultEta );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::rdrString() ).
	    setApplyDefaultValue( m_defaultRdr );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::npString() ).
	    setApplyDefaultValue( m_defaultNp );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::ndString() ).
	    setApplyDefaultValue( m_defaultNd );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::MString() ).
	    setApplyDefaultValue( m_defaultM );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::lamdacString() ).
	    setApplyDefaultValue( m_defaultLamdac );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::e0String() ).
	    setApplyDefaultValue( m_defaultE0 );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::ksiString() ).
	    setApplyDefaultValue( m_defaultKsi );
	this->getWrapper< array1d< real64 > >( viewKeyStruct::einString() ).
	    setApplyDefaultValue( m_defaultEin );

	this->getWrapper< array2d< integer > >( viewKeyStruct::initialCyclesString() ).
	    setApplyDefaultValue( -m_defaultInitialCycles );

}

void CycLiqCPSP::allocateConstitutiveData( dataRepository::Group & parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_initialCycles.resize( 0, numConstitutivePointsPerParentIndex );
  m_epsvir.resize( 0, numConstitutivePointsPerParentIndex );
  m_epsvre.resize( 0, numConstitutivePointsPerParentIndex );
  m_gammamono.resize( 0, numConstitutivePointsPerParentIndex );
  m_epsvc.resize( 0, numConstitutivePointsPerParentIndex );
  m_etam.resize( 0, numConstitutivePointsPerParentIndex );
  m_strain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_alpha.resize( 0, numConstitutivePointsPerParentIndex, 6 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CycLiqCPSP, string const &, Group * const )


} /* namespace constitutive */
} /* namespace geosx */
