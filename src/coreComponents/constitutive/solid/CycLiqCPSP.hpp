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
