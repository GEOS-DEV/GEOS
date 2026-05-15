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
 * @file Geomechanics.hpp
 * @brief To Do
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_
#define GEOS_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GeomechanicsUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GeomechanicsUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] b0 The value for the tangent elastic bulk modulus paramter 0
   * @param[in] b1 The value for the tangent elastic bulk modulus paramter 1
   * @param[in] b2 The value for the tangent elastic bulk modulus paramter 2
   * @param[in] b3 The value for the tangent elastic bulk modulus paramter 3
   * @param[in] b4 The value for the tangent elastic bulk modulus paramter 4
   * @param[in] dstrendh The value (constant) for hardened stren
   * @param[in] dfslopedh The value (constant) for hardened fslope
   * @param[in] dpeakI1dh The value (constant) for hardened peakI1
   * @param[in] dcrdh The value (constant) for hardened cr
   * @param[in] g0 The value for the tangent elastic shear modulus paramter 0
   * @param[in] g1 The value for the tangent elastic shear modulus paramter 1
   * @param[in] g2 The value for the tangent elastic shear modulus paramter 2
   * @param[in] g3 The value for the tangent elastic shear modulus paramter 3
   * @param[in] g4 The value for the tangent elastic shear modulus paramter 4
   * @param[in] p0 The value for the crush curve paramter 0
   * @param[in] p1 The value for the crush curve paramter 1
   * @param[in] p2 The value for the crush curve paramter 2
   * @param[in] p3 The value for the crush curve paramter 3
   * @param[in] p4 The value for the crush curve paramter 4
   * @param[in] peakI1 The value for the peak T1 shear limit parameter
   * @param[in] fSlope The value for the F slope shear limit parameter
   * @param[in] fSlopeFailed The value for the F slope shear limit parameter
   * @param[in] stren The value for the stren shear limit paramter
   * @param[in] ySlope The value for the Y Slope shear limit parameter
   * @param[in] beta The value for the nonassociativity parameter
   * @param[in] t1RateDependence The value for the rate dependence parameter 1
   * @param[in] t2RateDependence The value for the rate dependence parameter 2
   * @param[in] fractureEnergyReleaseRate The value for the fracture energy release rate
   * @param[in] fractureSofteningExponent controls shape of softening with damage
   * @param[in] fractureStress The root J2 value for the fracture stress 
   * @param[in] initialTemperature initial temperature
   * @param[in] Q activation energy
   * @param[in] brittleDuctileTransition The root J2 value for the fracture stress 
   * @param[in] damageEvolutionCriterion this is to trigger the damage for TXCo tests 0 or 1 for dilation vs pressure
   * @param[in] cr The value for the cap shape paramter
   * @param[in] fluidBulkModulus The value for the fluid bulk modulus
   * @param[in] fluidInitialPressure The value for the fluid initial pressure
   * @param[in] bulkModulus The ArrayView holding the bulk modulus for each element/particle.
   * @param[in] shearModulus The ArrayView holding the shear modulus for each element/particle.
   * @param[in] velocityGradient The ArrayView holding the velocity gradient for each element/particle.
   * @param[in] materialDirection The ArrayView holding the material direction for each element/particle.
   * @param[in] deformationGradient The ArrayView holding the deformation gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point.
   * @param[in] porosity The ArrayView holding the porosity for each quadrature point.
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] constitutiveUpdateFlag The ArrayView holding the constitutive update status for each quadrature point.
   * @param[in] temperature The ArrayView holding the temperature for each quardrature point.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   * @param[in] density The ArrayView holding the density data for each quadrature point.
   * @param[in] wavespeed The ArrayView holding the wavespeed data for each quadrature point.
   * @param[in] disableInelasticity Flag to disable plastic response/
   */
  GeomechanicsUpdates( real64 const & b0,
                       real64 const & b1,
                       real64 const & b2,
                       real64 const & b3,
                       real64 const & b4,
                       real64 const & dstrendh,
                       real64 const & dfslopedh,
                       real64 const & dpeakI1dh,
                       real64 const & dcrdh,
                       real64 const & g0,
                       real64 const & g1,
                       real64 const & g2,
                       real64 const & g3,
                       real64 const & g4,  
                       real64 const & p0,
                       real64 const & p1,
                       real64 const & p2,
                       real64 const & p3,
                       real64 const & p4,
                       real64 const & peakI1,
                       real64 const & fSlope,
                       real64 const & fSlopeFailed,
                       real64 const & stren,
                       real64 const & ySlope,
                       real64 const & beta,
                       real64 const & t1RateDependence,
                       real64 const & t2RateDependence,
                       real64 const & fractureEnergyReleaseRate,
                       real64 const & fractureSofteningExponent,
                       real64 const & fractureStress,
                       real64 const & initialTemperature,
                       real64 const & Q,
                       real64 const & brittleDuctileTransition,
                       int const & damageEvolutionCriterion,
                       real64 const & cr,
                       real64 const & fluidBulkModulus,
                       real64 const & fluidInitialPressure,
                       int const & enableBuckling,
                       real64 const & bucklingLength,
                       real64 const & bucklingAmplitude,
                       int const & enableCreep,
                       real64 const & creepC0,
                       real64 const & creepC1,
                       real64 const & creepC2,
                       real64 const & creepA,
                       real64 const & creepB,
                       real64 const & creepC,
                       real64 const & creepD,
                       real64 const & creepE,
                       real64 const & creepF,
                       real64 const & creepG,
                       real64 const & strainHardeningN,
                       real64 const & strainHardeningK,
                       real64 const & plasticStrainTolerance,
                       real64 const & stressReturnTolerance,
                       int const & maxAllowedSubcycles,
                       int const & failedStepResponse,
                       arrayView1d< real64 > const bulkModulus,
                       arrayView1d< real64 > const shearModulus,
                       arrayView3d< real64 const > const velocityGradient,
                       arrayView2d< real64 > const & materialDirection,
                       arrayView3d< real64 > const & deformationGradient,
                       arrayView3d< real64 > const plasticStrain,
                       arrayView2d< real64 > const porosity,
                       arrayView2d< real64 > const damage,
                       arrayView2d< integer > const constitutiveUpdateFlag,
                       arrayView1d< real64 > const temperature,
                       arrayView1d< real64 > const lengthScale,
                       arrayView1d< real64 > const strengthScale,
                       arrayView1d< real64 > const & thermalExpansionCoefficient,
                       arrayView3d< real64, solid::STRESS_USD > const & newStress,
                       arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                       arrayView2d< real64 > const & density,
                       arrayView2d< real64 > const & wavespeed,
                       bool const & disableInelasticity ):
     SolidBaseUpdates( newStress, 
                       oldStress,
                       density,
                       wavespeed,
                       thermalExpansionCoefficient, 
                       disableInelasticity ),
    m_b0( b0 ),
    m_b1( b1 ),
    m_b2( b2 ),
    m_b3( b3 ),
    m_b4( b4 ),
    m_dstrendh( dstrendh ),
    m_dfslopedh( dfslopedh ),
    m_dpeakI1dh( dpeakI1dh ),
    m_dcrdh( dcrdh ),
    m_g0( g0 ),
    m_g1( g1 ),
    m_g2( g2 ),
    m_g3( g3 ),
    m_g4( g4 ),
    m_p0( p0 ),
    m_p1( p1 ),
    m_p2( p2 ),
    m_p3( p3 ),
    m_p4( p4 ),
    m_peakI1( peakI1 ),
    m_fSlope( fSlope ),
    m_fSlopeFailed( fSlopeFailed ),
    m_stren( stren ),
    m_ySlope( ySlope ),
    m_beta( beta ),
    m_t1RateDependence( t1RateDependence ),
    m_t2RateDependence( t2RateDependence ),
    m_fractureEnergyReleaseRate( fractureEnergyReleaseRate ),
    m_fractureSofteningExponent( fractureSofteningExponent ),
    m_fractureStress( fractureStress ),
    m_initialTemperature( initialTemperature ),
    m_Q( Q ),
    m_brittleDuctileTransition( brittleDuctileTransition ),
    m_damageEvolutionCriterion ( damageEvolutionCriterion ),
    m_cr( cr ),
    m_fluidBulkModulus( fluidBulkModulus ),
    m_fluidInitialPressure( fluidInitialPressure ),
    m_enableBuckling( enableBuckling ),
    m_bucklingLength( bucklingLength ),
    m_bucklingAmplitude( bucklingAmplitude ),
    m_enableCreep( enableCreep ),
    m_creepC0( creepC0 ),
    m_creepC1( creepC1 ),
    m_creepC2( creepC2 ),
    m_creepA( creepA ),
    m_creepB( creepB ),
    m_creepC( creepC ),
    m_creepD( creepD ),
    m_creepE( creepE ),
    m_creepF( creepF ),
    m_creepG( creepG ),
    m_strainHardeningN( strainHardeningN ),
    m_strainHardeningK( strainHardeningK ),
    m_plasticStrainTolerance(plasticStrainTolerance),
    m_stressReturnTolerance(stressReturnTolerance),
    m_maxAllowedSubcycles(maxAllowedSubcycles),
    m_failedStepResponse(failedStepResponse),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus ),
    m_velocityGradient( velocityGradient ),
    m_materialDirection( materialDirection ),
    m_deformationGradient( deformationGradient ),    
    m_plasticStrain( plasticStrain ),
    m_porosity( porosity ),
    m_damage( damage ),
    m_constitutiveUpdateFlag( constitutiveUpdateFlag ),
    m_temperature( temperature ),
    m_lengthScale( lengthScale ),
    m_strengthScale( strengthScale ),
    m_geomechanicsDisableInelasticity( disableInelasticity )
  {
    GEOS_UNUSED_VAR( m_g3 );
    GEOS_UNUSED_VAR( m_g4 );  
    GEOS_UNUSED_VAR( m_p2 );
    GEOS_UNUSED_VAR( m_p4 );
    GEOS_UNUSED_VAR( m_t1RateDependence );
    GEOS_UNUSED_VAR( m_t2RateDependence );
    GEOS_UNUSED_VAR( m_fluidBulkModulus );
    GEOS_UNUSED_VAR( m_fluidInitialPressure );
  }

  /// Default copy constructor
  GeomechanicsUpdates( GeomechanicsUpdates const & ) = default;

  /// Default move constructor
  GeomechanicsUpdates( GeomechanicsUpdates && ) = default;

  /// Deleted default constructor
  GeomechanicsUpdates() = delete;

  /// Deleted copy assignment operator
  GeomechanicsUpdates & operator=( GeomechanicsUpdates const & ) = delete;

  /// Deleted move assignment operator
  GeomechanicsUpdates & operator=( GeomechanicsUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using SolidBaseUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const timeIncrement,
                                real64 const ( & beginningRotation )[3][3],
                                real64 const ( & endRotation )[3][3],
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  int computeStep( real64 const ( & D )[6],              
                   const real64 & Dt,
                   const real64 & lch,  
                   const real64 & strengthScale,  
                   const real64 & temperature,                    
                   const real64 & Zeta_n,                
                   const real64 & coher_n,               
                   const real64 & porosity_n,            
                   const real64 & buckling,               // scalar-valued buckling multiplier at start of step(t_n)
                   real64 const ( & sigma_n )[6],       
                   real64 const ( & ep_n )[6],          
                   real64 & Zeta_p,                      
                   real64 & coher_p,                     
                   real64 & porosity_p,                  
                   real64 ( & sigma_p )[6],   
                   real64 ( & ep_p )[6],
                   int & constitutiveUpdateFlag ) const;  

  GEOS_HOST_DEVICE
  void computeElasticProperties( real64 & bulk,
                                 real64 & shear ) const;

  GEOS_HOST_DEVICE
  void computeElasticProperties( real64 const ( &stress )[6],
                                 real64 const ( &ep )[6],
                                 const real64 & p3_crush_curve,
                                 const real64 & fluid_pressure_initial,
                                 const real64 & Km,     
                                 const real64 & Kf,     
                                 const real64 & C1,     
                                 const real64 & ev0,    
                                 const real64 & phi_i,
                                 real64 & bulk,
                                 real64 & shear ) const;
  GEOS_HOST_DEVICE
  void computeHypoelasticTrialStress( const real64 dt,                 
                                      const real64 bulk,               
                                      const real64 shear,              
                                      real64 const ( & stress_old )[6],
                                      real64 const ( & D )[6],         
                                      real64 ( & stress_new )[6] ) const;    

  GEOS_HOST_DEVICE
  void computeInvariants( real64 const ( & stress )[6],
                          real64 ( & S )[6],
                          real64 & I1,
                          real64 & J2,
                          real64 & rJ2 ) const;

  GEOS_HOST_DEVICE
  real64 derivative_f(real64 const & peakI1,
                      const real64 & cr,
                      const real64 & X,
                      const real64 & a1,
                      const real64 & a2,
                      const real64 & a3,
                      const real64 & a4,
                      const real64 & I1) const;


  GEOS_HOST_DEVICE
  real64 computeBD( const real64 & a1,
                    const real64 & a2,
                    const real64 & a3,
                    const real64 & a4) const;

  GEOS_HOST_DEVICE
  void computeCoher(const real64 & a1,
                    const real64 & a2,
                    const real64 & a3,
                    const real64 & a4,
                    const real64 & d_evp,         // increment in vol plastic strain
                    real64 const ( & d_ep )[6],    // increment in total plastic strain
                    const real64 & I1_trial,      // trial value of I1
                    const real64 & rJ2_trial,  // trial value of rootJ2
                    const real64 & I1_0,          // I1 value on yield surface
                    const real64 & rJ2_0,         // rJ2 value on yield surface
                    const real64 & lch,           // length scale
                    const real64 & strengthScale, // strength scale
                    const real64 & coher_old,     // old coherence = 1-damage
                          real64 & coher_new     // OUPUT: new value of coherconst real64 & lch,       // length scale
) const;    


  GEOS_HOST_DEVICE
  real64 computeX( const real64 & evp,
		               const real64 & phi_i,
		               const real64 & Km,
		               const real64 & Kf,
		               const real64 & C1,
		               const real64 & ev0,
		               const real64 & fluid_pressure_initial,
                   const real64 & buckling  // crush strength multiplier for cascading collapse.
                   ) const;

  GEOS_HOST_DEVICE
  int computeStepDivisions( const real64 & X,
		                        const real64 & Zeta,
		                        real64 const ( & ep )[6],
		                        real64 const ( & sigma_n )[6],
		                        real64 const ( & sigma_trial )[6],
		                        const real64 & fluid_B0,
		                        const real64 & fluid_pressure_initial,
		                        const real64 & Km,       
		                        const real64 & Kf,       
		                        const real64 & C1,       
		                        const real64 & ev0,      
		                        const real64 & phi_i,
                                int & constitutiveUpdateFlag ) const;   

  GEOS_HOST_DEVICE
  int computeSubstep( real64 const ( & D )[6],        
                      const real64 & dt,
                      const real64 & lch,
                      const real64 & strengthScale,              
                      real64 const ( & sigma_old )[6],
                      real64 const ( & ep_old )[6],   
                      const real64 & X_old,
                      const real64 & Zeta_old,
                      const real64 & coher_old,       
                      const real64 & buckling,         // scalar valued buckling crush-curve modifier
                      const real64 & fluid_pressure_initial,
                      const real64 & Km,
                      const real64 & Kf,
                      const real64 & C1,
                      const real64 & ev0,
                      const real64 & phi_i,                   
                      real64 ( & sigma_new )[6],           
                      real64 ( & ep_new )[6],              
                      real64 & X_new,                      
                      real64 & Zeta_new,                   
                      real64 & coher_new ) const;	                

  GEOS_HOST_DEVICE
  int nonHardeningReturn( const real64 & I1_trial,            
                          const real64 & rJ2_trial,          
                          real64 const ( & S_trial )[6],       
                          const real64 & I1_old,               
                          const real64 & rJ2_old,          
                          real64 const ( & S_old )[6],        
                          real64 const ( & d_e )[6],           
                          const real64 & X,                    
                          const real64 & Zeta,                 
                          const real64 & coher,     
                          const real64 & hardening,          
                          const real64 & strengthScale,
                          const real64 & buckling,
                          const real64 & bulk,                 
                          const real64 & shear,
                          real64 & I1_new,                     
                          real64 & rJ2_new,
                          real64 ( & S_new )[6],
                          real64 ( & d_ep_new )[6] ) const;          

  GEOS_HOST_DEVICE
  void transformedBisection( real64 & z_0,
                             real64 & r_0,
                             const real64 & z_trial,
                             const real64 & r_trial,
                             const real64 & X,
                             const real64 & Zeta,
                             const real64 & coher,
                             const real64 & hardening,
                             const real64 & strengthScale,
                             const real64 & buckling,
                             const real64 & a1,
                             const real64 & a2,
                             const real64 & a3,
                             const real64 & a4,
                             const real64 & r_to_rJ2 ) const;

  GEOS_HOST_DEVICE
  int transformedYieldFunction( const real64 & z,
                                const real64 & r,
                                const real64 & X,
                                const real64 & Zeta,
										            const real64 & coher,
                                const real64 & hardening,
                                const real64 & strengthScale,
                                const real64 & buckling,
										            const real64 & a1,
										            const real64 & a2,
										            const real64 & a3,
										            const real64 & a4,
                                const real64 & r_to_rJ2 ) const;

  GEOS_HOST_DEVICE
  int computeYieldFunction( const real64 & I1,
                            const real64 & rJ2,
                            const real64 & X,
                            const real64 & Zeta,
									          const real64 & coher,
                            const real64 & hardening,
                            const real64 & strengthScale,
                            const real64 & buckling,
									          const real64 & a1,
									          const real64 & a2,
									          const real64 & a3,
									          const real64 & a4 ) const;

  GEOS_HOST_DEVICE
  real64 computedZetadevp( real64 const & fluid_pressure_initial,
                           real64 const & Km,
                           real64 const & Kf,
                           real64 const & ev0,
                           real64 const & phi_i,
                           real64 Zeta,
                           real64 evp ) const;

  GEOS_HOST_DEVICE
  void computeLimitParameters( real64 & a1,
		                           real64 & a2,
		                           real64 & a3,
		                           real64 & a4,
		                           const real64 & coher,
                               const real64 & hardening,
                               const real64 & strengthScale,
                               const real64 & buckling ) const;

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;

private:
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int failStep( real64 const & Zeta_n,
                real64 const & coher_n,
                real64 const & porosity_n,
                real64 const ( & sigma_n )[6],
                real64 const ( & ep_n )[6],
                real64 & Zeta_p,
                real64 & coher_p,
                real64 & porosity_p,
                real64 ( & sigma_p )[6],
                real64 ( & ep_p )[6] ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int failSubstep( real64 const ( & sigma_old )[6],
                   real64 const ( & ep_old )[6],
                   real64 const & X_old,
                   real64 const & Zeta_old,
                   real64 const & coher_old,
                   real64 ( & sigma_new )[6],
                   real64 ( & ep_new )[6],
                   real64 & X_new,
                   real64 & Zeta_new,
                   real64 & coher_new ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void applyCreepRelaxation( real64 const & Dt,
                             real64 const & temperature,
                             real64 const & phi_i,
                             real64 const & Km,
                             real64 const & Kf,
                             real64 const & C1,
                             real64 const & ev0,
                             real64 const & fluid_pressure_initial,
                             real64 const & buckling,
                             real64 ( & sigma )[6],
                             real64 ( & ep )[6],
                             real64 & X ) const;

  real64 const & m_b0;
  real64 const & m_b1;
  real64 const & m_b2;
  real64 const & m_b3;
  real64 const & m_b4;
  real64 const & m_dstrendh;
  real64 const & m_dfslopedh;
  real64 const & m_dpeakI1dh;
  real64 const & m_dcrdh;
  real64 const & m_g0;
  real64 const & m_g1;
  real64 const & m_g2;
  real64 const & m_g3;
  real64 const & m_g4;  
  real64 const & m_p0;
  real64 const & m_p1;
  real64 const & m_p2;
  real64 const & m_p3;
  real64 const & m_p4;
  real64 const & m_peakI1;
  real64 const & m_fSlope;
  real64 const & m_fSlopeFailed;
  real64 const & m_stren;
  real64 const & m_ySlope;
  real64 const & m_beta;
  real64 const & m_t1RateDependence;
  real64 const & m_t2RateDependence;
  real64 const & m_fractureEnergyReleaseRate;
  real64 const & m_fractureSofteningExponent;
  real64 const & m_fractureStress;
  real64 const & m_initialTemperature;
  real64 const & m_Q;
  real64 const & m_brittleDuctileTransition;
  int const & m_damageEvolutionCriterion;
  real64 const & m_cr;
  real64 const & m_fluidBulkModulus;
  real64 const & m_fluidInitialPressure;
  int const & m_enableBuckling;
  real64 const & m_bucklingLength;
  real64 const & m_bucklingAmplitude;
  int const & m_enableCreep;
  real64 const & m_creepC0;
  real64 const & m_creepC1;
  real64 const & m_creepC2;  
  real64 const & m_creepA;
  real64 const & m_creepB;
  real64 const & m_creepC;
  real64 const & m_creepD;
  real64 const & m_creepE;
  real64 const & m_creepF;
  real64 const & m_creepG;
  real64 const & m_strainHardeningN;
  real64 const & m_strainHardeningK;
  real64 const & m_plasticStrainTolerance;
  real64 const & m_stressReturnTolerance;
  int const & m_maxAllowedSubcycles;
  int const & m_failedStepResponse;

  /// A reference to the ArrayView holding the bulk modulus for each element/particle.
  arrayView1d< real64 > const m_bulkModulus;
  
  /// A reference to the ArrayView holding the shear modulus for each element/particle.
  arrayView1d< real64 > const m_shearModulus;

  /// A reference to the ArrayView holding the velocity gradient for each element/particle.
  arrayView3d< real64 const > const m_velocityGradient;

  /// A reference to the ArrayView holding the material direction for each element/particle.
  arrayView2d< real64 > const m_materialDirection;

  /// A reference to the ArrayView holding the deformation gradient for each element/particle.
  arrayView3d< real64 > const m_deformationGradient;

  /// A reference to the ArrayView holding the plastic strain for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the porosity for each quadrature point.
  arrayView2d< real64 > const m_porosity;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the constitutive update status for each quadrature point.
  ///  0 = good update, 1 = warning/retried or capped subcycling, -1 = failed update/delete flag.
  arrayView2d< integer > const m_constitutiveUpdateFlag;

  /// A reference to the ArrayView holding the temperature for each quadrature point.
  arrayView1d< real64 > const m_temperature;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// A reference to the ArrayView holding the strength scale for each element/particle.
  arrayView1d< real64 > const m_strengthScale;

  /// Local copy of the inelasticity-disable flag, used by this update kernel.
  bool const & m_geomechanicsDisableInelasticity;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  GEOS_ERROR( "smallStrainUpdate not implemented for Geomechanics model" );
#endif
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( & stress )[6],
                                             DiscretizationOps & stiffness ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  GEOS_ERROR( "smallStrainUpdate not implemented for Geomechanics model" );
#endif
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( & strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for Geomechanics model" );
#endif
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( & beginningRotation )[3][3],
                                                    real64 const ( & endRotation )[3][3],
                                                    real64 const ( & strainIncrement )[6],
                                                    real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( strainIncrement );

  // call the constitutive model
  GeomechanicsUpdates::smallStrainUpdateHelper( k, 
                                                q, 
                                                timeIncrement,                
                                                beginningRotation,
                                                endRotation,
                                                stress );

  // GEOS_LOG_RANK( "After stress update | k: " << k << ", " << 
  //                 "porosity: " << m_porosity[k] << ", " << 
  //                 "damage: " << m_damage[k] << ", " << 
  //                 "stress: {" << stress[0] << ", " << 
  //                               stress[1] << ", " << 
  //                               stress[2] << ", " << 
  //                               stress[3] << ", " << 
  //                               stress[4] << ", " << 
  //                               stress[5] << "}" );

  // save new stress and return
  saveStress( k, q, stress );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdateHelper( localIndex const k,
                                                   localIndex const q,
                                                   real64 const timeIncrement,  
                                                   real64 const ( & beginningRotation )[3][3],
                                                   real64 const ( & endRotation )[3][3],
                                                   real64 ( & stress )[6] ) const
{   
  // GEOS_UNUSED_VAR( timeIncrement );
  // GEOS_UNUSED_VAR( endRotation );
  
  real64 beginningRotationTranspose[3][3] = { { 0.0 } };
  LvArray::tensorOps::transpose< 3, 3 >( beginningRotationTranspose, beginningRotation );

  //  temp matrix for computing rotations.
  real64 tempMat[ 3 ][ 3 ]= { { 0 } };
  
  // Unrotate velocity gradient:
  real64 unrotatedVelocityGradient[3][3] = { { 0 } };
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, m_velocityGradient[k] );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, tempMat, beginningRotation );
  real64 unrotatedVelocityGradientTranspose[3][3] = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

  // CC: Is there an LvArray operation to get the symmetric part of a matrix?
  // Symmetric part of unrotated velocity gradient, the strain "rate" for the step
  real64 denseD[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );
  real64 D[6] = { 0 };
  LvArray::tensorOps::denseToSymmetric<3>( D, denseD );

  // vvv-----------------------------------------------------------------------BUCKLING
  // The deformation gradient and material direction are used to compute a direction
  // strain for use with a micro-structure buckling model
  // make sure material direction is normalized.
  real64 materialDirection[3] = { 0 };
  LvArray::tensorOps::copy< 3 >( materialDirection, m_materialDirection[k] );
  LvArray::tensorOps::normalize< 3 >( materialDirection );


  // This form adjusts the buckling based on total strain, and then computes a return based on 
  // the bucklin-adjusted yield surface, but there will generally be cap hardening associated with
  // that return such that the returned state is never to a fully buckled value. (.e.g if buckling
  // amplitude = 0.2, the actual stress oscillations will be much less that 80% of)
  real64 buckling= 1.0;   // strength-scale multiplier from buckling of truss-like microstructure
  if( m_enableBuckling > 0)
  {

    // m_enableBuckling
    // m_bucklingLength
    // m_bucklingAmplitude

    // Number of unit cells per element, round up to nearest integer
    real64 beta = std::ceil( m_lengthScale[k] / m_bucklingLength );

    // Compute the length-scale dependent compaction-buckling scalar used to 
    // modify crush and shear strength:
    real64 J = 1.; // volumetric or directional stretch
    if (m_enableBuckling == 1)
    { // Isotropic:
      J = LvArray::tensorOps::determinant< 3 >( m_deformationGradient[k] );
    }
    else if (m_enableBuckling == 2)
    { // Anisotropic:
      real64 temp[3] = { 0 };
      LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( temp, m_deformationGradient[k], materialDirection  );
      J = LvArray::tensorOps::AiBi< 3 >( temp, materialDirection);
    }
    else
    {
      // Invalid buckling flags are rejected during input validation. Keep the
      // multiplier neutral here in case a device kernel reaches this branch.
      J = 1.0;
    }

    J = std::max( J, 1.0e-16 );

    real64 ev0 = m_p0 / ( 2.*m_b0 );
    real64 ev1 = ev0 - m_p3;
    real64 ev = log(J);     // volumetric or directional strain
    real64 normalizedStrain; // normalizedStrain

    if (ev >= ev0)
    {
      normalizedStrain = 0.;
    }
    else if (ev0 > ev && ev > ev1)
    {
      normalizedStrain = (ev-ev0)/(ev1-ev0);
    }
    else
    {
      normalizedStrain = 1.0;
    }
    // strength scale multiplier to crush curve (0: complete losss of strength, 1: no effect)
    real64 pi = 3.141592653589793;
    buckling = 1.0 - m_bucklingAmplitude*pow( sin( -1.0*beta*pi*normalizedStrain) , 2 );
    buckling = fmin(1.0,fmax(0.0,buckling));
  }
  // ^^^-----------------------------------------------------------------------BUCKLING




  // unrotated beginning-of-step stress:  sigma_old_unrotated = R_old^T*sigma_old*R_old
  real64 oldStress[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( oldStress, m_oldStress[k][q] ); 

  // Note stress gets un-rotated in solid utilities hypo update 2 
  
  // Unrotated beginning-of-step plastic strain. The stored state uses engineering
  // shear components, while Rij_eq_AikSymBklAjl expects tensor shear components.
  // Convert a local copy only; do not mutate particle state during the update.
  real64 oldPlasticStrain[6] = { 0.0 };
  real64 storedPlasticStrain[6] = { 0.0 };
  LvArray::tensorOps::copy< 6 >( storedPlasticStrain, m_plasticStrain[k][q] );
  storedPlasticStrain[3] *= 0.5;
  storedPlasticStrain[4] *= 0.5;
  storedPlasticStrain[5] *= 0.5;
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( oldPlasticStrain, beginningRotationTranspose, storedPlasticStrain );
  oldPlasticStrain[3] *= 2.0;
  oldPlasticStrain[4] *= 2.0;
  oldPlasticStrain[5] *= 2.0;

  real64 characteristicLength = m_lengthScale[k];
  
  real64 strengthScale = m_strengthScale[k];
  real64 temperature = m_temperature[k];

  m_constitutiveUpdateFlag[k][q] = 0;
  int constitutiveUpdateFlag = 0;

  real64 oldZeta = 0.0;
  real64 oldCoher = 1.0 - m_damage[k][q];
  real64 oldPorosity = m_porosity[k][q];
  real64 newZeta;
  real64 newCoher;
  real64 newPorosity;
  real64 newStress[6] = {0.};
  real64 newPlasticStrain[6] = {0.};

  // TODO: make elastic properties temperature dependent. Fluid effects are
  // intentionally disabled in this implementation until a persistent fluid
  // backstress state is added.
  computeElasticProperties( oldStress,
                            oldPlasticStrain,
                            0.0,
                            0.0,
                            0.0,             // Matrix bulk modulus, unused when Kf=0
                            0.0,             // Fluid bulk modulus disabled
                            0.0,             // Term to simplify the fluid model expressions
                            0.0,             // Zero fluid pressure volumetric strain
                            0.0,             // Initial porosity, unused when Kf=0
                            m_bulkModulus[k],
                            m_shearModulus[k] );

  // Use a conservative CFL stiffness: at least the current tangent stiffness and
  // at least the high-pressure dry-material limit.
  real64 cflBulk = 0.0;
  real64 cflShear = 0.0;
  computeElasticProperties( cflBulk, cflShear );
  cflBulk = std::max( cflBulk, m_bulkModulus[k] );
  cflShear = std::max( cflShear, m_shearModulus[k] );
  real64 const cflModulus = std::max( 0.0, cflBulk + ( 4.0 / 3.0 ) * cflShear );
  m_wavespeed[k][q] = m_density[k][q] > 0.0 ? sqrt( cflModulus / m_density[k][q] ) : 0.0;

  int errorFlag = 0;
  if( m_geomechanicsDisableInelasticity )
  {
    computeHypoelasticTrialStress( timeIncrement,
                                   m_bulkModulus[k],
                                   m_shearModulus[k],
                                   oldStress,
                                   D,
                                   newStress );
    LvArray::tensorOps::copy< 6 >( newPlasticStrain, oldPlasticStrain );
    newZeta = oldZeta;
    newCoher = oldCoher;
    newPorosity = oldPorosity;
  }
  else
  {
    errorFlag = computeStep( D,               // strain "rate"
                 timeIncrement,              // time step (s)
                 characteristicLength,        // length scale
                 strengthScale,               // scaler for strength
                 temperature,                 // scalar for temperature
                 oldZeta,                     // trace of isotropic backstress at start of step(t_n)
                 oldCoher,                    // scalar-valued coherence at start of step(t_n)
                 oldPorosity,                 // scalar-valued porosity at start of step(t_n)
                 buckling,                    // buckling strength multiplier.
                 oldStress,                   // unrotated stress at start of step(t_n)
                 oldPlasticStrain,            // plastic strain at start of step(t_n)
                 newZeta,                     // trace of isotropic backstress at end of step(t_n+1)
                 newCoher,                    // scalar-valued coherence at end of step(t_n+1)
                 newPorosity,                 // scalar-valued porosity at end of step(t_n+1)
                 newStress,                   // unrotated stress at end of step(t_n+1)
                 newPlasticStrain,            // plastic strain at end of step(t_n+1)
                 constitutiveUpdateFlag );     // update status for solver warning/delete handling
  }

  if ( errorFlag == 1 )
  {
    constitutiveUpdateFlag = -1;
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    if ( m_failedStepResponse == 1 || m_failedStepResponse == 2 || m_failedStepResponse == 3 )
    {
      GEOS_LOG_RANK( "Geomechanics model failed to converge for material point " << k <<
                     ", quadrature point " << q << ". Reverting to the beginning-of-step state and setting constitutiveUpdateFlag=-1." );
    }
#endif
    // The solver owns particle deletion. A failed constitutive update is exposed
    // through constitutiveUpdateFlag=-1 while the state is reverted by computeStep.
  }

  m_constitutiveUpdateFlag[k][q] = constitutiveUpdateFlag;

  // Note stress gets re-rotated in solid utilities hypo update 2 
  LvArray::tensorOps::copy< 6 >( stress, newStress ); // Copy to material data:

  // re-rotate end-of-step plastic strain: ep_new = R_new*ep_new*R_new^T
  newPlasticStrain[3] *= 0.5;
  newPlasticStrain[4] *= 0.5;
  newPlasticStrain[5] *= 0.5;
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( m_plasticStrain[k][q], endRotation, newPlasticStrain );
  m_plasticStrain[k][q][3] *= 2.0;
  m_plasticStrain[k][q][4] *= 2.0;
  m_plasticStrain[k][q][5] *= 2.0;

  m_damage[k][q] = 1. - newCoher; // Copy to material data:
  m_porosity[k][q] = newPorosity; // Copy to material data:

  //return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::failStep( real64 const & Zeta_n,
                                   real64 const & coher_n,
                                   real64 const & porosity_n,
                                   real64 const ( & sigma_n )[6],
                                   real64 const ( & ep_n )[6],
                                   real64 & Zeta_p,
                                   real64 & coher_p,
                                   real64 & porosity_p,
                                   real64 ( & sigma_p )[6],
                                   real64 ( & ep_p )[6] ) const
{
  Zeta_p = Zeta_n;
  coher_p = coher_n;
  porosity_p = porosity_n;

  LvArray::tensorOps::copy< 6 >( sigma_p, sigma_n );
  LvArray::tensorOps::copy< 6 >( ep_p, ep_n );

  return 1;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::failSubstep( real64 const ( & sigma_old )[6],
                                      real64 const ( & ep_old )[6],
                                      real64 const & X_old,
                                      real64 const & Zeta_old,
                                      real64 const & coher_old,
                                      real64 ( & sigma_new )[6],
                                      real64 ( & ep_new )[6],
                                      real64 & X_new,
                                      real64 & Zeta_new,
                                      real64 & coher_new ) const
{
  LvArray::tensorOps::copy< 6 >( sigma_new, sigma_old );
  LvArray::tensorOps::copy< 6 >( ep_new, ep_old );

  X_new = X_old;
  Zeta_new = Zeta_old;
  coher_new = coher_old;

  return 1;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::applyCreepRelaxation( real64 const & Dt,
                                                real64 const & temperature,
                                                real64 const & phi_i,
                                                real64 const & Km,
                                                real64 const & Kf,
                                                real64 const & C1,
                                                real64 const & ev0,
                                                real64 const & fluid_pressure_initial,
                                                real64 const & buckling,
                                                real64 ( & sigma )[6],
                                                real64 ( & ep )[6],
                                                real64 & X ) const
{
  if( m_enableCreep != 1 )
  {
    return;
  }

  // Creep is applied once over the full explicit step. The plasticity retry loop
  // then starts from this relaxed state for every subcycle attempt.
  real64 referenceTemperature = m_initialTemperature;
  real64 gasConstantR = 8.314;
  real64 creepActivationEnergy = m_Q;

  real64 creepRateTemperatureMultiplier = 1.0;
  if( temperature > 1.0e-10 && referenceTemperature > 1.0e-10 )
  {
    creepRateTemperatureMultiplier = exp( -creepActivationEnergy *
                                          ( 1.0 / ( gasConstantR * temperature ) -
                                            1.0 / ( gasConstantR * referenceTemperature ) ) );
  }

  real64 const equilibriumPorosityPressureExponent = m_creepD;
  real64 const equilibriumPorosityOffset = m_creepE;
  real64 const compactionRatePressureExponent = m_creepF;
  real64 const rootTwoThirds = 0.81649658092772603273242802490196;

  real64 stress_iso[6] = { 0.0 };
  real64 stress_dev[6] = { 0.0 };
  real64 iso_old = 0.0;
  real64 vonMisesStress_old = 0.0;
  twoInvariant::stressDecomposition( sigma,
                                     iso_old,
                                     vonMisesStress_old,
                                     stress_dev );
  LvArray::tensorOps::scale< 6 >( stress_dev, rootTwoThirds * vonMisesStress_old );
  stress_iso[0] = iso_old;
  stress_iso[1] = iso_old;
  stress_iso[2] = iso_old;

  real64 ep_iso[6] = { 0.0 };
  real64 ep_dev[6] = { 0.0 };
  real64 ep_iso_old = 0.0;
  real64 ep_rootJ2_old = 0.0;
  twoInvariant::stressDecomposition( ep,
                                     ep_iso_old,
                                     ep_rootJ2_old,
                                     ep_dev );
  LvArray::tensorOps::scale< 6 >( ep_dev, rootTwoThirds * ep_rootJ2_old );
  ep_iso[0] = ep_iso_old;
  ep_iso[1] = ep_iso_old;
  ep_iso[2] = ep_iso_old;

  real64 bulk = 0.0;
  real64 shear = 0.0;
  computeElasticProperties( sigma,
                            ep,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk,
                            shear );

  if( shear > 0.0 )
  {
    real64 const elasticVMShearStrain = vonMisesStress_old / ( 3.0 * shear );

    if( elasticVMShearStrain > 1.0e-12 && m_creepC2 > 1.0e-16 )
    {
      real64 const equilibriumShearStrainConstant = m_creepC0;
      real64 const equilibriumShearStrainExponent = m_creepC1;
      real64 const shearStrainRateConstant = creepRateTemperatureMultiplier * m_creepC2;

      real64 const plasticVMshearStrain = rootTwoThirds * ep_rootJ2_old;
      real64 const equilibriumVMShearStrain = equilibriumShearStrainConstant *
                                              std::pow( vonMisesStress_old,
                                                        equilibriumShearStrainExponent );
      real64 creepVMStrainIncrement = Dt * shearStrainRateConstant *
                                       std::max( equilibriumVMShearStrain - plasticVMshearStrain, 0.0 );
      creepVMStrainIncrement = std::min( elasticVMShearStrain, creepVMStrainIncrement );

      real64 const devStressCreepScale = std::max( elasticVMShearStrain - creepVMStrainIncrement, 0.0 ) /
                                         elasticVMShearStrain;

      real64 creepStrainIncrement[6] = { 0.0 };
      LvArray::tensorOps::copy< 6 >( creepStrainIncrement, stress_dev );
      LvArray::tensorOps::scale< 6 >( creepStrainIncrement,
                                      0.5 / shear * ( 1.0 - devStressCreepScale ) );

      LvArray::tensorOps::add< 6 >( ep_dev, creepStrainIncrement );
      LvArray::tensorOps::scale< 6 >( stress_dev, devStressCreepScale );
    }
  }

  // Volumetric creep. Pressure p is positive in compression.
  real64 const p = -iso_old;
  real64 const evp = ep[0] + ep[1] + ep[2];
  real64 const phi_p = std::max( 1.0e-10, 1.0 + exp( -evp ) * ( phi_i - 1.0 ) );

  real64 pressureTerm = 0.0;
  if( p > 1.0e-12 && m_creepB > 0.0 )
  {
    pressureTerm = std::pow( p / m_creepB, equilibriumPorosityPressureExponent );
  }

  real64 const phi_e = std::max( 1.0e-10,
                                 m_creepA * exp( -pressureTerm ) +
                                 equilibriumPorosityOffset +
                                 ( -3.0e-6 * std::pow( std::max( temperature - referenceTemperature, 0.0 ), 2.0 ) ) -
                                 ( 0.0002 * std::max( temperature - referenceTemperature, 0.0 ) ) );

  if( phi_p - phi_e > 1.0e-10 &&
      p > 1.0e-12 &&
      m_creepC > 1.0e-16 &&
      m_creepG > 0.0 &&
      evp + m_p3 > 1.0e-10 &&
      bulk > 0.0 )
  {
    real64 const dphidt = -creepRateTemperatureMultiplier *
                          std::pow( p / m_creepG, compactionRatePressureExponent ) *
                          m_creepC * ( phi_p - phi_e );

    real64 const phi_c = std::max( phi_e, phi_p + dphidt * Dt );
    real64 evp_c = log( ( phi_i - 1.0 ) / ( phi_c - 1.0 ) );
    evp_c = std::max( evp_c, -0.99999999 * m_p3 );
    real64 devp = evp_c - evp;

    real64 p_c = 0.0;
    if( devp < -p / bulk )
    {
      devp = -p / bulk;
      p_c = 0.0;
      evp_c = evp + devp;
    }
    else
    {
      p_c = std::max( 0.0, p + bulk * devp );
    }

    stress_iso[0] = -p_c;
    stress_iso[1] = -p_c;
    stress_iso[2] = -p_c;
    ep_iso[0] = evp_c / 3.0;
    ep_iso[1] = evp_c / 3.0;
    ep_iso[2] = evp_c / 3.0;

    X = computeX( evp_c,
                  phi_i,
                  Km,
                  Kf,
                  C1,
                  ev0,
                  fluid_pressure_initial,
                  buckling );
  }

  LvArray::tensorOps::copy< 6 >( sigma, stress_dev );
  LvArray::tensorOps::add< 6 >( sigma, stress_iso );

  LvArray::tensorOps::copy< 6 >( ep, ep_dev );
  LvArray::tensorOps::add< 6 >( ep, ep_iso );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeStep( real64 const ( & D )[6],               // strain "rate"
                                      const real64 & Dt,                     // time step (s)
                                      const real64 & lch,                    // length scale
                                      const real64 & strengthScale,          // strength scale
                                      const real64 & temperature,            // temperature
                                      const real64 & Zeta_n,                 // trace of isotropic backstress at start of step(t_n)
                                      const real64 & coher_n,                // scalar-valued coherence at start of step(t_n)
                                      const real64 & porosity_n,             // scalar-valued porosity at start of step(t_n)
                                      const real64 & buckling,               // scalar-valued buckling multiplier at start of step(t_n)
                                      real64 const ( & sigma_n )[6],         // unrotated stress at start of step(t_n)
                                      real64 const ( & ep_n )[6],            // plastic strain at start of step(t_n)
                                      real64 & Zeta_p,                       // trace of isotropic backstress at end of step(t_n+1)
                                      real64 & coher_p,                      // scalar-valued coherence at end of step(t_n+1)
                                      real64 & porosity_p,                   // scalar-valued porosity at end of step(t_n+1)
                                      real64 ( & sigma_p )[6],               // unrotated stress at end of step(t_n+1)
                                      real64 ( & ep_p )[6],                  // plastic strain at end of step(t_n+1)
                                      int & constitutiveUpdateFlag            // 0 good, 1 warning, -1 failed
) const
{
  constitutiveUpdateFlag = 0;

  // Fluid effects are intentionally disabled until a persistent fluid/backstress
  // state is added following the original Arensica/Homel implementation.
  real64 const fluid_B0 = 0.0;
  real64 const fluid_pressure_initial = 0.0;

  real64 const phi_i = 1.0 - exp( -m_p3 );
  real64 const Km = m_b0 + m_b1;
  real64 const Kf = fluid_B0;
  real64 const C1 = Kf * ( 1.0 - phi_i ) + Km * phi_i;
  real64 const ev0 = Kf > 0.0 ? C1 * fluid_pressure_initial / ( Kf * Km ) : 0.0;

  real64 const evp_n = LvArray::tensorOps::symTrace< 3 >( ep_n );
  real64 const X_n = computeX( evp_n,
                               phi_i,
                               Km,
                               Kf,
                               C1,
                               ev0,
                               fluid_pressure_initial,
                               buckling );

  // Conservative elastic properties are used for the initial trial stress and
  // substep-count estimate.
  real64 bulk = 0.0;
  real64 shear = 0.0;
  computeElasticProperties( bulk, shear );

  real64 sigma_trial[6] = { 0.0 };
  computeHypoelasticTrialStress( Dt,
                                 bulk,
                                 shear,
                                 sigma_n,
                                 D,
                                 sigma_trial );

  int nsub = computeStepDivisions( X_n,
                                   Zeta_n,
                                   ep_n,
                                   sigma_n,
                                   sigma_trial,
                                   fluid_B0,
                                   fluid_pressure_initial,
                                   Km,
                                   Kf,
                                   C1,
                                   ev0,
                                   phi_i,
                                   constitutiveUpdateFlag );

  if( nsub < 0 )
  {
    constitutiveUpdateFlag = -1;
    return failStep( Zeta_n,
                     coher_n,
                     porosity_n,
                     sigma_n,
                     ep_n,
                     Zeta_p,
                     coher_p,
                     porosity_p,
                     sigma_p,
                     ep_p );
  }

  nsub = std::max( 1, nsub );
  int const chimax = std::max( 1, m_maxAllowedSubcycles / nsub );

  real64 stepSigma[6] = { 0.0 };
  real64 stepEp[6] = { 0.0 };
  LvArray::tensorOps::copy< 6 >( stepSigma, sigma_n );
  LvArray::tensorOps::copy< 6 >( stepEp, ep_n );
  real64 stepX = X_n;

  applyCreepRelaxation( Dt,
                        temperature,
                        phi_i,
                        Km,
                        Kf,
                        C1,
                        ev0,
                        fluid_pressure_initial,
                        buckling,
                        stepSigma,
                        stepEp,
                        stepX );

  int chi = 1;
  while( chi <= chimax )
  {
    int const totalSubsteps = chi * nsub;
    real64 const dt = Dt / real64( totalSubsteps );

    real64 X_old = stepX;
    real64 X_new = stepX;
    real64 Zeta_old = Zeta_n;
    real64 Zeta_new = Zeta_n;
    real64 coher_old = coher_n;
    real64 coher_new = coher_n;
    real64 sigma_old[6] = { 0.0 };
    real64 sigma_new[6] = { 0.0 };
    real64 ep_old[6] = { 0.0 };
    real64 ep_new[6] = { 0.0 };

    LvArray::tensorOps::copy< 6 >( sigma_old, stepSigma );
    LvArray::tensorOps::copy< 6 >( ep_old, stepEp );

    bool attemptFailed = false;
    for( int substep = 0; substep < totalSubsteps; ++substep )
    {
      int const substepFlag = computeSubstep( D,
                                              dt,
                                              lch,
                                              strengthScale,
                                              sigma_old,
                                              ep_old,
                                              X_old,
                                              Zeta_old,
                                              coher_old,
                                              buckling,
                                              fluid_pressure_initial,
                                              Km,
                                              Kf,
                                              C1,
                                              ev0,
                                              phi_i,
                                              sigma_new,
                                              ep_new,
                                              X_new,
                                              Zeta_new,
                                              coher_new );

      if( substepFlag != 0 )
      {
        attemptFailed = true;
        break;
      }

      if( substep + 1 < totalSubsteps )
      {
        X_old = X_new;
        Zeta_old = Zeta_new;
        coher_old = coher_new;
        LvArray::tensorOps::copy< 6 >( sigma_old, sigma_new );
        LvArray::tensorOps::copy< 6 >( ep_old, ep_new );
      }
    }

    if( !attemptFailed )
    {
      Zeta_p = Zeta_new;
      coher_p = coher_new;
      LvArray::tensorOps::copy< 6 >( sigma_p, sigma_new );
      LvArray::tensorOps::copy< 6 >( ep_p, ep_new );

      real64 const evp_p = LvArray::tensorOps::symTrace< 3 >( ep_new );
      porosity_p = 1.0 + exp( -evp_p ) * ( phi_i - 1.0 );
      porosity_p = std::min( 1.0, std::max( 0.0, porosity_p ) );

      if( chi > 1 && constitutiveUpdateFlag == 0 )
      {
        constitutiveUpdateFlag = 1;
      }

      return 0;
    }

    if( chi >= chimax )
    {
      break;
    }

    int const nextChi = std::min( 2 * chi, chimax );
    if( nextChi <= chi )
    {
      break;
    }
    chi = nextChi;
  }

  constitutiveUpdateFlag = -1;
  return failStep( Zeta_n,
                   coher_n,
                   porosity_n,
                   sigma_n,
                   ep_n,
                   Zeta_p,
                   coher_p,
                   porosity_p,
                   sigma_p,
                   ep_p );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeHypoelasticTrialStress( const real64 dt,                    // time step
                                                         const real64 bulk,                  // bulk modulus
                                                         const real64 shear,                 // shear modulus
                                                         real64 const ( & stress_old )[6],   // stress at start of step.
                                                         real64 const ( & D )[6],            // D=sym(L)
                                                         real64 ( & stress_new )[6]         // stress at end of step
) const
{
  // Hypoelastic, isotropic, linear elasticity.
  real64 Diso[6] = { 0 }, 
         Ddev[6] = { 0 };

  // isotropicDeviatoricDecomposition( D, Diso, Ddev );
  LvArray::tensorOps::symAddIdentity< 3 >( Diso, LvArray::tensorOps::symTrace< 3 >( D )/3.0 );
  LvArray::tensorOps::copy< 6 >( Ddev, D );
  LvArray::tensorOps::subtract< 6 >( Ddev, Diso );

  // Cauchy trial stress:  sigma_n+1^trial = sigma_n + [ 3*K*iso(D) + 2*G*dev(D) ]*dt
  // Diso *= 3.*bulk*dt; // this is now the increment in isotropic stress
  // Ddev *= 2.*shear*dt; // this is now the increment in deviatoric stress
  // CC: TODO check if need to double off diagonal components for symmetric calculations
  LvArray::tensorOps::scale< 6 >( Diso, 3.*bulk*dt );
  LvArray::tensorOps::scale< 6 >( Ddev, 2.*shear*dt );

  // stress_new = stress_old;
  // stress_new += Diso;
  // stress_new += Ddev;
  LvArray::tensorOps::copy< 6 >( stress_new, stress_old );
  LvArray::tensorOps::add< 6 >( stress_new, Diso );
  LvArray::tensorOps::add< 6 >( stress_new, Ddev );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeElasticProperties( real64 & bulk,
                                                    real64 & shear 
) const
{
    // When computeElasticProperties() is called with two real64s as arguments, it
    // computes the high pressure limit tangent elastic shear and bulk modulus
    // This is used to estimate wave speeds and make conservative estimates of substepping.
    bulk  = m_b0 + m_b1;  // Bulk Modulus

    shear = m_g0;  // Default behavior is constant shear modulus   
    if( !isZero( m_g1 ) ) // Poisson ratio control.
    {
      real64 nu = m_g1 + m_g2; // high-pressure limit.

      // Force 0 <= nu < 0.5.
      nu = std::min( std::max( nu, 0.0 ), 0.499999 );
      shear = 1.5 * bulk * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
    }

    shear = fmax( shear, m_g0 );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeElasticProperties( real64 const ( &stress )[6],
                                                    real64 const ( &ep )[6],
		                                                const real64 & GEOS_UNUSED_PARAM( p3_crush_curve ),
		                                                const real64 & GEOS_UNUSED_PARAM( fluid_pressure_initial ),
		                                                const real64 & Km,             // Matrix bulk modulus
		                                                const real64 & Kf,             // Fluid bulk modulus
		                                                const real64 & C1,             // Term to simplify the fluid model expressions
		                                                const real64 & ev0,            // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                                const real64 & phi_i,			    // Initial porosity (inferred from crush curve, used for fluid model/
		                                                real64 & bulk,
		                                                real64 & shear
) const
{
  // Compute the nonlinear elastic tangent stiffness as a function of the pressure
  // plastic strain, and fluid parameters.
	real64 I1 = LvArray::tensorOps::symTrace< 3 >( stress ), //stress.Trace(),
		     evp = LvArray::tensorOps::symTrace< 3 >( ep ); //ep.Trace();

  // ..........................................................Undrained
  // The low pressure bulk and shear moduli are also used for the tensile response.
	bulk = m_b0;

	if( evp <= 0.0 )
    {
		if ( I1 < -1.e-12 )
    {
			real64 expb2byI1 = exp( m_b2 / I1 );
			bulk = bulk + m_b1 * expb2byI1;
		}

		// Elastic-plastic coupling
		if ( evp < -1.e-12 )
        {
            bulk = bulk - m_b3 * exp( m_b4 / evp );
        }
	}



    // In  compression, or with fluid effects if the strain is more compressive
    // than the zero fluid pressure volumetric strain:
	if ( evp <= ev0 && !isZero( Kf ) )
    {   // ..........................................................Undrained
		// Compute the porosity from the strain using Homel's simplified model, and
		// then use this in the Biot-Gassmann formula to compute the bulk modulus.

		// The dry bulk modulus, taken as the low pressure limit of the nonlinear
		// formulation:
		real64 Kd = m_b0;
		if ( evp < 0.0 )
        { 
            Kd = m_b0 - m_b3 * exp( m_b4 / evp );
        }

		// Current unloaded porosity (phi):
		real64 C2 = std::exp( evp * Km / C1 ) * phi_i;
		real64 phi = C2 / ( -std::exp( evp * Kf / C1 ) * ( phi_i - 1.0) + C2 );

		// Biot-Gassmann formula for the saturated bulk modulus, evaluated at the
		// current porosity.  This introduces some error since the Kd term is a
		// function of the initial porosity, but since large strains would also
		// modify the bulk modulus through damage
		real64 oneminusKdbyKm = 1.0 - Kd / Km;
		bulk = Kd + oneminusKdbyKm * oneminusKdbyKm / ( ( oneminusKdbyKm - phi ) / Km + ( 1.0 / Kf - 1.0 / Km ) * phi );
	}

  // don't allow elastic-plastic coupling to bring bulk modulus too low:
  bulk = fmax(bulk, m_b0);
  
  // To be thermodynamically consistent, the shear modulus in an isotropic model
	// must be constant, but the bulk modulus can depend on pressure.  However, this
	// leads to a Poisson's ratio that approaches 0.5 at high pressures, which is
	// inconsistent with experimental data for the Poisson's ratio, inferred from the
	// Young's modulus.  Induced anisotropy is likely the cause of the discrepency,
	// but it may be better to allow the shear modulus to vary so the Poisson's ratio
	// remains reasonable.
	//
	// If the user has specified a nonzero value of G1, the shear modulus will be defined
  // from a poisson's ratio nu = g1+g2*exp(g3/I1)
	// to vary with pressure so the drained Poisson's ratio transitions from G1+G2 to G1 as
	// the pressure increases relative to g3.  Treatment of fluid effects has not yet been developed.
  
  shear = m_g0;  // Default behavior is constant shear modulus
  
  if( !isZero( m_g1 ) ) // Poisson ratio control.
  {
    real64 nu = m_g1;
    if ( I1 < -1.e-12 ) // in compression scale the poisson ratio
    {
      nu += m_g2 * exp( m_g3 / I1 );
    }
    // Force 0 <= nu < 0.5.
    nu = std::min( std::max( nu, 0.0 ), 0.499999 );
		shear = 1.5 * bulk * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
	}

  shear = fmax(shear, m_g0);
}

// [nsub] = computeStepDivisions(X,Zeta,ep,sigma_n,sigma_trial)
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeStepDivisions( const real64 & X,
		                                           const real64 & GEOS_UNUSED_PARAM( Zeta ),
		                                           real64 const ( & ep )[6],
		                                           real64 const ( & sigma_n )[6],
		                                           real64 const ( & sigma_trial )[6],
		                                           const real64 & GEOS_UNUSED_PARAM( fluid_B0 ),
		                                           const real64 & fluid_pressure_initial,
		                                           const real64 & Km,          // Matrix bulk modulus
		                                           const real64 & Kf,          // Fluid bulk modulus
		                                           const real64 & C1,          // Term to simplify the fluid model expressions
		                                           const real64 & ev0,         // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                           const real64 & phi_i,        // Initial porosity (inferred from crush curve, used for fluid model/
                                                   int & constitutiveUpdateFlag // Set to 1 if requested subcycles are capped.
) const
{
  // compute the number of step divisions (substeps) based on a comparison
  // of the trial stress relative to the size of the yield surface, as well
  // as change in elastic properties between sigma_n and sigma_trial.
  int const nmax = std::max( 1, m_maxAllowedSubcycles );

  // Compute change in bulk modulus:
  real64 bulk_n,
         shear_n,
         bulk_trial,
         shear_trial;

  computeElasticProperties( sigma_n,
                            ep,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk_n,
                            shear_n );

  computeElasticProperties( sigma_trial,
                            ep,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk_trial,
                            shear_trial );
  int n_bulk = ceil( fabs( bulk_n - bulk_trial ) / bulk_n );

  // Compute trial stress increment relative to yield surface size:
  real64 d_sigma[6] = { 0 };
  // d_sigma = sigma_trial; //sigma_trial - sigma_n;
  // d_sigma -= sigma_n;
  LvArray::tensorOps::copy< 6 >( d_sigma, sigma_trial );
  LvArray::tensorOps::subtract< 6 >( d_sigma, sigma_n );

  real64 size = fabs( 0.5*( m_peakI1 - X ) );

  if( m_stren > 0.0 )
  {
	  size = std::min( size, m_stren );
  }

  real64 const sizeFloor = 1.0e-16 * std::max( 1.0, fabs( m_peakI1 ) + fabs( X ) + fabs( m_stren ) );
  size = std::max( size, sizeFloor );

  d_sigma[3] *= 2.0;
  d_sigma[4] *= 2.0;
  d_sigma[5] *= 2.0;
  int n_yield = int( std::ceil( 1.0e-4 * LvArray::tensorOps::l2Norm< 6 >( d_sigma ) / size ));
  d_sigma[3] *= 0.5;
  d_sigma[4] *= 0.5;
  d_sigma[5] *= 0.5;

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = std::max( n_bulk, n_yield);

  //if( nsub > nmax && ( m_failedStepResponse == 2 || m_failedStepResponse == 3 ) )
  if( nsub > nmax )
  {
    // Response 0/1 accepts a capped subcycle count and lets the return mapping
    // convergence decide whether the step is usable. Response 2/3 reports
    // failure to the caller. Logging is handled outside device-callable code.
    if( m_failedStepResponse <= 1 )
    {
      constitutiveUpdateFlag = 1;
      nsub = nmax;
    }
    else
    {
      nsub = -1;
    }
  }
  else
  {
    nsub = std::min( std::max( nsub, 1 ), nmax  );
  }
  return nsub;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeInvariants( real64 const ( & stress )[6],
                                             real64 ( & S )[6],
                                             real64 & I1,
                                             real64 & J2,
                                             real64 & rJ2 
) const
{
  // Compute the first invariants
//   I1 = LvArray::tensorOps::symTrace< 3 >( stress ); //stress.Trace();  //Pa

  // Compute the deviatoric part of the tensor
//   S = dev(stress); // stress - one_third*Identity*I1

//   real64 deviator[6] = { 0 };
  real64 vonMisesStress;
  twoInvariant::stressDecomposition( stress,
                                     I1,
                                     vonMisesStress, //this is actually von Mises
                                     S); //this gives a unit verctor
    I1 *= 3.0;
    LvArray::tensorOps::scale< 6 >( S, sqrt(2.0 / 3.0) * vonMisesStress ); //Stress tensor in voight notation.
    rJ2 = vonMisesStress/sqrt(3.);
    J2 = rJ2 * rJ2; 

  if( J2 < 1e-16*( I1 * I1 + J2 ) )
  {
    J2 = 0.0;
    rJ2 = 0.0;
  }
}

// Derivative of Ff*Fc with respect to I1, used to find the yield surface APEX for an 
// automatic estimate of the brittle-ductile transition
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GeomechanicsUpdates::derivative_f(const real64 & peakI1,
                                    const real64 & cr,
                                    const real64 & X,
									const real64 & a1,
									const real64 & a2,
									const real64 & a3,
                                    const real64 & a4,
									const real64 & I1
                                    ) const    
{
    real64 f_prime;
    const real64 Kappa = peakI1 + cr * (X - peakI1);
    //const real64 F =(a1 - a3) * (exp(a2 * I1)) - a4 * I1 * sqrt (1 - pow((-(I1) - (-Kappa)) / ((-X) - (-Kappa)), 2));
    // Compute intermediate values
    real64 exp_term = exp(a2 * I1);  // e^(a2 * I1)
    real64 numerator = (-I1 + Kappa);
    real64 denominator = (X - Kappa);
    real64 fraction = numerator / denominator;
    
    // Compute square root term
    real64 sqrt_term = sqrt(1 - std::pow(fraction, 2));
    
    // Compute g'(I1)
    real64 g_prime = -a3 * a2 * exp_term - a4;
    
    // Compute h'(I1)
    real64 h_prime = (-1.0 / denominator) / sqrt_term;
    
    // Compute g(I1)
    real64 g_I1 = a1 - a3 * exp_term - a4 * I1;

    // Final derivative calculation using the product rule
    f_prime = (g_prime * sqrt_term) + (g_I1 * h_prime);

    return f_prime;

}


//find the boundary between brittle ductile?
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE

real64 GeomechanicsUpdates::computeBD(const real64 & a1,
                                      const real64 & a2,
                                      const real64 & a3,
                                      const real64 & a4
) const    
{
   // TODO: This may need to be changed for general use if the 
   // value of m_peakI1 has been modified by coher or hardening
   // before these
   // values of a1,a2,a3,a4 were computed.  For now we just use this to decide
   // whether inelastic deformation should trigger damage evolution so 
   // we'll stick to 

   // we want initial value so setting X to p0
   real64 X = m_p0;

   const real64 kappa = m_peakI1 + m_cr * (X - m_peakI1);
   //const real64 F =(a1 - a3) * (exp(a2 * I1)) - a4 * I1 * sqrt (1 - pow((-(I1) - (-Kappa)) / ((-X) - (-Kappa)), 2));
    real64 x_0 = kappa;
    real64 x_1 = m_p0;
    real64 tolerance = 1.e-6;

    real64 x = 0.0;
    real64 check_val;

    while (std::abs(x_0 - x_1) > tolerance) {
        x = 0.5 * (x_0 + x_1);  // Compute the midpoint
        check_val = derivative_f(m_peakI1,m_cr,X,a1,a2,a3,a4,x);
        if ( check_val > 0) {
            x_1 = x;  // Update x1
        } else {
            x_0 = x;  // Update x0
        }
    }
 return x;
} 


// update the coherence (1-damage) variable based on dilational plastic work
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeCoher( const real64 & a1,
                                        const real64 & a2,
                                        const real64 & a3,
                                        const real64 & a4,
                                        const real64 & d_evp,         // increment in vol plastic strain
                                        real64 const ( & d_ep )[6],   // increment in total plastic strain
                                        const real64 & I1_trial,      // trial value of I1
                                        const real64 & rJ2_trial,  // trial value of rootJ2
			                                  const real64 & I1_0,          // I1 value on yield surface
                                        const real64 & rJ2_0,         // rJ2 value on yield surface
                                        const real64 & lch,           // length scale
                                        const real64 & strengthScale, // strength scale
                                        const real64 & coher_old,     // old coherence = 1-damage
                                              real64 & coher_new     // OUPUT: new value of coherconst real64 & lch,       // length scale
) const
{
	// material coherence (coher = 1 - damage) will remain constant unless damage criterion is met.
  coher_new = coher_old;

  real64 scaledFractureStress = m_fractureStress*strengthScale;

	if( m_fractureEnergyReleaseRate > 1.e-16 )  // energy release rate is 0 or negative to disable damage.
	{
		if ( m_damageEvolutionCriterion == 0 ) 
		{ // This approach is designed to increment damage when there is plastic loading with dilation above some
      // threshold stress.  The damage evolution is based on dilational plastic work.
      if (d_evp > 0 && ( I1_trial - I1_0 ) > 0 && ( (rJ2_0 > scaledFractureStress || coher_old < 1 ))   &&  m_damageEvolutionCriterion == 0)
      {
			// increment in work per unit volume.
			real64 d_dilationalPlasticWork = d_evp*0.5*(I1_trial - I1_0);
      // increment damage based on increment in plastic dilatational work relative to the fracture
      // energy release rate, normalized by the length scale, to be per-unit-volume     
			real64 d_damage = d_dilationalPlasticWork*lch / m_fractureEnergyReleaseRate; // forced to be positive since I1_trial>I1_0 and d_evp>0
      // force coher = 1-damage to be 0<=coher<=1
			coher_new = std::max( 0.0, coher_old - d_damage );
      }
		}
    else if ( m_damageEvolutionCriterion == 1 || m_damageEvolutionCriterion == 2 )
    {
      // Set the brittle-ductile transition pressure based on input value (1) or based on yield surface apex (2)
      real64 I_db = ( m_damageEvolutionCriterion == 1 ) ? -3.*m_brittleDuctileTransition : computeBD(a1, a2, a3, a4);

      // If pressure (p=-I1/3) is below the brittle ductile transition pressure, increment damage based on plastic work increment
      // relative to fracture energy regularized by length scale.
      if ( I1_0 > I_db && ( rJ2_0 > scaledFractureStress || coher_old < 1 ) )
      {
        // increment in work per unit volume.
        real64 d_dilationalPlasticWork = std::max(0., d_evp*0.5*(I1_trial - I1_0) );
        real64 d_shearPlasticWork = std::max(0., (rJ2_trial - rJ2_0)*std::sqrt( (2./3.)*( 
                                                        (d_ep[0]-d_ep[1])*(d_ep[0]-d_ep[1]) + 
                                                        (d_ep[0]-d_ep[2])*(d_ep[0]-d_ep[2]) + 
                                                        (d_ep[1]-d_ep[2])*(d_ep[1]-d_ep[2]) + 
                                                        d_ep[3]*d_ep[3] + d_ep[4]*d_ep[4] + d_ep[5]*d_ep[5] 
                                                        ) ) );
        // increment damage based on increment in plastic dilatational work relative to the fracture
        // energy release rate, normalized by the length scale, to be per-unit-volume     
        real64 d_damage = ( d_dilationalPlasticWork + d_shearPlasticWork ) * lch  / m_fractureEnergyReleaseRate;
        // force coher = 1-damage to be 0<=coher<=1
        coher_new = std::max( 0.0, coher_old - d_damage );
      }
		}
	}
}

// Compute state variable X, the Hydrostatic Compressive strength (cap position)
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GeomechanicsUpdates::computeX( const real64 & evp,
		                                  const real64 & phi_i, // Initial porosity (inferred from crush curve, used for fluid model/
		                                  const real64 & Km, // matrix bulk modulus
		                                  const real64 & Kf, // Fluid bulk modulus
		                                  const real64 & C1, // Term to simplify the fluid model expressions
		                                  const real64 & ev0, // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                  const real64 & fluid_pressure_initial,
                                      const real64 & buckling  // crush strength multiplier for cascading collapse.
) const
{
  // X is the value of (I1 - Zeta) at which the cap function crosses
  // the hydrostat. For the drained material in compression. X(evp)
  // is derived from the emprical Kayenta crush curve, but with p2 = 0.
  // In tension, M. Homel's piecewise formulation is used.

  // define and initialize some variables
  real64 X;
  real64 identity[6] = {1.0,1.0,1.0,0.0,0.0,0.0};
  // R2Tensor Identity = R2Tensor( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  if( evp <= -m_p3 ) { // ------------Plastic strain exceeds allowable limit--------------------------
    // The plastic strain for this iteration has exceed the allowable
    // value.  X is not defined in this region, so we set it to a large
    // negative number.
    //
    // The code should never have evp<-p3, but will have evp=-p3 if the
    // porosity approaches zero (within the specified tolerance).  By setting
    // X=1e12*p0, the material will respond as though there is no porosity.
    X = 1.0e12 * m_p0;
  }
  else { // ------------------Plastic strain is within allowable domain------------------------
    // We first compute the drained response.  If there are fluid effects, this value will
    // be used in determining the elastic volumetric strain to yield.
    if( evp <= 0.0 ){
      X = ( m_p0 * m_p1 + log( ( evp + m_p3 ) / m_p3 ) ) / m_p1;
    }
    else{
      X = m_p0 * std::pow( 1.0 + evp, 1.0 / ( m_p0 * m_p1 * m_p3 ) );
    }

    if( !isZero( Kf ) && evp <= ev0 ) { // ------------------------------------------- Fluid Effects
      // First we evaluate the elastic volumetric strain to yield from the
      // empirical crush curve (Xfit) and bulk modulus (Kfit) formula for
      // the drained material.  Xfit was computed as X above.

      // Kfit is the drained bulk modulus evaluated at evp, and for I1 = Xdry/2.
      real64 Kdry = m_b0 + m_b1 * exp( 2.0 * m_b2 / X );
      if ( evp < 0.0 )
      {
        Kdry = Kdry - m_b3 * exp( m_b4 / evp );
      }

      // Now we use our engineering model for the bulk modulus of the
      // saturated material (Keng) to compute the stress at our elastic strain to yield.
      // Since the stress and plastic strain tensors are not available in this scope, we call the
      // computeElasticProperties function with isotropic matrices that will have the
      // correct values of evp. (The saturated bulk modulus doesn't depend on I1).
      real64 Ksat, 
             Gsat;       // Not used, but needed to call computeElasticProperties()
      // This needs to be evaluated at the current value of pressure.

      real64 hydroStress[6] = { 0 }, 
             hydroStrain[6] = { 0 };

      // hydroStress = Identity; //(1./6.)*X*Identity;
      // hydroStress *= ( 1.0 / 6.0 ) * X;
      LvArray::tensorOps::copy< 6 >( hydroStress, identity );
      LvArray::tensorOps::scale< 6 >( hydroStress, ( 1.0 / 6.0) * X );

      // hydroStrain = Identity; //(1./3.)*evp*Identity;
      // hydroStrain *= (1./3.) * evp;
      LvArray::tensorOps::copy< 6 >( hydroStrain, identity );
      LvArray::tensorOps::scale< 6 >( hydroStrain, ( 1.0 / 3.0) * evp );

      computeElasticProperties( hydroStress,
			                          hydroStrain,
			                          m_p3,
			                          fluid_pressure_initial,
                                Km,
                                Kf,
                                C1,
                                ev0,
                                phi_i,
			                          Ksat, 
                                Gsat ); // OUTPUT: saturated bulk and shear modulus
			
      // Compute the stress to hydrostatic yield.
      // We are only in this loop if(evp <= ev0)
      X = X * Ksat / Kdry;   // This is X_sat = K_sat*eve = K_sat*(X_dry/K_dry)
    } // End fluid effects
  } // End plastic strain in allowable domain
  return buckling*X;
}


// Computes the updated stress state for a substep
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeSubstep( real64 const ( & D )[6],         // strain "rate" for substep D=sym(L)
                                         const real64 & dt,               // substep time interval
                                         const real64 & lch,              // length scale
                                         const real64 & strengthScale,    // strength scale
                                         real64 const ( & sigma_old )[6], // stress at start of substep
                                         real64 const ( & ep_old )[6],    // plastic strain at start of substep
                                         const real64 & X_old,            // hydrostatic compressive strength at start of substep
                                         const real64 & Zeta_old,         // trace of isotropic backstress at start of substep
                                         const real64 & coher_old,        // scalar valued coherence = 1-Damage
                                         const real64 & buckling,         // scalar valued buckling crush-curve modifier
                                         const real64 & fluid_pressure_initial,
                                         const real64 & Km,
                                         const real64 & Kf,
                                         const real64 & C1,
                                         const real64 & ev0,
                                         const real64 & phi_i,
                                         real64 ( & sigma_new )[6],       // stress at end of substep
                                         real64 ( & ep_new )[6],          // plastic strain at end of substep
                                         real64 & X_new,                  // hydrostatic compressive strength at end of substep
                                         real64 & Zeta_new,               // trace of isotropic backstress at end of substep
                                         real64 & coher_new               // coherence at the end of substep
) const
{
  // Computes the updated stress state for a substep that may be either elastic, plastic, or
  // partially elastic. Returns an integer flag 0/1 for a good/bad update.
  int returnFlag = 0;

  real64 const one_third = 1.0 / 3.0;
  real64 identity[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // (1) Compute the elastic properties based on the stress and plastic strain at
  // the start of the substep.
  real64 bulk = 0.0;
  real64 shear = 0.0;
  computeElasticProperties( sigma_old,
                            ep_old,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk,
                            shear );

  // (3) Compute the trial stress: [sigma_trial] = computeTrialStress(sigma_old,d_e,K,G)
  real64 sigma_trial[6] = { 0.0 };
  real64 S_trial[6] = { 0.0 };
  computeHypoelasticTrialStress( dt,
                                 bulk,
                                 shear,
                                 sigma_old,
                                 D,
                                 sigma_trial );

  real64 I1_trial = 0.0;
  real64 J2_trial = 0.0;
  real64 rJ2_trial = 0.0;
  computeInvariants( sigma_trial,
                     S_trial,
                     I1_trial,
                     J2_trial,
                     rJ2_trial );

  // (4) Evaluate the yield function at the trial stress.
  real64 hardening = 0.0;
  if( m_strainHardeningK > 0.0 )
  {
    real64 const ep_J2 = ep_old[0] * ep_old[0] + ep_old[1] * ep_old[1] + ep_old[2] * ep_old[2]
                         + 3.0 * ep_old[3] * ep_old[3] + 3.0 * ep_old[4] * ep_old[4] + 3.0 * ep_old[5] * ep_old[5]
                         - ep_old[1] * ep_old[2] - ep_old[0] * ep_old[1] - ep_old[0] * ep_old[2];

    if( ep_J2 > 0.0 )
    {
      real64 const equilivantPlasticStrain = ( 2.0 / 3.0 ) * sqrt( ep_J2 );
      hardening = m_strainHardeningK * ( 1.0 - exp( -m_strainHardeningN * equilivantPlasticStrain ) );
    }
  }

  real64 a1 = 0.0;
  real64 a2 = 0.0;
  real64 a3 = 0.0;
  real64 a4 = 0.0;
  computeLimitParameters( a1,
                          a2,
                          a3,
                          a4,
                          coher_old,
                          hardening,
                          strengthScale,
                          buckling );

  int YIELD = computeYieldFunction( I1_trial,
                                    rJ2_trial,
                                    X_old,
                                    Zeta_old,
                                    coher_old,
                                    hardening,
                                    strengthScale,
                                    buckling,
                                    a1,
                                    a2,
                                    a3,
                                    a4 );

  if( m_geomechanicsDisableInelasticity )
  {
    YIELD = -1;
  }

  if( YIELD <= 0 )
  {
    X_new = X_old;
    Zeta_new = Zeta_old;
    coher_new = coher_old;
    LvArray::tensorOps::copy< 6 >( sigma_new, sigma_trial );
    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    return 0;
  }

  if( YIELD != 1 )
  {
    return failSubstep( sigma_old,
                        ep_old,
                        X_old,
                        Zeta_old,
                        coher_old,
                        sigma_new,
                        ep_new,
                        X_new,
                        Zeta_new,
                        coher_new );
  }

  // (5) Compute non-hardening return to the initial yield surface.
  real64 I1_0 = 0.0;
  real64 rJ2_0 = 0.0;
  real64 const TOL = std::max( m_plasticStrainTolerance, 1.0e-16 );
  real64 S_0[6] = { 0.0 };
  real64 d_ep_0[6] = { 0.0 };
  real64 S_old[6] = { 0.0 };
  real64 I1_old = 0.0;
  real64 J2_old = 0.0;
  real64 rJ2_old = 0.0;
  real64 const evp_old = LvArray::tensorOps::symTrace< 3 >( ep_old );
  computeInvariants( sigma_old,
                     S_old,
                     I1_old,
                     J2_old,
                     rJ2_old );

  real64 d_e[6] = { 0.0 };
  LvArray::tensorOps::copy< 6 >( d_e, D );
  LvArray::tensorOps::scale< 6 >( d_e, dt );

  returnFlag = nonHardeningReturn( I1_trial,
                                   rJ2_trial,
                                   S_trial,
                                   I1_old,
                                   rJ2_old,
                                   S_old,
                                   d_e,
                                   X_old,
                                   Zeta_old,
                                   coher_old,
                                   hardening,
                                   strengthScale,
                                   buckling,
                                   bulk,
                                   shear,
                                   I1_0,
                                   rJ2_0,
                                   S_0,
                                   d_ep_0 );

  if( returnFlag != 0 )
  {
    return failSubstep( sigma_old,
                        ep_old,
                        X_old,
                        Zeta_old,
                        coher_old,
                        sigma_new,
                        ep_new,
                        X_new,
                        Zeta_new,
                        coher_new );
  }

  // If there is no porosity (p3=0) and no fluid effects (Kf=0), the non-hardening
  // return is the solution.
  if( isZero( m_p3 ) && isZero( Kf ) )
  {
    Zeta_new = Zeta_old;
    X_new = X_old;
    coher_new = coher_old;

    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    LvArray::tensorOps::add< 6 >( ep_new, d_ep_0 );

    LvArray::tensorOps::copy< 6 >( sigma_new, identity );
    LvArray::tensorOps::scale< 6 >( sigma_new, one_third * I1_0 );
    LvArray::tensorOps::add< 6 >( sigma_new, S_0 );

    return 0;
  }

  real64 const d_evp_0 = LvArray::tensorOps::symTrace< 3 >( d_ep_0 );

  // (6) Iterate to solve for plastic volumetric strain consistent with the updated
  // values for the cap (X) and isotropic backstress (Zeta). Use a bounded bisection
  // on the multiplier eta, where 0 < eta < 1.
  real64 eta_out = 1.0;
  real64 eta_in = 0.0;
  real64 const dZetadevp = computedZetadevp( fluid_pressure_initial, Km, Kf, ev0, phi_i, Zeta_old, evp_old );
  int const imax = std::max( 93, int( std::ceil( std::log( TOL ) / std::log( 0.5 ) ) ) + 2 );
  bool forcedLowerBracket = false;

  for( int i = 1; i <= imax + 1; ++i )
  {
    real64 const eta_mid = 0.5 * ( eta_out + eta_in );
    real64 const d_evp = eta_mid * d_evp_0;

    // Update X exactly.
    real64 evp_new = evp_old + d_evp;
    X_new = computeX( evp_new,
                      phi_i,
                      Km,
                      Kf,
                      C1,
                      ev0,
                      fluid_pressure_initial,
                      buckling );

    // Update the damage variable. This is untested and only active if Gf > 0.
    computeCoher( a1,
                  a2,
                  a3,
                  a4,
                  d_evp,
                  d_ep_0,
                  I1_trial,
                  rJ2_trial,
                  I1_0,
                  rJ2_0,
                  lch,
                  strengthScale,
                  coher_old,
                  coher_new );

    // Update zeta. min() eliminates tensile fluid pressure from explicit integration error.
    Zeta_new = std::min( Zeta_old + dZetadevp * d_evp, 0.0 );

    real64 a1_new = 0.0;
    real64 a2_new = 0.0;
    real64 a3_new = 0.0;
    real64 a4_new = 0.0;
    computeLimitParameters( a1_new,
                            a2_new,
                            a3_new,
                            a4_new,
                            coher_new,
                            hardening,
                            strengthScale,
                            buckling );

    // (8) Check if the updated yield surface encloses trial stress. If it does not,
    // there is too much plastic strain for this iteration, so adjust the bisection bracket.
    if( computeYieldFunction( I1_trial,
                              rJ2_trial,
                              X_new,
                              Zeta_new,
                              coher_new,
                              hardening,
                              strengthScale,
                              buckling,
                              a1_new,
                              a2_new,
                              a3_new,
                              a4_new ) != 1 )
    {
      if( i >= imax )
      {
        if( forcedLowerBracket )
        {
          return failSubstep( sigma_old,
                              ep_old,
                              X_old,
                              Zeta_old,
                              coher_old,
                              sigma_new,
                              ep_new,
                              X_new,
                              Zeta_new,
                              coher_new );
        }
        eta_out = eta_in;
        forcedLowerBracket = true;
      }
      else
      {
        eta_out = eta_mid;
      }
      continue;
    }

    // (9) Recompute the non-hardening return to the updated surface.
    real64 I1_new = 0.0;
    real64 rJ2_new = 0.0;
    real64 d_evp_new = 0.0;
    real64 S_new[6] = { 0.0 };
    real64 d_ep_new[6] = { 0.0 };
    returnFlag = nonHardeningReturn( I1_trial,
                                     rJ2_trial,
                                     S_trial,
                                     I1_old,
                                     rJ2_old,
                                     S_old,
                                     d_e,
                                     X_new,
                                     Zeta_new,
                                     coher_new,
                                     hardening,
                                     strengthScale,
                                     buckling,
                                     bulk,
                                     shear,
                                     I1_new,
                                     rJ2_new,
                                     S_new,
                                     d_ep_new );

    if( returnFlag != 0 )
    {
      return failSubstep( sigma_old,
                          ep_old,
                          X_old,
                          Zeta_old,
                          coher_old,
                          sigma_new,
                          ep_new,
                          X_new,
                          Zeta_new,
                          coher_new );
    }

    // (10) Check whether the isotropic component of the return has changed sign, as this
    // would indicate that the cap apex has moved past the trial stress.
    real64 const sgnI1tmn = ( ( I1_trial - I1_new ) < 0.0 ) ? ( -1.0 ) : ( 1.0 );
    real64 const sgnI1tm0 = ( ( I1_trial - I1_0 ) < 0.0 ) ? ( -1.0 ) : ( 1.0 );
    if( fabs( sgnI1tmn - sgnI1tm0 ) > 1.0e-12 )
    {
      if( i >= imax )
      {
        if( forcedLowerBracket )
        {
          return failSubstep( sigma_old,
                              ep_old,
                              X_old,
                              Zeta_old,
                              coher_old,
                              sigma_new,
                              ep_new,
                              X_new,
                              Zeta_new,
                              coher_new );
        }
        eta_out = eta_in;
        forcedLowerBracket = true;
      }
      else
      {
        eta_out = eta_mid;
      }
      continue;
    }

    // Compare magnitude of plastic strain with prior update.
    d_evp_new = LvArray::tensorOps::symTrace< 3 >( d_ep_new );
    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    LvArray::tensorOps::add< 6 >( ep_new, d_ep_new );

    if( fabs( eta_out - eta_in ) < TOL )
    {
      LvArray::tensorOps::copy< 6 >( sigma_new, identity );
      LvArray::tensorOps::scale< 6 >( sigma_new, one_third * I1_new );
      LvArray::tensorOps::add< 6 >( sigma_new, S_new );

      // If out of range, scale back isotropic plastic strain.
      if( LvArray::tensorOps::symTrace< 3 >( ep_new ) < -m_p3 )
      {
        real64 ep_new_iso[6] = { 0.0 };
        real64 ep_new_dev[6] = { 0.0 };
        LvArray::tensorOps::copy< 6 >( ep_new_iso, identity );
        LvArray::tensorOps::scale< 6 >( ep_new_iso, one_third * LvArray::tensorOps::symTrace< 3 >( ep_new ) );
        LvArray::tensorOps::copy< 6 >( ep_new_dev, ep_new );
        LvArray::tensorOps::subtract< 6 >( ep_new_dev, ep_new_iso );

        // Force value to be maximum limit evp=-p3.
        evp_new = -m_p3;
        d_evp_new = evp_new - LvArray::tensorOps::symTrace< 3 >( ep_old );

        LvArray::tensorOps::copy< 6 >( ep_new_iso, identity );
        LvArray::tensorOps::scale< 6 >( ep_new_iso, one_third * evp_new );
        LvArray::tensorOps::copy< 6 >( ep_new, ep_new_dev );
        LvArray::tensorOps::add< 6 >( ep_new, ep_new_iso );
      }

      // Update X exactly.
      real64 const evp = LvArray::tensorOps::symTrace< 3 >( ep_new );
      X_new = computeX( evp,
                        phi_i,
                        Km,
                        Kf,
                        C1,
                        ev0,
                        fluid_pressure_initial,
                        buckling );

      Zeta_new = std::min( Zeta_old + dZetadevp * d_evp_new, 0.0 );
      return 0;
    }

    if( i >= imax )
    {
      return failSubstep( sigma_old,
                          ep_old,
                          X_old,
                          Zeta_old,
                          coher_old,
                          sigma_new,
                          ep_new,
                          X_new,
                          Zeta_new,
                          coher_new );
    }

    // (11) Compare magnitude of the volumetric plastic strain and bisect on eta.
    if( fabs( d_evp_new ) > eta_mid * fabs( d_evp_0 ) )
    {
      eta_in = eta_mid;
    }
    else
    {
      eta_out = eta_mid;
    }
  }

  return failSubstep( sigma_old,
                      ep_old,
                      X_old,
                      Zeta_old,
                      coher_old,
                      sigma_new,
                      ep_new,
                      X_new,
                      Zeta_new,
                      coher_new );
}

// Compute nonhardening return from trial stress to some yield surface
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::nonHardeningReturn( const real64 & I1_trial,              // Trial Stress
                                             const real64 & rJ2_trial,          
                                             real64 const ( & S_trial )[6],          
                                             const real64 & I1_old,                // Stress at start of subtep
                                             const real64 & GEOS_UNUSED_PARAM( rJ2_old ),          
                                             real64 const ( & S_old )[6],          
                                             real64 const ( & d_e )[6],            // increment in total strain
                                             const real64 & X,                     // cap position
                                             const real64 & Zeta,                  // isotropic bacstress
                                             const real64 & coher,    
                                             const real64 & hardening,
                                             const real64 & strengthScale,      
                                             const real64 & buckling,    
                                             const real64 & bulk,                  // elastic bulk modulus
                                             const real64 & shear,                 // elastic shear modulus
                                             real64 & I1_new,                      // New stress state on yield surface
                                             real64 & rJ2_new,
                                             real64 ( & S_new )[6],
                                             real64 ( & d_ep_new )[6]            // increment in plastic strain for return
) const
{
  // Computes a non-hardening return to the yield surface in the meridional profile
  // (constant Lode angle) based on the current values of the internal state variables
  // and elastic properties.  Returns the updated stress and  the increment in plastic
  // strain corresponding to this return.
  //
  // NOTE: all values of r and z in this function are transformed!

	real64 one_sqrt_three = 0.5773502691896258;
	real64 sqrt_two = 1.414213562373095;
	real64 one_third = 0.3333333333333333;
	real64 sqrt_three = 1.732050807568877;
	real64 one_sixth = 0.1666666666666667;
	real64 one_ninth = 0.1111111111111111;

	real64 identity[6] = {1.0,1.0,1.0,0.0,0.0,0.0};

  const int nmax = 19;  // If this is changed, more entries may need to be added to sinV cosV.
  int n = 0,
	    returnFlag = 0,   // error flag = 0 for successful return.
      interior;

  // (1) Define an interior point, (nominally I1_0 = Zeta, also, J2_0 = 0 but no need to  create this variable.)
  real64 I1_0,
		     I1trialMinusZeta = I1_trial-Zeta;

  // The following are the input parameters, modified by hardening and or damage
  //  Any change to peakI1 should be copied in the non-hardening return and
  //  yield function updates, which have branch points based on the peakI1 value.
  real64 nonlinearCoher = std::pow(coher,m_fractureSofteningExponent);
  real64 peakI1_h = buckling*( m_peakI1 + hardening*m_dpeakI1dh )*nonlinearCoher*strengthScale; // dpeakI1dh has units of stress

  // It may be better to use an interior point at the center of the yield surface, rather than at zeta, in particular
  // when PEAKI1=0.  Picking the midpoint between PEAKI1 and X would be problematic when the user has specified
  // some no porosity condition (e.g. p0=-1e99)
  if( I1trialMinusZeta>= peakI1_h ) // Trial is past vertex
  { 
	  real64 lTrial = sqrt(I1trialMinusZeta * I1trialMinusZeta + rJ2_trial * rJ2_trial),
			     lYield = 0.5 * (peakI1_h - X);
	  I1_0 = Zeta + peakI1_h - std::min(lTrial, lYield);
  }
  else if( (I1trialMinusZeta < peakI1_h) && (I1trialMinusZeta > X) ){ // Trial is above yield surface
	  I1_0 = I1_trial;
  }
  else if( I1trialMinusZeta <= X ) // Trial is past X, use yield midpoint as interior point
  {
	  I1_0 = Zeta + 0.5 * (peakI1_h + X);
  }
  else
  { // Shouldn't get here
	  I1_0 = Zeta;
  }

  // (2) Transform the trial and interior points as follows where beta defines the degree
  //  of non-associativity.
  // multiplier to compute Lode R to sqrt(J2)
  real64 rJ2_to_r = sqrt_two * m_beta * sqrt(1.5 * bulk / shear);
  // multiplier to compute sqrt(J2) to Lode R
  real64 r_to_rJ2 = 1.0 / rJ2_to_r;
  real64 r_trial = rJ2_to_r * rJ2_trial,
         z_trial = I1_trial * one_sqrt_three,
         z_test,
         r_test,
         r_0 = 0.0,
         z_0 = I1_0 * one_sqrt_three;

  // Lookup tables for computing the sin() and cos() of th rotation angle.
  real64 sinV[]={0.7071067811865475,-0.5,0.3420201433256687,-0.2306158707424402,0.1545187928078405,
                  -0.1032426220806015,0.06889665647555759,-0.04595133277786571,0.03064021661344469,
                  -0.02042858745187096,0.01361958465478159,-0.009079879062402308,0.006053298918749807,
                  -0.004035546304539714,0.002690368259933135,-0.001793580042002626,0.001195720384163988,
                  -0.0007971470283055577,0.0005314313834717263,-0.00035428759824575,0.0002361917349088998};
  real64 cosV[]={0.7071067811865475,0.8660254037844386,0.9396926207859084,0.9730448705798238,
                  0.987989849476809,0.9946562024066014,0.9976238022052647,0.9989436796015769,
                  0.9995304783376449,0.9997913146325693,0.999907249155556,0.9999587770484402,
                  0.9999816786182636,0.999991857149859,0.9999963809527642,0.9999983915340229,
                  0.9999992851261259,0.9999996822782572,0.9999998587903324,0.9999999372401469,
                  0.9999999721067318};
  real64 sinTheta = sinV[0],
         cosTheta = cosV[0];

  // Compute the a1,a2,a3,a4 parameters from FSLOPE,YSLOPE,STREN and PEAKI1,
  // which are perturbed by variability according to coher.  These are then
  // passed down to the computeYieldFunction, to avoid the expense of computing a3
  real64 a1 = 0.0;
  real64 a2 = 0.0;
  real64 a3 = 0.0;
  real64 a4 = 0.0;         
  computeLimitParameters( a1,
                          a2,
                          a3,
                          a4,
                          coher,
                          hardening,
                          strengthScale,
                          buckling );

  // (3) Perform Bisection between in transformed space, to find the new point on the
  //  yield surface: [znew,rnew] = transformedBisection(z0,r0,z_trial,r_trial,X,Zeta,K,G)
  int k = 0;
  while ( (n < nmax)&&(k < 10*nmax) )
  {
    // transformed bisection to find a new interior point, just inside the boundary of the
    // yield surface.  This function overwrites the inputs for z_0 and r_0
    //  [z_0,r_0] = transformedBisection(z_0,r_0,z_trial,r_trial,X_Zeta,bulk,shear)
    transformedBisection( z_0,
                          r_0,
                          z_trial,
                          r_trial,
                          X,
                          Zeta,
                          coher,
                          hardening,
                          strengthScale,
                          buckling,
                          a1,
                          a2,
                          a3,
                          a4,
                          r_to_rJ2 );

    // (4) Perform a rotation of {z_new,r_new} about {z_trial,r_trial} until a new interior point
    // is found, set this as {z0,r0}
    interior = 0;
    n = std::max(n-4,0);  //
    // (5) Test for convergence:
    while ( (interior==0)&&(n < nmax) )
    {
		  k++;
      
      // To avoid the cost of computing pow() to get theta, and then sin(), cos(),
      // we use a lookup table defined above by sinV and cosV.
      //
      // theta = pi_fourth*Pow(-two_third,n);
      // z_test = z_trial + cos(theta)*(z_0-z_trial) - sin(theta)*(r_0-r_trial);
      // r_test = r_trial + sin(theta)*(z_0-z_trial) + cos(theta)*(r_0-r_trial);
      sinTheta = sinV[n];
      cosTheta = cosV[n];
      z_test = z_trial + cosTheta*(z_0-z_trial) - sinTheta*(r_0-r_trial);
      r_test = r_trial + sinTheta*(z_0-z_trial) + cosTheta*(r_0-r_trial);

      if ( transformedYieldFunction( z_test,
                                     r_test,
                                     X,
                                     Zeta,
                                     coher,
                                     hardening,
                                     strengthScale,
                                     buckling,
                                     a1,
                                     a2,
                                     a3,
                                     a4,
                                     r_to_rJ2) == -1 ) 
      { // new interior point
        interior = 1;
        z_0 = z_test;
        r_0 = r_test;
      }
      else { n++; }
    }
  }

  if (k>=10*nmax)
  {
	  returnFlag = 1;
  }

// (6) Solution Converged, Compute Untransformed Updated Stress:
  I1_new = sqrt_three*z_0;
  rJ2_new = r_to_rJ2*r_0;

  LvArray::tensorOps::copy< 6 >( S_new, S_trial );
  if ( !isZero( rJ2_trial ) )
  {
	  // S_new = S_trial; //S_trial*rJ2_new/rJ2_trial;
	  // S_new *= rJ2_new/rJ2_trial;
    LvArray::tensorOps::scale< 6 >( S_new, rJ2_new / rJ2_trial );
  }
  // else
  // {
	//   S_new = S_trial;
  // }

  real64 sigma_new[6],
         sigma_old[6],
         d_sigma[6];
  // sigma_new = Identity; // I1_new*one_third*Identity + S_new;
  // sigma_new *= I1_new*one_third;
  // sigma_new += S_new;
  LvArray::tensorOps::copy< 6 >( sigma_new, identity );
  LvArray::tensorOps::scale< 6 >( sigma_new, I1_new*one_third );
  LvArray::tensorOps::add< 6 >( sigma_new, S_new );

  // sigma_old = Identity; // I1_old*one_third*Identity + S_old;
  // sigma_old *= I1_old*one_third;
  // sigma_old += S_old;
  LvArray::tensorOps::copy< 6 >( sigma_old, identity );
  LvArray::tensorOps::scale< 6 >( sigma_old, I1_old*one_third );
  LvArray::tensorOps::add< 6 >( sigma_old, S_old );

  // d_sigma = sigma_new; // sigma_new - sigma_old;
  // d_sigma -= sigma_old;
  LvArray::tensorOps::copy< 6 >( d_sigma, sigma_new );
  LvArray::tensorOps::subtract< 6 >( d_sigma, sigma_old );

  // (7) Compute increment in plastic strain for return:
  //  d_ep0 = d_e - [C]^-1:(sigma_new-sigma_old)
  real64 d_ee[6], 
         d_ee_2[6];

  // [C]^-1:d_sigma =  0.5*d_sigma/shear + (one_ninth/bulk - one_sixth/shear)*d_sigma.Trace()*Identity;
  // d_ee = Identity; 
  // d_ee *= (one_ninth/bulk - one_sixth/shear)*d_sigma.Trace();
  LvArray::tensorOps::copy< 6 >( d_ee, identity );
  LvArray::tensorOps::scale< 6 >( d_ee, (one_ninth / bulk - one_sixth / shear) * LvArray::tensorOps::symTrace< 3 >( d_sigma ) );

  // d_ee_2 = d_sigma;
  // d_ee_2 *= (0.5/shear);
  // d_ee += d_ee_2;
  LvArray::tensorOps::copy< 6 >( d_ee_2, d_sigma );
  LvArray::tensorOps::scale< 6 >( d_ee_2, 0.5 / shear );
  LvArray::tensorOps::add< 6 >( d_ee, d_ee_2 );

  // d_ep_new = d_e; // d_e - d_ee;
  // d_ep_new -= d_ee;
  LvArray::tensorOps::copy< 6 >( d_ep_new, d_e );
  LvArray::tensorOps::subtract< 6 >( d_ep_new, d_ee );

  return returnFlag;
}

// Computes bisection between two points in transformed space
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::transformedBisection(real64 & z_0,
                                               real64 & r_0,
                                               const real64 & z_trial,
                                               const real64 & r_trial,
                                               const real64 & X,
                                               const real64 & Zeta,
									                             const real64 & coher,
                                               const real64 & hardening,
                                               const real64 & strengthScale,
                                               const real64 & buckling,
									                             const real64 & a1,
									                             const real64 & a2,
									                             const real64 & a3,
									                             const real64 & a4,
                                               const real64 & r_to_rJ2
) const
{
  // Computes a bisection in transformed stress space between point sigma_0 (interior to the
  // yield surface) and sigma_trial (exterior to the yield surface).  Returns this new point,
  // which will be just outside the yield surface, overwriting the input arguments for
  // z_0 and r_0.

  // After the first iteration of the non-hardening return, the subsequent bisections will likely
  // converge with eta << 1.  It may be faster to put in some logic to try to start bisection
  // with tighter bounds, and only expand them to 0<eta<1 if the first eta_mid is too large.


  // (1) initialize bisection
  real64 eta_out = 1.0,
         eta_in  = 0.0,
         eta_mid,
         TOL = std::max( m_stressReturnTolerance, 1.0e-16 ),
         r_test,
         z_test;

  // (2) Test for convergence
  while (eta_out-eta_in > TOL){

    // (3) Transformed test point
    eta_mid = 0.5*(eta_out+eta_in);
    z_test = z_0 + eta_mid*(z_trial-z_0);
    r_test = r_0 + eta_mid*(r_trial-r_0);
    // (4) Check if test point is within the yield surface:
    if ( transformedYieldFunction( z_test,
                                   r_test,
                                   X,
                                   Zeta,
                                   coher,
                                   hardening,
                                   strengthScale,
                                   buckling,
                                   a1,
                                   a2,
                                   a3,
                                   a4,
                                   r_to_rJ2 ) !=1 )
    {
    	eta_in = eta_mid;
    }
    else
    {
    	eta_out = eta_mid;
    }
  }
  
  // (5) Converged, return {z_new,r_new}={z_test,r_test}
  z_0 = z_0 + eta_out*(z_trial-z_0); //z_0 = z_test;
  r_0 = r_0 + eta_out*(r_trial-r_0); //r_0 = r_test;
}

// computeTransformedYieldFunction from transformed inputs
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::transformedYieldFunction( const real64 & z,
                                                   const real64 & r,
                                                   const real64 & X,
                                                   const real64 & Zeta,
										                               const real64 & coher,
                                                   const real64 & hardening,
                                                   const real64 & strengthScale,
                                                   const real64 & buckling,
										                               const real64 & a1,
										                               const real64 & a2,
										                               const real64 & a3,
										                               const real64 & a4,
                                                   const real64 & r_to_rJ2
) const
{
  // Evaluate the yield criteria and return:
  //  -1: elastic
  //   0: on yield surface within tolerance
  //   1: plastic
	real64 sqrt_three = 1.732050807568877;
  
  // Untransformed values:
  real64 I1  = sqrt_three * z,
		     rJ2 = r_to_rJ2 * r;
  int    YIELD = computeYieldFunction( I1,
                                       rJ2,
                                       X,
                                       Zeta,
                                       coher,
                                       hardening,
                                       strengthScale,
                                       buckling,
                                       a1,
                                       a2,
                                       a3,
                                       a4 );
  return YIELD;
}

// computeYieldFunction from untransformed inputs
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeYieldFunction( const real64 & I1,
                                               const real64 & rJ2,
                                               const real64 & X,
                                               const real64 & Zeta,
									                             const real64 & coher,
                                               const real64 & hardening,
                                               const real64 & strengthScale,
                                               const real64 & buckling,
									                             const real64 & a1,
									                             const real64 & a2,
									                             const real64 & a3,
									                             const real64 & a4
) const    
{
   //std::cout<<"I1 = "<<I1<<", rJ2 = "<<rJ2<<", Zeta = "<<Zeta<<", coher = "<<coher<<", hardening = "<<hardening<<", a1 = "<<a1<<", a2 = "<<a2<<", a3 = "<<a3<<", a4 = "<<a4<<std::endl;

	// Evaluate the yield criteria and return:
	//  -1: elastic
	//   0: on yield surface within tolerance (not used)
	//   1: plastic
	//
	//                        *** Developer Note ***
	// THIS FUNCTION IS DEEP WITHIN A NESTED LOOP AND IS CALLED THOUSANDS
	// OF TIMES PER TIMESTEP.  EVERYTHING IN THIS FUNCTION SHOULD BE
	// OPTIMIZED FOR SPEED.
	//
	int YIELD = -1;
	real64 I1mZ = I1 - Zeta;    // Shifted stress to evalue yield criteria

  // Parameters modified by damage and/or hardening
  real64 peakI1_h, cr_h; //
  real64 nonlinearCoher = std::pow(coher,m_fractureSofteningExponent);
  peakI1_h = buckling*( m_peakI1 + hardening*m_dpeakI1dh )*nonlinearCoher*strengthScale;

  // Branch point after hardening with limits (not this doesn't affect brittle-ductile transition)
  // point used to determine damamage evolution.
  cr_h = std::min(0.99999999999, std::max( 1.e-10, m_cr * ( 1 + hardening*m_dcrdh/m_g0 ) ) );

	// --------------------------------------------------------------------
	// *** SHEAR LIMIT FUNCTION (Ff) ***
	// --------------------------------------------------------------------
	// Read input parameters to specify strength model
	real64  Ff;
	Ff = a1 - a3*exp(a2*I1mZ) - a4*I1mZ;

	// --------------------------------------------------------------------
	// *** Branch Point (Kappa) ***
	// --------------------------------------------------------------------
	real64  Kappa  = peakI1_h - cr_h*(peakI1_h-X); // Branch Point

	// --------------------------------------------------------------------
	// *** COMPOSITE YIELD FUNCTION ***
	// --------------------------------------------------------------------
	// Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
	// The elseif statements have nested if statements, which is not equivalent
	// to them having a single elseif(A&&B&&C)
	if( I1mZ < X )
	{//---------------------------------------------------(I1<X)
		YIELD = 1;
	}
	else if(( I1mZ < Kappa )&&( I1mZ >= X )) 
  {// ---------------(X<I1<kappa)
		// p3 is the maximum achievable volumetric plastic strain in compresson
		// so if a value of 0 has been specified this indicates the user
		// wishes to run without porosity, and no cap function is used, i.e. fc=1

		// **Elliptical Cap Function: (fc)**
		// fc = sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
		// faster version: fc2 = fc^2
		real64 fc2 = 1.0 - ((Kappa-I1mZ)/(Kappa-X))*((Kappa-I1mZ)/(Kappa-X));
		if(rJ2*rJ2 > Ff*Ff*fc2 )
		{
			YIELD = 1;
		}
	}
	else if(( I1mZ <= peakI1_h )&&( I1mZ >= Kappa ))
  { // -----(kappa<I1<PEAKI1)
		if(rJ2 > Ff) {
			YIELD = 1;
		}
	}
	else if( I1mZ > peakI1_h )
	{// --------------------------------(peakI1<I1)
    YIELD = 1;
	};

  //std::cout<<"Ff = "<<Ff<<", kappa = "<<Kappa<<", YIELD = "<<YIELD<<std::endl;

  return YIELD;
} 

// Compute (dZeta/devp) Zeta and vol. plastic strain
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GeomechanicsUpdates::computedZetadevp( real64 const & fluid_pressure_initial,
                                              real64 const & Km,
                                              real64 const & Kf,
                                              real64 const & ev0,
                                              real64 const & phi_i,
                                              real64 Zeta,
                                              real64 evp
) const
{
  // Computes the partial derivative of the trace of the
  // isotropic backstress (Zeta) with respect to volumetric
  // plastic strain (evp).
  real64 dZetadevp = 0.0;           // Evolution rate of isotropic backstress

  if (evp <= ev0 && !isZero( Kf ) ) { // .................................... Fluid effects are active
    real64 pfi = fluid_pressure_initial; // initial fluid pressure

    // This is an expensive calculation, but fasterexp() seemed to cause errors.
    dZetadevp = (3.0*exp(evp)*Kf*Km)/(exp(evp)*(Kf + Km)
                                      + exp(Zeta/(3.0*Km))*Km*(-1.0 + phi_i)
                                      - exp((3.0*pfi + Zeta)/(3.0*Kf))*Kf*phi_i);
  }
  return dZetadevp;
} 

// // Compute compaction buckling multiplier to crush curve.
// GEOS_HOST_DEVICE
// GEOS_FORCE_INLINE
// void GeomechanicsUpdates::computeBuckling(real64 & buckling,          // this is overwritten
//                                           const real64 & ev,          // buckling strain, could be evp, directional strain, log(J), etc. 
// 		                                      const real64 & lengthScale  // element
// ) const 
// {   
//     // m_enableBuckling
//     // m_bucklingLength
//     // m_bucklingAmplitude

//     // // Number of unit cells per element, round up to nearest integer
//     real64 beta = std::ceil( lengthScale / m_bucklingLength );

//     // // Compute the length-scale dependent compaction-buckling scalar used to 
//     // // modify crush and shear strength:
//     // real64 J = 1.; // volumetric or directional stretch
//     // if (m_enableBuckling == 1)
//     // { // Isotropic:
//     //   J = LvArray::tensorOps::determinant< 3 >( deformationGradient );
//     // }
//     // else if (m_enableBuckling == 2)
//     // { // Anisotropic:
//     //   real64 temp[3] = { 0 };
//     //   LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( temp, deformationGradient, materialDirection  );
//     //   J = LvArray::tensorOps::AiBi< 3 >( temp, materialDirection);
//     // }
//     // else{
//     //   GEOS_LOG_RANK( "unsupported buckling type: " << m_enableBuckling );
//     // }
//     // real64 ev = log(J);     // volumetric or directional strain

//     real64 ev0 = m_p0 / ( 2.*m_b0 );
//     real64 ev1 = ev0 - m_p3;
//     real64 normalizedStrain; // normalizedStrain

//     if (ev >= ev0)
//     {
//       normalizedStrain = 0.;
//     }
//     else if (ev0 > ev && ev > ev1)
//     {
//       normalizedStrain = (ev-ev0)/(ev1-ev0);
//     }
//     else
//     {
//       normalizedStrain = 1.0;
//     }
//     // strength scale multiplier to crush curve (0: complete losss of strength, 1: no effect)
//     real64 pi = 3.141592653589793;
//     buckling = 1.0 - m_bucklingAmplitude*pow( sin( -1.0*beta*pi*normalizedStrain) , 2 );
//     buckling = fmin(1.0,fmax(0.0,buckling));
// }

// Compute (dZeta/devp) Zeta and vol. plastic strain
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeLimitParameters( real64 & a1,
		                                              real64 & a2,
		                                              real64 & a3,
		                                              real64 & a4,
		                                              const real64 & coher,
                                                  const real64 & hardening,
                                                  const real64 & strengthScale,
                                                  const real64 & buckling
) const 
{ // Value of I1 at strength=0 (Perturbed by variability)
  // The shear limit surface is defined in terms of the a1,a2,a3,a4 parameters, but
  // the user inputs are the more intuitive set of FSLOPE. YSLOPE, STREN, and PEAKI1.

  // This routine computes the a_i parameters from the user inputs.  The code was
  // originally written by R.M. Brannon, with modifications by M.S. Swan.
  // harden peakI1 and stren together to not change slope:
  // The following are the input parameters, modified by hardening and or damage
  //  Any change to peakI1 should be copied in the non-hardening return and
  //  yield function updates, which have branch points based on the peakI1 value.

  // Parameters modified by damage and/or hardening
  real64 stren_h, peakI1_h, fSlope_h, ySlope_h;  
  real64 nonlinearCoher = std::pow(coher,m_fractureSofteningExponent);
  
  stren_h = buckling*(m_stren + hardening*m_dstrendh); //  dstrendh has units of stress
  fSlope_h = nonlinearCoher*m_fSlope*( 1 + hardening*m_dfslopedh ) + ( 1. - nonlinearCoher )*m_fSlopeFailed; // dfslopedh is unitless
  peakI1_h = buckling*( m_peakI1 + hardening*m_dpeakI1dh )*nonlinearCoher*strengthScale; // dpeakI1dh has units of stress
  ySlope_h = std::min( 0.99999*fSlope_h, m_ySlope );

  if (fSlope_h > 0.0 && peakI1_h >= 0.0 && isZero( m_stren ) && isZero( ySlope_h) )
  {// ----------------------------------------------Linear Drucker-Prager
    a1 = peakI1_h * fSlope_h;
    a2 = 0.0;
    a3 = 0.0;
    a4 = fSlope_h;
  }
  else if ( isZero( fSlope_h ) && isZero( peakI1_h ) && stren_h > 0.0 && isZero( ySlope_h ) )
  { // ------------------------------------------------------- Von Mises
    a1 = stren_h*nonlinearCoher;
    a2 = 0.0;
    a3 = 0.0;
    a4 = 0.0;
  }
  else if ( fSlope_h > 0.0 && isZero( ySlope_h ) && stren_h > 0.0 && isZero( peakI1_h ) )
  { // ------------------------------------------------------- 0 PEAKI1 to vonMises
    a1 = stren_h;
    a2 = fSlope_h / stren_h;
    a3 = stren_h;
    a4 = 0.0;
  }
  else if (fSlope_h > ySlope_h && ySlope_h > 0.0 && stren_h > ySlope_h*peakI1_h && peakI1_h >= 0.0)
  { // ------------------------------------------------------- Nonlinear Drucker-Prager
    a1 = stren_h;
    a2 = (fSlope_h-ySlope_h )/(stren_h - ySlope_h*peakI1_h);
    a3 = (stren_h-ySlope_h*peakI1_h)*std::exp(-a2*peakI1_h);
    a4 = ySlope_h ;
  }
  else
  {
    // Invalid combinations are rejected during host-side input validation where
    // possible. Leave a zero surface here so the caller returns a failed step
    // instead of using host I/O from device-callable code.
    a1 = 0.0;
    a2 = 0.0;
    a3 = 0.0;
    a4 = 0.0;
  }

  //std::cout<<"m_peakI1 = "<<m_peakI1<<", m_stren = "<<m_stren<<", m_ySlope = "<<m_ySlope<<", m_fSlope = "<<m_fSlope<<std::endl;
  //std::cout<<"peakI1_h = "<<peakI1_h<<", stren_h = "<<stren_h<<", ySlope_h = "<<ySlope_h<<", fSlope_h = "<<fSlope_h<<std::endl;
  //std::cout<<"a1 = "<<a1<<", a2 = "<<a2<<", a3 = "<<a3<<", a4 = "<<a4<<std::endl;
}


/**
 * @class Geomechanics
 *
 * Geomechanics material model.
 */
class Geomechanics : public SolidBase
{
public:

  /// @typedef Alias for GeomechanicsUpdates
  using KernelWrapper = GeomechanicsUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Geomechanics( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Geomechanics() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Geomechanics";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for tangent elastic bulk modulus parameter 0 
    static constexpr char const * b0String() { return "b0"; }

    /// string/key for tangent elastic bulk modulus parameter 1 
    static constexpr char const * b1String() { return "b1"; }

    /// string/key for tangent elastic bulk modulus parameter 2 
    static constexpr char const * b2String() { return "b2"; }

    /// string/key for tangent elastic bulk modulus parameter 3 
    static constexpr char const * b3String() { return "b3"; }

    /// string/key for tangent elastic bulk modulus parameter 4 
    static constexpr char const * b4String() { return "b4"; }

    /// string/key for constant for tuning hardened stren
    static constexpr char const * dstrendhString() { return "dstrendh"; }

    /// string/key for constant for tuning  hardened fslope
    static constexpr char const * dfslopedhString() { return "dfslopedh"; }

    /// string/key for constant for tuning hardened peakI1
    static constexpr char const * dpeakI1dhString() { return "dpeakI1dh"; }

    /// string/key for constant for tuning hardened peakI1
    static constexpr char const * dcrdhString() { return "dcrdh"; }

    /// string/key for tangent elastic shear modulus parameter 0 
    static constexpr char const * g0String() { return "g0"; }

    /// string/key for tangent elastic shear modulus parameter 1 
    static constexpr char const * g1String() { return "g1"; }

    /// string/key for tangent elastic shear modulus parameter 2 
    static constexpr char const * g2String() { return "g2"; }

    /// string/key for tangent elastic shear modulus parameter 3 
    static constexpr char const * g3String() { return "g3"; }

    /// string/key for tangent elastic shear modulus parameter 4 
    static constexpr char const * g4String() { return "g4"; }

    /// string/key for crush curve parameter 0 
    static constexpr char const * p0String() { return "p0"; }

    /// string/key for crush curve parameter 1 
    static constexpr char const * p1String() { return "p1"; }

    /// string/key for crush curve parameter 2 
    static constexpr char const * p2String() { return "p2"; }

    /// string/key for crush curve parameter 3 
    static constexpr char const * p3String() { return "p3"; }

    /// string/key for crush curve parameter 4 
    static constexpr char const * p4String() { return "p4"; }

    /// string/key for cap shape parameter
    static constexpr char const * crString() { return "cr"; }

    /// string/key for fluid bulk modulus
    static constexpr char const * fluidBulkModulusString() { return "fluidBulkModulus"; }

    /// string/key for fluid initial pressure
    static constexpr char const * fluidInitialPressureString() { return "fluidInitialPressure"; }

    /// string/key for rate dependence parameter 1
    static constexpr char const * t1RateDependenceString() { return "t1RateDependence"; }

    /// string/key for rate dependence parameter 2
    static constexpr char const * t2RateDependenceString() { return "t2RateDependence"; }

    /// string/key for fracture energy release rate
    static constexpr char const * fractureEnergyReleaseRateString() { return "fractureEnergyReleaseRate"; }

    /// string/key for fracture softening shape parameter
    static constexpr char const * fractureSofteningExponentString() { return "fractureSofteningExponent"; }
   

    /// string/key for fracture stress
    static constexpr char const * fractureStressString() { return "fractureStress"; }

    /// string/key for initialDomainTemperature
    static constexpr char const * initialTemperatureString() { return "initialTemperature"; }

    /// string/key for activation energy
    static constexpr char const * QString() { return "Q"; }

    /// string/keay for brittleDuctileTransition
    static constexpr char const * brittleDuctileTransitionString() { return "brittleDuctileTransition"; }

    /// string/key for damage evolution criterion pressure
    static constexpr char const * damageEvolutionCriterionString() { return "damageEvolutionCriterion"; }

    /// string/key for peak t1 shear limit parameter
    static constexpr char const * peakI1String() { return "peakI1"; }

    /// string/key for F slope shear limit parameter
    static constexpr char const * fSlopeString() { return "fSlope"; }

    /// string/key for F slope shear limit parameter
    static constexpr char const * fSlopeFailedString() { return "fSlopeFailed"; }

    /// string/key for stren shear limit parameter
    static constexpr char const * strenString() { return "stren"; }

    /// string/key for Y slope shear limit parameter
    static constexpr char const * ySlopeString() { return "ySlope"; }

    /// string/key for nonassociativity parameter
    static constexpr char const * betaString() { return "beta"; }

    /// string/key for buckling flag
    static constexpr char const * enableBucklingString() { return "enableBuckling"; }

    /// string/key for buckling Length
    static constexpr char const * bucklingLengthString() { return "bucklingLength"; }

    /// string/key for buckling amplitude 
    static constexpr char const * bucklingAmplitudeString() { return "bucklingAmplitude"; }

    /// string/key for creep flag
    static constexpr char const * creepString() { return "enableCreep"; }

    /// string/key for creep C0 parameter
    static constexpr char const * creepC0String() { return "creepC0"; }

    /// string/key for creep C1 parameter
    static constexpr char const * creepC1String() { return "creepC1"; }

    /// string/key for creep C2 parameter
    static constexpr char const * creepC2String() { return "creepC2"; }

    /// string/key for creep A parameter
    static constexpr char const * creepAString() { return "creepA"; }

    /// string/key for creep B parameter
    static constexpr char const * creepBString() { return "creepB"; }

    /// string/key for creep C parameter
    static constexpr char const * creepCString() { return "creepC"; }

    /// string/key for creep D parameter
    static constexpr char const * creepDString() { return "creepD"; }

    /// string/key for creep E parameter
    static constexpr char const * creepEString() { return "creepE"; }

    /// string/key for creep F parameter
    static constexpr char const * creepFString() { return "creepF"; }

    /// string/key for creep G parameter
    static constexpr char const * creepGString() { return "creepG"; }

    /// string/key for strain-hardening N parameter
    static constexpr char const * strainHardeningNString() { return "strainHardeningN"; }

    /// string/key for strain-hardening K parameter
    static constexpr char const * strainHardeningKString() { return "strainHardeningK"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * plasticStrainToleranceString() { return "plasticStrainTolerance"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * stressReturnToleranceString() { return "stressReturnTolerance"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * maxAllowedSubcyclesString() { return "maxAllowedSubcycles"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * failedStepResponseString() { return "failedStepResponse"; }

    //string/key for element/particle bulk modulus value
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * shearModulusString() { return "shearModulus"; }

    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for material direction value
    static constexpr char const * materialDirectionString() { return "materialDirection"; }

    // string/key for element/particle deformation gradient value
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }

    /// string/key for quadrature point plasticStrain value 
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point porosity value 
    static constexpr char const * porosityString() { return "porosity"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point constitutive update status value
    static constexpr char const * constitutiveUpdateFlagString() { return "constitutiveUpdateFlag"; }

    /// string/key for quadrature point temperature value
    static constexpr char const * temperatureString() { return "temperature"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for element/particle strength scale
    static constexpr char const * strengthScaleString() { return "strengthScale"; }
  };

  /**
   * @brief Create a instantiation of the GeomechanicsUpdate class that refers to the data in this.
   * @return An instantiation of GeomechanicsUpdate.
   */
  GeomechanicsUpdates createKernelUpdates() const
  {
    return GeomechanicsUpdates( m_b0,
                                m_b1,
                                m_b2,
                                m_b3,
                                m_b4,
                                m_dstrendh,
                                m_dfslopedh,
                                m_dpeakI1dh,
                                m_dcrdh,
                                m_g0,
                                m_g1,
                                m_g2,
                                m_g3,
                                m_g4,
                                m_p0,
                                m_p1,
                                m_p2,
                                m_p3,
                                m_p4,
                                m_peakI1,
                                m_fSlope,
                                m_fSlopeFailed,
                                m_stren,
                                m_ySlope,
                                m_beta,
                                m_t1RateDependence,
                                m_t2RateDependence,
                                m_fractureEnergyReleaseRate,
                                m_fractureSofteningExponent,
                                m_fractureStress,
                                m_initialTemperature,
                                m_Q,
                                m_brittleDuctileTransition,
                                m_damageEvolutionCriterion,
                                m_cr,
                                m_fluidBulkModulus,
                                m_fluidInitialPressure,
                                m_enableBuckling,
                                m_bucklingLength,
                                m_bucklingAmplitude,
                                m_enableCreep,
                                m_creepC0,
                                m_creepC1,
                                m_creepC2,                                
                                m_creepA,
                                m_creepB,
                                m_creepC,
                                m_creepD,
                                m_creepE,
                                m_creepF,
                                m_creepG,
                                m_strainHardeningN,
                                m_strainHardeningK,
                                m_plasticStrainTolerance,
                                m_stressReturnTolerance,
                                m_maxAllowedSubcycles,
                                m_failedStepResponse,
                                m_bulkModulus,
                                m_shearModulus,
                                m_velocityGradient,
                                m_materialDirection,
                                m_deformationGradient,
                                m_plasticStrain,
                                m_porosity,
                                m_damage,
                                m_constitutiveUpdateFlag,
                                m_temperature,
                                m_lengthScale,
                                m_strengthScale,
                                m_thermalExpansionCoefficient,
                                m_newStress,
                                m_oldStress,
                                m_density,
                                m_wavespeed,
                                m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_b0,
                          m_b1,
                          m_b2,
                          m_b3,
                          m_b4,
                          m_dstrendh,
                          m_dfslopedh,
                          m_dpeakI1dh,
                          m_dcrdh,
                          m_g0,
                          m_g1,
                          m_g2,
                          m_g3,
                          m_g4,
                          m_p0,
                          m_p1,
                          m_p2,
                          m_p3,
                          m_p4,
                          m_peakI1,
                          m_fSlope,
                          m_fSlopeFailed,
                          m_stren,
                          m_ySlope,
                          m_beta,
                          m_t1RateDependence,
                          m_t2RateDependence,
                          m_fractureEnergyReleaseRate,
                          m_fractureSofteningExponent,
                          m_fractureStress,
                          m_initialTemperature,
                          m_Q,
                          m_brittleDuctileTransition,
                          m_damageEvolutionCriterion,
                          m_cr,
                          m_fluidBulkModulus,
                          m_fluidInitialPressure,
                          m_enableBuckling,
                          m_bucklingLength,
                          m_bucklingAmplitude,
                          m_enableCreep,
                          m_creepC0,
                          m_creepC1,
                          m_creepC2,
                          m_creepA,
                          m_creepB,
                          m_creepC,
                          m_creepD,
                          m_creepE,
                          m_creepF,
                          m_creepG,
                          m_strainHardeningN,
                          m_strainHardeningK,
                          m_plasticStrainTolerance,
                          m_stressReturnTolerance,
                          m_maxAllowedSubcycles,
                          m_failedStepResponse,
                          m_bulkModulus,
                          m_shearModulus,
                          m_velocityGradient,
                          m_materialDirection,
                          m_deformationGradient,
                          m_plasticStrain,
                          m_porosity,
                          m_damage,
                          m_constitutiveUpdateFlag,
                          m_temperature,
                          m_lengthScale,
                          m_strengthScale,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }

// CC: commmented out in case I need to use this for template of setters and getters for scalar values
//   /**
//    * @brief Getter for default transverse Young's modulus
//    * @return The value of the default transverse Young's modulus.
//    */
//   real64 getDefaultYoungModulusTransverse() const
//   {
//     return m_defaultYoungModulusTransverse;
//   }

//   /**
//    * @brief Setter for the default transverse Young's modulus.
//    * @param[in] input New value for the default transverse Young's modulus
//    */
//   void setDefaultYoungModulusTransverse( real64 const input )
//   {
//     m_defaultYoungModulusTransverse = input;
//   }

//   /**
//    * @brief Getter for default axial Young's modulus
//    * @return The value of the default axial Young's modulus.
//    */
//   real64 getDefaultYoungModulusAxial() const
//   {
//     return m_defaultYoungModulusAxial;
//   }

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         bulk modulus (at every element).
   */
  arrayView1d< real64 const > const bulkModulus() const { return m_bulkModulus; }

  /**
   * @brief Getter for bulk modulus.
   * @return reference to mutable bulk modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getBulkModulus() const override { return m_bulkModulus; }

  /**
   * @brief Accessor for shear modulus
   * @return A const reference to arrayView1d<real64> containing the shear
   *         modulus (at every element).
   */
  arrayView1d< real64 > const shearModulus() { return m_shearModulus; }

  /**
   * @brief Const accessor for shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         shear modulus (at every element).
   */
  arrayView1d< real64 const > const shearModulus() const { return m_shearModulus; }

  /**
   * @brief Getter for shear modulus.
   * @return reference to mutable shear modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getShearModulus() const override { return m_shearModulus; }

protected:
  virtual void postInputInitialization() override;

  // Tangent elastic bulk modulus parameters
  real64 m_b0;
  real64 m_b1;
  real64 m_b2;
  real64 m_b3;
  real64 m_b4;

  real64 m_dstrendh;
  real64 m_dfslopedh;
  real64 m_dpeakI1dh;
  real64 m_dcrdh;

  // Tangent elastic shear modulus parameters
  real64 m_g0;
  real64 m_g1;
  real64 m_g2;
  real64 m_g3;
  real64 m_g4;

  // Cruch curve parameters
  real64 m_p0;
  real64 m_p1;
  real64 m_p2;
  real64 m_p3;
  real64 m_p4;

  // Shear limit surface parameters
  real64 m_peakI1;
  real64 m_fSlope;
  real64 m_fSlopeFailed;
  real64 m_stren;
  real64 m_ySlope;

  // Nonassociativity parameter
  real64 m_beta;

  // Rate dependence paramters
  real64 m_t1RateDependence;
  real64 m_t2RateDependence;

  // Fracture energy release rate
  real64 m_fractureEnergyReleaseRate;
  real64 m_fractureSofteningExponent;  // shape parameter that controls softening with damage
  real64 m_fractureStress;

  real64 m_initialTemperature;
  real64 m_Q;

  real64 m_brittleDuctileTransition;
  int m_damageEvolutionCriterion;

  // Cap shape parameter
  real64 m_cr;

  // Fluid parameters
  real64 m_fluidBulkModulus;
  real64 m_fluidInitialPressure;

  // Buckling parameter
  int m_enableBuckling;
  real64 m_bucklingLength;
  real64 m_bucklingAmplitude;

  // Flag to enable creep
  int m_enableCreep;

  // deviatoric creep parameters
  real64 m_creepC0;
  real64 m_creepC1;
  real64 m_creepC2;  
  // compaction creep parameters
  real64 m_creepA;
  real64 m_creepB;
  real64 m_creepC;
  real64 m_creepD;
  real64 m_creepE;
  real64 m_creepF;
  real64 m_creepG;

  // strain-hardening parameters
  real64 m_strainHardeningN;
  real64 m_strainHardeningK;

  // Solution iteration control parameters
  real64 m_plasticStrainTolerance;
  real64 m_stressReturnTolerance;
  int m_maxAllowedSubcycles;
  int m_failedStepResponse;

  /// The bulk modulus for each element/particle
  array1d< real64 > m_bulkModulus;

  /// The shear modulus for each element/particle
  array1d< real64 > m_shearModulus;
 
  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  /// State variable: The material direction for each element/particle
  array2d< real64 > m_materialDirection;

  /// State variable: The deformation gradient values for each element/particle.
  array3d< real64 > m_deformationGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  // State variable: The porosity for each element/particle
  array2d< real64 > m_porosity;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: Constitutive update status for each quadrature point.
  ///  0 = good update, 1 = warning/retried or capped subcycling, -1 = failed update/delete flag.
  array2d< integer > m_constitutiveUpdateFlag;

  /// State variable: The temperature values for each quadrature point
  array1d< real64 > m_temperature;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_strengthScale;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_ */
