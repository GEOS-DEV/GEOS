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
 * @file Graphite.hpp
 * @brief Simple damage model for modeling material failure in brittle materials.
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. The model is really a hybrid plasticity/
 * damage model in the sense that we assume damaged material behaves like granular material
 * and hence follows a modified Mohr-Coulomb law. The modifications are that at low pressures,
 * the shape of the yield surface is modified to resemble a maximum principal stress criterion,
 * and at high pressures the shape converges on the von Mises yield surface. The damage
 * parameter results in softening of the deviatoric response i.e. causes the yield surface to
 * contract. Furthermore, damage is used to scale back tensile pressure: p = (1 - d) * pTrial.
 * pTrial is calculatd as pTrial = -k * log(J), where the Jacobian J of the material motion is
 * integrated and tracked by this model.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_

#include "ElasticTransverseIsotropicPressureDependent.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GraphiteUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GraphiteUpdates : public ElasticTransverseIsotropicPressureDependentUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] velocityGradient The ArrayView holding the velocity gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point.
   * @param[in] relaxation The ArrayView holding the relaxation for each quadrature point.
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] materialDirection The ArrayView holding the material direction for each element/particle.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] failureStrength The failure strength.
   * @param[in] crackSpeed The crack speed velocity.
   * @param[in] c11 The ArrayView holding the C11 data for each element.
   * @param[in] c13 The ArrayView holding the C13 data for each element.
   * @param[in] c33 The ArrayView holding the C33 data for each element.
   * @param[in] c44 The ArrayView holding the C44 data for each element.
   * @param[in] c66 The ArrayView holding the C66 data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  GraphiteUpdates( arrayView3d< real64 > const & velocityGradient,
                   arrayView3d< real64 > const & plasticStrain,
                   arrayView2d< real64 > const & relaxation,
                   arrayView2d< real64 > const & damage,
                   arrayView2d< real64 > const & jacobian,
                   arrayView1d< real64 > const & lengthScale,
                   real64 const & failureStrength,
                   real64 const & crackSpeed,
                   real64 const & damagedMaterialFrictionalSlope,
                   real64 const & distortionShearResponseX2,
                   real64 const & distortionShearResponseY1,
                   real64 const & distortionShearResponseY2,
                   real64 const & distortionShearResponseM1,
                   real64 const & inPlaneShearResponseX2,
                   real64 const & inPlaneShearResponseY1,
                   real64 const & inPlaneShearResponseY2,
                   real64 const & inPlaneShearResponseM1,
                   real64 const & coupledShearResponseX2,
                   real64 const & coupledShearResponseY1,
                   real64 const & coupledShearResponseY2,
                   real64 const & coupledShearResponseM1,
                   real64 const & maximumPlasticStrain,
                   real64 const & refC11,
                   real64 const & refC13,
                   real64 const & refC33,
                   real64 const & refC44,
                   real64 const & refC66,
                   real64 const & dc11dp,
                   real64 const & dc13dp,
                   real64 const & dc33dp,
                   real64 const & dc44dp,
                   real64 const & dc66dp,
                   arrayView1d< real64 > const & c11,
                   arrayView1d< real64 > const & c13,
                   arrayView1d< real64 > const & c33,
                   arrayView1d< real64 > const & c44,
                   arrayView1d< real64 > const & c66,
                   arrayView1d< real64 > const & effectiveBulkModulus,
                   arrayView1d< real64 > const & effectiveShearModulus,
                   arrayView2d< real64 > const & materialDirection,
                   arrayView1d< real64 > const & thermalExpansionCoefficient,
                   arrayView3d< real64, solid::STRESS_USD > const & newStress,
                   arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                   bool const & disableInelasticity ):
    ElasticTransverseIsotropicPressureDependentUpdates( refC11,
                                                        refC13,
                                                        refC33,
                                                        refC44,
                                                        refC66,
                                                        dc11dp,
                                                        dc13dp,
                                                        dc33dp,
                                                        dc44dp,
                                                        dc66dp,
                                                        c11, 
                                                        c13, 
                                                        c33, 
                                                        c44, 
                                                        c66, 
                                                        effectiveBulkModulus,
                                                        effectiveShearModulus,
                                                        materialDirection,
                                                        thermalExpansionCoefficient, 
                                                        newStress, 
                                                        oldStress, 
                                                        disableInelasticity ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain ),
    m_relaxation( relaxation ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_lengthScale( lengthScale ),
    m_failureStrength( failureStrength ),
    m_crackSpeed( crackSpeed ),
    m_damagedMaterialFrictionalSlope( damagedMaterialFrictionalSlope ),
    m_distortionShearResponseX2( distortionShearResponseX2 ),
    m_distortionShearResponseY1( distortionShearResponseY1 ),
    m_distortionShearResponseY2( distortionShearResponseY2 ),
    m_distortionShearResponseM1( distortionShearResponseM1 ),
    m_inPlaneShearResponseX2( inPlaneShearResponseX2 ),
    m_inPlaneShearResponseY1( inPlaneShearResponseY1 ),
    m_inPlaneShearResponseY2( inPlaneShearResponseY2 ),
    m_inPlaneShearResponseM1( inPlaneShearResponseM1 ),
    m_coupledShearResponseX2( coupledShearResponseX2 ),
    m_coupledShearResponseY1( coupledShearResponseY1 ),
    m_coupledShearResponseY2( coupledShearResponseY2 ),
    m_coupledShearResponseM1( coupledShearResponseM1 ),
    m_maximumPlasticStrain( maximumPlasticStrain )
  {}

  /// Default copy constructor
  GraphiteUpdates( GraphiteUpdates const & ) = default;

  /// Default move constructor
  GraphiteUpdates( GraphiteUpdates && ) = default;

  /// Deleted default constructor
  GraphiteUpdates() = delete;

  /// Deleted copy assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates const & ) = delete;

  /// Deleted move assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate;

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
  void computeTransverselyIsotropicTrialStress( const real64 timeIncrement,      
                                                const real64 Ez,                 
                                                const real64 Ep,                 
                                                const real64 nuzp,               
                                                const real64 nup,                
                                                const real64 Gzp,                
                                                real64 const (& materialDirection)[3], 
                                                real64 const (& oldStress)[6],
                                                real64 const (& D)[6],           
                                                real64 (& newStress) [6] ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicPlasticStrainIncrement( real64 const ( & velocityGradient )[3][3],
                                                           real64 const ( & oldStress )[6],
                                                           real64 const ( & newStress )[6],
                                                           const real64 Ez,             
                                                           const real64 Ep,             
                                                           const real64 nuzp,           
                                                           const real64 nup,            
                                                           const real64 Gzp,            
                                                           real64 const ( & materialDirection )[3],	
                                                           const real64 timeIncrement,
                                                           real64 ( & plasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB1( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB2( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB3( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB4( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB5( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;

  GEOS_HOST_DEVICE
  real64 delta( const int i,
                const int j)  const;

  GEOS_HOST_DEVICE
  real64 slopePoint0( const real64 x,
                      const real64 x1,
                      const real64 x2,
                      const real64 y1,
                      const real64 y2,
                      const real64 m1 ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticTransverseIsotropicPressureDependentUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the velocity gradient for each element/particle.
  arrayView3d< real64 > const m_velocityGradient;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_relaxation;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// The maximum theoretical strength
  real64 const m_failureStrength;

  // The time to failure
  real64 const m_crackSpeed;

  // The damaged material frictional slope
  real64 const m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope
  real64 const m_distortionShearResponseX2;
  real64 const m_distortionShearResponseY1;
  real64 const m_distortionShearResponseY2;
  real64 const m_distortionShearResponseM1;

  real64 const m_inPlaneShearResponseX2;
  real64 const m_inPlaneShearResponseY1;
  real64 const m_inPlaneShearResponseY2;
  real64 const m_inPlaneShearResponseM1;

  real64 const m_coupledShearResponseX2;
  real64 const m_coupledShearResponseY1;
  real64 const m_coupledShearResponseY2;
  real64 const m_coupledShearResponseM1;
   
  real64 const m_maximumPlasticStrain;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdate not implemented for Graphite model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, 
                     q, 
                     timeIncrement, 
                     strainIncrement, 
                     stress, 
                     stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for CeramicDamage model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( & beginningRotation )[3][3],
                                                    real64 const ( & endRotation )[3][3],
                                                    real64 const ( & strainIncrement )[6],
                                                    real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate_StressOnly( k,
                                                                                    q,
                                                                                    timeIncrement,
                                                                                    beginningRotation,
                                                                                    endRotation,
                                                                                    strainIncrement,
                                                                                    stress );

  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  GraphiteUpdates::smallStrainUpdateHelper( k, 
                                            q, 
                                            timeIncrement,                
                                            beginningRotation,
                                            endRotation,
                                            stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdateHelper( localIndex const k,
                                               localIndex const q,
                                               real64 const timeIncrement,  
                                               real64 const ( & beginningRotation )[3][3],
                                               real64 const ( & endRotation )[3][3],
                                               real64 ( & stress )[6] ) const
{
    GEOS_UNUSED_VAR( endRotation );

    real64 oldStress[6];
    LvArray::tensorOps::copy< 6 >(oldStress, stress);

    real64 rotationTranspose[3][3];
    LvArray::tensorOps::transpose< 3, 3 >(rotationTranspose, beginningRotation); 

    real64 unrotatedVelocityGradient[3][3];
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >(unrotatedVelocityGradient, rotationTranspose, m_velocityGradient[k]);
    
    real64 unrotatedVelocityGradientTranspose[3][3];
    LvArray::tensorOps::transpose< 3, 3 >(unrotatedVelocityGradientTranspose, unrotatedVelocityGradient);

    // CC: Is there an LvArray operation to get the symmetric part of a matrix?
    real64 denseD[3][3];
    LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
    LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
    LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );

    real64 D[6];
    LvArray::tensorOps::denseToSymmetric<3>(D, denseD);

    // make sure material direction is normalized.
    real64 materialDirection[3];
    LvArray::tensorOps::copy< 3 >( materialDirection, m_materialDirection[k] );
    LvArray::tensorOps::normalize< 3 >( materialDirection );

    // Use beginning of step normal stress to compute stress dependence of Ez
    // real64 temp[3];
    int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
    // LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, oldStress, materialDirection );
    // real64 oldPlaneNormalStress = LvArray::tensorOps::AiBi< 3 >( materialDirection, temp );  // CC: Unused?

    // Beginning of step pressure to compute pressure-dependence of elastic moduli
    // real64 oldPressure = (-1.0/3.0)*( oldStress[0] + oldStress[1] + oldStress[2] ); // CC: Unused?

    // This is a transversely isotropic material for graphite-like crystals having
    // some weak plane with plane-normal-stress- and pressure-dependent elastic properties:

    // CC: in old geos the elastic on two different lines of the same code, Mike used planeNormalStress in one but not the other
    // need to ask him about that

    // Elastic constants are updated based on pressure before in ElasticTransverseIsotropicPressureDependent update
    real64 c11 = m_c11[k];
    real64 c13 = m_c13[k];
    real64 c33 = m_c33[k];
    real64 c44 = m_c44[k];
    real64 c66 = m_c66[k];

    real64 Ez = c33 - c13 * c13 / ( c11 - c66 );
    real64 Ep = 4 * c66 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 );
    real64 Gzp = c44 / 2.0;
    real64 nuzp = c13 / ( 2 * ( c11 - c66 ) );
    real64 nup = Ep / c66 - 1; //4 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 ) - 1;

    // Hypoelastic trial stress update.
    // CC: this should already be done by ElasticTransverseIsotropicPressureDependent
    // computeTransverselyIsotropicTrialStress( timeIncrement,      // time step
    //                                          Ez,                 // Elastic modulus preferred direction
    //                                          Ep,                 // Elastic modulus transverse plane
    //                                          nuzp,               // Poisson ratio coupled
    //                                          nup,                // Poisson ratio transverse plane
    //                                          Gzp,                // Shear modulus coupled plane
    //                                          materialDirection,	 // preferred direction
    //                                          oldStress,          // stress at start of step
    //                                          D,                  // D=sym(L)
    //                                          stress );           // stress at end of step

    // Decompose stress tensor into pieces.
    real64 sigma1[6] = {0}; // axial
    real64 sigma2[6] = {0}; // in-plane normal
    real64 sigma4[6] = {0}; // in-plane total stress
    real64 sigma5[6] = {0}; // weak plane - shear
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int p=0; p<3; p++)
            {
                for(int w=0; w<3; w++)
                {
                    sigma1[voigtMap[i][j]] += transverselyIsotropicB1(materialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma2[voigtMap[i][j]] += transverselyIsotropicB2(materialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma4[voigtMap[i][j]] += transverselyIsotropicB4(materialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma5[voigtMap[i][j]] += transverselyIsotropicB5(materialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                }
            }
        }
    }

    // Trial pressure to compute pressure-dependence of strength
    real64 pressure = (-1.0/3.0)*( stress[0] + stress[1] + stress[2] );

    // Check for tensile failure in preferred direction
    real64 temp[3] = { 0 };
    LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, stress, materialDirection );
    real64 planeNormalStress = LvArray::tensorOps::AiBi< 3 >( materialDirection, temp );

    // increment damage, but enforce 0<=d<=1
    if (planeNormalStress > m_failureStrength)
    {
        real64 timeToFailure = m_lengthScale[k]/m_crackSpeed;
        m_damage[k][q] = std::min( m_damage[k][q] + timeIncrement / timeToFailure, 1.0 );
    }

    // strength scale factor, combining plastic softening and damage
    real64 fac = (1.0 - m_damage[k][q])*(1.0 - m_relaxation[k][q]);

    // Enforce damage, no tensile stress on plane, and frictional response to shear.
    if(m_damage[k][q] == 1.0)
    {
        if ( planeNormalStress > 0 )
        {
            LvArray::tensorOps::scale< 6 >(sigma1, 0.0);
        }
    }

    // find in-plane isotropic and deviatoric stress
    real64 inPlaneIso[6], inPlaneDev[6];
    LvArray::tensorOps::copy< 6 >(inPlaneIso, sigma2);
    LvArray::tensorOps::scale< 6 >(inPlaneIso, 0.5);
    LvArray::tensorOps::copy< 6 >(inPlaneDev, sigma4);
    LvArray::tensorOps::subtract< 6 >(inPlaneDev, inPlaneIso);

    // find distortion stress from in-plane isotropic and plane-normal stress
    real64 distortion[6], distortion_iso[6] = {0}, distortion_dev[6];
    LvArray::tensorOps::copy< 6 >(distortion, inPlaneIso);
    LvArray::tensorOps::add< 6 >(distortion, sigma1);

    real64 trialP;
    real64 trialQ; //  Check that this isn't a redeclaration
    twoInvariant::stressDecomposition( distortion,
                                       trialP,
                                       trialQ,
                                       distortion_dev );

    distortion_iso[0] = trialP;
    distortion_iso[1] = trialP;
    distortion_iso[2] = trialP;

    // Define 3 pressure-dependent shear strengths for different shear modes.
    real64 inPlaneShearStrength, coupledYieldStrength, totalShearStrength;
    real64 x1, x2, y1, y2, m1;

    // -------------------------------
    // "distortion" shear relating sigma normal and in-plane iso"
    x1 = 0;
    x2 = m_distortionShearResponseX2;
    y1 = m_distortionShearResponseY1;
    y2 = m_distortionShearResponseY2;
    m1 = m_distortionShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 = fac*y1;
    m1 = fac*m1 + (1. - fac)*std::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if(pressure<x1)
    {
        totalShearStrength=std::max(0.0,y1-(x2-pressure)*m1);
    }
    else if(pressure<x2)
    {
        totalShearStrength=slopePoint0(pressure,x1,x2,y1,y2,m1);
    }
    else
    {
        totalShearStrength=y2;
    }

    // -------------------------------------------
    // Coupled shear response (slip on weak plane)
    x1 = 0;
    x2 = m_coupledShearResponseX2;
    y1 = m_coupledShearResponseY1;
    y2 = m_coupledShearResponseY2;
    m1 = m_coupledShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 = fac*y1;
    m1 = fac*m1 + (1. - fac)*std::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if(pressure<x1)
    {
        coupledYieldStrength=std::max(0.0,y1-(x2-pressure)*m1);
    }
    else if(pressure<x2)
    {
        coupledYieldStrength=slopePoint0(pressure,x1,x2,y1,y2,m1);
    }
    else
    {
        coupledYieldStrength=y2;
    }

    // -------------------------------
    // In-plane shear response
    x1 = 0;
    x2 = m_inPlaneShearResponseX2;
    y1 = m_inPlaneShearResponseY1;
    y2 = m_inPlaneShearResponseY2;
    m1 = m_inPlaneShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 = fac*y1;
    m1 = fac*m1 + (1. - fac)*std::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if( pressure < x1 )
    {
        inPlaneShearStrength=std::max(0.0, y1  -( x2 - pressure ) * m1 );
    }
    else if( pressure < x2 )
    {
        inPlaneShearStrength=slopePoint0( pressure, x1, x2, y1, y2, m1);
    }
    else
    {
        inPlaneShearStrength = y2;
    }

    // Enforce strength ---------------------------------------------
    // flags indicating whether plastic strain needs to be updated.
    bool plastic = false;

    // Enforce distortion strain yield
    real64 totalShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( distortion_dev );
    // GEOS_LOG_RANK_0(k << ": totalShearStress " << totalShearStress);
    if ( totalShearStress > totalShearStrength )
    {
      LvArray::tensorOps::scale< 6 >( distortion_dev, totalShearStrength / totalShearStress );
      // distortion_dev *= totalShearStrength / totalShearStress;
      plastic = true;
    }

    // enforce in-plane yield
    real64 inPlaneShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( inPlaneDev ) ;
    // GEOS_LOG_RANK_0(k << ":inPlaneShearStress " << inPlaneShearStress);
    if ( inPlaneShearStress > inPlaneShearStrength )
    {
      LvArray::tensorOps::scale< 6 >( inPlaneDev, inPlaneShearStrength / inPlaneShearStress );
      // inPlaneDev *= inPlaneShearStrength / inPlaneShearStress;
      plastic = true;
    }

    // enforce coupled yield
    real64 coupledShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( sigma5 );
    // GEOS_LOG_RANK_0(k << ": coupledShearStress " << coupledShearStress);
    if ( coupledShearStress > coupledYieldStrength )
    {
      LvArray::tensorOps::scale< 6 >( sigma5, coupledYieldStrength / coupledShearStress );
      // sigma5 *= coupledYieldStrength / coupledShearStress;
      plastic = true;
    }

    // reassemble end-of-step stress
    real64 newStress[6] = {0};
    LvArray::tensorOps::add< 6 >(newStress, distortion_iso);
    LvArray::tensorOps::add< 6 >(newStress, distortion_dev);
    LvArray::tensorOps::add< 6 >(newStress, inPlaneDev);
    LvArray::tensorOps::add< 6 >(newStress, sigma5);

    //////////////////////////////////////////////////////
    // Evolve state variables.
    //
    // For inelastic response evolve damage if pressure is below the
    // brittle ductile transition pressure and evolve plastic strain
    // if pressure is above brittle ductile transition pressure.
    // plastic strain doesn't currently affect material response, but
    // it may be useful for plotting regions of high-pressure yield.
    if( plastic )
    {
        // increment plastic strain
        real64 plasticStrainIncrement[6] = {0};
        computeTransverselyIsotropicPlasticStrainIncrement( unrotatedVelocityGradient,          // Velocity gradient tensor
                                                            oldStress,                 // Stress at start of step
                                                            stress,                    // Stress at end of step
                                                            Ez,                        // Elastic modulus preferred direction
                                                            Ep,                        // Elastic modulus transverse plane
                                                            nuzp,                      // Poisson ratio coupled
                                                            nup,                       // Poisson ratio transverse plane
                                                            Gzp,                       // Shear modulus coupled plane
                                                            materialDirection,		     // material direction (unit vector, plane normal)
                                                            timeIncrement,             // timeStep
                                                            plasticStrainIncrement );
        LvArray::tensorOps::add< 6 >(m_plasticStrain[k][q], plasticStrainIncrement);

        // increment relaxation
        m_relaxation[k][q] += LvArray::tensorOps::l2Norm< 6 >( plasticStrainIncrement ) / m_maximumPlasticStrain;
        m_relaxation[k][q] = std::min(1.0, m_relaxation[k][q]);
    }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicTrialStress(const real64 timeIncrement,      // time step
                                                              const real64 Ez,                 // Elastic modulus preferred direction
                                                              const real64 Ep,                 // Elastic modulus transverse plane
                                                              const real64 nuzp,               // Poisson ratio coupled
                                                              const real64 nup,                // Poisson ratio transverse plane
                                                              const real64 Gzp,                // Shear modulus coupled plane
                                                              real64 const (& materialDirection)[3], // preferred direction
                                                              real64 const (& oldStress)[6],   // stress at start of step.
                                                              real64 const (& D)[6],           // D=sym(L)
                                                              real64 (& newStress) [6]         // stress at end of step
) const
{
  // These are the TI elastic stiffness coefficients using Brannon's TI basis tensors:
	real64 h1 = ( Ez*Ez*(-1 + nup) ) / ( Ez*(-1 + nup) + 2*Ep*nuzp*nuzp );
	real64 h2 = -( ( Ep*(Ez*nup + Ep*nuzp*nuzp) ) / ( (1 + nup)*(Ez*(-1 + nup) + 2*Ep*nuzp*nuzp ) ) );
	real64 h3 = ( Ep*Ez*nuzp ) / ( Ez - Ez*nup - 2*Ep*nuzp*nuzp );
	real64 h4 = Ep/( 1 + nup );
	real64 h5 = 2*Gzp;

	LvArray::tensorOps::copy< 6 >(newStress, oldStress);
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int p=0; p<3; p++)
			{
				for(int w=0; w<3; w++)
				{
					newStress[voigtMap[i][j]] += (h1*transverselyIsotropicB1(materialDirection,i,j,p,w) +
                                        h2*transverselyIsotropicB2(materialDirection,i,j,p,w) +
                                        h3*transverselyIsotropicB3(materialDirection,i,j,p,w) +
                                        h4*transverselyIsotropicB4(materialDirection,i,j,p,w) +
                                        h5*transverselyIsotropicB5(materialDirection,i,j,p,w))*D[voigtMap[p][w]]*timeIncrement;
				}
			}
		}
	}

}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicPlasticStrainIncrement(real64 const ( & velocityGradient )[3][3],          // Velocity gradient tensor
                                                                         real64 const ( & oldStress )[6], // Stress at start of step
                                                                         real64 const ( & newStress )[6], // Stress at end of step
                                                                         const real64 Ez,             // Elastic modulus preferred direction
                                                                         const real64 Ep,             // Elastic modulus transverse plane
                                                                         const real64 nuzp,           // Poisson ratio coupled
                                                                         const real64 nup,            // Poisson ratio transverse plane
                                                                         const real64 Gzp,            // Shear modulus coupled plane
                                                                         real64 const ( & materialDirection )[3],			// material direction (unit vector)
                                                                         const real64 timeIncrement,
                                                                         real64 ( & plasticStrainIncrement )[6]
) const
{  // For hypo-elastic transversely isotropic models we compute the increment in plastic strain assuming
	// for some increment in total strain and stress and elastic properties.

	// New stress minus old stress
    real64 stressIncrement[6];
    LvArray::tensorOps::copy< 6 >( stressIncrement, newStress );
    LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );

	// These are the TI elastic compliance coefficients using Brannon's TI basis tensors:
	real64 s1 = 1.0/Ez;
	real64 s2 = -nup/Ep;
	real64 s3 = -nuzp/Ez;
	real64 s4 = (1+nup)/Ep;
	real64 s5 = 1./(2.*Gzp);

	real64 elasticStrainIncrement[6];
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int p=0; p<3; p++)
			{
				for(int w=0; w<3; w++)
				{
					elasticStrainIncrement[voigtMap[i][j]] += ( s1*transverselyIsotropicB1(materialDirection,i,j,p,w) +
                                                      s2*transverselyIsotropicB2(materialDirection,i,j,p,w) +
                                                      s3*transverselyIsotropicB3(materialDirection,i,j,p,w) +
                                                      s4*transverselyIsotropicB4(materialDirection,i,j,p,w) +
                                                      s5*transverselyIsotropicB5(materialDirection,i,j,p,w) )*stressIncrement[voigtMap[p][w]];
				}
			}
		}
	}

	// plastic strain increment = Total strain increment - elastic strain increment
	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
            //Old GEOS code did not index plasticStrainIncrement, how did this not through an error
            //does the voigt notation include the factore of 2 for off axis entries of the tensor?
			plasticStrainIncrement[voigtMap[i][j]] += 0.5 * ( velocityGradient[i][j] + velocityGradient[j][i] ) * timeIncrement - elasticStrainIncrement[voigtMap[i][j]];
		}
	}

}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB1( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B1_ijkl component, of the transversely isotropic basis tensor B1
  // described in Brannon's rotation/tensor book.

	real64 B1 = materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B1;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB2( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B2_ijkl component, of the transversely isotropic basis tensor B2
  // described in Brannon's rotation/tensor book.

	real64 B2 = delta(i,j)*delta(p,w) - 
               ( materialDirection[i]*materialDirection[j]*delta(p,w) + 
                 delta(i,j)*materialDirection[p]*materialDirection[w] ) + 
                 materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B2;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB3( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B3_ijkl component, of the transversely isotropic basis tensor B3
  // described in Brannon's rotation/tensor book.

	real64 B3 = materialDirection[i]*materialDirection[j]*delta(p,w) + 
                delta(i,j)*materialDirection[p]*materialDirection[w] - 
                2.0*materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B3;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB4( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B4_ijkl component, of the transversely isotropic basis tensor B4
  // described in Brannon's rotation/tensor book.

	real64 B4 = 0.5*( delta(i,p)*delta(j,w) + delta(i,w)*delta(j,p) ) -
              0.5*( delta(i,w)*materialDirection[j]*materialDirection[p] + materialDirection[i]*delta(j,p)*materialDirection[w] + 
                    materialDirection[i]*delta(j,w)*materialDirection[p] + delta(i,p)*materialDirection[j]*materialDirection[w] ) + 
                    materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B4;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB5( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B5_ijkl component, of the transversely isotropic basis tensor B5
  // described in Brannon's rotation/tensor book.

	real64 B5 = 0.5*( delta(i,p)*materialDirection[j]*materialDirection[w] + 
                      delta(i,w)*materialDirection[j]*materialDirection[p] + 
                      materialDirection[i]*delta(j,w)*materialDirection[p] + 
                      materialDirection[i]*delta(j,p)*materialDirection[w] ) - 
                      2.0*materialDirection[i]*materialDirection[p]*materialDirection[j]*materialDirection[w];
	return B5;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::delta( const int i,
		                       const int j ) const
{ // Dirac delta function that returns a real
	return i == j ? 1.0 : 0.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::slopePoint0( const real64 x,
                                     const real64 x1,
                                     const real64 x2,
                                     const real64 y1,
                                     const real64 y2,
                                     const real64 m1 ) const
{
  // returns a function h(x), for x1 <= x <= x2, with h(x1)=y1, h(x2)=y2, h'(x1)=m1, x'(x2)=0
  //realT beta = (m1*x1 - m1*x2 - y1 + y2)/(y1 - y2);
  //return -(((m1*(x1 - x2) + m1*std::pow((x - x2)/(x1 - x2),beta)*(-x + x2) + (m1*(-x1 + x2)*y1)/(y1 - y2))*(y1 - y2))/(m1*(x1 - x2)));
  // This should be the same, but simplified:
  return y2 + (y1-y2)*std::pow( (x-x2)/(x1-x2) , m1*(x1-x2)/(y1-y2) );
}


/**
 * @class Graphite
 *
 * Graphite material model.
 */
class Graphite : public ElasticTransverseIsotropicPressureDependent
{
public:

  /// @typedef Alias for GraphiteUpdates
  using KernelWrapper = GraphiteUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Graphite( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Graphite() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Graphite";

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
    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for quadrature point plasticStrain value 
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point relaxation value
    static constexpr char const * relaxationString() { return "relaxation"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for maximum strength
    static constexpr char const * failureStrengthString() { return "failureStrength"; }

    /// string/key for crack speed
    static constexpr char const * crackSpeedString() { return "crackSpeed"; }

    /// string/key for damaged material frictional slope
    static constexpr char const * damagedMaterialFrictionalSlopeString() { return "damagedMaterialFrictionalSlope"; }
    
    // string/key for distortion shear response parameter x2
    static constexpr char const * distortionShearResponseX2String() { return "distortionShearResponseX2"; }
    
    // string/key for distortion shear response parameter y1
    static constexpr char const * distortionShearResponseY1String() { return "distortionShearResponseY1"; }
    
    // string/key for distortion shear response parameter y2
    static constexpr char const * distortionShearResponseY2String() { return "distortionShearResponseY2"; }
    
    // string/key for distortion shear response parameter m1
    static constexpr char const * distortionShearResponseM1String() { return "distortionShearResponseM1"; }
    
    // string/key for in plane shear response x2
    static constexpr char const * inPlaneShearResponseX2String() { return "inPlaneShearResponseX2"; }
    
    // string/key for in plane shear response y1
    static constexpr char const * inPlaneShearResponseY1String() { return "inPlaneShearResponseY1"; }
    
    // string/key for in plane shear response y2
    static constexpr char const * inPlaneShearResponseY2String() { return "inPlaneShearResponseY2"; }
    
    // string/key for in plane shear response m1
    static constexpr char const * inPlaneShearResponseM1String() { return "inPlaneShearResponseM1"; }

    // string/key for coupled shear response x2
    static constexpr char const * coupledShearResponseX2String() { return "coupledShearResponseX2"; }
    
    // string/key for coupled shear response y1 
    static constexpr char const * coupledShearResponseY1String() { return "coupledShearResponseY1"; }
    
    // string/key for coupled shear response y2 
    static constexpr char const * coupledShearResponseY2String() { return "coupledShearResponseY2"; }
    
    // string/key for coupled shear response m1 
    static constexpr char const * coupledShearResponseM1String() { return "coupledShearResponseM1"; }
    
    // string/key for maximum plastic strain
    static constexpr char const * maximumPlasticStrainString() { return "maximumPlasticStrain"; }
  };

  /**
   * @brief Create a instantiation of the GraphiteUpdate class that refers to the data in this.
   * @return An instantiation of GraphiteUpdate.
   */
  GraphiteUpdates createKernelUpdates() const
  {
    return GraphiteUpdates( m_velocityGradient,
                            m_plasticStrain,
                            m_relaxation,
                            m_damage,
                            m_jacobian,
                            m_lengthScale,
                            m_failureStrength,
                            m_crackSpeed,
                            m_damagedMaterialFrictionalSlope,
                            m_distortionShearResponseX2,
                            m_distortionShearResponseY1,
                            m_distortionShearResponseY2,
                            m_distortionShearResponseM1,
                            m_inPlaneShearResponseX2,
                            m_inPlaneShearResponseY1,
                            m_inPlaneShearResponseY2,
                            m_inPlaneShearResponseM1,
                            m_coupledShearResponseX2,
                            m_coupledShearResponseY1,
                            m_coupledShearResponseY2,
                            m_coupledShearResponseM1,
                            m_maximumPlasticStrain,
                            m_refC11,
                            m_refC13,
                            m_refC33,
                            m_refC44,
                            m_refC66,
                            m_dc11dp,
                            m_dc13dp,
                            m_dc33dp,
                            m_dc44dp,
                            m_dc66dp,
                            m_c11,
                            m_c13,
                            m_c33,
                            m_c44,
                            m_c66,
                            m_effectiveBulkModulus,
                            m_effectiveShearModulus,
                            m_materialDirection,
                            m_thermalExpansionCoefficient,
                            m_newStress,
                            m_oldStress,
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
                          m_velocityGradient,
                          m_plasticStrain,
                          m_relaxation,
                          m_damage,
                          m_jacobian,
                          m_lengthScale,
                          m_failureStrength,
                          m_crackSpeed,
                          m_damagedMaterialFrictionalSlope,
                          m_distortionShearResponseX2,
                          m_distortionShearResponseY1,
                          m_distortionShearResponseY2,
                          m_distortionShearResponseM1,
                          m_inPlaneShearResponseX2,
                          m_inPlaneShearResponseY1,
                          m_inPlaneShearResponseY2,
                          m_inPlaneShearResponseM1,
                          m_coupledShearResponseX2,
                          m_coupledShearResponseY1,
                          m_coupledShearResponseY2,
                          m_coupledShearResponseM1,
                          m_maximumPlasticStrain,
                          m_refC11,
                          m_refC13,
                          m_refC33,
                          m_refC44,
                          m_refC66,
                          m_dc11dp,
                          m_dc13dp,
                          m_dc33dp,
                          m_dc44dp,
                          m_dc66dp,
                          m_c11,
                          m_c13,
                          m_c33,
                          m_c44,
                          m_c66,
                          m_effectiveBulkModulus,
                          m_effectiveShearModulus,
                          m_materialDirection,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postProcessInput() override;
  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The relaxation values for each quadrature point
  array2d< real64 > m_relaxation;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Material parameter: The value of the failure strength
  real64 m_failureStrength;

  /// Material parameter: The value of crack speed
  real64 m_crackSpeed;

  /// Material parameter: The value of the damaged material frictional slope
  real64 m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope
  real64 m_distortionShearResponseX2;
  real64 m_distortionShearResponseY1;
  real64 m_distortionShearResponseY2;
  real64 m_distortionShearResponseM1;

  real64 m_inPlaneShearResponseX2;
  real64 m_inPlaneShearResponseY1;
  real64 m_inPlaneShearResponseY2;
  real64 m_inPlaneShearResponseM1;

  real64 m_coupledShearResponseX2;
  real64 m_coupledShearResponseY1;
  real64 m_coupledShearResponseY2;
  real64 m_coupledShearResponseM1;
  
  real64 m_maximumPlasticStrain;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP_ */
