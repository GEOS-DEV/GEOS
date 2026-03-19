/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BartonBandisPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_BARTONBANDISPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_BARTONBANDISPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class BartonBandisPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  BartonBandisPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure,
                                    //arrayView4d< real64 > const & dPerm_dDispJump,
                                    bool const updateTransversalComponent,
                                    real64 const aperture0,
                                       real64 const biot,
                                       real64 const poisson,
                                       real64 const normalStiffness,
                                       real64 const referencePressure,
                                       R1Tensor const &referenceTotalStress )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    //m_dPerm_dDispJump( dPerm_dDispJump ),
    m_numDimensionsToUpdate( 3 ),
    m_aperture0( aperture0 ),
    m_biot( biot ),
    m_poisson( poisson ),
    m_normalStiffness( normalStiffness ), // Kni
    m_referencePressure( referencePressure ),
    m_referenceTotalStress( referenceTotalStress )
  {
    m_numDimensionsToUpdate = updateTransversalComponent ? 3 : 2;

    m_referenceEffectiveStress[0] = m_referenceTotalStress[0] - m_biot*m_referencePressure; 
    m_referenceEffectiveStress[1] = m_referenceTotalStress[1] - m_biot*m_referencePressure; 
    m_referenceEffectiveStress[2] = m_referenceTotalStress[2] - m_biot*m_referencePressure; 
  }

  GEOS_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                //real64 const & dHydraulicAperture_dNormalJump,
                arraySlice1d< real64 > const & permeability,
                //arraySlice2d< real64 > const & dPerm_dDispJump,
                real64 & dPerm_dHydraulicAperture  ) const
  {
    GEOS_UNUSED_VAR( oldHydraulicAperture );

    real64 const perm  = newHydraulicAperture*newHydraulicAperture / 12.0;
    dPerm_dHydraulicAperture = newHydraulicAperture / 6.0;

    for( int dim=0; dim < m_numDimensionsToUpdate; dim++ )
    {
      permeability[dim]        = perm;
      /*dPerm_dDispJump[dim][0]  = dPerm_dHydraulicAperture * dHydraulicAperture_dNormalJump;
      dPerm_dDispJump[dim][1]  = 0.0;
      dPerm_dDispJump[dim][2]  = 0.0;*/
    }
  }

  /*GEOS_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & oldHydraulicAperture,
                                   real64 const & newHydraulicAperture,
                                   real64 const & dHydraulicAperture_dNormalJump ) const override final
  {
    GEOS_UNUSED_VAR( q );
    compute( oldHydraulicAperture,
             newHydraulicAperture,
             dHydraulicAperture_dNormalJump,
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0] );
  }*/

  GEOS_HOST_DEVICE
  virtual void updateFromPressureApertureAndNormal( localIndex const k,
                                   localIndex const q,
                                   real64 const & pressure,
                                   real64 const & oldHydraulicAperture,
                                   real64 & newHydraulicAperture,
                                   array1d< real64 > const & normal,
                                   real64 const & dHydraulicAperture_dNormalJump ) const override final
  {
    GEOS_UNUSED_VAR( q, dHydraulicAperture_dNormalJump);
    // compute effective normal stress on the fracture
    real64 dStress_dPressure = -1.0;
    real64 const fractureStress = computeFractureStress(pressure, normal, dStress_dPressure);
    // compute new aperture using Barton Bandis model
    real64 dAperture_dStress = -1.0;
    real64 hydraulicAperture = computeHydraulicAperture(pressure, fractureStress, normal, dAperture_dStress, k);
    //newHydraulicAperture = hydraulicAperture;
    
    real64 dPerm_dHydraulicAperture = -1.0;
    compute( oldHydraulicAperture,
             hydraulicAperture, //newHydraulicAperture,
             //dHydraulicAperture_dNormalJump,
             m_permeability[k][0],
             //m_dPerm_dDispJump[k][0],
             dPerm_dHydraulicAperture );

    real64 const dPerm_dPressure = dPerm_dHydraulicAperture * dAperture_dStress * dStress_dPressure;
    //if ( k == 0 || k == 10 )
    {
      GEOS_LOG_RANK_0( GEOS_FMT("  pressure                 {:.15e}", pressure) );
      GEOS_LOG_RANK_0( GEOS_FMT("  dStress_dPressure        {:.15e} fractureStress {:.15e}", dStress_dPressure, fractureStress ) );
      GEOS_LOG_RANK_0( GEOS_FMT("  dAperture_dStress        {:.15e} apertureBB {:.15e}", dAperture_dStress , hydraulicAperture) );
      GEOS_LOG_RANK_0( GEOS_FMT("  dPerm_dHydraulicAperture {:.15e} m_permeability[5][0][0] {:.15e}", dPerm_dHydraulicAperture, m_permeability[k][0][0] ) );
      GEOS_LOG_RANK_0( GEOS_FMT("  dPerm_dPressure          {:.15e} apertureField {:.15e}", dPerm_dPressure, newHydraulicAperture) );
    
    
      for( localIndex i=0; i < m_permeability[k][0].size(); i++ ) // size = 3
      {
        m_dPerm_dPressure[k][0][i] = dPerm_dPressure; 
        GEOS_LOG_RANK_0( GEOS_FMT("  m_dPerm_dPressure[{}][0][{}]   {:.15e} ", k, i, m_dPerm_dPressure[k][0][i]) );
      }
    
    }
  }

  /*GEOS_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & dHydraulicAperture_dNormalJump,
                                                       real64 const & pressure,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const override final
  {
    GEOS_UNUSED_VAR( dispJump, traction, pressure );

    updateFromAperture( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump );
  }*/

private:

  //arrayView4d< real64 > m_dPerm_dDispJump;
  int m_numDimensionsToUpdate;


  real64 m_aperture0;
  
  real64 m_biot;
  real64 m_poisson;
  real64 m_normalStiffness; // Kni
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
  R1Tensor m_referenceEffectiveStress; // sigma_0

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeFractureStress( real64 const pressure, array1d< real64 > const & normal, real64 & dStress_dPressure ) const
  {  
    // TODO: remove this lambda expression
    auto matmul = [](real64 const (&u)[3], array1d< real64 > const &v, array1d< real64 > &r) -> void
    {
      r[0] = u[0]*v[0];
      r[1] = u[1]*v[1];
      r[2] = u[2]*v[2];
    }; 
    
    real64 const deltaSigmaZ = m_biot * (pressure - m_referencePressure);
    real64 const poisson_deltaSigma = deltaSigmaZ * m_poisson/(1.0 - m_poisson);
    // sigma: matrix diagonal
    real64 effectiveStress[3] = { m_referenceEffectiveStress[0] - poisson_deltaSigma,
                                  m_referenceEffectiveStress[1] - poisson_deltaSigma,
                                  m_referenceEffectiveStress[2] - deltaSigmaZ };
    array1d< real64 > effectiveStressOnFracture(3); // sigma_c
    matmul(effectiveStress, normal, effectiveStressOnFracture);
    real64 normalComponentOfStressOnFracture = dot(effectiveStressOnFracture, normal); // sigmaN_N
    // derivative
    dStress_dPressure = -m_biot;

    return normalComponentOfStressOnFracture;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeHydraulicAperture( real64 const pressure, real64 const normalComponentOfStressOnFracture, 
                                   array1d< real64 > const & normal, real64 & dAperture_dStress, int k ) const
  {
    real64 const biot_pressure = m_biot * m_referencePressure; // biot is alpha in the equations

    // Computation of maximum fracture closure (Barton-Bandis parameter)
    // Fracture traction via Terzaghi's Principle
    real64 const sigma_c0[3] = { m_referenceTotalStress[0] * normal[0] - biot_pressure * normal[0],
                                m_referenceTotalStress[1] * normal[1] - biot_pressure * normal[1],
                                m_referenceTotalStress[2] * normal[2] - biot_pressure * normal[2] };
    real64 const sigma_n0 = sigma_c0[0]*normal[0] + 
                            sigma_c0[1]*normal[1] + 
                            sigma_c0[2]*normal[2];
    real64 const g0 = (-m_normalStiffness*m_aperture0 + 
                        std::sqrt((m_normalStiffness*m_aperture0)*
                        (m_normalStiffness*m_aperture0) + 
                        4.0*m_normalStiffness*sigma_n0*m_aperture0)) / (2.0*m_normalStiffness);
    real64 const maximumFractureClosure = g0 + m_aperture0; // Vm -> a_m -> aperture at free-stress state

    // Normal effective stress on the fracture
    real64 const fractureClosure = normalComponentOfStressOnFracture*maximumFractureClosure/(m_normalStiffness*maximumFractureClosure + normalComponentOfStressOnFracture); // gn_BB

    // Compute the new aperture which is equal to the aperture at the free-stress state 
    // minus the closure from the free-stress state to the current state
    real64 const newHydraulicAperture = maximumFractureClosure - fractureClosure;

    real64 const Kni_apert = m_normalStiffness*maximumFractureClosure;
    real64 const Kni_aper_stress = Kni_apert + normalComponentOfStressOnFracture;
    dAperture_dStress = -(Kni_apert*maximumFractureClosure) / (Kni_aper_stress * Kni_aper_stress);

    //if ( k == 0 || k ==10 )
    {
      GEOS_LOG_RANK_0( GEOS_FMT("k {}     m_aperture0 {}", k, m_aperture0) );
      GEOS_LOG_RANK_0( GEOS_FMT("  sigmaN_N               {:.15e}", normalComponentOfStressOnFracture) );
      GEOS_LOG_RANK_0( GEOS_FMT("  maximumFractureClosure {:.15e}", maximumFractureClosure) );
      GEOS_LOG_RANK_0( GEOS_FMT("  fractureClosure        {:.15e}", fractureClosure) );
      GEOS_LOG_RANK_0( GEOS_FMT("  g0   {:.15e} sigma_n0 {:.15e}", g0, sigma_n0) );
    }


    return newHydraulicAperture;
  }
};


class BartonBandisPermeability : public PermeabilityBase
{
public:

  BartonBandisPermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "BartonBandisPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /*virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;*/

  virtual void initializeState() const override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BartonBandisPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          //m_dPerm_dDispJump,
                          m_updateTransversalComponent, 
                          m_aperture0, 
                          m_biot, 
                          m_poisson, 
                          m_normalStiffness, 
                          m_referencePressure, 
                          m_referenceTotalStress);
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * transversalPermeabilityString() { return "transversalPermeability"; }

    /// string/key for aperture under zero normal stress
    static constexpr char const * apertureZeroString() { return "referenceAperture"; }
    static constexpr char const * biotString()                  { return "biot"; }
    static constexpr char const * poissonString()               { return "poisson"; }
    static constexpr char const * normalStiffnessString()       { return "normalStiffness"; }
    static constexpr char const * referencePressureString()     { return "referencePressure"; }
    static constexpr char const * referenceTotalStressString()  { return "referenceTotalStress"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  array3d< real64 > m_dPerm_dAperture;

  /// Derivative of fracture permeability w.r.t. displacement jump
  //array4d< real64 > m_dPerm_dDispJump;

  real64 m_transversalPermeability;

  bool m_updateTransversalComponent;

  /// Reference hydraulic aperture. Aperture at zero normal stress
  real64 m_aperture0;  /// TODO: this will replace what is currently called defaultAperture.
  real64 m_biot;
  real64 m_poisson;
  real64 m_normalStiffness; // Kni  
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_BARTONBANDISPERMEABILITY_HPP_
