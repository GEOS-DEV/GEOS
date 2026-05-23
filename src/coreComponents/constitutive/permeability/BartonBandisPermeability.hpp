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
#include "LvArray/src/genericTensorOps.hpp"


namespace geos
{
namespace constitutive
{

class BartonBandisPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  BartonBandisPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                  arrayView3d< real64 > const & dPerm_dPressure,
                                  bool const updateTransversalComponent,
                                  real64 const aperture0,
                                  real64 const biot,
                                  real64 const poisson,
                                  real64 const normalStiffness,
                                  real64 const referencePressure,
                                  R1Tensor const &referenceTotalStress )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_numDimensionsToUpdate( 3 ),
    m_aperture0( aperture0 ),
    m_biotCoefficient( biot ),
    m_poissonRatio( poisson ),
    m_normalStiffness( normalStiffness ), // Kni
    m_referencePressure( referencePressure ),
    m_referenceTotalStress( referenceTotalStress )
  {
    m_numDimensionsToUpdate = updateTransversalComponent ? 3 : 2;

    m_referenceEffectiveStress[0] = m_referenceTotalStress[0] - m_biotCoefficient*m_referencePressure; 
    m_referenceEffectiveStress[1] = m_referenceTotalStress[1] - m_biotCoefficient*m_referencePressure; 
    m_referenceEffectiveStress[2] = m_referenceTotalStress[2] - m_biotCoefficient*m_referencePressure; 
  }

  GEOS_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                arraySlice1d< real64 > const & permeability,
                real64 & dPerm_dHydraulicAperture  ) const
  {
    GEOS_UNUSED_VAR( oldHydraulicAperture );

    real64 const perm  = newHydraulicAperture*newHydraulicAperture / 12.0;
    dPerm_dHydraulicAperture = newHydraulicAperture / 6.0;

    for( int dim=0; dim < m_numDimensionsToUpdate; dim++ )
    {
      permeability[dim] = perm;
    }
  }

  GEOS_HOST_DEVICE
  virtual void updateFromPressureApertureAndNormal( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & pressure,
                                                    real64 const & oldHydraulicAperture,
                                                    real64 const & newHydraulicAperture,
                                                    arraySlice1d< real64 const > const & normal,
                                                    real64 const & dHydraulicAperture_dNormalJump ) const override final
  {
    GEOS_UNUSED_VAR( q, dHydraulicAperture_dNormalJump);
    // compute effective normal stress on the fracture
    real64 dStress_dPressure = -1.0;
    real64 const fractureStress = computeFractureStress(pressure, normal, dStress_dPressure);
    // compute new aperture using Barton Bandis model
    real64 dAperture_dStress = -1.0;
    real64 hydraulicAperture = computeHydraulicAperture(pressure, fractureStress, normal, dAperture_dStress, k);
    
    real64 dPerm_dHydraulicAperture = -1.0;
    compute( oldHydraulicAperture,
             hydraulicAperture, 
             m_permeability[k][0],
             dPerm_dHydraulicAperture );

    real64 const dPerm_dPressure = dPerm_dHydraulicAperture * dAperture_dStress * dStress_dPressure;
    for( localIndex i=0; i < m_permeability[k][0].size(); i++ ) // size = 3
    {
      m_dPerm_dPressure[k][0][i] = dPerm_dPressure; 
    }
  }



  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeHydraulicAperture( real64 const pressure, real64 const normalComponentOfStressOnFracture, 
                                   arraySlice1d< real64 const > const & normal, real64 & dAperture_dStress, int k ) const
  {
    real64 const referenceTotalStress[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3 (m_referenceTotalStress); 
    real64 const biot_pressure = m_biotCoefficient * m_referencePressure; // biot is alpha in the equations

    // Computation of maximum fracture closure (Barton-Bandis parameter)
    // Fracture traction via Terzaghi's Principle
    real64 sigma_c0[3] = {0.0};
    LvArray::tensorOps::hadamardProduct< 3 >( sigma_c0, referenceTotalStress, normal );
    LvArray::tensorOps::scaledAdd< 3 >(sigma_c0, normal, -biot_pressure);
    
    real64 const sigma_n0 = LvArray::tensorOps::AiBi< 3 >( sigma_c0, normal );

    // \dfrac{-K_{ni}a_0 + \sqrt{(K_{ni}a_0)^2 + 4K_{ni}a_0\sigma_{n0}}}{2K_{ni}}.
    real64 const normalStiffApertureProduct = m_normalStiffness * m_aperture0;
    real64 const sqrtTerm = normalStiffApertureProduct * (normalStiffApertureProduct + 4.0*sigma_n0) ;
    real64 const g0 = ( -normalStiffApertureProduct + std::sqrt(sqrtTerm) ) / (2.0 * m_normalStiffness);

    real64 const maximumFractureClosure = g0 + m_aperture0; // V_m or a_m -> aperture at stress-free state

    // Normal effective stress on the fracture
    // g_n(\sigma_n) = \dfrac{\sigma_n * V_m}{K_{ni} * V_m + \sigma_n}
    real64 const fractureClosure = (normalComponentOfStressOnFracture*maximumFractureClosure) / (m_normalStiffness*maximumFractureClosure + normalComponentOfStressOnFracture); 

    real64 const newHydraulicAperture = maximumFractureClosure - fractureClosure;

    // derivative
    // \frac{da}{d\sigma}(\sigma) = -\dfrac{K_{ni} V_m^2}{(K_{ni} V_m + \sigma)^2}
    real64 const normalStiffMaxClosureProduct = m_normalStiffness*maximumFractureClosure;
    real64 const denom = normalStiffMaxClosureProduct + normalComponentOfStressOnFracture;
    dAperture_dStress = -(normalStiffMaxClosureProduct*maximumFractureClosure) / (denom * denom);

    return newHydraulicAperture;
  }

private:

  int m_numDimensionsToUpdate;


  real64 m_aperture0;
  
  real64 m_biotCoefficient;
  real64 m_poissonRatio;
  real64 m_normalStiffness; // Kni
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
  R1Tensor m_referenceEffectiveStress; // sigma_0

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeFractureStress( real64 const pressure, arraySlice1d< real64 const > const & normal, real64 & dStress_dPressure ) const
  {  
    //real64 const normal[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3 (normal_);
    
    real64 const deltaSigmaZ = m_biotCoefficient * (pressure - m_referencePressure);
    real64 const poisson_deltaSigma = deltaSigmaZ * m_poissonRatio/(1.0 - m_poissonRatio);
    // sigma: matrix diagonal
    real64 effectiveStress[3] = { m_referenceEffectiveStress[0] - poisson_deltaSigma,
                                  m_referenceEffectiveStress[1] - poisson_deltaSigma,
                                  m_referenceEffectiveStress[2] - deltaSigmaZ };
    real64 effectiveStressOnFracture[3] = {0.0}; // sigma_c
    LvArray::tensorOps::hadamardProduct< 3 >( effectiveStressOnFracture, normal, effectiveStress );
    real64 const normalComponentOfStressOnFracture = LvArray::tensorOps::AiBi< 3 >(effectiveStressOnFracture, normal); // sigmaN_N
    
    // derivative 
    dStress_dPressure = -m_biotCoefficient;

    return normalComponentOfStressOnFracture;
  }

};


class BartonBandisPermeability : public PermeabilityBase
{
public:

  BartonBandisPermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "BartonBandisPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

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
                          m_updateTransversalComponent, 
                          m_aperture0, 
                          m_biotCoefficient, 
                          m_poissonRatio, 
                          m_normalStiffness, 
                          m_referencePressure, 
                          m_referenceTotalStress);
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * transversalPermeabilityString() { return "transversalPermeability"; }

    /// string/key for aperture under zero normal stress
    static constexpr char const * apertureZeroString() { return "referenceAperture"; }
    static constexpr char const * biotCoefficientString()       { return "biotCoefficient"; }
    static constexpr char const * poissonRatioString()          { return "poissonRatio"; }
    static constexpr char const * normalStiffnessString()       { return "normalStiffness"; }
    static constexpr char const * referencePressureString()     { return "referencePressure"; }
    static constexpr char const * referenceTotalStressString()  { return "referenceTotalStress"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  real64 m_transversalPermeability;

  bool m_updateTransversalComponent;

  /// Reference hydraulic aperture. Aperture at zero normal stress
  real64 m_aperture0;  /// TODO: this will replace what is currently called defaultAperture.
  real64 m_biotCoefficient;
  real64 m_poissonRatio;
  real64 m_normalStiffness; // Kni  
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_BARTONBANDISPERMEABILITY_HPP_
