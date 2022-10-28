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
 * @file ExponentialDecayPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_EXPONENTIALDECAYPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_EXPONENTIALDECAYPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ExponentialDecayPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ExponentialDecayPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                      arrayView3d< real64 > const & dPerm_dPressure,
                                      arrayView4d< real64 > const & dPerm_dTraction,
                                      real64 const empiricalConstant,
                                      R1Tensor const & initialPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dTraction( dPerm_dTraction ),
    m_empiricalConstant( empiricalConstant ),
    m_initialPermeability( initialPermeability )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const ( &traction )[3],
                R1Tensor const & initialPermeability,
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dTraction ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & pressure,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const override
  {
    GEOSX_UNUSED_VAR( q, oldHydraulicAperture, newHydraulicAperture, dispJump, pressure );

    compute( traction,
             m_initialPermeability,
             m_permeability[k][0],
             m_dPerm_dTraction[k][0] );
  }

private:

  /// Derivative of fracture permeability to effective normal stress applied on the fracture surfaces
  arrayView4d< real64 > m_dPerm_dTraction;

  /// An empirical constant
  real64 m_empiricalConstant;

  /// Initial permeability tensor
  R1Tensor m_initialPermeability;

};


class ExponentialDecayPermeability : public PermeabilityBase
{
public:

  ExponentialDecayPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ExponentialDecayPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ExponentialDecayPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dTraction,
                          m_empiricalConstant,
                          m_initialPermeability );
  }

private:

  /// Derivative of fracture permeability to traction acting on fracture surfaces
  array4d< real64 > m_dPerm_dTraction;

  /// An empirical constant
  real64 m_empiricalConstant;

  /// Initial permeability tensor
  R1Tensor m_initialPermeability;

  struct viewKeyStruct
  {
    static constexpr char const * empiricalConstantString() { return "empiricalConstant"; }
    static constexpr char const * initialPermeabilityString() { return "initialPermeability"; }
  };

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ExponentialDecayPermeabilityUpdate::compute( real64 const ( &traction )[3],
                                                  R1Tensor const & initialPermeability,
                                                  arraySlice1d< real64 > const & permeability,
                                                  arraySlice2d< real64 > const & dPerm_dTraction ) const
{
  real64 const effNormalStress = -traction[0];

  real64 const permMultiplier = std::exp( -m_empiricalConstant * effNormalStress );

  real64 const dpermMultiplier_dTraction = std::exp( -m_empiricalConstant * effNormalStress ) * m_empiricalConstant;


  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permMultiplier * initialPermeability[i];
    dPerm_dTraction[i][0] = initialPermeability[i] * dpermMultiplier_dTraction;
    dPerm_dTraction[i][1] = 0.0;
    dPerm_dTraction[i][2] = 0.0;
  }
}



} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
