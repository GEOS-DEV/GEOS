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
 * @file ConstantPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_CONSTANTPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_CONSTANTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class ConstantPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ConstantPermeabilityUpdate( real64 const pressureDependenceConstant,
                              arrayView1d< real64 const > const & referencePressure,
                              arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & initialPermeability,
                              arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_pressureDependenceConstant( pressureDependenceConstant ),
    m_referencePressure( referencePressure ),
    m_initialPermeability( initialPermeability )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & deltaPressure,
                real64 const pressureDependenceConstant,
                real64 const (&initialPermeability)[3],
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressure( localIndex const k,
                                   localIndex const q,
                                   real64 const & pressure_n,
                                   real64 const & pressure ) const override
  {
    GEOS_UNUSED_VAR( q, pressure );

    real64 const deltaPressure = pressure_n - m_referencePressure[k];

    real64 initialPermeability[3];

    initialPermeability[0] = m_initialPermeability[k][0][0];
    initialPermeability[1] = m_initialPermeability[k][0][1];
    initialPermeability[2] = m_initialPermeability[k][0][2];

    compute( deltaPressure,
             m_pressureDependenceConstant,
             initialPermeability,
             m_permeability[k][0],
             m_dPerm_dPressure[k][0] );
  }

private:

  /// Pressure dependence constant
  real64 m_pressureDependenceConstant;

  /// Reference pressure
  arrayView1d< real64 const > const m_referencePressure;

  arrayView3d< real64 > m_initialPermeability;

};


class ConstantPermeability : public PermeabilityBase
{
public:

  ConstantPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void initializeState() const override;

  static string catalogName() { return "ConstantPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_pressureDependenceConstant,
                          m_referencePressure,
                          m_permeability,
                          m_initialPermeability,
                          m_dPerm_dPressure );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * permeabilityComponentsString() { return "permeabilityComponents"; }
    static constexpr char const * pressureDependenceConstantString() { return "pressureDependenceConstant"; }
    static constexpr char const * defaultReferencePressureString() { return "defaultReferencePressure"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * initialPermeabilityString() { return "initialPermeability"; }
  } viewKeys;

  virtual void initializeState() const override final;

protected:

  virtual void postProcessInput() override;

private:

  R1Tensor m_permeabilityComponents;

  real64 m_pressureDependenceConstant;

  real64 m_defaultReferencePressure;

  array1d< real64 > m_referencePressure;

  array3d< real64 > m_initialPermeability;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ConstantPermeabilityUpdate::compute( real64 const & deltaPressure,
                                          real64 const pressureDependenceConstant,
                                          real64 const (&initialPermeability)[3],
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const perm = initialPermeability[i] * std::exp( pressureDependenceConstant * deltaPressure );

    permeability[i] = perm;
    dPerm_dPressure[i] = 0;
  }
}

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
