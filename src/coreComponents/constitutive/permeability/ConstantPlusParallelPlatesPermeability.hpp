/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ConstantPlusParallelPlatesPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_CONSTANTPLUSPARALLELPLATESPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_CONSTANTPLUSPARALLELPLATESPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ConstantPlusParallelPlatesPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ConstantPlusParallelPlatesPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                                arrayView3d< real64 > const & dPerm_dPressure,
                                                arrayView3d< real64 > const & dPerm_dAper,
                                                real64 const defaultConductivity )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dAperture( dPerm_dAper ),
    m_defaultConductivity( defaultConductivity )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dAperture ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & oldHydraulicAperture,
                                   real64 const & newHydraulicAperture ) const override final
  {
    GEOSX_UNUSED_VAR( q );
 
    compute( oldHydraulicAperture,
             newHydraulicAperture,
             m_permeability[k][0],
             m_dPerm_dAperture[k][0] );
  }

private:

  arrayView3d< real64 > m_dPerm_dAperture;

  /// Default conductivity
  real64 m_defaultConductivity;

};


class ConstantPlusParallelPlatesPermeability : public PermeabilityBase
{
public:

  ConstantPlusParallelPlatesPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ConstantPlusParallelPlatesPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantPlusParallelPlatesPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dAperture,
                          m_defaultConductivity );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * defaultConductivityString() { return "defaultConductivity"; }
  } viewKeys;

private:

  array3d< real64 > m_dPerm_dAperture;

  /// Default conductivity
  real64 m_defaultConductivity;

};

GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ConstantPlusParallelPlatesPermeabilityUpdate::compute( real64 const & oldHydraulicAperture,
                                                            real64 const & newHydraulicAperture,
                                                            arraySlice1d< real64 > const & permeability,
                                                            arraySlice1d< real64 > const & dPerm_dAperture ) const
{

  GEOSX_UNUSED_VAR( oldHydraulicAperture );

  real64 const perm  = m_defaultConductivity + newHydraulicAperture*newHydraulicAperture*newHydraulicAperture / 12.0;
  real64 const dPerm = newHydraulicAperture*newHydraulicAperture / 4.0;

  for( int dim=0; dim < 3; dim++ )
  {
    permeability[dim]     = perm;
    dPerm_dAperture[dim]  = dPerm;
  }
}

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_CONSTANTPLUSPARALLELPLATESPERMEABILITY_HPP_
