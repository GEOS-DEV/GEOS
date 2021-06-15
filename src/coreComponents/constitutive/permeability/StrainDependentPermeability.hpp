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
 * @file StrainDependentPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class StrainDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  StrainDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                     arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure )
  {}

  /// Default copy constructor
  StrainDependentPermeabilityUpdate( StrainDependentPermeabilityUpdate const & ) = default;

  /// Default move constructor
  StrainDependentPermeabilityUpdate( StrainDependentPermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  StrainDependentPermeabilityUpdate & operator=( StrainDependentPermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  StrainDependentPermeabilityUpdate & operator=( StrainDependentPermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void updatePressureStrain( localIndex const k,
                                     localIndex const q,
                                     real64 const & pressure,
                                     real64 const & volStrain,
                                     real64 ( & dPerm_dVolStrain )[3] ) const override
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( pressure );
    GEOSX_UNUSED_VAR( volStrain );
    GEOSX_UNUSED_VAR( dPerm_dVolStrain );
  }

private:

};


class StrainDependentPermeability : public PermeabilityBase
{
public:
  StrainDependentPermeability( string const & name, Group * const parent );

  virtual ~StrainDependentPermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "StrainDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = StrainDependentPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure );
  }



protected:
  virtual void postProcessInput() override;

private:

};



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
