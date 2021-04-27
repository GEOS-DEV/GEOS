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
                                  arrayView3d< real64 > const & dPerm_dPorosity,
                                  real64 const particleDiameter,
                                  real64 const sphericity )
    : PermeabilityBaseUpdate( permeability ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_particleDiameter( particleDiameter ),
    m_sphericity( sphericity )
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
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & pressure,
                       real64 const & volStrain,
                       real64 ( dPerm_dVolStrain )[3] ) override
  {
    compute( porosity,
             m_permeability[k][q],
             m_dPerm_dPorosity[k][q] );
  }

private:

  /// dPermeability_dPorosity
  arrayView3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;

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
                          m_dPerm_dPorosity,
                          m_particleDiameter,
                          m_sphericity );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dPorosityString() { return "dPerm_dPorosity"; }
    static constexpr char const * particleDiameterString() { return "particleDiameter"; }
    static constexpr char const * sphericityString() { return "sphericity"; }
  } viewKeys;

protected:
  virtual void postProcessInput() override;

private:
  /// dPermeability_dPorosity
  array3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;
};




}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
