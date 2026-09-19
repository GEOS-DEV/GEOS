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
 * @file StrainDependentPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"
#include "LvArray/src/tensorOps.hpp"


namespace geos
{
namespace constitutive
{

class StrainDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  StrainDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                     arrayView3d< real64 > const & dPerm_dPressure,
                                     arrayView3d< real64 > const & dPerm_dPorosity,
                                     arrayView3d< real64 > const & referencePermeability,
                                     R1Tensor const strainDependenceConstants )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_referencePermeability( referencePermeability ),
    m_strainDependenceConstants( strainDependenceConstants )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & referencePorosity,
                real64 const (&referencePermeability)[3],
                R1Tensor const strainDependenceConstants,
                real64 const & currentPorosity,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPorosityAndStrain( localIndex const k,
                                            real64 const & currentPorosity,
                                            real64 const & referencePorosity ) const override
  {
    real64 referencePermeability[3];

    referencePermeability[0] = m_referencePermeability[k][0][0];
    referencePermeability[1] = m_referencePermeability[k][0][1];
    referencePermeability[2] = m_referencePermeability[k][0][2];

    compute( referencePorosity,
             referencePermeability,
             m_strainDependenceConstants,
             currentPorosity,
             m_permeability[k][0],
             m_dPerm_dPorosity[k][0] );
  }

private:

  /// dPermeability_dPorosity
  arrayView3d< real64 > m_dPerm_dPorosity;

  /// Reference permeability
  arrayView3d< real64 > m_referencePermeability;

  /// Volumetric strain dependence coefficients for each permeability component
  R1Tensor m_strainDependenceConstants;
};


class StrainDependentPermeability : public PermeabilityBase
{
public:

  StrainDependentPermeability( string const & name, dataRepository::Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    dataRepository::Group * const parent ) const override;

  static string catalogName() { return "StrainDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = StrainDependentPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dPorosity,
                          m_referencePermeability,
                          m_strainDependenceConstants );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * dPerm_dPorosityString() { return "dPerm_dPorosity"; }
    static constexpr char const * strainDependenceConstantsString() { return "strainDependenceConstants"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
  };

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the reference pressure
  R1Tensor m_referencePermeabilityComponents;

  /// dPermeability_dPorosity
  array3d< real64 > m_dPerm_dPorosity;

  /// Reference permeability
  array3d< real64 > m_referencePermeability;

  /// Volumetric strain dependence coefficients for each permeability component
  R1Tensor m_strainDependenceConstants;
};


GEOS_HOST_DEVICE
inline
void StrainDependentPermeabilityUpdate::compute( real64 const & referencePorosity,
                                                 real64 const (&referencePermeability)[3],
                                                 R1Tensor const strainDependenceConstants,
                                                 real64 const & currentPorosity,
                                                 arraySlice1d< real64 > const & permeability,
                                                 arraySlice1d< real64 > const & dPerm_dPorosity ) const
{
  real64 const porosityRatio = currentPorosity / referencePorosity;

  for( localIndex i = 0; i < permeability.size(); ++i )
  {
    permeability[i] = referencePermeability[i] * pow( porosityRatio, strainDependenceConstants[i] );
    dPerm_dPorosity[i] = strainDependenceConstants[i] * permeability[i] / currentPorosity;
  }
}


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
