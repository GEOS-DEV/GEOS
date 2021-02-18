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
 * @file CarmanKozenyPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_

#include "constitutive/permeability/CarmanKozenyPermeability.hpp"


namespace geosx
{
namespace constitutive
{

class CarmanKozenyPermeabilityUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_permeability.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_permeability.size( 1 ); }


  CarmanKozenyPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                  arrayView3d< real64 > const & dPerm_dPorosity )
    : m_permeability( permeability ),
      m_dPerm_dPorosity( dPerm_dPorosity )
  {}

  /// Default copy constructor
  CarmanKozenyPermeabilityUpdate( FracturePermeabilityUpdate const & ) = default;

  /// Default move constructor
  CarmanKozenyPermeabilityUpdate( FracturePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  CarmanKozenyPermeabilityUpdate & operator=( FracturePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  CarmanKozenyPermeabilityUpdate & operator=( FracturePermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( real64 const & effectiveAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void update( localIndex const k,
               localIndex const q,
               real64 const & porosity ) const
  {
    compute( porosity,
             m_permeability[k][q],
             m_dPerm_dPorosity[k][q] );
  }

private:

  arrayView3d< real64 > m_permeability;
  arrayView3d< real64 > m_dPerm_dPorosity;

};


class CarmanKozenyPermeability : public PermeabilityBase
{
public:
  CarmanKozenyPermeability( string const & name, Group * const parent );

  virtual ~CarmanKozenyPermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "FracturePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FracturePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPorosity );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

protected:
  virtual void postProcessInput() override;

private:

  array3d< real64 > m_dPerm_dPorosity;
  real64 m_Dp;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CarmanKozenyPermeabilityUpdate::compute( real64 const & porosity,
                                              arraySlice1d< real64 > const & permeability,
                                              arraySlice1d< real64 > const & dPerm_dPorosity ) const
{
  for ( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = porosity; // carmanKozeny eq. here;
    dPerm_dPorosity[i] = 1;
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
