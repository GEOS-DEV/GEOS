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
 * @file FracturePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"


namespace geosx
{
namespace constitutive
{

class FracturePermeabilityBaseUpdate
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

protected:

  FracturePermeabilityBaseUpdate( arrayView3d< real64 > const & permeability,
                                  arrayView3d< real64 > const & dPerm_dAperture )
    : m_permeability( permeability ),
      m_dPerm_dAperture( dPerm_dAperture )
  {}

  /// Default copy constructor
  FracturePermeabilityBaseUpdate( FracturePermeabilityBaseUpdate const & ) = default;

  /// Default move constructor
  FracturePermeabilityBaseUpdate( FracturePermeabilityBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  FracturePermeabilityBaseUpdate & operator=( FracturePermeabilityBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  FracturePermeabilityBaseUpdate & operator=( FracturePermeabilityBaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const & effectiveAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dAperture ) const = 0;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
               localIndex const q,
               real64 const & effectiveAperture ) const = 0;
private:

  arrayView3d< real64 > m_permeability;
  arrayView3d< real64 > m_dPerm_dAperture;

};


class FracturePermeabilityBase : public ConstitutiveBase
{
public:
  FracturePermeabilityBase( string const & name, Group * const parent );

  virtual ~FracturePermeabilityBase() override;

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
                          m_dPerm_dAperture );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

protected:
  virtual void postProcessInput() override;

private:

  array3d< real64 > m_permeability;
  array3d< real64 > m_dPerm_dAperture;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
