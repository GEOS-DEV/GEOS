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
 * @file PorousSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE
 */
template< typename SOLID_TYPE >
class PorousSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity >
{
public:
  /**
   * @brief Constructor
   */
  PorousSolidUpdates( SOLID_TYPE * solidModel,
                      BiotPorosity * porosityModel ):
  CoupledSolidUpdates< SOLID_TYPE, BiotPorosity >( solidModel, porosityModel )
  {}

   /// Deleted default constructor
  PorousSolidUpdates() = delete;

  /// Default copy constructor
  PorousSolidUpdates( PorousSolidUpdates const & ) = default;

  /// Default move constructor
  PorousSolidUpdates( PorousSolidUpdates && ) = default;

  /// Deleted copy assignment operator
  PorousSolidUpdates & operator=( PorousSolidUpdates const & ) = delete;

  /// Deleted move assignment operator
  PorousSolidUpdates & operator=( PorousSolidUpdates && ) =  delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & pressure,
                          real64 const & deltaPressure,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 & dPorosity_dPressure,
                          real64 & dPorosity_dVolStrainIncrement,
                          real64 & dTotalStress_dPressure ) const
  {
    SOLID_TYPE::KernelWrapper::DiscretizationOps stiffness;

    m_solidUpdate.smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

    m_porosityUpdate.updatePorosity( k,
                                     q,
                                     pressure,
                                     deltaPressure,
                                     strainIncrement,
                                     dPorosity_dPressure,
                                     dPorosity_dVolStrainIncrement,
                                     dTotalStress_dPressure );
  }

};


class PorousSolidBase : public ConstitutiveBase
{};

/**
 * @brief Class to represent a coupled solid model
 */
template< typename SOLID_TYPE >
class PorousSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity >
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Poro" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

};

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
