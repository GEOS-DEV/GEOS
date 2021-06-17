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
 * @file CompressibleSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/NullModel.hpp"

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
template< typename PORO_TYPE >
class CompressibleSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE >
{
public:

  using CoupledSolidUpdates< NullModel, PORO_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE >::m_porosityUpdate;

  /**
   * @brief Constructor
   */
  CompressibleSolidUpdates( NullModel * solidModel,
                            PORO_TYPE * porosityModel ):
    CoupledSolid< NullModel, PORO_TYPE >( solidModel, porosityModel )
  {}

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void pressureUpdate( localIndex const k,
                       localIndex const q,
                       real64 const & pressure,
                       real64 const & deltaPressure ) const
  {
    m_porosityUpdate.updatePorosity( k, q, pressure, deltaPressure );
  }

};

template< typename PORO_TYPE >
class CompressibleSolid : public CoupledSolid< NullModel, PORO_TYPE >
{
public:

  using CoupledSolid< NullModel, PORO_TYPE >::m_solidModel;
  using CoupledSolid< NullModel, PORO_TYPE >::m_porosityModel;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CompressibleSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CompressibleSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "CompressibleSolid" ); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


  /**
   * @brief Create a instantiation of the CompressibleSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of CompressibleSolidUpdates.
   */
  CompressibleSolidUpdates< PORO_TYPE > createKernelUpdates() const
  {

    return CompressibleSolidUpdates< PORO_TYPE >( m_solidModel,
                                                  m_porosityModel );
  }

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
