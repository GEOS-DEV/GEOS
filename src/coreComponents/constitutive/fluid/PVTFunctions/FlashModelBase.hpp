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
 * @file FlashModelBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class FlashModelBaseUpdate
{
public:

  FlashModelBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    :
    m_componentMolarWeight( componentMolarWeight )
  {}

  /// Default virtual destructor
  virtual ~FlashModelBaseUpdate() = default;

  /// Default copy constructor
  FlashModelBaseUpdate( FlashModelBaseUpdate const & ) = default;

  /// Default move constructor
  FlashModelBaseUpdate( FlashModelBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  FlashModelBaseUpdate & operator=( FlashModelBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  FlashModelBaseUpdate & operator=( FlashModelBaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & compFraction,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64 > const & dPhaseFraction_dCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFraction_dCompFraction ) const = 0;

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  virtual void move( LvArray::MemorySpace const space, bool const touch )
  {
    m_componentMolarWeight.move( space, touch );
  }

protected:

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

};


class FlashModelBase
{
public:

  FlashModelBase( string const & name,
                  string_array const & componentNames,
                  array1d< real64 > const & componentMolarWeight ):
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~FlashModelBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< FlashModelBase, string_array const &,
                                                             string_array const &,
                                                             string_array const &,
                                                             array1d< real64 > const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & flashModelName() const { return m_modelName; }

protected:

  /// Name the solubility model
  string m_modelName;

  /// Array storing the name of the components
  string_array m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_
