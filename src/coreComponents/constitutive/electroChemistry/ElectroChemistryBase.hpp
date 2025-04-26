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
 * @file ElectroChemistryBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_ELECTROCHEMISTRYBASE_HPP
#define GEOS_CONSTITUTIVE_ELECTROCHEMISTRYBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class for electrochemistry material
 */
class ElectroChemistryBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the electro conductivity updates
   * @param conductivity the array of cell-wise conductivities in the subregion
   */
  ElectroChemistryBaseUpdate(arrayView1d<real64> const& conductivity)
  :
  m_conductivity(conductivity)
  {}

  /**
   * @brief Number of elements storing conductivity data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_conductivity.size();
  }

  /**
   * @brief Get conductivity
   * @param[in] k index of the cell in the subRegion
   * @return the conductivity of element k
   */
  GEOS_HOST_DEVICE
  inline
  real64 getConductivity(localIndex const k) const
  {
    return m_conductivity[k];
  }

protected:

  /// View on the cell-wise conductivities
  arrayView1d<real64> m_conductivity;
};

/**
 * @brief The abstract base class for electro conductivity
 */
class ElectroChemistryBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  ElectroChemistryBase( string const & name, dataRepository::Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElectroChemistryBase() override;

  static string catalogName() { return "ElectroChemistryBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const* conductivityString() {return "conductivity";}
    static constexpr char const* defaultConductivityString() {return "defaultConductivity";}
  };

  // virtual void allocateConstitutiveData( dataRepository::Group & parent,
  //                                        localIndex const numConstitutivePointsPerParentIndex ) override;

  using KernelWrapper = ElectroChemistryBaseUpdate;
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper(m_conductivity);
  }

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  array1d<real64> m_conductivity;

  real64 m_defaultConductivity = 1.0;
};
} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_ELECTROCHEMISTRYBASE_HPP