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
 * @file ButlerVolmerReaction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_BUTLERVOLMERREACTION_HPP
#define GEOS_CONSTITUTIVE_BUTLERVOLMERREACTION_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class for reaction coefficient update
 */
class ButlerVolmerReactionUpdate
{
public:

  /**
   * @brief Constructor for the class performing the reactivity coefficient updates
   * @param reactivityCoeff the array of face-wise reactivity coefficient in the subregion
   */
  ButlerVolmerReactionUpdate( arrayView1d< real64 > const & reactivityCoeff )
    :
    m_reactivityCoeff( reactivityCoeff )
  {}

  /**
   * @brief Number of elements storing reactivity coefficient data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_reactivityCoeff.size();
  }

  /**
   * @brief Get reactivity coefficient
   * @param[in] k index of the element in the subRegion
   * @return the reactivity coefficient of element k
   */
  GEOS_HOST_DEVICE
  inline
  real64 getReactCoeff( localIndex const k ) const
  {
    return m_reactivityCoeff[k];
  }

protected:

  /// View on the face-wise reactivity coefficient
  /// Also commonly referred to as k_rxn in electrochemistry
  arrayView1d< real64 > m_reactivityCoeff;
};

/**
 * @brief The class for interface Butler-Volmer kinetics
 */
class ButlerVolmerInterface : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  ButlerVolmerInterface( string const & name, dataRepository::Group * const parent );

  /**
   * Destructor
   */
  virtual ~ButlerVolmerInterface() override;

  static string catalogName() { return "ButlerVolmerInterface"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * reactivityCoefficientString() {return "reactivityCoefficient";}
    static constexpr char const * defaultReactivityCoefficientString() {return "defaultReactivityCoefficient";}
  };

  using KernelWrapper = ButlerVolmerReactionUpdate;
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_reactivityCoeff );
  }

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  array1d< real64 > m_reactivityCoeff;

  real64 m_defaultKrxn = 1.0;
};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_BUTLERVOLMERREACTION_HPP
