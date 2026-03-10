/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CohesiveZoneBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_COHESIVEZONEBASE_HPP_
#define GEOS_CONSTITUTIVE_COHESIVEZONEBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "common/logger/Logger.hpp"
#include "common/DataTypes.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Base class for all continuum constitutive kernel wrapper classes.
 *
 * The responsibility of this base is to:
 *
 * 1) Contain views to state and parameter data for solid models.
 * 2) Specify an interface for state update functions.
 *
 * In general, the ArrayView data in the wrapper is specified to be of type
 * "arrayView<T> const" or "arrayView<T const> const". The "const-ness"
 * of the data indicates whether it is a parameter" or a state variable,
 * with the parameters being "T const" and state variables being "T".
 *
 * @note
 * If an allocation occurs on the underlying Array after a KernelWrapper is created,
 * then the ArrayView members of that KernelWrapper are silently invalid.
 */
class CohesiveZoneBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] newTraction The new traction data from the cohesive zone model class.
   * @param[in] oldTraction The old traction data from the cohesive zone model class.
   */
  CohesiveZoneBaseUpdates( arrayView2d< real64 > const & newNormalStress,
                           arrayView2d< real64 > const & newShearStress,
                           arrayView2d< real64 > const & oldNormalStress,
                           arrayView2d< real64 > const & oldShearStress ):
    m_newNormalStress( newNormalStress ),
    m_newShearStress( newShearStress ),
    m_oldNormalStress( oldNormalStress ),
    m_oldShearStress( oldShearStress )
  {}

  /// Deleted default constructor
  CohesiveZoneBaseUpdates() = delete;

  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  CohesiveZoneBaseUpdates( CohesiveZoneBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  CohesiveZoneBaseUpdates( CohesiveZoneBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  CohesiveZoneBaseUpdates & operator=( CohesiveZoneBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  CohesiveZoneBaseUpdates & operator=( CohesiveZoneBaseUpdates && ) =  delete;

public:

  /// A reference the current material traction at quadrature points.
  arrayView2d< real64 > const m_newNormalStress;

  /// A reference the current material traction at quadrature points.
  arrayView2d< real64 > const m_newShearStress;

  /// A reference the previous material traction at quadrature points.
  arrayView2d< real64 > const m_oldNormalStress;

    /// A reference the previous material traction at quadrature points.
  arrayView2d< real64 > const m_oldShearStress;

  /**
   * @name Update Interfaces: Traction
   *
   * This group of interfaces returns traction base on cohesive zone jump displacements
   */
  ///@{

  /**
   * @brief Jump displacement update.
   *
   * @param[in] k Element index.
   * @param[in] normalDisplaceent The displacement normal to the cohesive interface
   * @param[in] tangentialDisplacement The dispalcement tangential to the cohesive interface
   * @param[out] traction The cohesive traction
   */
  GEOS_HOST_DEVICE
  /**
   * this function is not virtual to avoid a compilation warning with nvcc.
   */
  void jumpDisplacementUpdate( localIndex const k,
                               real64 const & normalDisplacement,
                               real64 const & tangentialDisplacement,
                               real64 & normalStress,
                               real64 & shearStress ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( normalDisplacement );
    GEOS_UNUSED_VAR( tangentialDisplacement );
    GEOS_UNUSED_VAR( normalStress );
    GEOS_UNUSED_VAR( shearStress );
    GEOS_ERROR( "jumpDisplacementUpdate() not implemented for this model" );
  }

  /**
   * @brief Save converged state data at index (k,q)
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   */
  GEOS_HOST_DEVICE
  inline
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const
  {
    m_oldNormalStress[k][q] = m_newNormalStress[k][q];
    m_oldShearStress[k][q] = m_newShearStress[k][q];
  }

  ///@}
};


/**
 * @class CohesiveZoneBase
 * This class serves as the base class for cohesive zone constitutive models.
 */
class CohesiveZoneBase : public constitutive::ConstitutiveBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the CohesiveZoneBase object in the repository.
   * @param parent The parent group of the CohesiveZoneBase object.
   */
  CohesiveZoneBase( string const & name,
                    Group * const parent );

  /**
   * Destructor
   */
  virtual ~CohesiveZoneBase() override = default;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CohesiveZoneBase";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * newNormalStressString() { return "normalStress"; }
    static constexpr char const * newShearStressString() { return "shearStress"; }
    static constexpr char const * oldNormalStressString() { return "oldNormalStress"; }
    static constexpr char const * oldShearStressString() { return "oldShearStress"; }
  };

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Accessors
   */
  ///@{

  /**
   * @brief Non-const/mutable accessor for normal stress
   * @return Accessor
   */
  arrayView2d< real64 > const getNormalStress()
  {
    return m_newNormalStress;
  }

  /**
   * @brief Const/non-mutable accessor for normal stress
   * @return Accessor
   */
  arrayView2d< real64 const > const getNormalStress() const
  {
    return m_newNormalStress;
  }

  /**
   * @brief Non-const/mutable accessor for old normal stress
   * @return Accessor
   */
  arrayView2d< real64 > const getOldNormalStress()
  {
    return m_oldNormalStress;
  }

  /**
   * @brief Const/non-mutable accessor for old normal stress
   * @return Accessor
   */
  arrayView2d< real64 const > const getOldNormalStress() const
  {
    return m_oldNormalStress;
  }

  /**
   * @brief Non-const/mutable accessor for shear stress
   * @return Accessor
   */
  arrayView2d< real64 > const getShearStress()
  {
    return m_newShearStress;
  }

  /**
   * @brief Const/non-mutable accessor for shear stress
   * @return Accessor
   */
  arrayView2d< real64 const > const getSheearStress() const
  {
    return m_newShearStress;
  }

  /**
   * @brief Non-const/mutable accessor for old shear stress
   * @return Accessor
   */
  arrayView2d< real64 > const getOldShearStress()
  {
    return m_oldShearStress;
  }

  /**
   * @brief Const/non-mutable accessor for old shear stress
   * @return Accessor
   */
  arrayView2d< real64 const > const getOldShearStress() const
  {
    return m_oldShearStress;
  }

  ///@}
  //

protected:
  /// Post-process XML input
  virtual void postInputInitialization() override;

  /// The current normal stress at a quadrature point (i.e. at timestep n, global newton iteration k)
  array2d< real64 > m_newNormalStress; // CC TODO: Check permutation of array

  /// The previous normal stress at a quadrature point (i.e. at timestep (n-1))
  array2d< real64 > m_oldNormalStress; // CC TODO: Check permutation of array

  /// The current shear stress at a quadrature point (i.e. at timestep n, global newton iteration k)
  array2d< real64 > m_newShearStress; // CC TODO: Check permutation of array

  /// The previous shear stress at a quadrature point (i.e. at timestep (n-1))
  array2d< real64 > m_oldShearStress; // CC TODO: Check permutation of array
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_COHESIVEZONEBASE_HPP_ */
