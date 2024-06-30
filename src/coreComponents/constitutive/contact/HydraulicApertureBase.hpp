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
 * @file HydraulicApertureBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTUREBASE_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTUREBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class HydraulicApertureBaseUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class HydraulicApertureBaseUpdates
{
public:

  HydraulicApertureBaseUpdates( TableFunction const & apertureTable )
    : m_apertureTable( apertureTable.createKernelWrapper() )
  {}

  /// Default copy constructor
  HydraulicApertureBaseUpdates( HydraulicApertureBaseUpdates const & ) = default;

  /// Default move constructor
  HydraulicApertureBaseUpdates( HydraulicApertureBaseUpdates && ) = default;

  /// Deleted default constructor
  HydraulicApertureBaseUpdates() = default;

  /// Deleted copy assignment operator
  HydraulicApertureBaseUpdates & operator=( HydraulicApertureBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  HydraulicApertureBaseUpdates & operator=( HydraulicApertureBaseUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  virtual real64 computeHydraulicAperture( real64 const aperture,
                                           real64 & dHydraulicAperture_dAperture ) const;

protected:

  /// The aperture table function wrapper
  TableFunction::KernelWrapper m_apertureTable;
};


/**
 * @class HydraulicApertureBase
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class HydraulicApertureBase : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  HydraulicApertureBase( string const & name,
                         Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~HydraulicApertureBase() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;



  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = HydraulicApertureBaseUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    /// string/key for aperture tolerance
    static constexpr char const * apertureToleranceString() { return "apertureTolerance"; }

    /// string/key for aperture table name
    static constexpr char const * apertureTableNameString() { return "apertureTableName"; }
  };

protected:
  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Validate the values provided in the aperture table
   * @param[in] apertureTable the effective aperture vs aperture table
   */
  void validateApertureTable( TableFunction const & apertureTable ) const;

  /// The aperture tolerance to avoid floating point errors in expressions involving aperture
  real64 m_apertureTolerance;

  /// The name of the aperture table, if any
  string m_apertureTableName;

  /// Pointer to the function that limits the model aperture to a physically admissible value.
  TableFunction const * m_apertureTable;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 HydraulicApertureBaseUpdates::computeHydraulicAperture( real64 const aperture,
                                                               real64 & dHydraulicAperture_dAperture ) const
{
  return m_apertureTable.compute( &aperture, &dHydraulicAperture_dAperture );
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_ */
