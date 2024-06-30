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
 * @file ContactBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_CONTACT_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_CONTACT_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class Contact
 *
 * This class is used for in-kernel contact relation updates
 */
template< typename FRICTION_TYPE, typename APERTURE_TYPE >
class ContactUpdates
{
public:

  ContactUpdates( real64 const & displacementJumpThreshold )
    :
  {}

  /// Default copy constructor
  ContactUpdates( ContactBaseUpdates const & ) = default;

  /// Default move constructor
  ContactUpdates( ContactBaseUpdates && ) = default;

  /// Deleted default constructor
  ContactUpdates() = default;

  /// Deleted copy assignment operator
  ContactUpdates & operator=( ContactBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  ContactUpdates & operator=( ContactBaseUpdates && ) =  delete;

  /**
   * @brief Update the traction with the pressure term
   * @param[in] pressure the pressure term
   * @param[in] isOpen a flag specifying whether the fracture is open or closed
   * @param[inout] traction the current tractionVector
   * @param[out] dTraction_dPressure the derivative of the fist component of traction wrt pressure
   * @return the updated traction
   */
  GEOS_HOST_DEVICE
  virtual void addPressureToTraction( real64 const & pressure,
                                      arraySlice1d< real64 >const & tractionVector,
                                      real64 & dTraction_dPressure ) const;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFractures object through the object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FRICTION_TYPE, CoulombFriction > )
    {
      return "Coulomb";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFractures";
    }
  }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


protected:

  typename FRICTION_TYPE::KernelWrapper const m_frictionUpdate;

  typename APERTURE_TYPE::KernelWrapper const m_apertureUpdate;
};


/**
 * @class Contact
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
template< typename FRICTION_TYPE, typename APERTURE_TYPE >
class Contact : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  Contact( string const & name,
           Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ContContactactBase() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ContactBaseUpdates;

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
    /// string/key for friction model
    static constexpr char const * frictionModelModelNameString() { return "frictionModelName"; }
    /// string/key for aperture model
    static constexpr char const * apertureModelNameString() { return "apertureModelName"; }
  };


  FRICTION_TYPE const & getFrictionModel() const
  { return this->getParent().template getGroup< FRICTION_TYPE >( m_solidModelName ); }

  APERTURE_TYPE const & getApertureModel() const
  { return this->getParent().template getGroup< APERTURE_TYPE >( m_porosityModelName ); }

private:

  /// the name of the friction model
  string m_frictionModelName;

  /// the name of the hydraulic aperture model
  string m_apertureModelName;

  /**
   * @brief get a Friction base constant reference to the friction model
   * return a constant FrictionBase reference to the friction model
   */
  FrictionBase const & getBaseFrictionModel() const
  { return this->getParent().template getGroup< PermeabilityBase >( m_frictionModelName ); }

  /**
   *@brief get a ApertureBase constant reference to the aperture model
   * return a constant ApertureBase reference to the aperture model
   */
  ApertureBase const & getBaseApertureModel() const
  { return this->getParent().template getGroup< SolidBase >( m_apertureModelName ); }

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ContactUpdates::addPressureToTraction( real64 const & pressure,
                                            arraySlice1d< real64 > const & tractionVector,
                                            real64 & dTraction_dPressure ) const
{
  tractionVector[0] -= pressure;
  dTraction_dPressure = -1.0;
}


} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_CONTACT_HPP_ */
