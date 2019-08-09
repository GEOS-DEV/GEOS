/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ContactRelationBase.hpp
 */

#ifndef GEOSX_SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_
#define GEOSX_SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "managers/Functions/FunctionBase.hpp"

namespace geosx
{


namespace constitutive
{

/**
 * @class ContactRelationBase
 *
 * This class will serve as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 *
 * Currently this class actually implements something. It should be converted to a base class when
 * it grows up.
 */
class ContactRelationBase : public ConstitutiveBase
{
public:
  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent ManagedGroup that holds this relation object.
   */
  ContactRelationBase( string const & name,
                       ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ContactRelationBase() override;

  /**
   * @brief Name that is used to register this a type of "ContactRelationBase" in the object catalog
   * @return See description
   */
  static string CatalogName() { return "Contact"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override {}

  virtual ManagedGroup * CreateChild( string const & catalogKey,
                                      string const & name ) override;

  virtual void InitializePreSubGroups( ManagedGroup * const ) override;


  /// accessor for penalty stiffness
  inline real64 stiffness() const { return m_penaltyStiffness; }

  /**
   * @brief evaluation of effective aperture given in input model aperture/gap
   * @param[in] aperture the model aperture/gap
   * @return And effective physical aperture that is always > 0
   */
  inline real64 effectiveAperture( real64 const aperture ) const { return m_apertureFunction->Evaluate( & aperture ); }

  /**
   * @brief evaluation of the derivative of the effective physical aperture
   * @param[in] aperture the model aperture/gap
   * @return
   */
  inline real64 dEffectiveAperture_dAperture( real64 const aperture ) const
  {
    real64 aperPlus = aperture ;
    real64 aperMinus = aperture - 1.0e-6;
    real64 slope = (m_apertureFunction->Evaluate( &aperPlus ) - m_apertureFunction->Evaluate( &aperMinus ) ) / 1.0e-6;
    return slope;
  }


  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct: public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto penaltyStiffnessString  = "penaltyStiffness";
  };

private:

  /// the value of penalty to penetration
  real64 m_penaltyStiffness;

  /// pointer to the function that limits the model aperture to a physically admissible value.
  FunctionBase * m_apertureFunction;

};

}
} /* namespace geosx */

#endif /* GEOSX_SRC_CORECOMPONENTS_CONSTITUTIVE_CONTACTRELATIONBASE_HPP_ */
