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
 * @file ConstitutiveBase.hpp
 */

#ifndef CONSTITUTIVEBASE_HPP_
#define CONSTITUTIVEBASE_HPP_

#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"
#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

namespace constitutive
{


class ConstitutiveBase : public dataRepository::ManagedGroup
{
public:

  /**
   * Single point of reference for generated constitutive field naming convention
   *
   * @param prefix name prefix (e.g. constitutive model name)
   * @param name actual field name
   * @return prefixed field name that is used to access data
   */
  inline static string makeFieldName(string const & prefix, string const & name) { return prefix + "_" + name; }


  ConstitutiveBase( string const & name,
                    ManagedGroup * const parent );

  virtual ~ConstitutiveBase() override;

  /**
   * @brief create a clone of this constitutive model
   * @param[in]  name   The name of the clone in the repository
   * @param[in]  parent A pointer to the group that contains the instance of the new clone
   * @param[out] clone  A reference to a unique_ptr  that will hold the clone.
   */
  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const = 0;


  virtual void StateUpdatePointPressure( real64 const & pres,
                                         localIndex const k,
                                         localIndex const q ) {}

  /**
   * @brief function to resize the fields in this constitutive model
   * @param[in] newSize the new size of the fields
   */
  virtual void resize( localIndex newSize ) override;


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  /**
   * @brief function to return the catalog name of the derived class
   * @return a string that contains the catalog name of the derived class
   */
  virtual string GetCatalogName() = 0;

  /**
   * @brief Allocate constitutive data and make views to data on parent objects
   * @param[in] parent pointer to the group that holds the constitutive relation
   * @param[in] numConstitutivePointsPerParentIndex number of quadrature points
   *
   * This function does 2 things:
   *   1) Allocate data according to the size of parent and numConstitutivePointsPerParentIndex
   *   2) Create wrappers to the constitutive data in the parent for easier access
   */
  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex );

  struct viewKeyStruct
  {
    static constexpr auto poreVolumeMultiplierString  = "poreVolumeMultiplier";
    static constexpr auto dPVMult_dPresString  = "dPVMult_dDensity";

  };

  struct groupKeyStruct
  {};


  localIndex numQuadraturePoints() const { return m_numQuadraturePoints; }

protected:

private:
  localIndex m_numQuadraturePoints;
  ManagedGroup * m_constitutiveDataGroup = nullptr;

  ConstitutiveBase( ConstitutiveBase const & ) = delete;
  ConstitutiveBase( ConstitutiveBase && ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase const & ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase && ) = delete;

};



}
}
#endif
