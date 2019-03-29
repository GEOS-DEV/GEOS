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

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const = 0;


  virtual void StateUpdatePointPressure( real64 const & pres,
                                         localIndex const k,
                                         localIndex const q ) {}

  virtual void resize( localIndex ) override;


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  virtual string GetCatalogName() = 0;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex );

  struct viewKeyStruct
  {
    static constexpr auto poreVolumeMultiplierString  = "poreVolumeMultiplier";
    static constexpr auto dPVMult_dPresString  = "dPVMult_dDensity";

  } m_ConstitutiveBaseViewKeys;

  struct groupKeyStruct
  {} m_ConstitutiveBaseGroupKeys;

  virtual viewKeyStruct       & viewKeys()        { return m_ConstitutiveBaseViewKeys; }
  virtual viewKeyStruct const & viewKeys() const  { return m_ConstitutiveBaseViewKeys; }

  virtual groupKeyStruct       & groupKeys()       { return m_ConstitutiveBaseGroupKeys; }
  virtual groupKeyStruct const & groupKeys() const { return m_ConstitutiveBaseGroupKeys; }

protected:

private:
  ManagedGroup * m_constitutiveDataGroup = nullptr;

  ConstitutiveBase( ConstitutiveBase const & ) = delete;
  ConstitutiveBase( ConstitutiveBase && ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase const & ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase && ) = delete;

};



}
}
#endif
