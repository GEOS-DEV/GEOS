/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file SimpleWell.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP

#include "WellBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto simpleWell = "SimpleWell";
}
}

class DomainPartition;

class SimpleWell : public WellBase
{
public:
  explicit SimpleWell( string const & name, dataRepository::ManagedGroup * const parent );
  ~SimpleWell() override;

  SimpleWell() = delete;
  SimpleWell( SimpleWell const &) = delete;
  SimpleWell( SimpleWell && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return dataRepository::keys::simpleWell; }

  void FillDocumentationNode() override;

  virtual void InitializationOrder(string_array & order) override;

  virtual void FinalInitialization( ManagedGroup * const group ) override;

  // update each connection pressure from bhp and hydrostatic head
  void StateUpdate( DomainPartition const * domain, localIndex fluidIndex );

  real64 GetTotalFlowRate();

  struct viewKeyStruct : public WellBase::viewKeyStruct
  {

    static constexpr auto pressureString = "pressure";
    static constexpr auto flowRateString = "flowRate";
    static constexpr auto bhpString      = "bhp";

    dataRepository::ViewKey pressure = { pressureString };
    dataRepository::ViewKey flowRate = { flowRateString };
    dataRepository::ViewKey bhp      = { bhpString };

  } viewKeysSimpleWell;

  struct groupKeyStruct : public WellBase::groupKeyStruct
  {

  } groupKeysSimpleWell;

private:

  array1d<real64> m_bhp;

};

} //namespace geosx


#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP
