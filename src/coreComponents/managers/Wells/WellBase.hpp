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
 * @file WellBase.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class WellBase : public ObjectManagerBase
{
public:

  explicit WellBase( string const & name, dataRepository::ManagedGroup * const parent );
  ~WellBase() override;

  WellBase() = delete;
  WellBase( WellBase const &) = delete;
  WellBase( WellBase && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return "WellBase"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< WellBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual void FillDocumentationNode() override;

  struct viewKeyStruct
  {

    static constexpr auto referenceDepthString = "referenceDepth";
    static constexpr auto typeString = "type";

    dataRepository::ViewKey referenceDepth = { referenceDepthString };
    dataRepository::ViewKey type = { typeString };

  } viewKeys;

private:

  real64 m_referenceDepth;
  string m_type;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_
