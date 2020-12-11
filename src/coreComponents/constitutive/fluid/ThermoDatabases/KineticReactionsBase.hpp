
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
 * @file KineticReactionsBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_KINETICREACTIONSBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_KINETICREACTIONSBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

namespace constitutive
{

struct KineticReaction
{
  string name;
  real64 MW;
  real64 density;
  array1d< real64 > stochs;
  array1d< localIndex > basisSpeciesIndices;
  real64 logK;
  real64 E;
  real64 rateConst;
};


class KineticReactionsBase
{
public:

  KineticReactionsBase( const Path & fileName ):
    m_fileName( fileName )
  {}

  virtual ~KineticReactionsBase(){}

  using CatalogInterface = dataRepository::CatalogInterface< KineticReactionsBase, Path const &, string_array const & >;

  static typename CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string GetCatalogName() = 0;

  virtual const array1d< KineticReaction > & GetKineticReactions() const
  {
    return m_kineticReactions;
  }

  virtual localIndex numReaction() const
  {
    return m_kineticReactions.size();
  }

protected:
  Path m_fileName;
  array1d< KineticReaction > m_kineticReactions;

};

typedef std::unique_ptr< KineticReactionsBase > KineticReactions;

}

}

#endif
