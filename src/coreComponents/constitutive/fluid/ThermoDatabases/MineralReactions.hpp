
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
 * @file MineralReactions.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MINERALREACTIONS_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MINERALREACTIONS_HPP

#include "KineticReactionsBase.hpp"

namespace geosx
{

namespace constitutive
{

class MineralReactions : public KineticReactionsBase
{
public:

  MineralReactions( const Path & fileName,
                    const string_array & basisSpeciesNames );

  ~MineralReactions() override
  {}

  static constexpr auto m_catalogName = "MineralReactions";

  static string CatalogName()                    { return m_catalogName; }

  virtual string GetCatalogName() override final { return CatalogName(); }


private:

  /** read database and generate species properties and stoch matrix */

  void ReadMineralReactions( const string_array & basisSpeciesNames );


};

}

}

#endif
