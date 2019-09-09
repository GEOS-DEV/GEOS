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

/*
 * @file SimpleGeometricObjectBase.hpp
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

//#include "common/Common.h"
#include "dataRepository/Group.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "ObjectCatalog.hpp"

class Function;

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const geometricObjects("GeometricObjects");
}
}


class SimpleGeometricObjectBase : public dataRepository::Group
{
public:

  explicit SimpleGeometricObjectBase( std::string const & name,
                                      Group * const parent );

  virtual ~SimpleGeometricObjectBase();

  static string CatalogName() { return "SimpleGeometricObjectBase"; }

  virtual bool IsCoordInObject( const R1Tensor& coord ) const = 0;

  using CatalogInterface = cxx_utilities::CatalogInterface< SimpleGeometricObjectBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};


}
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
