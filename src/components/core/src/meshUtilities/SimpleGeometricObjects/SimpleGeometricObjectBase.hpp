// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * SimpleGeometricObjects.h
 *
 *  Created on: Dec 4, 2012
 *      Author: settgast1
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

//#include "common/Common.h"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ManagedGroup.hpp"
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


class SimpleGeometricObjectBase : public dataRepository::ManagedGroup
{
public:

  explicit SimpleGeometricObjectBase( std::string const & name,
                                      ManagedGroup * const parent );

  virtual ~SimpleGeometricObjectBase();

  static string CatalogName() { return "SimpleGeometricObjectBase"; }

  SimpleGeometricObjectBase() = default;
  SimpleGeometricObjectBase( SimpleGeometricObjectBase const & ) = default;
  SimpleGeometricObjectBase( SimpleGeometricObjectBase &&) = default;
  SimpleGeometricObjectBase& operator=( SimpleGeometricObjectBase const & ) = default;
  SimpleGeometricObjectBase& operator=( SimpleGeometricObjectBase&& ) = default;

  virtual bool IsCoordInObject( const R1Tensor& coord ) const = 0;

  using CatalogInterface = cxx_utilities::CatalogInterface< SimpleGeometricObjectBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};


}
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
