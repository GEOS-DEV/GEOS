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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef QUADRATURE_H
#define QUADRATURE_H

/**
 * @file Quadrature.h
 * @author white230
 */

#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"

#include <cassert>
#include "fileIO/xmlWrapper.hpp"

/*
 * Pure virtual base class representing a generic quadrature object.
 */
namespace geosx
{

class QuadratureBase
{
public:

  using CatalogInterface = cxx_utilities::CatalogInterface< QuadratureBase >;

  static CatalogInterface::CatalogType& GetCatalog();


  QuadratureBase() = default;
  virtual ~QuadratureBase();

  virtual int size() const = 0;
  virtual R1Tensor integration_point( const int index ) const = 0;
  virtual double integration_weight( const int index ) const = 0;

  virtual void ReadXML( xmlWrapper::xmlNode const & xmlNode ) = 0;

private:

};
}

#endif
