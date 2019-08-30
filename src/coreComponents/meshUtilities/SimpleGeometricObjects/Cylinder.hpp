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
 * @file Cylinder.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Cylinder : public SimpleGeometricObjectBase
{
public:
  Cylinder( const std::string& name,
       ManagedGroup * const parent );

  virtual ~Cylinder() override;

  static string CatalogName() { return "Cylinder"; }

  bool IsCoordInObject( const R1Tensor& coord ) const override final;


private:

  R1Tensor m_point1;
  R1Tensor m_point2;
  realT m_radius = 0.0;

  struct viewKeyStruct
  {
    static constexpr auto point1String = "point1";
    static constexpr auto point2String = "point2";
    static constexpr auto radiusString = "radius";

  };


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
        */
