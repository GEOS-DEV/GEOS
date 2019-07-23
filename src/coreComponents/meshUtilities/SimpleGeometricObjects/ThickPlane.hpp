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
 * @file ThickPlane.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class ThickPlane : public SimpleGeometricObjectBase
{
public:
  ThickPlane( const std::string& name,
              ManagedGroup * const parent );

  virtual ~ThickPlane() override;

  static string CatalogName() { return "ThickPlane"; }

  bool IsCoordInObject( const R1Tensor& coord ) const override final;

protected:
  virtual void PostProcessInput() override final;

private:

  R1Tensor m_origin;
  R1Tensor m_normal;
  real64   m_thickness;

  struct viewKeyStruct
  {
    static constexpr auto originString = "origin";
    static constexpr auto normalString = "normal";
    static constexpr auto thicknessString = "thickness";
  };


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
        */

