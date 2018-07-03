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

/*
 * Box.hpp
 *
 *  Created on: Aug 4, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Box : public SimpleGeometricObjectBase
{
public:
  Box( const std::string& name,
       ManagedGroup * const parent );

  virtual ~Box() override;

  static string CatalogName() { return "Box"; }

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override final;

  bool IsCoordInObject( const R1Tensor& coord ) const override final;
private:
  R1Tensor m_min;
  R1Tensor m_max;
  realT m_strikeAngle=0.0;
  R1Tensor m_boxCenter={0.0,0.0,0.0};
  realT m_cosStrike=0.0, m_sinStrike=0.0;

  struct viewKeyStruct
  {
    dataRepository::ViewKey xmin = { "xMin" };
    dataRepository::ViewKey xmax = { "xMax" };
  } viewKeys;


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
        */
