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
 * FiniteElementSpace.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#include "../dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"

namespace geosx
{

class NodeManager;
class CellBlockManager;

namespace dataRepository
{
namespace keys
{
string const finiteElementSpace = "finiteElementSpace";
string const basis = "basis";
string const quadrature = "quadrature";
string const dNdX = "dNdX";
string const detJ = "detJ";
}
}

class BasisBase;
class QuadratureBase;
class FiniteElementBase;

class FiniteElementSpace : public dataRepository::ManagedGroup
{
public:

  FiniteElementSpace() = delete;

  explicit FiniteElementSpace( std::string const & name, ManagedGroup * const parent );

  ~FiniteElementSpace() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static string CatalogName() { return dataRepository::keys::finiteElementSpace; }

  ///@}

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const parent ) override;

  void FillDocumentationNode() override;

  void ApplySpaceToTargetCells( dataRepository::ManagedGroup * const group ) const;

  void ReadXML_PostProcess() override final;

  virtual void InitializePreSubGroups( ManagedGroup * const group ) override;

  void CalculateShapeFunctionGradients( r1_array const & X,
                                        dataRepository::ManagedGroup * const cellBlock ) const;

public:

  BasisBase const *    m_basis    = nullptr;
  QuadratureBase const * m_quadrature = nullptr;
  FiniteElementBase * m_finiteElement = nullptr;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_ */
