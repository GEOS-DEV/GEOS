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
 * TableFunction.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef TABLEFUNCTION_HPP_
#define TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

namespace geosx
{


class TableFunction : public FunctionBase
{
public:
  TableFunction( const std::string& name,
                 dataRepository::ManagedGroup * const parent );

  virtual ~TableFunction() override;
  static string CatalogName() { return "TableFunction"; }
  virtual void FillDocumentationNode() override final;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override final;

  virtual void InitializeFunction() override;

  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         lSet const & sets,
                         real64_array & result ) const override final;

  virtual real64 Evaluate( real64 const * const input) const override final;

private:
  array<real64_array> m_coordinates;
  real64_array m_values;
  static localIndex constexpr m_maxDimensions = 4;
  localIndex m_dimensions;
  localIndex_array m_size;
  localIndex_array m_indexIncrement;

  // m_corners should be of size m_maxDimensions x (2^m_maxDimensions)
  localIndex m_corners[m_maxDimensions][16];
  localIndex m_numCorners;
};


} /* namespace geosx */

#endif /* TABLEFUNCTION_HPP_ */
