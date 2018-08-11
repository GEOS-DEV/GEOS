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

/**
 * @file TableFunction.hpp
 */

#ifndef TABLEFUNCTION_HPP_
#define TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

namespace geosx
{

/**
 * @class TableFunction
 *
 * An interface for a dense table-based function
 */
class TableFunction : public FunctionBase
{
public:
  /// Main constructor
  TableFunction( const std::string& name,
                 dataRepository::ManagedGroup * const parent );

  /// Destructor
  virtual ~TableFunction() override;

  /// Catalog name interface
  static string CatalogName() { return "TableFunction"; }

  /// Documentation assignment
  virtual void FillDocumentationNode() override final;
  
  /// Initialize the function
  virtual void InitializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         set<localIndex> const & sets,
                         real64_array & result ) const override final;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   */
  virtual real64 Evaluate( real64 const * const input) const override final;

private:
  /// An array of table axes
  array<real64_array> m_coordinates;

  /// Table values (in fortran order)
  real64_array m_values;

  /// Table size indicators
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
