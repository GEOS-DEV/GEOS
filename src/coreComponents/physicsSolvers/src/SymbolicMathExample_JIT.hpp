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
 * SymbolicMathExample_JIT.hpp
 *
 *  Created on: Jun 13, 2017
 *      Author: sherman
 */

#ifndef SYMBOLICMATHEXAMPLE_JIT_HPP_
#define SYMBOLICMATHEXAMPLE_JIT_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class DomainPartition;

class SymbolicMathExample_JIT : public SolverBase
{
public:
  SymbolicMathExample_JIT( const std::string& name,
                           ManagedGroup * const parent );


  virtual ~SymbolicMathExample_JIT() override;

  static string CatalogName() { return "SymbolicMathExample_JIT"; }

  virtual void FillDocumentationNode() override;

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void Initialize( ManagedGroup * const problemManager ) override final;

  virtual real64 SolverStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain ) override;

};

} /* namespace geosx */

#endif /* SYMBOLICMATHEXAMPLE_JIT_HPP_ */
