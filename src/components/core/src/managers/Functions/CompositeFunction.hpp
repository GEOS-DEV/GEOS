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
 * CompositeFunction.hpp
 *
 *  Created on: August 18, 2017
 *      Author: sherman
 */

#ifndef COMPOSITEFUNCTION_HPP_
#define COMPOSITEFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{


class CompositeFunction : public FunctionBase
{
public:
  CompositeFunction( const std::string& name,
                     dataRepository::ManagedGroup * const parent );

  virtual ~CompositeFunction() override;
  static string CatalogName() { return "CompositeFunction"; }
  virtual void FillDocumentationNode() override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void InitializeFunction() override;

  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         lSet const & sets,
                         real64_array & result ) const override final;

  virtual real64 Evaluate( real64 const * const input) const override final;

private:
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  localIndex m_numSubFunctions;
  static constexpr localIndex m_maxNumSubFunctions = 10;
  std::vector<FunctionBase*> m_subFunctions;

};


} /* namespace geosx */

#endif /* COMPOSITEFUNCTION_HPP_ */
