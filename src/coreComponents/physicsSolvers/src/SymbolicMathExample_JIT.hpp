/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2019 Total, S.A
* Copyright (c) 2019-     GEOSX Contributors
* All right reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
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
class Group;
}
class DomainPartition;

class SymbolicMathExample_JIT : public SolverBase
{
public:
SymbolicMathExample_JIT( const std::string& name,
Group * const parent );


virtual ~SymbolicMathExample_JIT() override;

static string CatalogName() { return "SymbolicMathExample_JIT"; }

virtual void FillDocumentationNode() override;

virtual void BuildDataStructure( dataRepository::Group * const domain ) override;

virtual void Initialize( Group * const problemManager ) override final;

virtual real64 SolverStep( real64 const& time_n,
real64 const& dt,
integer const cycleNumber,
DomainPartition * domain ) override;

};

} /* namespace geosx */

#endif /* SYMBOLICMATHEXAMPLE_JIT_HPP_ */
