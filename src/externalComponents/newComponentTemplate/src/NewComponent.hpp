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
* NewComponent.hpp
*
*  Created on: Jun 8, 2016
*      Author: settgast
*/

#ifndef COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#define COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{
namespace dataRepository
{
class Group;
}
class DomainPartition;

class NewComponent : public SolverBase
{
public:
NewComponent( std::string const & name,
Group * const parent);
virtual ~NewComponent() override;

static std::string CatalogName() { return "NewComponent"; }


virtual real64 SolverStep( real64 const& time_n,
real64 const& dt,
integer const cycleNumber,
DomainPartition * domain ) override;

private:
NewComponent() = delete;
NewComponent(const NewComponent&) = delete;
NewComponent(const NewComponent&&) = delete;
NewComponent& operator=(const NewComponent&) = delete;
NewComponent& operator=(const NewComponent&&) = delete;
};

} /* namespace geosx */

#endif /* COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_ */
