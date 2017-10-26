/*
 * KeyNames.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYNAMES_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYNAMES_HPP_

#include <string>

namespace geosx
{
namespace dataRepository
{
namespace keys
{

//static constexpr auto ReferencePosition = "ReferencePosition";
static constexpr auto referencePositionString = "ReferencePosition";

static constexpr auto TotalDisplacement = "TotalDisplacement";
static constexpr auto IncrementalDisplacement = "IncrementalDisplacement";
static constexpr auto Velocity = "Velocity";
static constexpr auto Acceleration = "Acceleration";
static constexpr auto Mass = "Mass";
static constexpr auto Force = "Force";
static constexpr auto Strain = "Strain";
static constexpr auto Name = "name";
static constexpr auto Size = "size";
static constexpr auto ProblemManager = "ProblemManager";
static constexpr auto ConstitutiveManager = "ConstitutiveManager";
static constexpr auto ConstitutiveBase = "ConstitutiveBase";
static constexpr auto solverNames = "solverNames";

static constexpr auto schema = "schema";
static constexpr auto beginTime = "time";
static constexpr auto endTime = "endTime";
static constexpr auto dt = "dt";
static constexpr auto cycle = "cycle";

static constexpr auto domain  = "domain";
static constexpr auto solvers = "solvers";
static constexpr auto simulationParameterMap = "simulationParameterMap";
static constexpr auto FE_Space    = "FE_Space";
//static constexpr auto FEM_Nodes    = "FEM_Nodes";
//static constexpr auto FEM_Edges    = "FEM_Edges";
//static constexpr auto FEM_Faces    = "FEM_Faces";
//static constexpr auto FEM_Elements = "FEM_Elements";
static constexpr auto cellManager = "cellMananger";
static constexpr auto functionManager = "FunctionManager";

}
}
}
#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYNAMES_HPP_ */
