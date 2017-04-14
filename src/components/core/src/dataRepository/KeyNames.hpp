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

std::string const ReferencePosition = "ReferencePosition";
std::string const TotalDisplacement = "TotalDisplacement";
std::string const IncrementalDisplacement = "IncrementalDisplacement";
std::string const Velocity = "Velocity";
std::string const Acceleration = "Acceleration";
std::string const Mass = "Mass";
std::string const Force = "Force";
std::string const Strain = "Strain";
std::string const Name = "name";
std::string const Size = "size";
std::string const ProblemManager = "ProblemManager";
std::string const ConstitutiveManager = "ConstitutiveManager";
std::string const ConstitutiveBase = "ConstitutiveBase";
std::string const solverNames = "solverNames";

std::string const schema = "schema";
std::string const beginTime = "beginTime";
std::string const endTime = "endTime";
std::string const dt = "dt";

std::string const domain  = "domain";
std::string const solvers = "solvers";
std::string const simulationParameterMap = "simulationParameterMap";
std::string const FE_Space    = "FE_Space";
std::string const FEM_Nodes    = "FEM_Nodes";
std::string const FEM_Edges    = "FEM_Edges";
std::string const FEM_Faces    = "FEM_Faces";
std::string const FEM_Elements = "FEM_Elements";

}
}
}
#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYNAMES_HPP_ */
