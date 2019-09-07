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
* SymbolicMathExample_JIT.cpp
*
*  Created on: June 9, 2017
*      Author: sherman
*/

#include "SymbolicMathExample_JIT.hpp"

#include <vector>
#include <math.h>

#include "common/DataTypes.hpp"
#include "mesh/NodeManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "managers/Functions/FunctionBase.hpp"
#include <numeric>

#include "dataRepository/Group.hpp"
namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const FunctionName = "functionName";
std::string const TargetObject = "targetObject";
std::string const TargetName = "targetName";
}
}

using namespace dataRepository;


SymbolicMathExample_JIT::SymbolicMathExample_JIT( const std::string& name,
Group * const parent ):
SolverBase( name, parent )
{}



SymbolicMathExample_JIT::~SymbolicMathExample_JIT()
{
// TODO Auto-generated destructor stub
}


void SymbolicMathExample_JIT::FillDocumentationNode()
{
cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
SolverBase::FillDocumentationNode();

docNode->setName(this->CatalogName());
docNode->setSchemaType("Node");
docNode->setShortDescription("An example solid mechanics solver");

docNode->AllocateChildNode( keys::FunctionName,
keys::FunctionName,
-1,
"string",
"string",
"Name of the symbolic math function",
"Name of the symbolic math function",
"default",
"",
1,
1,
0 );

docNode->AllocateChildNode( keys::TargetObject,
keys::TargetObject,
-1,
"string",
"string",
"Target object that holds field that we are trying to set",
"Target object that holds field that we are trying to set",
"default",
"",
1,
1,
0 );

docNode->AllocateChildNode( keys::TargetName,
keys::TargetName,
-1,
"string",
"string",
"Name of the field to hold results of symbolic expression",
"Name of the field to hold results of symbolic expression",
"default",
"",
1,
1,
0 );

}


void SymbolicMathExample_JIT::BuildDataStructure( Group * const domain )
{
SolverBase::BuildDataStructure( domain );
}


void SymbolicMathExample_JIT::Initialize( Group * const problemManager )
{
// Check to see if targets are valid
Group * domain = problemManager->GetGroup(keys::domain);
std::string targetObjectStr = getData<std::string>(keys::TargetObject);
std::string targetName = getData<std::string>(keys::TargetName);

if (domain->hasGroup(targetObjectStr))
{
dataRepository::Group * targetObject = domain->GetGroup<Group>(targetObjectStr);
if (targetObject->hasView(targetName))
{
GEOS_LOG("Symbolic expression is setting " << targetObjectStr << "/" << targetName);
}
else
{
GEOS_LOG("Target view wrapper is not present..  Initializing " << targetObjectStr << "/" << targetName);
targetObject->registerWrapper<real64_array>(targetName);
}
}
else
{
GEOS_ERROR("No matching group for symbolic solver: " + targetObjectStr);
}
}


real64 SymbolicMathExample_JIT::SolverStep( real64 const& time_n,
real64 const& dt,
const int cycleNumber,
DomainPartition * domain )
{
// Get target objects
dataRepository::Group * targetObject = domain->GetGroup<Group>(getData<std::string>(keys::TargetObject));
real64_array & targetField = targetObject->getReference<real64_array>(getData<std::string>(keys::TargetName));
view_rtype_const<r1_array> X = targetObject->getData<r1_array>(keys::referencePositionString);

// Evaluate symbolic math
NewFunctionManager * functionManager = NewFunctionManager::Instance();
FunctionBase * function = functionManager->GetGroup<FunctionBase>(getData<std::string>(keys::FunctionName));

// Test automatic evaluation for the entire mesh
set<localIndex> set;
for(localIndex ii = 0 ; ii < targetObject->size() ; ++ii)
{
set.insert(ii);
}
function->Evaluate(targetObject, time_n, set, targetField);

// Test manual serial evaluation
double data[4];
data[3] = time_n;
real64 error = 0.0;
for (localIndex ii = 0 ; ii < targetObject->size() ; ++ii)
{
data[0] = X[ii][0];
data[1] = X[ii][1];
data[2] = X[ii][2];
real64 res = function->Evaluate(data);
error += fabs(res - targetField[ii]);
}

// Throw an error if the function evaluation methods return different results
GEOS_ERROR_IF(error > 1e-6, "Serial and vector function evaluation do not match!");
return dt;
}



REGISTER_CATALOG_ENTRY( SolverBase, SymbolicMathExample_JIT, std::string const &, Group * const )
} /* namespace ANST */
