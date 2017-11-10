//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file UpdateFieldWithFunction.cpp
 * @author walsh24
 * @date June 19, 2011
 */

#include "SolverFactory.h"
#include "UpdateFieldWithFunction.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

UpdateFieldWithFunction::UpdateFieldWithFunction(  const std::string& name,
                                                   ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

UpdateFieldWithFunction::~UpdateFieldWithFunction()
{
  // TODO Auto-generated destructor stub
}

/*
 * <UpdateFieldWithFunction name="CaUpdate"       * name of the solver
 *            object="Face"                    * location of field to set (Default=Node)
 *            fieldtype="Scalar"                    * type of field to set    (Default=Scalar)
 *            fieldname="Ca"                       * name of field to update
 *            function="aFunction"                * name of the function 
 *            variables="Ca CO2"                  * function variables
 *            variableTypes="Scalar Scalar" />   * variableTypes (assumed scalar for all if omitted) 
 */
void UpdateFieldWithFunction::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );

  m_functionName = hdn->GetAttributeString("function");
  m_fieldName = hdn->GetAttributeString("fieldname");
  m_component = hdn->GetAttributeOrDefault("component",0);

  m_setNames = hdn->GetStringVector("setnames");
  if(m_setNames.empty()) m_setNames = hdn->GetStringVector("setname");
  

  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementNodeManagerStr() );
  {
    std::string oldStr = hdn->GetAttributeStringOrDefault("objecttype", "" ); // throw error if using old syntax
    if(!oldStr.empty()) {
      throw GPException("UpdateFieldWithFunction: Attempting to set objecttype - use 'object' instead.");
    }
  }
  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//=======
//  std::string objectTypeStr = hdn->GetAttributeString("objecttype");
//  if(objectTypeStr.empty())
//    throw GPException("Must define an object type for updating a field with a function");
//  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//>>>>>>> .r1185
  
  std::string fieldTypeStr = hdn->GetAttributeStringOrDefault("fieldtype",FieldInfo::RealStr);
  m_fieldType = fromString<FieldType>(fieldTypeStr);
  
  m_regionName = hdn->GetAttributeStringOrDefault("regionname","");
  
  std::string varsStr = hdn->GetAttributeString("variables");
  m_variables = Tokenize(varsStr," ");
  std::string varTypesStr = hdn->GetAttributeStringOrDefault("variabletypes","");
  
  if(varTypesStr == ""){
    m_variable_types = array<FieldType>(m_variables.size(),FieldInfo::realField);
  } else {
    array<string> vTypesVect = Tokenize(varTypesStr," ");
    m_variable_types.resize(vTypesVect.size());   
    
    if(m_variable_types.size() != m_variables.size()) 
      throw GPException("UpdateFieldWithFunction: Number of variable types not equal to number of variables.");
  
    for(size_t i =0; i < vTypesVect.size();++i)
      m_variable_types[i] = fromString<FieldType>(vTypesVect[i]);
    
  }
  
    
}

void UpdateFieldWithFunction::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  objectManager.AddKeylessDataField( m_fieldType,  m_fieldName, true, true );  
  // register variables
  for(size_t i =0; i < m_variables.size(); ++i){
    if(!streq(m_variables[i],"dt") && !streq(m_variables[i],"time") ){
  	  objectManager.AddKeylessDataField( m_variable_types[i],  m_variables[i], true, true );
    }
  }
}



/**
 * 
 * 
**/

double UpdateFieldWithFunction::TimeStep( const realT& time,
                                        const realT& dt,
                                        const int cycleNumber,
                                        PhysicalDomainT& domain,
                                        const array<string>& namesOfSolverRegions ,
                                        SpatialPartition& partition ,
                                        FractunatorBase* const fractunator )
{  
  realT dt_return = dt;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);
  if(m_setNames.empty()){
    objectManager.SetFieldEqualToFunction(m_fieldType,  m_fieldName, m_functionName, m_variables, m_variable_types, m_component,time,dt);
  } else {
    for( array<string>::size_type i =0; i < m_setNames.size(); ++i){
      lSet& set = objectManager.GetSet(m_setNames[i]);
      objectManager.SetFieldEqualToFunction(m_fieldType,  m_fieldName, m_functionName, m_variables, m_variable_types,set, m_component,time,dt);
    }
  }
  return dt_return;
}


REGISTER_SOLVER( UpdateFieldWithFunction )
