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
 * @file InitialCondition.cpp
 * @author walsh24
 * @date Feb 28, 2011
 */

#include <algorithm>
#include "InitialConditions.h"
#include "ObjectManagers/PhysicalDomainT.h"
//#include "ObjectManagers/ProblemManagerT.h"
#include "Utilities/Utilities.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"
#include "ObjectManagers/TableManager.h"

////////////////////////////
//
/// Initial condition base class
/**
 * @author walsh24
 * @brief Initial condition base class - should not be instantiated directly.
 * 
 **/
InitialConditionBase::InitialConditionBase(  TICPP::HierarchicalDataNode*  hdn,
                                            const ProblemManagerT* const pm  )
{
  ReadXML(hdn); 
}

/// explicit constructor (in most cases the hdn constructor should be invoked). 
InitialConditionBase::InitialConditionBase( const PhysicalDomainT::ObjectDataStructureKeys& objectType, const std::string& regionName,
                      const std::string& fieldName,const  std::vector<std::string>& setNames,  const FieldType& fieldType):
  objectType_(objectType), 
  regionName_(regionName), 
  fieldName_(fieldName), 
  setNames_(setNames), 
  fieldType_(fieldType)
{
                        	
}

InitialConditionBase::~InitialConditionBase()
{
  // TODO Auto-generated destructor stub
}


/**
 * @author walsh24
 * @brief Records object and field data.
 * 
 * <GENERIC_INITIAL_CONDITION  ... 
 *      object="Element"               * Object type, Node,Face,Element etc (Default = Node)
 *      regionname="Region1"           * Apply to all elements in Region1 (for Element fields only - no effect otherwise)
 *      setname="someFieldSubset"      * Subset to apply initial condition to - leave out or empty to apply to whole field 
 *      fieldname="Concentration"      * Name of field
 *      fieldtype="Scalar"/>           * Type of field (Default = Scalar)
 * 
 * 
 **/
void InitialConditionBase::ReadXML( TICPP::HierarchicalDataNode*  hdn)
{
  {
    std::string objectTypeStr = hdn->GetAttributeString("object");
    if(objectTypeStr.empty())
      throw GPException("Cannot specify an initial condition without an object type");
    objectType_ = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
  }
  
  std::string fieldTypeStr = hdn->GetAttributeStringOrDefault("fieldtype", FieldInfo::RealStr);
  fieldType_ = fromString<FieldType>(fieldTypeStr);
  
  fieldName_ = hdn->GetAttributeString("fieldname");
  regionName_ = hdn->GetAttributeStringOrDefault("toregions","");

  std::size_t foundSpace = regionName_.find(" ");
  std::size_t foundComma = regionName_.find(",");
  if (foundComma != std::string::npos || foundSpace != std::string::npos )
    throw GPException("Error.  You are trying to apply an initial condition to multiple regions at once.  This is illegal.");

  if (objectType_ == PhysicalDomainT::FiniteElementElementRegion && regionName_.empty())
  {
    throw GPException("Error.  You are applying an initial condition to Element but didn't provide a region name via toregions=");
  }

  {
    array<string> tempSetName;
    tempSetName = hdn->GetStringVector("setname");
    if (!tempSetName.empty())
      throw GPException("ERROR!!! 'setname' is no longer supported for initial conditions.  Use 'setnames' instead.");
  }
  setNames_ = hdn->GetStringVector("setnames");
  m_additive = hdn->GetAttributeOrDefault<bool>("additive", false);
  // This flag allows superpose a field onto another one.
  // It has been implemented in constant and tabled initial conditions.  I am skipping the function-based one for now.
  if (m_additive)
    m_scale = hdn->GetAttributeOrDefault<realT>("scale", 1.0);

}


bool InitialConditionBase::IsThisInitialStress()
{
  return (fieldName_.substr(0,6) == "sigma_");
}

////////////////////////////
//
/// ReadInitialConditionFromFile

/**
 * @author walsh24
 * @brief Read initial values from an ascii file.
 * 
 */
ReadInitialConditionFromFile::ReadInitialConditionFromFile(  TICPP::HierarchicalDataNode*  hdn, const ProblemManagerT* const pm):
  InitialConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

/**
 *     <ReadInitialConditionFromFile filename="TTF.conc"  
 *                                 object="node"
 *                                 setname="someFieldSubset"      * Leave out or empty to apply to whole field 
 *                                 fieldname="Concentration" 
 *                                 fieldtype="Scalar"/>
 **/
void ReadInitialConditionFromFile::ReadXML(TICPP::HierarchicalDataNode*  hdn)
{  
  InitialConditionBase::ReadXML(hdn);
  filename_ = hdn->GetAttributeString("filename");
  isIndexedFile_ = hdn->GetAttributeOrDefault<bool>("isIndexedFile",false);
}

void ReadInitialConditionFromFile::RegisterFields( PhysicalDomainT& domain ){
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
  objectManager.AddKeylessDataField(fieldType_,  fieldName_, true, true ); 
}

void ReadInitialConditionFromFile::Apply( PhysicalDomainT& domain )
{
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);

  if( setNames_.empty() ){
	if(isIndexedFile_){
	  objectManager.ReadIndexedAsciiFieldData(fieldType_, fieldName_, filename_);
	} else {
      objectManager.ReadAsciiFieldData(fieldType_, fieldName_, filename_);
	}
  } else {
	if(isIndexedFile_){
	      throw GPException("Error ReadInitialConditionFromFile: Specifying an subset with an indexed file is currently not supported.");
	} else {
  	  for( array<string>::size_type i =0; i < setNames_.size(); ++i){
        lSet& set = objectManager.GetSet(setNames_[i]);
        objectManager.ReadAsciiFieldData(fieldType_, fieldName_, filename_, set);
  	  }
	}
  }
}

REGISTER_InitialCondition( ReadInitialConditionFromFile )

////////////////////////////
//
// ConstantInitialCondition

/**
 * @author walsh24
 * @brief Sets the intial condition for a field or subset of a field to a constant value. 
 * 
 */
ConstantInitialCondition::ConstantInitialCondition(  TICPP::HierarchicalDataNode*  hdn, const ProblemManagerT* const pm):
  InitialConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

/**    <ConstantInitialCondition fieldname="ConcentrationFlux" 
 *                              fieldtype="Vector"  
 *                              object="Face"
 *                              setname="TopBoundaryFaces"   * Leave out or empty to apply to whole field
 *                              value="0 0 0"/>
 */
void ConstantInitialCondition::ReadXML( TICPP::HierarchicalDataNode* hdn)
{
  valueStr_ = hdn->GetAttributeString("value");
  if (valueStr_.empty())
    throw GPException("Error. Value of at least one constant initial condition was not properly set.");
}

void ConstantInitialCondition::RegisterFields( PhysicalDomainT& domain ){
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
  objectManager.AddKeylessDataField( fieldType_,  fieldName_, true, true ); 
}

void ConstantInitialCondition::Apply( PhysicalDomainT& domain )
{
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
  if( setNames_.empty() ){
    objectManager.SetFieldToConstantFromString( fieldType_,  fieldName_, valueStr_, m_additive);
  } else {
  	for( array<string>::size_type i =0; i < setNames_.size(); ++i){
      lSet& set = objectManager.GetSet(setNames_[i]);
      objectManager.SetFieldToConstantFromString( fieldType_,  fieldName_, valueStr_, set, m_additive);
  	}
  }
}

REGISTER_InitialCondition( ConstantInitialCondition )

////////////////////////////
//
//  InitialConditionTable
//

/**
 * @author johnson346
 * @brief Sets the initial condition for a field or subset of a field using a named table.
 *
 */
InitialConditionTable::InitialConditionTable(  TICPP::HierarchicalDataNode*  hdn,
                                              const ProblemManagerT* const pm): InitialConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

/**
 * <InitialConditionTable fieldname="Velocity"
 *                        fieldtype="Vector"
 *                        component="0"             * Component of tensor fields - default is 0
 *                        object="Node"
 *                        table="my_table"/>        * name of TableManager table
 *
 */
void InitialConditionTable::ReadXML( TICPP::HierarchicalDataNode*  hdn)
{

  m_tableName = hdn->GetAttributeString("table");
  if (m_tableName.empty())
    throw GPException("Must declare a table in InitialConditionTable");

  m_component = hdn->GetAttributeOrDefault("component", -1);
}


void InitialConditionTable::RegisterFields(PhysicalDomainT& domain)
{
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_, regionName_);
  objectManager.AddKeylessDataField(fieldType_, fieldName_, true, true);
}

void InitialConditionTable::Apply(PhysicalDomainT& domain)
{
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_, regionName_);
  objectManager.AddKeylessDataField(fieldType_, fieldName_, true, true);

  //get temporary pointers to fields to set
  array<R1Tensor>* r1ptr = fieldType_ != FieldInfo::realField ? objectManager.GetFieldDataPointer<R1Tensor>(fieldName_) : 0;
  array<real64>* r0ptr =           fieldType_ != FieldInfo::realField ? 0 : objectManager.GetFieldDataPointer<realT>(fieldName_);

  //get whether this is a finite element region
  const bool isFE = objectManager.GetObjectType() == ObjectDataStructureBaseT::ElementRegion;
  const array<R1Tensor>* pos = !isFE ? objectManager.GetFieldDataPointer<FieldInfo::referencePosition>() : 0;
  {
    if ((!isFE) && (!pos))
      throw GPException(
        "InitialConditionTable::Apply - when trying to apply you must either apply to an ElementRegionT or to an ObjectDataStructureBaseT where FieldInfo::referencePosition has been added as a data field");
  }

  //get the 3D table
  const Table3D* t3dp = stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_tableName);

  //if this is a 3D scalar table and the field it references is not a scalar field,
  //make sure to default m_component to 0
  if(t3dp && (fieldType_ != FieldInfo::realField) && m_component < 0)
    m_component = 0;

  //get the 3D vector field
  const VectorField3D* vf3dp = t3dp ? 0 : stlMapLookupPointer(TableManager::Instance().VectorFields<3>(), m_tableName);

  //determine which to use: scalar or vectorfield (scalar takes precedent)
  //then if vector field, make sure to check whether the field is a vector
  //and if so make sure that m_component is not specified
  if(vf3dp && fieldType_ != FieldInfo::realField && m_component >= 0)
    m_component = -1;

  //make sure something was found; otherwise, throw an exception
  if(!t3dp && !vf3dp)
      throw GPException("InitialConditionTable::Apply : unrecognized table or vector field");

  //fill the sets
  array<lSet> sets(setNames_.size());
  for (array<string>::const_iterator ss = setNames_.begin(); ss != setNames_.end(); ++ss)
    sets.push_back(objectManager.GetSet(*ss));

  //now that we know what's going on, let's call the right function template ... we need to know at compile
  //time, so you need to explicitly declare all cases
  if (!m_additive)
  {
    if(m_component >= 0 && !SetField(t3dp,  m_component, pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr))
    {
      //if you have a R1 tensor field in which you are assigning by component, the only rational way to do this is
      //with a scalar table ... otherwise, throw an exception
      throw GPException("InitialConditionTable::Apply - unhandled case for component-wise cases");
    }
    else
    {
      //with no component, you may set a scalar with a scalar
      bool ok = SetField(t3dp,  pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r0ptr);
      // ... you may set a vector with a scalar (BUT BE CAREFUL!!)
      if(!ok)
        ok =    SetField(t3dp,  pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr);
      // ... you may set a vector with a vector
      if(!ok)
        ok =    SetField(vf3dp, pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr);
      // ... or this is unhandled!!
      if(!ok)
        throw GPException("InitialConditionTable::Apply - unhandled case");
    }
  }
  else
  {
    if(m_component >= 0 && !AddToField(t3dp,  m_component, pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr, m_scale))
    {
      //if you have a R1 tensor field in which you are assigning by component, the only rational way to do this is
      //with a scalar table ... otherwise, throw an exception
      throw GPException("InitialConditionTable::Apply - unhandled case for component-wise cases");
    }
    else
    {
      //with no component, you may set a scalar with a scalar
      bool ok = AddToField(t3dp,  pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r0ptr, m_scale);
      // ... you may set a vector with a scalar (BUT BE CAREFUL!!)
      if(!ok)
        ok =    AddToField(t3dp,  pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr, m_scale);
      // ... you may set a vector with a vector
      if(!ok)
        ok =    AddToField(vf3dp, pos, domain.m_feNodeManager, (ElementRegionT&) objectManager, sets, r1ptr, m_scale);
      // ... or this is unhandled!!
      if(!ok)
        throw GPException("InitialConditionTable::Apply - unhandled case");
    }
  }

}

REGISTER_InitialCondition( InitialConditionTable )

////////////////////////////
//
//  InitialConditionFunction
//

/** 
 * @author walsh24
 * @brief Sets the initial condition for a field or subset of a field using a named function.
 * 
 */
InitialConditionFunction::InitialConditionFunction(  TICPP::HierarchicalDataNode*  hdn, const ProblemManagerT* const pm):
  InitialConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

/**
 * <InitialConditionFunction fieldname="Concentration"        
 *                           fieldtype="Scalar"          
 *                           component="0"                    * Component of tensor fields - default is 0   
 *                           object="Node"                          
 *                           function="gaussianDistribution"  * function name
 *                           variables="ReferencePosition"    * space separated list of variables in order used in function
 *                           variableTypes="Vector"/>         * space separated list of variable types
 * 
 */
void InitialConditionFunction::ReadXML( TICPP::HierarchicalDataNode*  hdn){
	
  functionName_ = hdn->GetAttributeString("function");
  component_ = hdn->GetAttributeOrDefault("component",0);
  
  std::string varsStr = hdn->GetAttributeStringOrDefault("variables","");
  if(varsStr.empty()){
  	variableNames_.resize(0);
  } else {
    variableNames_ = Tokenize(varsStr," ");
  }
  
  std::string varTypesStr = hdn->GetAttributeStringOrDefault("variableTypes","");
  if( varTypesStr.empty() ){
    variableTypes_ = array<FieldType>(variableNames_.size(),FieldInfo::realField);
  } else {
    array<string> vTypesVect = Tokenize(varTypesStr," ");
    variableTypes_.resize(vTypesVect.size()); 
      
    if(variableTypes_.size() != variableNames_.size()) 
      throw GPException("Error InitialConditionFunction: Number of variable types not equal to number of variables.");
  
    for( array<string>::size_type i=0 ; i < vTypesVect.size() ; ++i )
      variableTypes_[i] = fromString<FieldType>(vTypesVect[i]);
    
  }
}


void InitialConditionFunction::RegisterFields( PhysicalDomainT& domain ){
  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
  objectManager.AddKeylessDataField( fieldType_,  fieldName_, true, true ); 
}

void InitialConditionFunction::Apply( PhysicalDomainT& domain ){

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
  
  objectManager.AddKeylessDataField( fieldType_,  fieldName_, true, true );  
  if( setNames_.empty() ){
    objectManager.SetFieldEqualToFunction(fieldType_,  fieldName_, functionName_, variableNames_, variableTypes_,component_);
  } else {
  	for( array<string>::size_type i =0; i < setNames_.size(); ++i){
      lSet& set = objectManager.GetSet(setNames_[i]);
      objectManager.SetFieldEqualToFunction(fieldType_,  fieldName_, functionName_, variableNames_, variableTypes_,set,component_);
  	}
  }
 
}

REGISTER_InitialCondition( InitialConditionFunction )


//////////////////////////////
////
////  InitialConditionFrom3DTable
////
//
///**
// * @author fu4
// * @brief Sets the initial condition for a field or subset of a field from a 3D table. Currently only implemented for scalar fields.
// *
// */
//InitialConditionFrom3DTable::InitialConditionFrom3DTable( const TICPP::HierarchicalDataNode* const hdn, const ProblemManagerT* const pm):
//  InitialConditionBase(hdn,pm)
//{
//  ReadXML(hdn);
//}
//
///**
// * <InitialConditionFunction fieldname="sigma_x"
// *                           fieldtype="Scalar"
// *                           object="Element"
// *                           tablename="stress" />s
// *
// */
//void InitialConditionFrom3DTable::ReadXML(const TICPP::HierarchicalDataNode* const hdn){
//
//  m_tableName= hdn->GetAttributeString("tablename");
//}
//
//
//void InitialConditionFrom3DTable::RegisterFields( PhysicalDomainT& domain ){
//  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
//  objectManager.AddKeylessDataField( fieldType_,  fieldName_, true, true );
//}
//
//void InitialConditionFrom3DTable::Apply( PhysicalDomainT& domain ){
//
//  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(objectType_,regionName_);
//  const TableManager& tableManager = TableManager::Instance();
//  FieldTypeMultiPtr theFieldPtr;
//
//  theFieldPtr.SetFieldPtr(objectManager, fieldType_, fieldName_);
//
//  const std::map<std::string,Table3D >::const_iterator it3 = tableManager.Tables<3>().find(m_tableName);
//  if(it3 == tableManager.Tables<3>().end())
//    throw GPException("ConstitutivePropertiesTable::Apply : Cannot find requested table in the table manager: " + m_tableName);
//  const Table3D& t3d = it3->second;
//
//  for (localIndex i=0; i<objectManager.DataLengths(); ++i)
//  {
//    R1Tensor cpos(0);
//    for (localIndex j=0; j<((ElementRegionT&)objectManager).m_numNodesPerElem; j++)
//    {
//      cpos += (*domain.m_feNodeManager.m_refposition)[((ElementRegionT&)objectManager).m_toNodesRelation[i][j]];
//    }
//    cpos /= ((ElementRegionT&)objectManager).m_numNodesPerElem;
//    theFieldPtr.SetValue(i, 0, t3d.Lookup(cpos, true));
//  }
//
//}
//
//REGISTER_InitialCondition( InitialConditionFrom3DTable )
//

////////////////////////////
//
//  CalculateFaceCenters
//

/**
 * @author walsh24
 * @brief Assigns a face center field to each face. 
 * 
 */
CalculateFaceCenters::CalculateFaceCenters( TICPP::HierarchicalDataNode* hdn,
                                            const ProblemManagerT* const pm  ):
  InitialConditionBase(PhysicalDomainT::FiniteElementFaceManager,
                       "",
                       PS_STR::FaceCenterStr,
                       std::vector<std::string>(),
                       FieldInfo::R1TensorField)
{
  setNames_ = hdn->GetStringVector("setnames");
}

void CalculateFaceCenters::RegisterFields( PhysicalDomainT& domain ){
  domain.m_feFaceManager.AddKeylessDataField( FieldInfo::R1TensorField,  PS_STR::FaceCenterStr, true, true ); 
}

/*
 * <InitialConditions>
 *   <CalculateFaceCenters/>
 * </InitialConditions>
 */
 
void CalculateFaceCenters::Apply( PhysicalDomainT& domain )
{
	
  array<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( PS_STR::FaceCenterStr );
  
  if( setNames_.empty() ){
    
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
        domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, kf, faceCenter[kf] );
    
  } else {
  	
  	for( array<string>::size_type i =0; i < setNames_.size(); ++i){
      lSet& subset = domain.m_feFaceManager.GetSet(setNames_[i]);
  	  for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si ){
        localIndex kf = *si;
        domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, kf, faceCenter[kf] );
	  }
  	}
  }    

}

REGISTER_InitialCondition( CalculateFaceCenters )

////////////////////////////
//
//  CalculateElementCenters
//

/**
 * @author annavarapusr1
 * @brief Assigns an element center field to each element.
 *
 */
CalculateElementCenters::CalculateElementCenters( TICPP::HierarchicalDataNode* hdn,
                                                  const ProblemManagerT* const pm  ):
  InitialConditionBase(PhysicalDomainT::FiniteElementElementManager,
                       hdn->GetAttributeStringOrDefault("toregions",""),
                       PS_STR::ElementCenterStr,
                       std::vector<std::string>(),
                       FieldInfo::R1TensorField)
{
  setNames_ = hdn->GetStringVector("setnames");
}

void CalculateElementCenters::RegisterFields( PhysicalDomainT& domain )
{
  for (std::map<std::string, ElementRegionT>::iterator elementRegionIter=domain.m_feElementManager.m_ElementRegions.begin() ;
      elementRegionIter!= domain.m_feElementManager.m_ElementRegions.end(); ++elementRegionIter)
  {
    //const std::string& elementRegionName = elementRegionIter->first;
    ElementRegionT& elementRegion = elementRegionIter->second;
    elementRegion.AddKeylessDataField( FieldInfo::R1TensorField,  PS_STR::ElementCenterStr, true, true );
  }
}

/*
 * <InitialConditions>
 *   <CalculateElementCenters/>
 * </InitialConditions>
 */

void CalculateElementCenters::Apply( PhysicalDomainT& domain)
{
  for (std::map<std::string, ElementRegionT>::iterator elementRegionIter=domain.m_feElementManager.m_ElementRegions.begin() ;
      elementRegionIter!= domain.m_feElementManager.m_ElementRegions.end(); ++elementRegionIter)
  {
    //const std::string& elementRegionName = elementRegionIter->first;
    ElementRegionT& elementRegion = elementRegionIter->second;

    array<R1Tensor>& elementCenter = elementRegion.GetFieldData<R1Tensor>( PS_STR::ElementCenterStr );

    if( setNames_.empty() )
    {
      for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
      {
        elementCenter[k] = elementRegion.GetElementCenter( k, domain.m_feNodeManager);
      };
    }
    else
    {
      for( array<string>::size_type i =0; i < setNames_.size(); ++i)
      {
        lSet& subset = elementRegion.GetSet(setNames_[i]);
        for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si )
        {
          localIndex k = *si;
          elementCenter[k] = elementRegion.GetElementCenter( k, domain.m_feNodeManager);
        }
      }
    }
  }
}

REGISTER_InitialCondition( CalculateElementCenters )

////////////////////////////
//
//  CalculateFaceNormals
//

/**
 * @author walsh24
 * @brief Assigns a face normal field to each face. 
 * 
 */
CalculateFaceNormals::CalculateFaceNormals( TICPP::HierarchicalDataNode* hdn,
                                            const ProblemManagerT* const pm  ):
  InitialConditionBase(PhysicalDomainT::FiniteElementFaceManager, "",
                       PS_STR::FaceNormalStr, std::vector<std::string>(),  FieldInfo::R1TensorField)
{
  setNames_ = hdn->GetStringVector("setnames");
}

void CalculateFaceNormals::RegisterFields( PhysicalDomainT& domain ){
  domain.m_feFaceManager.AddKeylessDataField( FieldInfo::R1TensorField,  PS_STR::FaceNormalStr, true, true ); 
}

/*
 * <InitialConditions>
 *   <CalculateFaceNormals/>
 * </InitialConditions>
 */
 
void CalculateFaceNormals::Apply( PhysicalDomainT& domain )
{
	
  array<R1Tensor>& faceNormal = domain.m_feFaceManager.GetFieldData<R1Tensor>( PS_STR::FaceNormalStr );
  
  if( setNames_.empty() ){
    
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
        domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, kf, faceNormal[kf] );
    
  } else {
  	
  	for( array<string>::size_type i =0; i < setNames_.size(); ++i){
      lSet& subset = domain.m_feFaceManager.GetSet(setNames_[i]);
  	  for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si ){
        localIndex kf = *si;
          domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, kf, faceNormal[kf] );
	  }
  	}
  }    

}

REGISTER_InitialCondition( CalculateFaceNormals )

////////////////////////////
//
//  CalculateAperture
//

/**
 * @author walsh24
 * @brief Assigns an aperture to each face calculated from the common plane contact 
 * 
 */
CalculateAperture::CalculateAperture( TICPP::HierarchicalDataNode* hdn,
                                      const ProblemManagerT* const pm  ):
  InitialConditionBase(PhysicalDomainT::FiniteElementFaceManager, "",
                       "aperture", std::vector<std::string>(),  FieldInfo::realField)
{
  setNames_ = hdn->GetStringVector("setnames");
}

void CalculateAperture::RegisterFields( PhysicalDomainT& domain ){
	
  // fixme these fields are or may become redundant in the future as the common-plane face manager is cleaned up
  domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  "Aperture", true, true ); 
  domain.m_externalFaces.AddKeylessDataField( FieldInfo::realField,  "aperture", true, true ); 
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::contactForce>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::velocity>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::acceleration>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::force>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::mass>();
  domain.m_externalFaces.AddKeylessDataField( FieldInfo::R1TensorField,  "faceVelocity", true, true );
  domain.m_externalFaces.AddKeylessDataField( FieldInfo::realField,   "frictionSlope", true, true );
}

/*
 * <InitialConditions>
 *   <CalculateAperture/>
 * </InitialConditions>
 */
void CalculateAperture::Apply( PhysicalDomainT& domain )
{
	  
//  array<R1Tensor>& faceCenter = domain.m_faceManager.GetFieldData<R1Tensor>( PS_STR::FaceCenterStr );
  //std::cout << "Calculating Aperture" << std::endl;
  const array<integer>& isExternal = domain.m_feFaceManager.m_isExternal;
  const lArray1d& externalFaceIndex = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");
  //array<real64>& external_aperture = domain.m_externalFaces.GetFieldData<realT>("aperture");
  const array<real64>& normal_approach = domain.m_externalFaces.GetFieldData<realT>("normalApproach");

  array<real64>& face_aperture = domain.m_feFaceManager.GetFieldData<realT>("Aperture");
  
  ///////////////////////////////////////
  
  array<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce> ();
  contactForce = 0.0;
  
  //update face geometry and sort faces if necessary
  realT dt = 1.0;
  domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
                                                 domain.m_discreteElementSurfaceNodes,
                                                 domain.m_discreteElementManager,
                                                 dt, true);
  
  //if a resort has been triggered, then you also need to update the contact manager
  //if(resort)
    domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);


  domain.m_externalFaces.UpdateGeometricContactProperties(dt, domain);
  
  //////////////////////////////////
  
  if( setNames_.empty() ){
    
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if(isExternal[kf])
      {
        face_aperture[kf] = std::max(-normal_approach[externalFaceIndex[kf]],0.0);
      }
    }
    
  } else {
  	
  	for( array<string>::size_type i =0; i < setNames_.size(); ++i){
      lSet& subset = domain.m_feFaceManager.GetSet(setNames_[i]);
  	  for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si ){
        localIndex kf = *si;
        if(isExternal[kf])
        {
          face_aperture[kf] = std::max( -normal_approach[externalFaceIndex[kf]],0.0);
        }
	  }
  	}
  	
  }    
  //std::cout << "Done Calculating Aperture" << std::endl;

}

REGISTER_InitialCondition( CalculateAperture )


////////////////////////////
//
//  JoinFaces
//

/**
 * @author walsh24
 * @brief make parent and child pairs of opposing faces and nodes.
 *
 */
LinkFractureFaces::LinkFractureFaces(  TICPP::HierarchicalDataNode*  hdn, const ProblemManagerT* const pm):
  InitialConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

/**
 * <JoinFaces setnames="FractureRoof FractureFloor" />
 */
void LinkFractureFaces::ReadXML( TICPP::HierarchicalDataNode*  hdn )
{
  /* empty*/
}

void LinkFractureFaces::RegisterFields( PhysicalDomainT& domain ){
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::referencePosition>();
  domain.m_feFaceManager.AddKeylessDataField( FieldInfo::R1TensorField,  PS_STR::FaceCenterStr, true, true );
}

void LinkFractureFaces::Apply( PhysicalDomainT& domain )
{
  if(setNames_.size() != 2){
    throw GPException("LinkFractureFaces:: two sets are required");
  }

  const array<R1Tensor>& u = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement> ();
  const array<R1Tensor>& X = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition> ();

  array<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>("FaceCenter");

  //nodes
  {
    std::map< std::string, lSet >::const_iterator setMapA = domain.m_feNodeManager.m_Sets.find( setNames_[0] );
    std::map< std::string, lSet >::const_iterator setMapB = domain.m_feNodeManager.m_Sets.find( setNames_[1] );

    if( setMapA != domain.m_feNodeManager.m_Sets.end() &&  setMapB != domain.m_feNodeManager.m_Sets.end() ){

      const lSet& setA = setMapA->second;
      const lSet& setB = setMapB->second;

      // fixme - brute force search for nearest neighbors.
      for( lSet::const_iterator ndA=setA.begin() ; ndA!=setA.end() ; ++ndA ) {
        realT minSqrdDist = std::numeric_limits<realT>::max();
        localIndex nbr = 0;
        R1Tensor posA =  u[*ndA] + X[*ndA];
        for( lSet::const_iterator ndB=setB.begin(); ndB!=setB.end() ; ++ndB ) {
          R1Tensor l =  u[*ndB] + X[*ndB] - posA;

          realT ll = Dot(l,l);
          if(ll < minSqrdDist){
            minSqrdDist = ll;
            nbr = *ndB;
          }
        }

        if(minSqrdDist < std::numeric_limits<realT>::max()){
          //m_nearestNodeNeighborMap[*ndA] = nbr;
          domain.m_feNodeManager.m_childIndices[*ndA] = lArray1d(1,nbr);
          domain.m_feNodeManager.m_parentIndex[nbr] = *ndA;

          if( domain.m_feNodeManager.m_childIndices[nbr].size() > 0 ){
            localIndex oldNdA = domain.m_feNodeManager.m_childIndices[nbr][0];
            R1Tensor lb =  u[nbr] + X[nbr] - u[oldNdA]-X[oldNdA];
            if(minSqrdDist < Dot(lb,lb) ){
              //m_nearestNodeNeighborMap[nbr] = *ndA;
              domain.m_feNodeManager.m_childIndices[nbr] = lArray1d(1,*ndA);
              domain.m_feNodeManager.m_parentIndex[*ndA] = nbr;
            }
          }else {
            //m_nearestNodeNeighborMap[nbr] = *ndA;
            domain.m_feNodeManager.m_childIndices[nbr] = lArray1d(1,*ndA);
            domain.m_feNodeManager.m_parentIndex[*ndA] = nbr;
          }
        }
      }
    }
  }

  //faces
  {
    std::map< std::string, lSet >::const_iterator setMapA = domain.m_feFaceManager.m_Sets.find( setNames_[0] );
    std::map< std::string, lSet >::const_iterator setMapB = domain.m_feFaceManager.m_Sets.find( setNames_[1] );

    if( setMapA != domain.m_feFaceManager.m_Sets.end() &&  setMapB != domain.m_feFaceManager.m_Sets.end() ){

      const lSet& setA = setMapA->second;
      const lSet& setB = setMapB->second;

      // update face centers
      for( lSet::const_iterator fcA=setA.begin() ; fcA!=setA.end() ; ++fcA ) {
        domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager,*fcA, faceCenter(*fcA));
      }
      for( lSet::const_iterator fcB=setB.begin() ; fcB!=setB.end() ; ++fcB ) {
        domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager,*fcB, faceCenter(*fcB));
      }

      // brute force search for nearest neighbors.
      for( lSet::const_iterator fcA=setA.begin() ; fcA!=setA.end() ; ++fcA ) {
        realT minSqrdDist = std::numeric_limits<realT>::max();
        localIndex nbr = 0;
        R1Tensor posA = faceCenter(*fcA);
        for( lSet::const_iterator fcB=setB.begin(); fcB!=setB.end() ; ++fcB ) {
          R1Tensor l =  faceCenter(*fcB) - posA;

          realT ll = Dot(l,l);
          if(ll < minSqrdDist){
            minSqrdDist = ll;
            nbr = *fcB;
          }
        }

        if(minSqrdDist < std::numeric_limits<realT>::max()){
          //m_nearestFaceNeighborMap[*fcA] = nbr;
          domain.m_feFaceManager.m_childIndices[*fcA] = lArray1d(1,nbr);
          domain.m_feFaceManager.m_parentIndex[nbr] = *fcA;

          if( domain.m_feFaceManager.m_childIndices[nbr].size() > 0 ){
            localIndex oldFcA = domain.m_feFaceManager.m_childIndices[nbr][0];
            R1Tensor lb =  faceCenter(nbr) - faceCenter(oldFcA);
            if(minSqrdDist < Dot(lb,lb) ){
              //m_nearestFaceNeighborMap[nbr] = *fcA;
              domain.m_feFaceManager.m_childIndices[nbr] = lArray1d(1, *fcA);
              domain.m_feFaceManager.m_parentIndex[*fcA] = nbr;
            }

          } else {
            //m_nearestFaceNeighborMap[nbr] = *fcA;
            domain.m_feFaceManager.m_childIndices[nbr] = lArray1d(1, *fcA);
            domain.m_feFaceManager.m_parentIndex[*fcA] = nbr;
          }
        }
      }
    }
  }
}

REGISTER_InitialCondition( LinkFractureFaces )



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////
//
// Initial condition factory

typedef std::map<std::string, InitialConditionInitializer*> InitialConditionCatalogueType; 

InitialConditionCatalogueType & getInitialConditionCatalogue(){
  static InitialConditionCatalogueType theCatalogue ;
  return theCatalogue;
}

void getInitialConditionNames( std::vector<std::string>& nameList){
  for(InitialConditionCatalogueType::const_iterator it = getInitialConditionCatalogue().begin(); 
      it != getInitialConditionCatalogue().end(); ++it){
        nameList.push_back(it->first);
  }
}

InitialConditionBase* newInitialCondition(const std::string& InitialConditionName , TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm) 
{
  
  InitialConditionInitializer* InitialConditionInitializer = getInitialConditionCatalogue()[InitialConditionName];
  InitialConditionBase *theNewInitialCondition = NULL;
  
  if(!InitialConditionInitializer)
      throw GPException("Could not create unrecognized InitialCondition"+ InitialConditionName);

  theNewInitialCondition = InitialConditionInitializer->initializeInitialCondition( hdn,pm );
  
  return theNewInitialCondition;
}
