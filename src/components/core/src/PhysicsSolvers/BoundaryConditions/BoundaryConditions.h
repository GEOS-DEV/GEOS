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
 * @file BoundaryConditions.h
 * @author walsh24
 * @date December 5, 2011 
 */

#ifndef BoundaryConditions_H_
#define BoundaryConditions_H_

#include "DataStructures/Tables/Table.h"
#include "Common/Common.h"
#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/FieldTypeMultiPtr.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "DataStructures/Tables/Table.h"

#include<map>
#include<string>
#include<vector>
#include "../../legacy/Common/GPException.h.old"

class ProblemManagerT;
class Function;

/// Base class for all initial conditions
class BoundaryConditionBase
{
public:
  BoundaryConditionBase(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                        const ProblemManagerT* const problemManager);

  virtual ~BoundaryConditionBase();

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  virtual void RegisterFields(PhysicalDomainT& domain)
  {/*empty*/
  }


  virtual const char* GetBoundaryConditionName() = 0;

  virtual realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si,
                         realT time) = 0;

  // This function call should be used rather than upcasting directly to enable switch and other conditional boundary conditions
  // eg.
  // TractionBoundaryCondition* trbc = bc->UpcastActiveBCPointer<TractionBoundaryCondition*>();
  template<typename Type>
  Type* UpcastActiveBCPointer(realT time)
  {
    return dynamic_cast<Type*>(this->GetActiveBCPointer(time));
  }

  virtual BoundaryConditionBase* GetActiveBCPointer(realT time)
  {
    return this;
  } // enables switch and other conditional boundary conditions

  virtual const std::string& GetFieldName(realT time)
  {
    return m_fieldName;
  }

  virtual const FieldType& GetFieldType(realT time)
  {
    return m_fieldType;
  }
  virtual int GetComponent(realT time)
  {
    return m_component;
  }
  virtual const R1Tensor& GetDirection(realT time)
  {
    return m_direction;
  }


  virtual void SetUpdateTime(realT time){
	  m_time = time;
  }

  virtual realT GetLastUpdateTime(){
	  return m_time;
  }

  realT GetTimeFactor() const {return m_timeFactor;}

  realT GetStartTime() const { return m_startTime; }
  void SetStartTime( realT time ) {  m_startTime = time; }

  realT GetEndTime() const { return m_endTime; }
  void SetEndTime( realT time ) {  m_endTime = time; }

  // class data
  PhysicalDomainT::ObjectDataStructureKeys m_objectKey; // object the boundary condition is applied to: eg. face, node, edge.
  sArray1d m_setNames; // sets the boundary condition is applied to
  std::string m_regionName;

  bool m_isConstantInSpace;
  bool m_isConstantInTime;

  int getOption() const {return m_option;}
protected:

  std::string m_fieldName; // the name of the field the boundary condition is applied to or a description of the boundary condition.
  FieldType m_fieldType; // the variable type stored in the boundary condition field.
  int m_component; // the component the boundary condition acts on (-ve indicates that direction should be used).
  R1Tensor m_direction; // the direction the boundary condition acts in.

  realT m_time; // last time bc value was calculated
  realT m_timeFactor; // a general time value for use by the BC.

  realT m_startTime = -std::numeric_limits<realT>::max();
  realT m_endTime = std::numeric_limits<realT>::max();
  int m_option = 0;

};

/// Class to replicate behavior of original boundary condition data class
class SimpleBoundaryCondition: public BoundaryConditionBase
{
public:
  SimpleBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                          const ProblemManagerT* const problemManager);
  virtual ~SimpleBoundaryCondition()
  {
  }


  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  static const char* BoundaryConditionName()
  {
    return "BoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  std::string m_timeTableName;
  realT m_scale;
  R1Tensor m_position;

  Function* m_function;
  std::string m_functionName;

  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si,
                 realT time);

protected:
  realT m_value;

};


/// Set boundary condition using a function based on fields described on the object.
class BoundaryConditionFunction: public BoundaryConditionBase
{
public:
	BoundaryConditionFunction( TICPP::HierarchicalDataNode* BoundaryConditionNode, const ProblemManagerT* const problemManager);
  virtual ~BoundaryConditionFunction(){};

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* BoundaryConditionName(){return "BoundaryConditionFunction";};
  virtual const char* GetBoundaryConditionName(){return BoundaryConditionName();};

  realT GetValue(const ObjectDataStructureBaseT& object,
                 const lSet::const_iterator& si,realT time);

  protected:
    realT m_value;
    realT m_scale;
    realT m_offset;

  private:
    localIndex m_nVars; // number of variables
    //localIndex m_xLength; // length of input vector = sum (variables*number of components)


    std::vector<FieldTypeMultiPtr> m_fieldPtrs;
    Array1dT<FieldType> m_variableTypes;
    sArray1d m_variableNames;
    std::vector<realT> m_x;  // function input vector

    Function* m_function;
    std::string m_functionName;
};


/// Traction Boundary condition
class TractionBoundaryCondition: public SimpleBoundaryCondition
{
public:
  TractionBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                            const ProblemManagerT* const problemManager);
  virtual ~TractionBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "TractionBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  virtual bool IsNormalTraction()
  {
    return m_useNormalFlag;
  }

  virtual R1Tensor GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,
                                     realT& time);

protected:
  bool m_useNormalFlag;

};


/// Set traction boundary condition using a function based on fields described on the object.
class TractionBoundaryConditionFunction: public TractionBoundaryCondition
{
public:
  TractionBoundaryConditionFunction( TICPP::HierarchicalDataNode* BoundaryConditionNode, const ProblemManagerT* const problemManager);
  virtual ~TractionBoundaryConditionFunction(){};

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* BoundaryConditionName(){return "TractionBoundaryConditionFunction";};
  virtual const char* GetBoundaryConditionName(){return BoundaryConditionName();};

  virtual R1Tensor GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,
                                     realT& time);

  private:
    localIndex m_nVars; // number of variables
    //localIndex m_xLength; // length of input vector = sum (variables*number of components)

    std::vector<FieldTypeMultiPtr> m_fieldPtrs;
    Array1dT<FieldType> m_variableTypes;
    sArray1d m_variableNames;
    std::vector<realT> m_x;  // function input vector
    realT m_scale;
    realT m_offset;

};


/// Hydraulic Pressure Boundary condition
class HydraulicPressureBoundaryCondition: public TractionBoundaryCondition
{
public:
  HydraulicPressureBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                     const ProblemManagerT* const problemManager);
  virtual ~HydraulicPressureBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "HydraulicPressureBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  virtual R1Tensor GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,
                                     realT& time);
protected:
  bool m_useNormalFlag;

};

/// Uniform Pressure Boundary condition
class UniformPressureBoundaryCondition: public TractionBoundaryCondition
{
public:
  UniformPressureBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                   const ProblemManagerT* const problemManager);
  virtual ~UniformPressureBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "UniformPressureBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  virtual R1Tensor GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,
                                     realT& time);
protected:
  bool m_useNormalFlag;

};

class MultiVarBoundaryConditionBase: public BoundaryConditionBase
{
public:

  MultiVarBoundaryConditionBase(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                const ProblemManagerT* const problemManager);
  virtual ~MultiVarBoundaryConditionBase()
  {
  }


  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  static const char* BoundaryConditionName()
  {
    return "MultiVarBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

  sArray1d m_varName;
protected:

  //realT m_time; // last time bc value was calculated
  rArray1d m_value;
  

};

class MultiVarDirichletBoundaryCondition: public MultiVarBoundaryConditionBase
{
public:

  MultiVarDirichletBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                     const ProblemManagerT* const problemManager);
  virtual ~MultiVarDirichletBoundaryCondition()
  {
  }


  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  static const char* BoundaryConditionName()
  {
    return "MultiVarDirichletBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si, realT time)
  {
    return 0.0;
  }

  const rArray1d& GetValues(realT time);

  void CheckVars(const sArray1d &nvarName);

  bool isClamped() const
  {
    return m_isClamped;
  }

protected:

  bool m_isClamped;
  sArray1d m_tables;

};

class MultiVarSrcFluxBoundaryCondition: public MultiVarBoundaryConditionBase
{
public:

  MultiVarSrcFluxBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                   const ProblemManagerT* const problemManager);
  virtual ~MultiVarSrcFluxBoundaryCondition()
  {
  }


  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  static const char* BoundaryConditionName()
  {
    return "MultiVarSrcFluxBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }


  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si, realT time)
  {
    return 0.0;
  }


  const rArray1d& GetValues(realT time);

  rArray1d &allocFac()
  {
    return m_allocFac;
  }

  const rArray1d &allocFac() const
  {
    return m_allocFac;
  }

  bool isAllocedByWeight() const
  {
    return m_allocedByWeight;
  }

  void CheckVars(const sArray1d &varNames);

protected:

  bool m_allocedByWeight;
  rArray1d m_allocFac;
  sArray1d m_tables;

};



/// Wall boundary condition
class WallBoundaryCondition: public SimpleBoundaryCondition
{
public:
  WallBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                        const ProblemManagerT* const problemManager);
  virtual ~WallBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "WallBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

};

/// Outflow boundary condition
class OutflowBoundaryCondition: public SimpleBoundaryCondition
{
public:
  OutflowBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                           const ProblemManagerT* const problemManager);
  virtual ~OutflowBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "OutflowBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

};

/// Inflow boundary condition
class InflowBoundaryCondition: public SimpleBoundaryCondition
{
public:
  InflowBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                          const ProblemManagerT* const problemManager);
  virtual ~InflowBoundaryCondition()
  {
  }


  static const char* BoundaryConditionName()
  {
    return "InflowBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

};

/// Non-penetrating boundary condition

class NonPenetratingBoundaryCondition: public BoundaryConditionBase
{
public:
  NonPenetratingBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                  const ProblemManagerT* const problemManager);
  virtual ~NonPenetratingBoundaryCondition()
  {
  }

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si,
                 realT time);

  static const char* BoundaryConditionName()
  {
    return "NonPenetratingBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

  virtual void RegisterFields(PhysicalDomainT& domain);

  void UpdateNearestNeighborMaps(PhysicalDomainT& domain);

  std::map<localIndex, localIndex> m_nearestNodeNeighborMap;
  std::map<localIndex, localIndex> m_nearestFaceNeighborMap;
  lSet m_contactingNodes;

  realT m_bulkModulus;

  bool m_updatePressure;
};

/// Single partition periodic boundary condition
/// Warning - this is a temporary measure until periodic boundary conditions are fully supported.
class SinglePartitionPeriodicBoundaryCondition: public BoundaryConditionBase
{
public:
  SinglePartitionPeriodicBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                           const ProblemManagerT* const problemManager);
  virtual ~SinglePartitionPeriodicBoundaryCondition()
  {
  }

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si,
                 realT time);

  static const char* BoundaryConditionName()
  {
    return "SinglePartitionPeriodicBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

  virtual void RegisterFields(PhysicalDomainT& domain);

  void SetNeighborMaps(PhysicalDomainT& domain);

  std::map<localIndex, localIndex> m_nodeNeighborMapA;
  std::map<localIndex, localIndex> m_nodeNeighborMapB;
  std::map<localIndex, localIndex> m_faceNeighborMapA;
  std::map<localIndex, localIndex> m_faceNeighborMapB;
  std::map<localIndex, localIndex> m_edgeNeighborMapA;
  std::map<localIndex, localIndex> m_edgeNeighborMapB;

  int m_dimension;
};

/// Alternate between multiple boundary conditions
class SwitchBoundaryConditions: public BoundaryConditionBase
{
public:
  SwitchBoundaryConditions(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                           const ProblemManagerT* const problemManager);
  virtual ~SwitchBoundaryConditions()
  {
  }

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);

  realT GetValue(const ObjectDataStructureBaseT& object, const lSet::const_iterator& si,
                 realT time);

  static const char* BoundaryConditionName()
  {
    return "SwitchBoundaryConditions";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

  const std::string& GetFieldName(realT time);

  const FieldType& GetFieldType(realT time);
  int GetComponent(realT time);
  const R1Tensor& GetDirection(realT time);

protected:

  virtual BoundaryConditionBase* GetActiveBCPointer(realT time);

  std::vector<BoundaryConditionBase*> m_children;
  Function* m_function;
  std::string m_functionName;
  unsigned m_default_index;
  const ProblemManagerT* m_pmPtr;
};

//////////////////////////

// Boundary Condition Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newBoundaryCondition"
//   * A base class to derive the functions to generate BoundaryCondition pointers: "BoundaryConditionInitializer"
//   * A String-to-BoundaryCondition-Intializer map hidden behind the getBoundaryConditionCatalogue function
//   * A template to create BoundaryCondition initializers: "BoundaryConditionRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_BoundaryCondition"
// 
// Most boundary conditions will only need one or two of the parts:
//   * To register a new BoundaryCondition in the factory: REGISTER_BoundaryCondition( BoundaryConditionClassName )
//   * To load a BoundaryCondition pointer from the factory:       BoundaryConditionBase* aBoundaryConditionPtr = newBoundaryCondition(BoundaryConditionString, args );

/// The BoundaryCondition Factory.
BoundaryConditionBase* newBoundaryCondition(const std::string& BoundaryConditionName,
                                            TICPP::HierarchicalDataNode* hdn,
                                            const ProblemManagerT* const pm);

/// Base class to generate new BoundaryCondition pointers
class BoundaryConditionInitializer
{
public:
  virtual BoundaryConditionBase* initializeBoundaryCondition(TICPP::HierarchicalDataNode* hdn,
                                                             const ProblemManagerT* const pm) = 0;
  virtual ~BoundaryConditionInitializer() = 0;
};

inline BoundaryConditionInitializer::~BoundaryConditionInitializer()
{
}

/// Interface to the BoundaryCondition name -> BoundaryCondition initializer map
std::map<std::string, BoundaryConditionInitializer*> & getBoundaryConditionCatalogue();

/// Return a list of supported BoundaryCondition names
void getBoundaryConditionNames(std::vector<std::string>& nameList);

/// Template for creating classes derived from BoundaryConditionInitializer
template<class BoundaryConditionType>
class BoundaryConditionRegistrator: public BoundaryConditionInitializer
{

public:
  BoundaryConditionRegistrator(void)
  {
    std::string BoundaryConditionName = std::string(BoundaryConditionType::BoundaryConditionName());
    getBoundaryConditionCatalogue()[BoundaryConditionName] = this;
  }

  BoundaryConditionBase* initializeBoundaryCondition(TICPP::HierarchicalDataNode* hdn,
                                                     const ProblemManagerT* const pm)
  {
    return new BoundaryConditionType(hdn, pm);
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_BoundaryCondition( ClassName ) namespace{ BoundaryConditionRegistrator<ClassName> reg_##ClassName; }

#endif
