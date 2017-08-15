

#ifndef BOUNDARYCONDITIONBASE_H
#define BOUNDARYCONDITIONBASE_H

#include "common/DataTypes.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "managers/TableManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"

namespace geosx
{
class Function;

namespace dataRepository
{
namespace keys
{
string const objectPath = "objectPath";
string const setNames = "setNames";
string const fieldName = "fieldName";
string const dataType = "dataType";
string const component("component");
string const direction("direction");
string const timeFunctionName("timeFunction");
string const spaceFunctionName("spaceFunction");
string const bcApplicationTableName("bcApplicationTableName");
string const scale("scale");
string const functionName("functionName");
}
}

class BoundaryConditionBase : public dataRepository::ManagedGroup
{
public:

  using CatalogInterface = cxx_utilities::CatalogInterface< BoundaryConditionBase, string const &, dataRepository::ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  BoundaryConditionBase( string const & name, dataRepository::ManagedGroup *const parent );

  virtual ~BoundaryConditionBase();



  void FillDocumentationNode( dataRepository::ManagedGroup * const ) override;

  void ReadXML_PostProcess() override final;




  real64 GetValue( realT time ) const;


  template< typename T >
  void ApplyBounaryConditionDefaultMethod( lSet const & set,
                                           real64 const time,
                                           array<R1Tensor> const & X,
                                           array<T> & field );

  void ApplyBounaryConditionDefaultMethod( lSet const & set,
                                           real64 const time,
                                           array<R1Tensor> const & X,
                                           array<R1Tensor> & field );

  void ApplyBounaryConditionDefaultMethod( lSet const & set,
                                           real64 const time,
                                           array<R1Tensor> const & X,
                                           dataRepository::ManagedGroup & dataGroup,
                                           string const & fieldname );


  virtual const string& GetFieldName()
  {
    return m_fieldName;
  }

  virtual int GetComponent() const
  {
    return m_component;
  }

  virtual const R1Tensor& GetDirection(realT time)
  {
    return m_direction;
  }

  real64 GetStartTime()
  {
    return -1;
  }

  real64 GetEndTime()
  {
    return 1.0e9;
  }

  string_array const & GetSetNames() const
  {
    return m_setNames;
  }

  int initialCondidion() const
  {
    return m_initialCondition;
  }

protected:

  string_array m_setNames; // sets the boundary condition is applied to

  string m_fieldName;    // the name of the field the boundary condition is applied to or a description of the boundary condition.

  string m_dataType;
  // TODO get rid of components. Replace with direction only.

  int m_component;       // the component the boundary condition acts on (-ve indicates that direction should be used).
  R1Tensor m_direction;  // the direction the boundary condition acts in.

  int m_initialCondition;

  string m_timeFunctionName;
  string m_spaceFunctionName;
  string m_bcApplicationFunctionName;

  real64 m_scale;




};


template< typename T >
void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const & set,
                                                                real64 const time,
                                                                array<R1Tensor> const & X,
                                                                array<T> & field )
{

  string const spaceFunctionName = getData<string>(dataRepository::keys::spaceFunctionName);
  string const timeFunctionName = getData<string>(dataRepository::keys::timeFunctionName);
  NewFunctionManager & functionManager = NewFunctionManager::Instance();

  FunctionBase const * const spaceFunction = functionManager.GetGroupPtr<FunctionBase>(spaceFunctionName);
  FunctionBase const * const timeFunction  = functionManager.GetGroupPtr<FunctionBase>(timeFunctionName);

  if( timeFunction!=nullptr && spaceFunction!=nullptr )
  {
//    real64 const tfactor = m_scale * ( timeFunction->Evaluate( &time ) );
//    for( auto a : set )
//    {
//      field[a][component] = tfactor * ( spaceFunction->Evaluate( &(X[a]) ) );
//    }
  }
  else if( timeFunction!=nullptr )
  {
    real64 const factor = m_scale * ( timeFunction->Evaluate( &time ) );
    for( auto a : set )
    {
      field[a] = factor;
    }
  }
  else if( spaceFunction!=nullptr )
  {
    for( auto a : set )
    {
      field[a] = m_scale * ( spaceFunction->Evaluate( X[a].Data() ) );
    }

  }
  else
  {
    for( auto a : set )
    {
      field[a] = m_scale;
    }
  }
}





}
#endif
