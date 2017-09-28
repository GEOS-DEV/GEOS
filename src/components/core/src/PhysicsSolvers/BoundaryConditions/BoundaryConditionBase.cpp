#include "BoundaryConditionBase.hpp"

namespace geosx
{
using namespace dataRepository;

BoundaryConditionBase::BoundaryConditionBase( string const & name, ManagedGroup * parent ):
  ManagedGroup(name,parent)
{}


BoundaryConditionBase::~BoundaryConditionBase()
{}

BoundaryConditionBase::CatalogInterface::CatalogType& BoundaryConditionBase::GetCatalog()
{
  static BoundaryConditionBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void BoundaryConditionBase::FillDocumentationNode( dataRepository::ManagedGroup * UNUSED )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( keys::initialCondition,
                              keys::initialCondition,
                              -1,
                              "int32",
                              "int32",
                              "BC is applied as an initial condition.",
                              "",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::setNames,
                              keys::setNames,
                              -1,
                              "string_array",
                              "string_array",
                              "Name of sets that boundary condition is applied to.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );


  docNode->AllocateChildNode( keys::elementRegionName,
                              keys::elementRegionName,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::fieldName,
                              keys::fieldName,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::dataType,
                              keys::dataType,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::component,
                              keys::component,
                              -1,
                              "int32",
                              "int32",
                              "Component of field (if tensor) to apply boundary condition to",
                              "",
                              "-1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::direction,
                              keys::direction,
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Direction to apply boundary condition to",
                              "",
                              "{0,0,0}",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::functionName,
                              keys::functionName,
                              -1,
                              "string",
                              "string",
                              "Name of table that specifies variation of BC with time.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::bcApplicationTableName,
                              keys::bcApplicationTableName,
                              -1,
                              "string",
                              "string",
                              "Name of table that specifies the on/off application of the bc.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::scale,
                              keys::scale,
                              -1,
                              "real64",
                              "real64",
                              "Component of field (if tensor) to apply boundary condition to",
                              "",
                              "-1",
                              "",
                              0,
                              1,
                              0 );
}

void BoundaryConditionBase::ReadXML_PostProcess()
{
  m_setNames = this->getReference<string_array>( keys::setNames );
  m_elementRegionName = this->getReference<string>( keys::elementRegionName );
  m_fieldName = this->getReference<string>( keys::fieldName );
  m_component = this->getReference<int32>( keys::component );
  m_direction = this->getReference<R1Tensor>( keys::direction );
  m_functionName          = this->getReference<string>( keys::functionName );
  m_bcApplicationFunctionName = this->getReference<string>( keys::bcApplicationTableName );
  m_scale                  = this->getReference<real64>( keys::scale );
  m_initialCondition = this->getReference<int32>(keys::initialCondition);
}

//real64 BoundaryConditionBase::GetValue( realT time ) const
//{
//
//  real64 rval = m_scale;
//  if (!(m_functionName.empty()))
//  {
//    rArray1d t(1);
//    t[0] = time;
//    real64 const tableval = TableManager::Instance().LookupTable<1>(m_functionName, t);
//    rval = m_scale * tableval;
//  }
////  else if (!(m_spaceFunctionName.empty()))
////  {
//////    rval = m_scale * (*m_function)(time);
////  }
//  return rval;
//}
//
//void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const & set,
//                                                                real64 const time,
//                                                                array<R1Tensor> const & X,
//                                                                array<R1Tensor> & field )
//{
//
//  int32 const component = GetComponent();
//  string const functionName = getData<string>(dataRepository::keys::functionName);
//  NewFunctionManager * functionManager = NewFunctionManager::Instance();
//
//  FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
//
//  if( function!=nullptr )
//  {
//    real64 const factor = m_scale * ( timeFunction->Evaluate( &time ) );
//    for( auto a : set )
//    {
//      field[a][component] = factor;
//    }
//  }
//  else
//  {
//    for( auto a : set )
//    {
//      field[a][component] = m_scale;
//    }
//  }
//}

void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const & set,
                                                                real64 const time,
                                                                ManagedGroup * dataGroup,
                                                                string const & fieldName ) const
{

  int32 const component = GetComponent();
  string const functionName = getData<string>(dataRepository::keys::functionName);
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index(vw->get_typeid());

  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID(typeIndex) , [&]( auto type ) -> void
  {
    using fieldType = decltype(type);
    ViewWrapper<fieldType> & view = dynamic_cast< ViewWrapper<fieldType> & >(*vw);
    view_rtype<fieldType> field = view.data();
    if( functionName.empty() )
    {
      for( auto a : set )
      {
//        field[a] = m_scale;
        rtTypes::equate( field[a], component, m_scale );
      }
    }
    else
    {
      FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
      if( function!=nullptr)
      {
        if( function->isFunctionOfTime()==2 )
        {
          real64 value = m_scale * function->Evaluate( &time );
          for( auto a : set )
          {
            rtTypes::equate( field[a], component, value );
          }
        }
        else
        {
          real64_array result(set.size());
          function->Evaluate( dataGroup, time, set, result );
          int32 count=0;
          for( auto a : set )
          {
            rtTypes::equate( field[a], component, result[count] );
            ++count;
          }
        }
      }
    }
  });
}




}
