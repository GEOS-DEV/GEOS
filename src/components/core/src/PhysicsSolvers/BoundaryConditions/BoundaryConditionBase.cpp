#include "BoundaryConditionBase.hpp"

namespace geosx
{
using namespace dataRepository;

BoundaryConditionBase::BoundaryConditionBase( string const & name, ManagedGroup *const parent ):
  ManagedGroup(name,parent)
{}


BoundaryConditionBase::~BoundaryConditionBase()
{}

BoundaryConditionBase::CatalogInterface::CatalogType& BoundaryConditionBase::GetCatalog()
{
  static BoundaryConditionBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void BoundaryConditionBase::FillDocumentationNode( dataRepository::ManagedGroup * const )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

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
                              "",
                              "Direction to apply boundary condition to",
                              "",
                              "{0,0,0}",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::timeTableName,
                              keys::timeTableName,
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
}


real64 BoundaryConditionBase::GetValue( realT time )
{
  real64 rval = m_scale;
  if (!(m_timeTableName.empty()))
  {
    rArray1d t(1);
    t[0] = time;
    real32 const tableval = TableManager::Instance().LookupTable<1>(m_timeTableName, t);
    rval = m_scale * tableval;
  }
  else if (!(m_functionName.empty()))
  {
//    rval = m_scale * (*m_function)(time);
  }
  return rval;
}





}
