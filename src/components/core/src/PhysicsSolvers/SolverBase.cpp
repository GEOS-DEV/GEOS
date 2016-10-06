/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::ReadXML( pugi::xml_node const & solverNode )
{
  *(this->getData<real64>(keys::courant)) = solverNode.attribute("courant").as_double(0.5);
}

void SolverBase::BuildDataStructure( dataRepository::ManagedGroup * const /*domain*/ )
{
  this->RegisterViewWrapper<real64>(keys::maxDt);
  this->RegisterViewWrapper<real64>(keys::courant);
}

void SolverBase::FillDocumentationNode( dataRepository::ManagedGroup * const /*group*/ )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed groups, this could be done automatically
  docNode->setSchemaType("Node");

  cxx_utilities::DocumentationNode docVar;
  
  // Document cfl number
  docVar.m_name               = "cfl";
  docVar.m_stringKey          = "cfl";
  docVar.m_intKey             = -1;
  docVar.m_dataType           = "real64";
  docVar.m_schemaType         = "double";
  docVar.m_shortDescription   = "Courant–Friedrichs–Lewy (CFL) factor";
  docVar.m_longDescription    = "Courant–Friedrichs–Lewy (CFL) factor is multiplied with CFL condition to reduce/increase "
                                "allowable timestep.";
  docVar.m_default            = "1.0";
//  docVar.m_path               = "";
  docVar.m_level              = m_docNode->m_level + 1;
  docVar.m_isInput            = 1;

  docNode->m_child.insert( { docVar.m_name, docVar } );




}


void SolverBase::Initialize( dataRepository::ManagedGroup& /*domain*/ )
{
  *(this->getData<real64>(keys::courant)) = std::numeric_limits<real64>::max();
}

} /* namespace ANST */
