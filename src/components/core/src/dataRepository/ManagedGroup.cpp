/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "dataRepository/SidreWrapper.hpp"

namespace geosx
{
namespace dataRepository
{

ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent ) :
  ManagedGroup( name, parent, nullptr )
{}


ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent,
                            cxx_utilities::DocumentationNode * docNode ) :
  m_docNode(docNode),
  m_keyLookup(),
  m_wrappers(),
  m_parent(parent),
  m_subGroups(),
  m_sidreGroup(nullptr)
{

  // SIDRE interaction
  asctoolkit::sidre::DataGroup * sidreParent = nullptr;
  if( m_parent==nullptr )
  {
    sidreParent = SidreWrapper::dataStore().getRoot();
  }
  else
  {
    sidreParent = parent->m_sidreGroup;
  }

  if( sidreParent->hasGroup(name) )
  {
    m_sidreGroup = sidreParent->getGroup(name);
  }
  else
  {
    m_sidreGroup = sidreParent->createGroup(name);
  }



  // Setup DocumentationNode
  if( parent != nullptr )
  {
    if( parent->m_docNode != nullptr )
    {
      m_docNode = parent->m_docNode->AllocateChildNode( name,
                                                        name,
                                                        0,
                                                        "ManagedGroup",
                                                        "",
                                                        "ManagedGroup",
                                                        "ManagedGroup",
                                                        "",
                                                        parent->getName(),
                                                        0,
                                                        0 ) ;

    }
  }


  m_docNode->AllocateChildNode( "size",
                                "size",
                                -1,
                                "int32",
                                "",
                                "size of group",
                                "Number of entries in this group.",
                                "0",
                                "",
                                0,
                                0 );

  m_docNode->AllocateChildNode( "name",
                                "name",
                                -1,
                                "string",
                                "string",
                                "name of group",
                                "name of group.",
                                name,
                                "",
                                0,
                                0 );

  *(RegisterViewWrapper<int32>( "size" ).data()) = 0;
  RegisterViewWrapper<std::string>( "name" ).reference() = name;
  RegisterViewWrapper<std::string>( "path" );

}

ManagedGroup::~ManagedGroup()
{
  // TODO Auto-generated destructor stub
}

//DataObjectManager::DataObjectManager( DataObjectManager const & source ):
//    m_size( source.m_size ),
//    m_name( source.m_name ),
//    m_path( source.m_path ),
//    m_keyLookup( source.m_keyLookup ),
//    m_dataObjects( source.m_dataObjects ),
//    m_parent( source.m_parent )
//{}

ManagedGroup::ManagedGroup( ManagedGroup&& source ) :
  m_keyLookup( std::move(source.m_keyLookup) ),
  m_wrappers( std::move(source.m_wrappers) ),
  m_parent( std::move(source.m_parent) )
{}

ManagedGroup::CatalogInterface::CatalogType& ManagedGroup::GetCatalog()
{
  static ManagedGroup::CatalogInterface::CatalogType catalog;
  return catalog;
}

ViewWrapperBase& ManagedGroup::RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type )
{
  return *( rtTypes::ApplyTypeLambda( type,
                                      [this, &name]( auto a ) -> ViewWrapperBase*
      {
        return &( this->RegisterViewWrapper<decltype(a)>(name) );
      } ) );
}

void ManagedGroup::resize( int32 const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    i->resize(newsize);
  }
  *(this->getWrapper<int32>( keys::Size ).data())=newsize;
}



void ManagedGroup::RegisterDocumentationNodes()
{
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
    if( subNode.second.getDataType() != "DocumentationNode" )
    {
      RegisterViewWrapper( subNode.second.getStringKey(),
                           rtTypes::typeID(subNode.second.getDataType() ) );
    }
  }
}

void ManagedGroup::BuildDataStructure( dataRepository::ManagedGroup * const rootGroup )
{
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->BuildDataStructure( rootGroup );
  }
}

void ManagedGroup::FillDocumentationNode( dataRepository::ManagedGroup * const  )
{

}


void ManagedGroup::SetDocumentationNodes( dataRepository::ManagedGroup * const group )
{
  FillDocumentationNode(group);
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetDocumentationNodes(group);
  }
}


void ManagedGroup::ReadXML( pugi::xml_node const & targetNode )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  for( auto const & subDocEntry : docNode->m_child )
  {
    cxx_utilities::DocumentationNode subDocNode = subDocEntry.second;
    
    if (subDocNode.getIsInput() == 1)
    {
      std::string childType = subDocNode.getSchemaType();    

      switch (rtTypes::typeID(childType))
      {
        case rtTypes::TypeIDs::real64_id:
        {
          real64 defVal = atof(subDocNode.getDefault().c_str());
          real64 xmlVal = targetNode.attribute(subDocNode.getStringKey().c_str()).as_double(defVal);
          *(this->getData<real64>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::int32_id:
        {
          int32 defVal = atoi(subDocNode.getDefault().c_str());
          int32 xmlVal = targetNode.attribute(subDocNode.getStringKey().c_str()).as_int(defVal);
          *(this->getData<int32>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::uint32_id:
        {
          uint32 defVal = atol(subDocNode.getDefault().c_str());
          uint32 xmlVal = targetNode.attribute(subDocNode.getStringKey().c_str()).as_uint(defVal);
          *(this->getData<uint64>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::string_id:
        {
          string defVal = subDocNode.getDefault();
          string xmlVal = targetNode.attribute(subDocNode.getStringKey().c_str()).value();
          this->getData<string>(subDocNode.getStringKey()) = xmlVal.empty() ? defVal : xmlVal;
          break;
        }
        // TODO: Define the assignment operator for std::vector<type> in viewWrapper
        //       Alternatively, pass the array reference into the load_type_array methods directly
        //       (This would require the push_back method to work on the viewWrappers)
        case rtTypes::TypeIDs::real64_array_id:
        {
          string defVal = subDocNode.getDefault();
          std::vector<real64> xmlVal;
          targetNode.attribute(subDocNode.getStringKey().c_str()).load_double_array(xmlVal, defVal);
          // *(this->getData<real64_array>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::int32_array_id:
        {
          string defVal = subDocNode.getDefault();
          std::vector<int32> xmlVal;
          targetNode.attribute(subDocNode.getStringKey().c_str()).load_int_array(xmlVal, defVal);
          // *(this->getData<int32_array>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::uint32_array_id:
        {
          string defVal = subDocNode.getDefault();
          std::vector<uint32> xmlVal;
          targetNode.attribute(subDocNode.getStringKey().c_str()).load_uint_array(xmlVal, defVal);
          // *(this->getData<uint32_array>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        case rtTypes::TypeIDs::string_array_id:
        {
          string defVal = subDocNode.getDefault();
          std::vector<string> xmlVal;
          targetNode.attribute(subDocNode.getStringKey().c_str()).load_string_array(xmlVal, defVal);
          // *(this->getData<string_array>(subDocNode.getStringKey())) = xmlVal;
          break;
        }
        default:
        {
          // TODO: Define missing cases
          throw std::invalid_argument("XML auto read method not defined");
        }
      }
    }
  }
}



}
} /* namespace ODS */
