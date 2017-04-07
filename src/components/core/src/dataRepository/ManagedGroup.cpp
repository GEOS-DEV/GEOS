/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "dataRepository/SidreWrapper.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{
namespace dataRepository
{

ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent ) :
  m_docNode(nullptr),
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
    if( parent->m_docNode != nullptr  )
    {
      if( this->m_docNode == nullptr )
      {
        m_docNode = parent->m_docNode->AllocateChildNode( name,
                                                          name,
                                                          0,
                                                          "",
                                                          "Node",
                                                          "",
                                                          "",
                                                          "",
                                                          parent->getName(),
                                                          0,
                                                          0 ) ;
      }
    }
    else
    {
      m_docNode = new cxx_utilities::DocumentationNode( name,
                                                        name,
                                                        -1,
                                                        "Node",
                                                        "Node",
                                                        "The Root DocumentationNode for " + name,
                                                        "",
                                                        "",
                                                        "",
                                                        0,
                                                        0,
                                                        0,
                                                        nullptr );
    }
  }
  else
  {
    m_docNode = new cxx_utilities::DocumentationNode( name,
                                                      name,
                                                      -1,
                                                      "Node",
                                                      "Node",
                                                      "The Root DocumentationNode for " + name,
                                                      "",
                                                      "",
                                                      "",
                                                      0,
                                                      0,
                                                      0,
                                                      nullptr );
  }


  m_docNode->AllocateChildNode( "size",
                                "size",
                                -1,
                                "int32",
                                "int32",
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
    if( parent->m_docNode != nullptr && this->m_docNode != nullptr )
    {
      m_docNode = parent->m_docNode->AllocateChildNode( name,
                                                        name,
                                                        0,
                                                        "ManagedGroup",
                                                        "Node",
                                                        "ManagedGroup",
                                                        "ManagedGroup",
                                                        "",
                                                        parent->getName(),
                                                        0,
                                                        0 ) ;

    }
    else
    {

    }
  }


  m_docNode->AllocateChildNode( "size",
                                "size",
                                -1,
                                "int32",
                                "int32",
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
  delete m_docNode;
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
  return *( rtTypes::ApplyTypeLambda1( type,
                                       [this, &name]( auto a ) -> ViewWrapperBase*
      {
        return &( this->RegisterViewWrapper<decltype(a)>(name) );
      } ) );
}

void ManagedGroup::resize( int32 const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    if( i->sizedFromParent() == 1 )
    {
      i->resize(newsize);
    }
  }
  *(this->getWrapper<int32>( keys::Size ).data())=newsize;
}



void ManagedGroup::RegisterDocumentationNodes()
{
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
//    std::cout<<subNode.first<<", "<<subNode.second.getName()<<std::endl;
    if( ( subNode.second.getSchemaType() != "DocumentationNode" ) &&
        ( subNode.second.getSchemaType() != "Node" ) )
    {
      std::cout<<"Register "<<subNode.second.getStringKey()<<" of type "<<subNode.second.getDataType()<<std::endl;
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
  
  ReadXMLsub( targetNode );

  for( auto const & subDocEntry : docNode->m_child )
  {
    cxx_utilities::DocumentationNode subDocNode = subDocEntry.second;
    
    if (subDocNode.getIsInput() == 1)
    {
      std::string childType = subDocNode.getSchemaType();    

#if 1
      rtTypes::TypeIDs const typeID = rtTypes::typeID(childType);
      rtTypes::ApplyIntrinsicTypeLambda2 ( typeID,
                                           [this, typeID, &targetNode, &subDocNode]( auto a, auto b ) -> void
      {
        string defVal = subDocNode.getDefault();

        pugi::xml_attribute xmlatt = targetNode.attribute(subDocNode.getStringKey().c_str());
        ViewWrapper<decltype(a)>& dataView = this->getWrapper<decltype(a)>(subDocNode.getStringKey());
        std::vector<decltype(b)> xmlVal;
        typename ViewWrapper<decltype(a)>::rtype data = dataView.data();

        if( !xmlatt.empty() )
        {
          xmlatt.as_type(xmlVal, defVal);
        }
        else
        {
          if( defVal == "REQUIRED")
          {
            string message = "variable " + subDocNode.getName() + " is required in " + targetNode.path();
            SLIC_ERROR( message );
          }
          else
          {
            stringutilities::StringToType( xmlVal, defVal );
          }


        }
        localIndex const size = xmlVal.size();
        dataView.resize( size );
        cxx_utilities::equateStlVector(data,xmlVal);
      });


#else
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
          (this->getData<real64_array>(subDocNode.getStringKey())) = xmlVal;
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
          
          std::cout << "Resizing string_array " << subDocNode.getStringKey().c_str() << ".  Warning: this currently requires that the entire managed group to be resized." << std::endl;
          this->resize(xmlVal.size());
          ViewWrapper<string_array>::rtype saVal = this->getData<string_array>(subDocNode.getStringKey());
          for (uint jj=0; jj<xmlVal.size(); ++jj)
          {
            saVal[static_cast<int>(jj)] = xmlVal[jj];
          }

          break;
        }
        default:
        {
          // TODO: Define missing cases
          throw std::invalid_argument("XML auto read method not defined");
        }
      }
#endif
    }
  }

  ReadXML_PostProcess();
}


void ManagedGroup::PrintDataHierarchy()
{
  for( auto& view : this->m_wrappers )
  {
    std::cout<<view->getName()<<", "<<view->get_typeid().name()<<std::endl;
  }

  for( auto& group : this->m_subGroups )
  {
    std::cout<<group.first<<std::endl;
    group.second->PrintDataHierarchy();
  }
}


}
} /* namespace ODS */
