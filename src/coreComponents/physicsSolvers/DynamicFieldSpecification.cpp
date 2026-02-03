#include "DynamicFieldSpecification.hpp"
#include "events/tasks/TasksManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{
using namespace dataRepository;

DynamicFieldSpecification::DynamicFieldSpecification( const string & name,
                                                      Group * const parent ):
  TaskBase( name, parent ),
  m_fieldSpecificationNames()
{

  registerWrapper( viewKeyStruct::fieldSpecificationNamesString(), &m_fieldSpecificationNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Array containing the field specifications to apply" );

}

DynamicFieldSpecification::~DynamicFieldSpecification() = default;

void DynamicFieldSpecification::postInputInitialization()
{
  TaskBase::postInputInitialization();
}


bool
DynamicFieldSpecification::
execute( real64 const time_n,
         real64 const dt,
         integer const cycleNumber,
         integer const eventCounter,
         real64 const eventProgress,
         DomainPartition & domain )
{

  FieldSpecificationManager & fsm = FieldSpecificationManager::getInstance();

  for( string const & fsName : m_fieldSpecificationNames )
  {
    GEOS_THROW_IF( !fsm.hasGroup( fsName ),
                   GEOS_FMT( "{}: FieldSpecification named {} not found",
                             getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ),
                             fsName ),
                   InputError, getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ) );

    FieldSpecificationBase const & fs = fsm.getGroup< FieldSpecificationBase >( fsName );

    for( auto & meshBodyPair : domain.getMeshBodies().getSubGroups() )
    {
      if( MeshBody * const meshBody = dynamic_cast< MeshBody * >( meshBodyPair.second ) )
      {
        for( auto & meshLevelPair : meshBody->getMeshLevels().getSubGroups() )
        {
          if( MeshLevel * const meshLevel = dynamic_cast< MeshLevel * >( meshLevelPair.second ) )
          {
            fs.apply< dataRepository::Group >( *meshLevel,
                                               [&]( FieldSpecificationBase const & bc,
                                                    string const &,
                                                    SortedArrayView< localIndex const > const & targetObject,
                                                    Group & targetGroup,
                                                    string const fieldName )
              {
                bc.applyFieldValue< FieldSpecificationEqual >( targetObject, 0.0, targetGroup, fieldName );
              } );
          }
        }
      }
    }
  }


  return false;
}


REGISTER_CATALOG_ENTRY( TaskBase, DynamicFieldSpecification, string const &, Group * const )


} // namespace geos
