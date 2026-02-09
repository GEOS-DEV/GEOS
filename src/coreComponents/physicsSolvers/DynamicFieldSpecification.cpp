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
  GEOS_UNUSED_VAR( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

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
              string const targetFieldName = getTargetFieldName( fieldName );
              bc.applyFieldValue< FieldSpecificationEqual >( targetObject, 0.0, targetGroup, targetFieldName );
            } );
          }
        }
      }
    }
  }


  return false;
}

string DynamicFieldSpecification::getTargetFieldName( string const & fieldName ) const
{
  // HydrostaticEquilibrium is defined programatically a field specification, but it isn’t actually a field itself.
  // Instead, it updates the pressure field, which is the target field being modified.
  if( fieldName == "HydrostaticEquilibrium" )
  {
    return "pressure";
  }
  return fieldName;
}


REGISTER_CATALOG_ENTRY( TaskBase, DynamicFieldSpecification, string const &, Group * const )


} // namespace geos
