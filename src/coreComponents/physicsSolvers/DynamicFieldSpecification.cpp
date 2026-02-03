

#include "DynamicFieldSpecification.hpp"

namespace geos
{
DynamicFieldSpecification::DynamicFieldSpecification( const string & name,
                                                      Group * const parent ):
  TaskBase( name, parent ),
  m_fieldSpecificationNames()
{

  registerWrapper( viewKeyStruct::fieldSpecificationNamesString(), &m_fieldSpecificationNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Array containing the field specifications to apply" );

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

    fs.apply< dataRepository::Group >( mesh,
                                       [&]( FieldSpecificationBase const & bc,
                                            string const &,
                                            SortedArrayView< localIndex const > const & targetObject,
                                            Group & targetGroup,
                                            string const fieldName )
      {
        bc.applyFieldValue< FieldSpecificationEqual >( targetObject, 0.0, targetGroup, fieldName );
      } );
  }


  return false;
}


REGISTER_CATALOG_ENTRY( TaskBase, DynamicFieldSpecification, string const &, Group * const )


} // namespace geos