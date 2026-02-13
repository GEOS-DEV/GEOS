#include "DynamicFieldSpecification.hpp"
#include "events/tasks/TasksManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/ElementSubRegionBase.hpp"

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
  GEOS_UNUSED_VAR( time_n, dt, cycleNumber, eventCounter, eventProgress );

  FieldSpecificationManager & fsm = FieldSpecificationManager::getInstance();
  FunctionManager & functionManager = FunctionManager::getInstance();

  for( string const & fsName : m_fieldSpecificationNames )
  {
    GEOS_THROW_IF( !fsm.hasGroup( fsName ),
                   GEOS_FMT( "{}: FieldSpecification named {} not found",
                             getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ),
                             fsName ),
                   InputError, getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ) );

    FieldSpecificationBase const & fs = fsm.getGroup< FieldSpecificationBase >( fsName );

    // Check if this is an EquilibriumInitialCondition (HydrostaticEquilibrium)
    EquilibriumInitialCondition const * equilIC = dynamic_cast< EquilibriumInitialCondition const * >( &fs );

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
                                                    SortedArrayView< localIndex const > const & targetSet,
                                                    Group & targetGroup,
                                                    string const fieldName )
            {
              // For HydrostaticEquilibrium, we need special handling since:
              // 1. Pressure requires complex hydrostatic computation with phase-dependent interpolation
              // 2. Temperature and globalCompFraction can be applied directly from elevation tables
              if( equilIC != nullptr )
              {
                applyEquilibriumInitialConditionFields( *equilIC, targetSet, targetGroup, functionManager );
              }
              else
              {
                // For regular field specifications, apply the field value normally
                string const targetFieldName = getTargetFieldName( fieldName );
                bc.applyFieldValue< FieldSpecificationEqual >( targetSet, 0.0, targetGroup, targetFieldName );
              }
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
  // HydrostaticEquilibrium is defined programatically a field specification, but it isn't actually a field itself.
  // Instead, it updates the pressure field, which is the target field being modified.
  if( fieldName == "HydrostaticEquilibrium" )
  {
    return "pressure";
  }
  return fieldName;
}


void DynamicFieldSpecification::applyEquilibriumInitialConditionFields( EquilibriumInitialCondition const & equilIC,
                                                                        SortedArrayView< localIndex const > const & targetSet,
                                                                        dataRepository::Group & targetGroup,
                                                                        FunctionManager & functionManager )
{
  // For HydrostaticEquilibrium, we need to set pressure, temperature, and globalCompFraction fields
  // using the elevation-based tables.

  // Get element center coordinates for elevation lookup
  if( !targetGroup.hasWrapper( ElementSubRegionBase::viewKeyStruct::elementCenterString() ) )
  {
    // If the target group doesn't have element centers, we cannot apply elevation-based fields
    return;
  }

  arrayView2d< real64 const > const elemCenter =
    targetGroup.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

  // Apply pressure field
  // For single-phase flow, a pressure table {equilIC.getName()}_{subRegion.getName()}_table is created
  // For compositional flow, no such table exists, so we use hydrostatic pressure calculation
  if( targetGroup.hasWrapper( "pressure" ) )
  {
    string const presTableName = equilIC.getName() + "_" + targetGroup.getName() + "_table";
    if( functionManager.hasGroup< TableFunction >( presTableName ) )
    {
      // Single-phase flow: use the pre-computed pressure table
      TableFunction const & presTable = functionManager.getGroup< TableFunction >( presTableName );
      TableFunction::KernelWrapper presTableWrapper = presTable.createKernelWrapper();

      arrayView1d< real64 > const pres = targetGroup.getReference< array1d< real64 > >( "pressure" );

      forAll< parallelHostPolicy >( targetSet.size(), [=] ( localIndex const i )
      {
        localIndex const k = targetSet[i];
        real64 const elevation = elemCenter[k][2];
        pres[k] = presTableWrapper.compute( &elevation );
      } );
    }
    else
    {
      // Compositional flow: compute hydrostatic pressure from datum
      // P(z) = P_datum - rho * g * (z - z_datum)
      // Note: This is a simplified approach. For more accurate results with compositional flow,
      // the full computeHydrostaticEquilibrium logic would be needed.
      real64 const datumPressure = equilIC.getDatumPressure();
      real64 const datumElevation = equilIC.getDatumElevation();

      arrayView1d< real64 > const pres = targetGroup.getReference< array1d< real64 > >( "pressure" );

      // Use a characteristic density for hydrostatic pressure calculation
      // For a more accurate approach, this would need access to the fluid model
      constexpr real64 defaultDensity = 1000.0;  // kg/m^3, water-like density
      constexpr real64 defaultGravity = 9.81;    // m/s^2

      forAll< parallelHostPolicy >( targetSet.size(), [=] ( localIndex const i )
      {
        localIndex const k = targetSet[i];
        real64 const elevation = elemCenter[k][2];
        // Hydrostatic pressure: P = P_datum + rho * g * (z_datum - z)
        pres[k] = datumPressure + defaultDensity * defaultGravity * ( datumElevation - elevation );
      } );
    }
  }

  // Apply temperature field if temperature table is specified
  string const & tempTableName = equilIC.getTemperatureVsElevationTableName();
  if( !tempTableName.empty() && targetGroup.hasWrapper( "temperature" ) )
  {
    TableFunction const & tempTable = functionManager.getGroup< TableFunction >( tempTableName );
    TableFunction::KernelWrapper tempTableWrapper = tempTable.createKernelWrapper();

    arrayView1d< real64 > const temp = targetGroup.getReference< array1d< real64 > >( "temperature" );

    forAll< parallelHostPolicy >( targetSet.size(), [=] ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elemCenter[k][2];
      temp[k] = tempTableWrapper.compute( &elevation );
    } );
  }

  // Apply global component fractions if component fraction tables are specified
  string_array const & compFracTableNames = equilIC.getComponentFractionVsElevationTableNames();
  integer const numComps = compFracTableNames.size();

  if( numComps > 0 && targetGroup.hasWrapper( "globalCompFraction" ) )
  {
    arrayView2d< real64 > const compFrac =
      targetGroup.getReference< array2d< real64 > >( "globalCompFraction" );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( compFracTableNames[ic] );
      TableFunction::KernelWrapper compFracTableWrapper = compFracTable.createKernelWrapper();

      forAll< parallelHostPolicy >( targetSet.size(), [=] ( localIndex const i )
      {
        localIndex const k = targetSet[i];
        real64 const elevation = elemCenter[k][2];
        compFrac[k][ic] = compFracTableWrapper.compute( &elevation );
      } );
    }
  }
}


REGISTER_CATALOG_ENTRY( TaskBase, DynamicFieldSpecification, string const &, Group * const )


} // namespace geos
