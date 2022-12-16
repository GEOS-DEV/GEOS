/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "RESQMLWriterInterface.hpp"

#include "common/Logger.hpp"
#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "fileIO/Outputs/OutputUtilities.hpp"

#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

#include "fesapi/resqml2_0_1/PropertyKind.h"
#include "fesapi/eml2/AbstractHdfProxy.h"
#include "fesapi/eml2/TimeSeries.h"
#include "fesapi/tools/TimeTools.h"

// System includes
#include <unordered_set>
#include <random>
#include <sstream>


namespace geosx
{


/// @brief Non standard UUID generation function
/// placeholder to avoid linking to boost
namespace uuid {
    static std::random_device              rd;
    static std::mt19937                    gen(rd());
    static std::uniform_int_distribution<> dis(0, 15);
    static std::uniform_int_distribution<> dis2(8, 11);

    string generate_uuid_v4() {
        std::stringstream ss;
        int i;
        ss << std::hex;
        for (i = 0; i < 8; i++) {
            ss << dis(gen);
        }
        ss << "-";
        for (i = 0; i < 4; i++) {
            ss << dis(gen);
        }
        ss << "-4";
        for (i = 0; i < 3; i++) {
            ss << dis(gen);
        }
        ss << "-";
        ss << dis2(gen);
        for (i = 0; i < 3; i++) {
            ss << dis(gen);
        }
        ss << "-";
        for (i = 0; i < 12; i++) {
            ss << dis(gen);
        };
        return ss.str();
    }
}

using namespace dataRepository;

/**
 * @brief Build/expand a list of strings used as default component labels on demand.
 * @param size number of labels requested
 * @return a span over range of strings (stored permanently in memory)
 */
static Span< string const >
getDefaultLabels( localIndex const size )
{
  static std::vector< string > labels;
  localIndex oldSize = LvArray::integerConversion< localIndex >( labels.size() );
  std::generate_n( std::back_inserter( labels ), size - oldSize, [&] { return std::to_string( oldSize++ ); } );
  return { labels.begin(), labels.begin() + size };
}


/**
 * @brief Checks consistency of user-provided per-dimension labels (must match the size of array).
 * @param wrapper the array wrapper
 * @param dim dimension index to check
 */
template< typename T, int NDIM, typename PERM >
void checkLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper, int const dim )
{
  GEOSX_ERROR_IF_NE_MSG( LvArray::integerConversion< localIndex >( wrapper.getDimLabels( dim ).size() ),
                         wrapper.reference().size( dim ),
                         "VTK writer: component names are set, but don't match the array size.\n"
                         "This is likely a bug in physics module (solver or constitutive model)." );
}


/**
 * @brief Get a list of component labels for a particular dimension of an array.
 * @param wrapper array wrapper
 * @param dim dimension index
 * @return a span over range of strings representing labels
 */
template< typename T, int NDIM, typename PERM >
static Span< string const >
getDimLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
              int const dim )
{
  Span< string const > const labels = wrapper.getDimLabels( dim );
  if( labels.empty() )
  {
    return getDefaultLabels( wrapper.reference().size( dim ) );
  }
  checkLabels( wrapper, dim );
  return labels;
}

/**
 * @brief Build a multidimensional component name out of distinct dimension-wise labels.
 * @tparam Ts types of indices
 * @tparam Is dummy template argument required for positional expansion of the pack
 * @param dimLabels per-dimension component labels
 * @param indices per-dimension component indices
 * @return combined component name
 */
template< typename ... Ts, integer ... Is >
static string
makeComponentName( std::vector< string >(&dimLabels)[sizeof...( Ts )],
                   std::integer_sequence< integer, Is... >,
                   Ts const & ... indices )
{
  return stringutilities::concat( '/', dimLabels[Is][indices] ... );
}


/**
 * @brief Specialized component metadata handler for 1D arrays.
 * @param data VTK typed data array
 */
template< typename T, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, 1, PERM > > const &,
                      vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( 1 );
}

/**
 * @brief Specialized component metadata handler for 2D arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 *
 * This exists because we want to keep default VTK handling for unlabeled components
 * (i.e. X/Y/Z for 1-3 components, numeric indices for higher) for the time being.
 * This function can be removed if we force each physics package to always set its labels.
 */
template< typename T, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, 2, PERM > > const & wrapper,
                      vtkAOSDataArrayTemplate< T > * data )
{
  auto const view = wrapper.referenceAsView();
  data->SetNumberOfComponents( view.size( 1 ) );

  Span< string const > const labels = wrapper.getDimLabels( 1 );
  if( !labels.empty() )
  {
    checkLabels( wrapper, 1 );
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      data->SetComponentName( i, labels[i].c_str() );
    }
  }
}

/**
 * @brief Produces a temporary array slice from a view that can be looped over.
 * @param view the source view
 * @return a fake slice that does not point to real data but has correct dims/strides.
 * @note The slice is only valid as long as the @p view is in scope.
 *       Values in the slice may be uninitialized and should not be used.
 */
template< typename T, int NDIM, int USD >
static ArraySlice< T const, NDIM - 1, USD - 1 >
makeTemporarySlice( ArrayView< T const, NDIM, USD > const & view )
{
  // The following works in all compilers, but technically invokes undefined behavior:
  // return ArraySlice< T, NDIM - 1, USD - 1 >( nullptr, view.dims() + 1, view.strides() + 1 );
  localIndex const numComp = LvArray::indexing::multiplyAll< NDIM - 1 >( view.dims() + 1 );
  static array1d< T > arr;
  arr.template resizeWithoutInitializationOrDestruction( numComp );
  return ArraySlice< T const, NDIM - 1, USD - 1 >( arr.data(), view.dims() + 1, view.strides() + 1 );
}

/**
 * @brief Generic component metadata handler for multidimensional arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 */
template< typename T, int NDIM, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
                      vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( wrapper.numArrayComp() );

  std::vector< string > labels[NDIM-1];
  for( integer dim = 1; dim < NDIM; ++dim )
  {
    Span< string const > dimLabels = getDimLabels( wrapper, dim );
    labels[dim-1].assign( dimLabels.begin(), dimLabels.end() );
  }

  auto const view = wrapper.referenceAsView();
  auto const slice = view.size( 0 ) > 0 ? view[0] : makeTemporarySlice( view );

  integer compIndex = 0;
  LvArray::forValuesInSliceWithIndices( slice, [&]( T const &, auto const ... indices )
  {
    using idx_seq = std::make_integer_sequence< integer, sizeof...(indices) >;
    data->SetComponentName( compIndex++, makeComponentName( labels, idx_seq{}, indices ... ).c_str() );
  } );
}

static void
writeField( WrapperBase const & wrapper,
            localIndex const offset,
            vtkDataArray * data )
{
  types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > * typedData = vtkAOSDataArrayTemplate< T >::FastDownCast( data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    forAll< parallelHostPolicy >( sourceArray.size( 0 ), [sourceArray, offset, typedData]( localIndex const i )
    {
      LvArray::forValuesInSliceWithIndices( sourceArray[i], [&, compIndex = 0]( T const & value, auto const ... indices  ) mutable
      {
        typedData->SetTypedComponent( offset + i, compIndex++, value );
        std::cout <<  LvArray::indexing::getIndexString(indices ...) << std::endl;

      } );
    } );
  } );
}

void
RESQMLWriterInterface::writeElementField( Group const & subRegions,
                   string const & field,
                   COMMON_NS::DataObjectRepository* outputRepository,
                   RESQML2_0_1_NS::UnstructuredGridRepresentation* parent,
                   EML2_NS::TimeSeries * timeSeries,
                   time_t timestamp) 
{
  // instantiate vtk array of the correct type
  vtkSmartPointer< vtkDataArray > data;
  localIndex numElements = 0;
  bool first = true;
  int numDims = 0;
  subRegions.forSubGroups< Group >( [&]( Group const & subRegion )
  {

    numElements += subRegion.size();
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    if( first )
    {
      types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
      {
        using ArrayType = decltype( array );
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data.TakeReference( typedData );
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), typedData );
      } );
      first = false;
      numDims = wrapper.numArrayDims();

      RESQML2_0_1_NS::ContinuousProperty* contProp1;
      auto search = m_property_uuid.find(field);
      if(search != m_property_uuid.end())
      {
        contProp1 = search->second;
      }
      else
      {
        string property = uuid::generate_uuid_v4();

        GEOSX_LOG_RANK_VAR( field );  
        GEOSX_LOG_RANK_VAR( property );

        MpiWrapper::broadcast( property, 0 );
          
        contProp1 = outputRepository->createContinuousProperty(parent, property, field, data->GetNumberOfComponents(),
                gsoap_eml2_3::resqml22__IndexableElement::cells, gsoap_resqml2_0_1::resqml20__ResqmlUom::m, gsoap_resqml2_0_1::resqml20__ResqmlPropertyKind::length);
      
        contProp1->setTimeSeries(timeSeries);
        contProp1->setSingleTimestamp(timestamp);

        if(numDims == 1)
          contProp1->pushBackHdf5Array1dOfValues(COMMON_NS::AbstractObject::numericalDatatypeEnum::DOUBLE, numElements);
        else if (numDims == 2)
          contProp1->pushBackHdf5Array2dOfValues(COMMON_NS::AbstractObject::numericalDatatypeEnum::DOUBLE, numElements, numElements);
        else if (numDims == 3)
          contProp1->pushBackHdf5Array3dOfValues(COMMON_NS::AbstractObject::numericalDatatypeEnum::DOUBLE, numElements, numElements, numElements);

        m_property_uuid[field] = contProp1;  
      }



    }
    else
    {
      // Sanity check
      GEOSX_ERROR_IF_NE_MSG( wrapper.numArrayDims(), numDims,
                             "VTK writer: sanity check failed for " << field << " (inconsistent array dimensions)" );
      GEOSX_ERROR_IF_NE_MSG( wrapper.numArrayComp(), data->GetNumberOfComponents(),
                             "VTK writer: sanity check failed for " << field << " (inconsistent array sizes)" );
    }
  } );

  data->SetNumberOfTuples( numElements );
  data->SetName( field.c_str() );

  // std::cout << field << std::endl;
  // std::cout << numElements << std::endl;


  // write each subregion in turn, keeping track of element offset
  localIndex offset = 0;
  subRegions.forSubGroups< Group >( [&]( Group const & subRegion )
  {
    ElementSubRegionBase const & subr = dynamic_cast<ElementSubRegionBase const &>(subRegion);
    

        localIndex const numElems = subr.size();
        //arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( wrapperIter.first );
        arrayView1d< globalIndex const > const & localToGlobal = subr.localToGlobalMap();

        std::cout << field <<  " numElements: " << numElems << std::endl;
        for( localIndex k = 0; k < numElems; ++k )
        {
          std::cout <<  localToGlobal[k] << std::endl;
        }
        std::cout << std::endl;


    globalIndex startIdx = localToGlobal[0];
    localIndex valueCount = subr.size();



    std::cout << "offset: " << offset << std::endl;
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    writeField( wrapper, offset, data.GetPointer() );    
    offset += subRegion.size();    

    if(numDims == 1)
      m_property_uuid[field]->setValuesOfDoubleHdf5Array1dOfValues(static_cast<double*>(data->GetVoidPointer(0)), valueCount, startIdx);
  
  } );



  std::cout << MpiWrapper::commRank() << std::endl;
  
  std::cout << "total offset: " << offset << std::endl;

    
  // std::cout <<" plop " << std::endl;
  // data->Print(std::cout);
  // std::cout << std::endl;
  

  std::cout << std::endl << std::endl;
  
}



RESQMLWriterInterface::RESQMLWriterInterface( string name ):
    VTKPolyDataWriterInterface(name)
  , m_outputRepository(new COMMON_NS::DataObjectRepository())
  , m_propertyKind( nullptr )
{
  // if( MpiWrapper::commRank() == 0) 
  string propertyKind = uuid::generate_uuid_v4();  
  MpiWrapper::broadcast( propertyKind, 0 );  
  
  m_propertyKind = m_outputRepository->createPropertyKind(propertyKind, "propType1", "F2I", gsoap_resqml2_0_1::resqml20__ResqmlUom::m, gsoap_resqml2_0_1::resqml20__ResqmlPropertyKind::length);
  
  string timeSeries = uuid::generate_uuid_v4();
  MpiWrapper::broadcast( timeSeries, 0 );
    
  m_timeSeries = m_outputRepository->createTimeSeries(timeSeries, "Testing time series");
  
  //Create default MPI Proxy
  string hdfProxy = uuid::generate_uuid_v4();
  MpiWrapper::broadcast( hdfProxy, 0 );
    
  EML2_NS::AbstractHdfProxy* m_hdfProxy = m_outputRepository->createHdfProxy(hdfProxy, "Parallel Hdf Proxy", m_outputDir, m_outputName + ".h5", COMMON_NS::DataObjectRepository::openingMode::OVERWRITE);
	m_outputRepository->setDefaultHdfProxy(m_hdfProxy);

    //TODO need local3dCrs ?
    // local3dCrs = repo.createLocalDepth3dCrs("", "Default local CRS", .0, .0, .0, .0, gsoap_resqml2_0_1::eml20__LengthUom::m, 23031, gsoap_resqml2_0_1::eml20__LengthUom::m, "Unknown", false);
    // repo.setDefaultCrs(local3dCrs);
}

void RESQMLWriterInterface::generateOutput()
{
  string outputFilename = joinPath( m_outputDir, m_outputName ) + ".epc" ;
  GEOSX_LOG_RANK_0( GEOSX_FMT( "Creating: {}", outputFilename) );
  COMMON_NS::EpcDocument pck(outputFilename);
  GEOSX_LOG_RANK_0( GEOSX_FMT("Start serialization of {} in {}", pck.getName(), pck.getStorageDirectory().empty() ? "working directory." : pck.getStorageDirectory()));
  pck.serializeFrom(*m_outputRepository);
}

void RESQMLWriterInterface::write( real64 const time,
                                        integer const GEOSX_UNUSED_PARAM( cycle ),
                                        DomainPartition const & domain )
{
  m_property_uuid.clear();
  auto as_duration = std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<real64>(time));
  // duration seconds from start of epoch
  std::chrono::system_clock::time_point time_point(as_duration);    
  time_t timestamp = std::chrono::system_clock::to_time_t(time_point);
  
  m_timeSeries->pushBackTimestamp(timestamp);

  // loop over all mesh levels and mesh bodies
  domain.forMeshBodies( [&]( MeshBody const & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
    {


  //     if( meshLevel.isShallowCopy() )
  //     {
  //       return;
  //     }

      ElementRegionManager const & elemManager = meshLevel.getElemManager();
      NodeManager const & nodeManager = meshLevel.getNodeManager();
      //string const & meshLevelName = meshLevel.getName();
      //string const & meshBodyName = meshBody.getName();


  //     if( m_outputRegionType == VTKRegionTypes::CELL || m_outputRegionType == VTKRegionTypes::ALL )
  //     {
      writeCellElementRegions( time, elemManager, nodeManager, timestamp );
  //     }
    } );
  } );


}

void RESQMLWriterInterface::writeCellElementRegions( real64 const time,
                                                     ElementRegionManager const & elemManager,
                                                     NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                    time_t timestamp ) 
{
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & region )
  {
    
    // CellData VTKCells = getVtkCells( region, nodeManager.size() );
    // vtkSmartPointer< vtkPoints > const VTKPoints = getVtkPoints( nodeManager, VTKCells.nodes );

    // auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    // ug->SetCells( VTKCells.cellTypes.data(), VTKCells.cells );
    // ug->SetPoints( VTKPoints );

    // // writeTimestamp( ug.GetPointer(), time );
    writeElementFields( region, time, timestamp );
    // writeNodeFields( nodeManager, VTKCells.nodes, ug->GetPointData() );

  } );
}


void RESQMLWriterInterface::writeElementFields( ElementRegionBase const & region,
                                                real64 const GEOSX_UNUSED_PARAM( time ),
                                                time_t timestamp ) 
{
  // Collect a list of regular fields (filter out material field wrappers)
  // TODO: this can be removed if we stop hanging constitutive wrappers on the mesh
  std::unordered_set< string > regularFields;
  region.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
  {
    for( auto const & wrapperIter : subRegion.wrappers() )
    {
      if( isFieldPlotEnabled( *wrapperIter.second ))
      {
        regularFields.insert( wrapperIter.first );


        // localIndex const numElems = subRegion.size();
        // //arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( wrapperIter.first );
        // arrayView1d< globalIndex const > const & localToGlobal = subRegion.localToGlobalMap();

        // std::cout << wrapperIter.first <<  " numElements: " << numElems << std::endl;
        // for( localIndex k = 0; k < numElems; ++k )
        // {
        //   std::cout <<  localToGlobal[k] << std::endl;
        // }
        // std::cout << std::endl;

      }
    }
  } );

  // Write regular fields
  for( string const & field : regularFields )
  {  
    //arrayView1d<const globalIndex >  localToGlobal = region.localToGlobalMap();


    // for( std::ptrdiff_t i = 0; i < localToGlobal.size( 0 ); ++i )
    // {
    //     std::cout << localToGlobal[i] << " " ;
    // }
    // std::cout << std::endl;

    std::cout << " region has field " << field << " ? " << region.hasSubRegion(field) << std::endl;

    writeElementField( region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() ), field, m_outputRepository, m_parent, m_timeSeries, timestamp);
  }
}


}