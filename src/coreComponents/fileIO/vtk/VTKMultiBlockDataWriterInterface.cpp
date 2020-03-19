/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "VTKMultiBlockDataWriterInterface.hpp"

#include "dataRepository/WrapperBase.hpp"
#include "dataRepository/Group.hpp"


#include "Array2DVTKDataArray.hpp"
#include "CellVTKDataArray.hpp"

#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>

#include <unordered_set>



namespace geosx
{
namespace vtk
{
  template< typename LAMBDA >
  static auto ApplyVTKTypeLambda( const rtTypes::TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
      case ( rtTypes::TypeIDs::integer_id ):
      {
        return lambda( integer( 1 ) );
      }
      case ( rtTypes::TypeIDs::real32_id ):
      {
        return lambda( real32( 1 ) );
      }
      case ( rtTypes::TypeIDs::real64_id ):
      {
        return lambda( real64( 1 ) );
      }
      case ( rtTypes::TypeIDs::localIndex_id ):
      {
        return lambda( localIndex( 1 ) );
      }
      case ( rtTypes::TypeIDs::globalIndex_id ):
      {
        return lambda( globalIndex( 1 ) );
      }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
      }
    }
  };
  VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ) :
    m_outputFolder( outputName ),
    m_pvd( outputName )
  {
    mode_t mode = 0733;
    mkdir( outputName.c_str(), mode );
  }

  vtkSmartPointer< vtkPoints >  VTKPolyDataWriterInterface::GetVTKPoints( NodeManager  const * const nodeManager ) const
  {
    vtkSmartPointer< vtkPoints > points = vtkPoints::New();
    points->SetNumberOfPoints( nodeManager->size() );
    auto connectivity = nodeManager->referencePosition();
    for( localIndex v = 0 ; v < nodeManager->size(); v++)
    {
      points->SetPoint( v, connectivity[v][0], connectivity[v][1], connectivity[v][2] );
    }
    return points;
  }

  std::tuple< vtkSmartPointer< vtkPoints >,  vtkSmartPointer< vtkCellArray > >VTKPolyDataWriterInterface::GetWell( WellElementSubRegion  const * const esr, NodeManager const * const nodeManager ) const
  {
    vtkSmartPointer< vtkPoints > points = vtkPoints::New();
    points->SetNumberOfPoints( esr->size() + 1 );
    vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
    cellsArray->SetNumberOfCells( esr->size() );
    for( localIndex edge = 0 ; edge < esr->size(); edge++)
    {
      localIndex firstPoint = esr->nodeList()[edge][0];
      auto point = nodeManager->referencePosition()[firstPoint];
      points->SetPoint( edge, point[0], point[1], point[2] );
      std::vector< vtkIdType > connectivity(esr->numNodesPerElement() );
      connectivity[0] = edge;
      connectivity[1] = edge+1;
      cellsArray->InsertNextCell( 2, connectivity.data() );
    }
    if( esr->size() > 0 )
    {
      localIndex lastPoint = esr->nodeList()[ esr->size() -1  ][1];
      auto point = nodeManager->referencePosition()[lastPoint];
      points->SetPoint( esr->size(), point[0], point[1], point[2] );
    }
    return std::make_tuple( points, cellsArray );
  }

  std::tuple< vtkSmartPointer< vtkCellArray >,  std::vector< int> > VTKPolyDataWriterInterface::GetVTKCells( CellElementRegion const * const er ) const
  {
    vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
    cellsArray->SetNumberOfCells( er->getNumberOfElements< CellElementRegion >() );
    std::vector< int > cellType;
    cellType.reserve( er->getNumberOfElements< CellElementRegion >() );
    er->forElementSubRegions<CellElementSubRegion>([&](CellElementSubRegion const * const esr )->void
    {
      for( localIndex c = 0; c < esr->size(); c++ )
      {
        std::vector< vtkIdType > connectivity(esr->numNodesPerElement() );
        if( esr->GetElementTypeString() == "C3D8" )
        {
          connectivity[0] = esr->nodeList(c,0);
          connectivity[1] = esr->nodeList(c,1);
          connectivity[2] = esr->nodeList(c,3);
          connectivity[3] = esr->nodeList(c,2);
          connectivity[4] = esr->nodeList(c,4);
          connectivity[5] = esr->nodeList(c,5);
          connectivity[6] = esr->nodeList(c,7);
          connectivity[7] = esr->nodeList(c,6);
        }
        cellType.push_back( VTK_HEXAHEDRON );
        cellsArray->InsertNextCell( esr->numNodesPerElement(), connectivity.data() );
      }
    });
    return std::make_tuple(cellsArray, cellType);
  }

  void VTKPolyDataWriterInterface::SetNodeFields( vtkSmartPointer< vtkPointData > const  pointdata, NodeManager const * const nodeManager) const
  {
    for( auto const & wrapperIter : nodeManager->wrappers() )
    {
      auto const * const wrapper = wrapperIter.second;
      if( wrapper->getPlotLevel() <= m_plotLevel )
      {
        vtkSmartPointer < vtkAOSDataArrayTemplate<real64 > > data =
          vtkAOSDataArrayTemplate< real64 >::New();
        data->SetNumberOfValues( nodeManager->size() );
        data->SetName( wrapper->getName().c_str() );;
        std::type_info const & typeID = wrapper->get_typeid();

        if( typeID==typeid(array1d<real64>) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array1d<real64>> const & >( *wrapper );
          arrayView1d< real64 const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 1 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            data->InsertValue(i, array[i]);
          }
        }
        if( typeID==typeid(array2d<real64>) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array2d<real64>> const & >( *wrapper );
          arrayView2d< real64 const > const & array = wrapperT.referenceAsView();
          integer dim = array.dims()[1];
          data->vtkAbstractArray::SetNumberOfComponents( dim );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            for( localIndex j = 0; j < dim; j++ )
            {
              data->InsertValue( dim * i + j ,array[i][j]);
            }
          }
        }
        if( typeID==typeid(array2d<real64,RAJA::PERM_JI>) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array2d<real64,RAJA::PERM_JI>> const & >( *wrapper );
          arrayView2d< real64 const, LvArray::getStrideOneDimension( RAJA::PERM_JI {} ) > const &
            array = wrapperT.referenceAsView();
          integer dim = array.dims()[1];
          data->vtkAbstractArray::SetNumberOfComponents( dim );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            for( localIndex j = 0; j < dim; j++ )
            {
              data->InsertValue( dim * i + j ,array[i][j]);
            }
          }
        }
        if( typeID==typeid(array3d<real64>) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array3d<real64>> const & >( *wrapper );
          arrayView3d< real64 const > const & array = wrapperT.referenceAsView();
          integer dim1 = array.dims()[1];
          integer dim2 = array.dims()[2];
          data->vtkAbstractArray::SetNumberOfComponents( dim1 * dim2 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            for( localIndex j = 0; j < dim1; j++ )
            {
              for( localIndex k = 0; k < dim2; k++ )
              {
                data->InsertValue( dim1*dim2 * i +dim2*j + k ,array[i][j][k]);
              }
            }
          }
        }
        if( typeID==typeid(r1_array) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<r1_array> const & >( *wrapper );
          arrayView1d< R1Tensor const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 3 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            for( localIndex j = 0; j < 3; j++ )
            {
              data->InsertValue( 3 * i + j ,array[i][j]);
            }
          }
        }
        if( typeID==typeid(integer_array) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<integer_array> const & >( *wrapper );
          arrayView1d< integer const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 1 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            data->InsertValue(i, static_cast< real64 > ( array[i] ) );
          }
        }
        if( typeID==typeid(localIndex_array) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<localIndex_array> const & >( *wrapper );
          arrayView1d< localIndex const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 1 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            data->InsertValue( i, static_cast< real64 > ( array[i] ) );
          }
        }
        if( typeID==typeid(globalIndex_array) )
        {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<globalIndex_array> const & >( *wrapper );
          arrayView1d< globalIndex const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 1 );
          for( localIndex i = 0; i < nodeManager->size(); i++ )
          {
            data->InsertValue(i, static_cast< real64 > ( array[i] ) );
          }
        }
        pointdata->AddArray( data );
      }
    }
  }

  void VTKPolyDataWriterInterface::SetCellFields( vtkSmartPointer< vtkCellData > const celldata, CellElementRegion const * const er ) const
  {
    std::unordered_set< string > allFields;
    er->forElementSubRegions<CellElementSubRegion>([&]( auto const * const esr )
    {
      for( auto const & wrapperIter : esr->wrappers() )
      {
        auto const * const wrapper = wrapperIter.second;
        if( wrapper->getPlotLevel() <= m_plotLevel )
        {
          allFields.insert(wrapperIter.first);
        }
      }
    });

    for( auto const & field : allFields )
    {
      vtkSmartPointer < vtkAOSDataArrayTemplate<real64 > > data =
        vtkAOSDataArrayTemplate< real64 >::New();
      data->SetNumberOfValues( er->getNumberOfElements<CellElementSubRegion> () );
      data->SetName( field.c_str() );

      localIndex count = 0;
      er->forElementSubRegions<CellElementSubRegion>([&]( auto const * const esr )
      {
          auto const * const wrapper = esr->getWrapperBase( field ) ;
          std::type_info const & typeID = wrapper->get_typeid();

          if( typeID==typeid(array1d<real64>) )
          {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array1d<real64>> const & >( *wrapper );
          arrayView1d< real64 const > const & array = wrapperT.referenceAsView();
          data->vtkAbstractArray::SetNumberOfComponents( 1 );
          for( localIndex i = 0; i < esr->size(); i++ )
          {
          data->InsertValue(count++, array[i]);
          }
          }
          if( typeID==typeid(array2d<real64>) )
          {
          auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array2d<real64>> const & >( *wrapper );
          arrayView2d< real64 const > const & array = wrapperT.referenceAsView();
          integer dim = array.dims()[1];
          data->vtkAbstractArray::SetNumberOfComponents( dim );
          for( localIndex i = 0; i < esr->size(); i++ )
          {
            for( localIndex j = 0; j < dim; j++ )
            {
              data->InsertValue( dim * count + j ,array[i][j]);
            }
            count++;
          }
          }
          if( typeID==typeid(array2d<real64,RAJA::PERM_JI>) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array2d<real64,RAJA::PERM_JI>> const & >( *wrapper );
            arrayView2d< real64 const, LvArray::getStrideOneDimension( RAJA::PERM_JI {} ) > const &
              array = wrapperT.referenceAsView();
            integer dim = array.dims()[1];
            data->vtkAbstractArray::SetNumberOfComponents( dim );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              for( localIndex j = 0; j < dim; j++ )
              {
                data->InsertValue( dim * count + j ,array[i][j]);
              }
              count++;
            }
          }
          if( typeID==typeid(array3d<real64>) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<array3d<real64>> const & >( *wrapper );
            arrayView3d< real64 const > const & array = wrapperT.referenceAsView();
            integer dim1 = array.dims()[1];
            integer dim2 = array.dims()[2];
            data->vtkAbstractArray::SetNumberOfComponents( dim1 * dim2 );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              for( localIndex j = 0; j < dim1; j++ )
              {
                for( localIndex k = 0; k < dim2; k++ )
                {
                  data->InsertValue( dim1*dim2 * count +dim2*j + k ,array[i][j][k]);
                }
              }
              count++;
            }
          }
          if( typeID==typeid(r1_array) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<r1_array> const & >( *wrapper );
            arrayView1d< R1Tensor const > const & array = wrapperT.referenceAsView();
            data->vtkAbstractArray::SetNumberOfComponents( 3 );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              for( localIndex j = 0; j < 3; j++ )
              {
                data->InsertValue( 3 * count + j ,array[i][j]);
              }
              count++;
            }
          }
          if( typeID==typeid(integer_array) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<integer_array> const & >( *wrapper );
            arrayView1d< integer const > const & array = wrapperT.referenceAsView();
            data->vtkAbstractArray::SetNumberOfComponents( 1 );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              data->InsertValue(count++, static_cast< real64 > ( array[i] ) );
            }
          }
          if( typeID==typeid(localIndex_array) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<localIndex_array> const & >( *wrapper );
            arrayView1d< localIndex const > const & array = wrapperT.referenceAsView();
            data->vtkAbstractArray::SetNumberOfComponents( 1 );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              data->InsertValue( count++, static_cast< real64 > ( array[i] ) );
            }
          }
          if( typeID==typeid(globalIndex_array) )
          {
            auto const & wrapperT = dynamic_cast< dataRepository::Wrapper<globalIndex_array> const & >( *wrapper );
            arrayView1d< globalIndex const > const & array = wrapperT.referenceAsView();
            data->vtkAbstractArray::SetNumberOfComponents( 1 );
            for( localIndex i = 0; i < esr->size(); i++ )
            {
              data->InsertValue(count++, static_cast< real64 > ( array[i] ) );
            }
          }
          });
      celldata->AddArray( data );
    }
  }
  void VTKPolyDataWriterInterface::WriteMeshFiles( double time, DomainPartition const * const domain ) const
  {
    string timeStepSubFolder = VTKPolyDataWriterInterface::GetTimeStepSubFolder( time );
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    elemManager->forElementRegions<CellElementRegion>([&](CellElementRegion const * const er)->void
        {
        vtkSmartPointer< vtkUnstructuredGrid > vtkUg = vtkUnstructuredGrid::New();
        vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtkUgWriter = vtkXMLUnstructuredGridWriter::New();
        auto VTKPoints = GetVTKPoints( domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
        vtkUg->SetPoints( VTKPoints );
        auto VTKCells = GetVTKCells( er );
        vtkUg->SetCells( std::get<1>(VTKCells).data(), std::get<0>(VTKCells) );
        SetCellFields( vtkUg->GetCellData(), er );
        SetNodeFields( vtkUg->GetPointData(),domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
        string vtuFilePath = timeStepSubFolder + "/" + std::to_string( mpiRank) +"_" + er->getName() + ".vtu";
        vtkUgWriter->SetFileName( vtuFilePath.c_str() );
        vtkUgWriter->SetInputData( vtkUg );
        vtkUgWriter->SetDataModeToAscii();
        vtkUgWriter->Write();
        });
  }

  void VTKPolyDataWriterInterface::WriteWellFiles( double time, DomainPartition const * const domain ) const
  {
    string timeStepSubFolder = VTKPolyDataWriterInterface::GetTimeStepSubFolder( time );
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    elemManager->forElementRegions<WellElementRegion>([&](WellElementRegion const * const er)->void
        {
        auto esr = er->GetSubRegion(0)->group_cast<WellElementSubRegion const *>();
        auto VTKWell = GetWell( esr, domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
        vtkSmartPointer<vtkUnstructuredGrid> ug = vtkUnstructuredGrid::New();
        ug->SetPoints(std::get<0>(VTKWell));
        ug->SetCells(VTK_LINE, std::get<1>(VTKWell));
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter =vtkXMLUnstructuredGridWriter::New();
        vtuWriter->SetInputData( ug );
        string vtuFilePath = timeStepSubFolder + "/" + std::to_string( mpiRank) +"_" + er->getName() + ".vtu";
        vtuWriter->SetFileName( vtuFilePath.c_str() );
        vtuWriter->SetDataModeToAscii();
        vtuWriter->Write();
        });
  }

  void VTKPolyDataWriterInterface::WriteVTMFile( double time, DomainPartition const * const domain, VTKVTMWriter const& vtmWriter ) const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    int const mpiSize = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
    if( mpiRank == 0 )
    {
      // Cells
      vtmWriter.AddBlock( CellElementRegion::CatalogName() );
      ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
      elemManager->forElementRegions<CellElementRegion>([&](CellElementRegion const * const er)->void
          {

          vtmWriter.AddSubBlock( CellElementRegion::CatalogName(), er->getName() );
          for(int i = 0; i < mpiSize; i++ )
          {
          vtmWriter.AddDataToSubBlock(CellElementRegion::CatalogName(), er->getName(), std::to_string(time) + "/" + std::to_string( i ) + "_" + er->getName() + ".vtu", i );
          }
          });

      // Wells
      vtmWriter.AddBlock( WellElementRegion::CatalogName() );
      elemManager->forElementRegions<WellElementRegion>([&](WellElementRegion const * const er)->void
          {
          vtmWriter.AddSubBlock( WellElementRegion::CatalogName(), er->getName() );
          for(int i = 0; i < mpiSize; i++ )
          {
          vtmWriter.AddDataToSubBlock( WellElementRegion::CatalogName(), er->getName(), std::to_string(time) + "/" + std::to_string( i ) +"_" + er->getName() + ".vtu", i );
          }
          });

      vtmWriter.Save();
    }
  }

  void VTKPolyDataWriterInterface::CreateTimeStepSubFolder( double time ) const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    if( mpiRank == 0 )
    {
      mode_t mode = 0773;
      string timeStepSubFolder = m_outputFolder + "/" + std::to_string( time );
      mkdir( timeStepSubFolder.c_str(), mode );
    }
    MpiWrapper::Barrier();
  }

  string VTKPolyDataWriterInterface::GetTimeStepSubFolder( double time ) const
  {
    return  m_outputFolder + "/" + std::to_string( time );
  }

  void VTKPolyDataWriterInterface::Write( real64 time, DomainPartition const * const domain  ) const
  {
    CreateTimeStepSubFolder( time );
    string vtmPath = m_outputFolder + "/" + std::to_string( time ) + ".vtm";
    VTKVTMWriter vtmWriter( vtmPath );
    WriteMeshFiles( time, domain);
    WriteWellFiles( time, domain);
    WriteVTMFile( time, domain, vtmWriter );
    m_pvd.AddData( time, vtmPath );
    m_pvd.Save();
  }
}
}
