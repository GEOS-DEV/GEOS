/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/DomainPartition.hpp"

using namespace geosx;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;


void RegisterAndApplyField( DomainPartition * domain,
                            string const & fieldName,
                            string const & objectPath,
                            real64 value )
{
  auto fieldSpecificationManager = FieldSpecificationManager::get();

  auto fieldSpec = fieldSpecificationManager->RegisterGroup<FieldSpecificationBase>(fieldName);
  fieldSpec->SetFieldName(fieldName);
  fieldSpec->SetObjectPath(objectPath);
  fieldSpec->SetScale(value);
  fieldSpec->InitialCondition(true);
  fieldSpec->AddSetName("all");

  fieldSpecificationManager->Apply( 0., domain, "", "",
                                    [&] ( FieldSpecificationBase const * const bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * const targetGroup,
                                          string const fieldNamee )
                                    {
                                      bc->ApplyFieldValue<FieldSpecificationEqual>( targetSet, 0.0, targetGroup, fieldNamee );
                                    });
}

TEST(FieldSpecification, Recursive)
{

  // Mesh Definitions
  localIndex nbTetReg0 = 30;
  localIndex nbHexReg0 = 60;
  localIndex nbTetReg1 = 40;
  localIndex nbHexReg1 = 50;
  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", nullptr ) );
  auto meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( "body" );
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>(std::string("Level0"));
  NodeManager * nodeManager = meshLevel0->getNodeManager();

  CellBlockManager * cellBlockManager = domain->GetGroup<CellBlockManager>( keys::cellManager );

  CellBlock * reg0Hex = cellBlockManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("reg0hex");
  reg0Hex->SetElementType("C3D8");
  reg0Hex->resize(nbHexReg0);
  auto & cellToVertexreg0Hex = reg0Hex->nodeList();
  cellToVertexreg0Hex.resize( nbHexReg0, 8);

  CellBlock * reg0Tet= cellBlockManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("reg0tet");
  reg0Tet->SetElementType("C3D4");
  reg0Tet->resize(nbTetReg0);
  auto & cellToVertexreg0Tet = reg0Tet->nodeList();
  cellToVertexreg0Tet.resize( nbTetReg0, 4);

  CellBlock * reg1Hex = cellBlockManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("reg1hex");
  reg1Hex->SetElementType("C3D8");
  reg1Hex->resize(nbHexReg1);
  auto & cellToVertexreg1Hex = reg1Hex->nodeList();
  cellToVertexreg1Hex.resize( nbHexReg1, 8);

  CellBlock * reg1Tet= cellBlockManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("reg1tet");
  reg1Tet->SetElementType("C3D4");
  reg1Tet->resize(nbTetReg1);
  auto & cellToVertexreg1Tet = reg1Tet->nodeList();
  cellToVertexreg1Tet.resize( nbTetReg1, 4);

  ElementRegionManager * elemManager = meshLevel0->getElemManager();
  ElementRegion* reg0 = elemManager->CreateChild("ElementRegion","reg0")->group_cast<ElementRegion*>();
  reg0->AddCellBlockName(reg0Hex->getName());
  reg0->AddCellBlockName(reg0Tet->getName());
  ElementRegion* reg1 = elemManager->CreateChild("ElementRegion","reg1")->group_cast<ElementRegion*>();
  reg1->AddCellBlockName(reg1Hex->getName());
  reg1->AddCellBlockName(reg1Tet->getName());
  reg0->GenerateMesh(cellBlockManager->GetGroup(keys::cellBlocks) );
  reg1->GenerateMesh(cellBlockManager->GetGroup(keys::cellBlocks) );


  /// Field Definition
  auto fieldSpecificationManager = FieldSpecificationManager::get();

  reg0->GetSubRegion("reg0hex")->RegisterViewWrapper< array1d<real64> >( "field0" );
  reg0->GetSubRegion("reg0tet")->RegisterViewWrapper< array1d<real64> >( "field0" );
  reg1->GetSubRegion("reg1tet")->RegisterViewWrapper< array1d<real64> >( "field0" );
  reg1->GetSubRegion("reg1hex")->RegisterViewWrapper< array1d<real64> >( "field0" );

  reg0->GetSubRegion("reg0hex")->RegisterViewWrapper< array1d<real64> >( "field1" );
  reg0->GetSubRegion("reg0tet")->RegisterViewWrapper< array1d<real64> >( "field1" );

  reg0->GetSubRegion("reg0hex")->RegisterViewWrapper< array1d<real64> >( "field2" );

  reg1->GetSubRegion("reg1tet")->RegisterViewWrapper< array1d<real64> >( "field3" );

  auto set0hex = reg0->GetSubRegion("reg0hex")->GetGroup("sets")->RegisterViewWrapper<localIndex_set>( std::string("all") );
  set0hex->resize(nbHexReg0);
  for(localIndex i = 0; i < set0hex->size() ; i++)
  {
    set0hex->dataPtr()[i]=i;
  }
  auto set0tet = reg0->GetSubRegion("reg0tet")->GetGroup("sets")->RegisterViewWrapper<localIndex_set>( std::string("all") );
  set0tet->resize(nbTetReg0);
  for(localIndex i = 0; i < set0tet->size() ; i++)
  {
    set0tet->dataPtr()[i] = i;
  }
  auto set1hex = reg1->GetSubRegion("reg1hex")->GetGroup("sets")->RegisterViewWrapper<localIndex_set>( std::string("all") );
  set1hex->resize(nbHexReg1);
  for(localIndex i = 0; i < set1hex->size() ; i++)
  {
    set1hex->dataPtr()[i] = i;
  }
  auto set1tet = reg1->GetSubRegion("reg1tet")->GetGroup("sets")->RegisterViewWrapper<localIndex_set>( std::string("all") );
  set1tet->resize(nbTetReg1);
  for(localIndex i = 0; i < set1tet->size() ; i++)
  {
    set1tet->dataPtr()[i] = i;
  }

  RegisterAndApplyField(domain.get(), "field0", "ElementRegions", 1.);
  RegisterAndApplyField(domain.get(), "field1", "ElementRegions/reg0", 2.);
  RegisterAndApplyField(domain.get(), "field2", "ElementRegions/reg0/elementSubRegions/reg0hex", 3.);
  RegisterAndApplyField(domain.get(), "field3", "ElementRegions/reg1/elementSubRegions/reg1tet", 4.);

  /// Check if the values are well set

  auto field0 = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "field0" );
  auto field1 = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "field1" );
  auto field2 = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "field2" );
  auto field3 = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "field3" );
  elemManager->forElementSubRegionsComplete( [&] ( localIndex er, localIndex esr,
                                                   ElementRegion * const region,
                                                   ElementSubRegionBase const * const subRegion )
  {
    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      GEOS_ERROR_IF(field0[er][esr][ei] < 1. || field0[er][esr][ei] > 1., "Recursive fields are not set");
    });
  });

  reg0->forElementSubRegionsIndex( [&](localIndex const esr, ElementSubRegionBase * const subRegion )->void
  {
    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      GEOS_ERROR_IF(field1[0][esr][ei] < 2. || field1[0][esr][ei] > 2., "Recursive fields are not set");
    });
  });

  forall_in_range<elemPolicy>( 0, reg0Hex->size(), GEOSX_LAMBDA ( localIndex ei )
  {
    GEOS_ERROR_IF(field2[0][0][ei] < 1. || field2[0][0][ei] > 3., "Recursive fields are not set");
  });

  forall_in_range<elemPolicy>( 0, reg1Tet->size(), GEOSX_LAMBDA ( localIndex ei )
  {
    GEOS_ERROR_IF(field3[1][1][ei] < 4. || field3[1][1][ei] > 4., "Recursive fields are not set");
  });


}


int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI

  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
