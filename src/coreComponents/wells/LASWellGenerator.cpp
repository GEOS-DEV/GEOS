/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "LASWellGenerator.hpp"


namespace geosx
{
using namespace dataRepository;

LASWellGenerator::LASWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  registerWrapper(keys::fileName, &m_fileName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("Path to the las file");

  registerWrapper(keys::geometryLogIndexInFile, &m_logIndexToTakeForGeometry, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDefaultValue( -1 )->
    setSizedFromParent(0)->
    setDescription("Position of the log to take if there are several log sections defined in the LAS file ");
}

LASWellGenerator::~LASWellGenerator()
{
}

void LASWellGenerator::PostProcessInput()
{
}

Group * LASWellGenerator::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}

void LASWellGenerator::GenerateMesh( DomainPartition * const domain )
{
  LASFile lasFile;
  lasFile.Load(  m_fileName );
  if( lasFile.HasLog( "X" ) && lasFile.HasLog( "Y" ) )
  {
    GeneratePolyLineFromXYZ( lasFile );
  }
  else
  {
    GeneratePolyLineFromDepth( lasFile );
  }
}

void LASWellGenerator::GeneratePolyLineFromXYZ( LASFile const & lasFile )
{
  auto X = lasFile.GetLog( "X" );
  auto Y = lasFile.GetLog( "Y" );
  auto Z = lasFile.GetLog( "TVDSS" );
  localIndex nbLogEntries = lasFile.LogSize( "X" );
  GEOS_ASSERT( nbLogEntries == lasFile.LogSize( "Y" ) );
  GEOS_ASSERT( nbLogEntries == lasFile.LogSize( "TVDSS" ) );
  m_polyNodeCoords.resize( nbLogEntries );
  for( localIndex i = 0; i < nbLogEntries; i++)
  {
    m_polyNodeCoords[i] = { X[i], Y[i], Z[i] };
  }
}

void LASWellGenerator::GeneratePolyLineFromDepth( LASFile const & lasFile )
{
  auto topDepths = lasFile.GetLASLines< LASWellInformationSection >("ELEV");
  auto Xs = lasFile.GetLASLines< LASWellInformationSection >("XCOORD");
  auto Ys = lasFile.GetLASLines< LASWellInformationSection >("YCOORD");
  GEOS_ASSERT( topDepths.size() == Xs.size() );
  GEOS_ASSERT( topDepths.size() == Ys.size() );
  localIndex wellSectionIndex = 0;
  if( topDepths.size() > 1 && m_logIndexToTakeForGeometry == -1 )
  {
    GEOS_LOG_RANK_0("Warning : " << this->getName() << " corresponding LAS file has more than 1 log section defined "
                    << "please specify the index of the log section you want to take into account to write the well "
                    << "into the GEOSX data structure. You have to use the keyword " << keys::geometryLogIndexInFile
                    << " . Taking the first one by default.");
  }
  else if( topDepths.size() > 1 && m_logIndexToTakeForGeometry != -1 )
  {
    wellSectionIndex = m_logIndexToTakeForGeometry;
  }

  real64 topDepth = topDepths[wellSectionIndex]->GetDataAsReal64();
  real64 topX = Xs[wellSectionIndex]->GetDataAsReal64();
  real64 topY = Ys[wellSectionIndex]->GetDataAsReal64();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, LASWellGenerator, std::string const &, Group * const )

} // namespace
