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

/**
 * @file VTKOutput.cpp
 */

#include "VTKOutput.hpp"
#include "fileIO/vtk/VTKFile.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

VTKOutput::VTKOutput( std::string const & name,
                        ManagedGroup * const parent ):
  OutputBase( name, parent),
  m_plotFileRoot(),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_vtkFile(name)
{
  RegisterViewWrapper(viewKeysStruct::plotFileRoot, &m_plotFileRoot, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  RegisterViewWrapper(viewKeysStruct::writeFEMFaces, &m_writeFaceMesh, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  RegisterViewWrapper(viewKeysStruct::plotLevel, &m_plotLevel, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  RegisterViewWrapper(viewKeysStruct::binaryString, &m_writeBinaryData, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Output the data in binary format");

  m_vtkFile.SetPlotLevel( m_plotLevel );
  m_vtkFile.SetBinaryMode( m_writeBinaryData) ;

}

VTKOutput::~VTKOutput()
{}



void VTKOutput::Execute(real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const eventProgress,
                         ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  m_vtkFile.Write( time_n, *domainPartition);
}


REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
