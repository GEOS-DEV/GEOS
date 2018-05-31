// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#include "BinStream.h"

// ****************************************************************************
// ***** BINSTREAM CLASS DEFINITIONS ******************************************
// ****************************************************************************
BinStream::BinStream(void)
{}

BinStream::~BinStream(void)
{}

// ****************************************************************************
// ***** oBINSTREAM CLASS DEFINITIONS *****************************************
// ****************************************************************************
oBinStream::oBinStream(void)
{}

oBinStream::~oBinStream(void)
{}

void oBinStream::open( const char* filename,  const bool truncate )
{
  if( truncate )
  {
    output.open( filename, std::ios::binary | std::ios::out | std::ios::trunc );
  }
  else
  {
    output.open( filename, std::ios::binary | std::ios::out );
  }
}

void oBinStream::close( void )
{ output.close(); }



// ****************************************************************************
// ***** iBINSTREAM CLASS DEFINITIONS *****************************************
// ****************************************************************************
iBinStream::iBinStream(void)
{}

iBinStream::~iBinStream(void)
{}

void iBinStream::open( const char* filename, const bool )
{ input.open( filename, std::ios::binary | std::ios::in ); }

void iBinStream::close( void )
{ input.close(); }



namespace BinStreamUtils
{}
