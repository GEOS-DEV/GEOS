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
 * @file PartitionBase.cpp
 */

#include "PartitionBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/NodeManager.hpp"

#include <limits.h>

namespace geosx
{

using namespace dataRepository;

PartitionBase::PartitionBase( ):
  m_domain( nullptr ),
  m_hasLocalGhosts( false )
{
  //maxComm
}


PartitionBase::PartitionBase( const unsigned int numPartitions, const unsigned int thisPartiton ):
  m_size( numPartitions ),
  m_rank( thisPartiton ),
  m_color( 0 ),
  m_numColors( 1 ),
  m_domain( nullptr ),
  m_t1( 0.0 ),
  m_t2( 0.0 ),
  m_t3( 0.0 ),
  m_t4( 0.0 ),
  m_hasLocalGhosts( false )
{
  //maxComm
}


PartitionBase::~PartitionBase()
{}

/**
 * @brief Call SetDomain on each neighbor
 */
void PartitionBase::SetDomain( DomainPartition * domain )
{
  // set the const pointer "m_domain" by casting away the the const
  // on the address of the pointer, and modifying what the address of
  // the pointer.
  DomainPartition** temp = const_cast<DomainPartition**>(&m_domain);
  *temp = domain;
}

}

