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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "LinearAlgebra/TrilinosDoFHandler.h"

namespace Trilinos
{

void DoFHandler :: initialize(PhysicalDomainT & domain,
                              SpatialPartition& partition,
                              unsigned dofs_per_node)
{
  MPI_Comm_size(MPI_COMM_WORLD,&n_mpi_processes);
  MPI_Comm_rank(MPI_COMM_WORLD,&this_mpi_process);

  sizing_info.n_dofs_per_node = dofs_per_node;
  sizing_info.n_local_dofs    = dofs_per_node*domain.m_feNodeManager.DataLengths();

  std::vector<int> scatter(n_mpi_processes,0);
  scatter[this_mpi_process] = sizing_info.n_local_dofs;

  std::vector<int> gather(n_mpi_processes);

  MPI_Allreduce(&scatter.front(),
                &gather.front(),
                n_mpi_processes,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD);

  sizing_info.first_local_dof = 0;
  sizing_info.n_global_dofs = 0;

  for(int p=0 ; p<n_mpi_processes ; ++p)
  {
    sizing_info.n_global_dofs += gather[p];
    if(p<this_mpi_process)
      sizing_info.first_local_dof += gather[p];
  }
}


unsigned DoFHandler :: get_n_global_dofs   () { return sizing_info.n_global_dofs; }
unsigned DoFHandler :: get_n_local_dofs    () { return sizing_info.n_local_dofs; }
unsigned DoFHandler :: get_first_local_dof () { return sizing_info.first_local_dof; }
}
