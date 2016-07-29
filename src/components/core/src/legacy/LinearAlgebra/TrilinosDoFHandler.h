#ifndef TRILINOS_DOF_HANDLER_H
#define TRILINOS_DOF_HANDLER_H

/**
 * @file TrilinosDoFHandler.h
 * @author white230
 */

#include "Common/Common.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "MPI_Communications/SpatialPartition.h"

/** \namespace Trilinos */

namespace Trilinos
{

  class DoFHandler
  {
    public:

      void initialize(PhysicalDomainT&  domain,
                      SpatialPartition& partition, 
                      unsigned dofs_per_node);

      unsigned get_n_global_dofs();
      unsigned get_n_local_dofs();
      unsigned get_first_local_dof();

    private:

      int this_mpi_process;
      int n_mpi_processes;

      struct
      {
        unsigned n_dofs_per_node;
        unsigned n_local_dofs;
        unsigned n_global_dofs;
        unsigned first_local_dof;
      } 
      sizing_info;
      
  };

}

#endif

