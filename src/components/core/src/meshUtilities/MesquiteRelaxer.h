/**
 * @file MesquireRelaxer.h
 * @author herbold, settgast
 */
#ifndef MESQUITE_RELAXER_H_
#define MESQUITE_RELAXER_H_

#include "Mesquite_all_headers.hpp"
#include "ObjectManagers/PhysicalDomainT.h"
#include <iostream>
#include <sstream>

class MesquiteRelaxer
{
public:

  enum class ElementTypes
  {
    QUADS = 9,
    HEXES = 12,
    TETS  = 11,
    TRIS  = 8
  };

  enum class SolverTypes
  {
    SMOOTH_OFF         = 1, // do nothing if smooting is turned off
    LAPLACE            = 2,
    SMART_LAPLACE      = 3,
    CONJUGATE_GRADIENT = 4,
    FEASIBLE_NEWTON    = 5,
    STEEPEST_DESCENT   = 6,
    BOUNDARY           = 7
  };

  enum class PatchType
  {
    PATCH_GLOBAL       = 1,
    PATCH_NASH         = 2,
    PATCH_BLOCK        = 3,
    PATCH_CULL         = 4,
    PATCH_JACOBI       = 5
  };

  void RelaxMesh( PhysicalDomainT& domain );

  /*
     void RelaxMesh(int *dimen,
                 unsigned long *num_nod,
                 double *coords,
                 const int *fixed,
                 int *num_elem,
                 const unsigned long int *quads,
                 int *num_iter,
                 const ElementTypes topo,
                 const SolverTypes solver,
                 const PatchType patch);
   */

  // EBH: you may want to declare a default constructor/destructor for
  // MesquiteRelaxer
  //MesquiteRelaxer();
  //~MesquiteRelaxer();

private:

  double OF_value = 1.e-5;

  MESQUITE_NS::ArrayMesh mesh;

  void RunLaplace(const ElementTypes topo, const int num_iter);
  void RunSmartLaplace(const ElementTypes topo, const int num_iter);
  void RunConjugateGradient(const ElementTypes topo,
                            const int num_iter,
                            const PatchType patch);
  void RunFeasibleNewton(const ElementTypes topo,
                         const int num_iter,
                         const PatchType patch);
  void RunSteepestDescent(const ElementTypes topo,
                          const int num_iter,
                          const PatchType patch);
  void ConstructMesh(unsigned long num_nod,
                     double *coords,
                     const int *fixed,
                     unsigned long num_elem,
                     const unsigned long *quads,
                     const ElementTypes topo);
};

#endif
