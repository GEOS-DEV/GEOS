/**
 * @file MesquireRelaxer.cpp
 * @author herbold, settgast
 */
#include "MesquiteRelaxer.h"
#include <iostream>
#include <sstream>


using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ostringstream;

#include <memory>
using std::auto_ptr;

#include <ctype.h>

// Libraries for sorting
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
//#include <iostream>

using namespace MESQUITE_NS;

typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;

void MesquiteRelaxer::RelaxMesh( PhysicalDomainT& domain )
{


  // TODO: EBH: Have the user select the solver type for the mesh relaxation.
  // For now, just hard code something that works
  //            most of the time.
  SolverTypes solver = SolverTypes::SMART_LAPLACE;
  int num_iter = 10;
  PatchType patch = PatchType::PATCH_GLOBAL;


  array<R1Tensor>& Xref = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
  const array<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  const array<integer>& isNodeExternal = domain.m_feNodeManager.m_isExternal;
  R1Tensor nodePosition;
  int xyzIdx, nodeIdx;
  ElementTypes topo;

  unsigned long int numNodes = domain.m_feNodeManager.m_numNodes;
  array<real64> coords(3*numNodes);
  // Interleave x/y/z coordinates into "coords"
  for (localIndex a = 0 ; a < numNodes ; ++a)
  {
    xyzIdx = 3*a;
    coords[xyzIdx  ] = Xref[a][0]+disp[a][0];
    coords[xyzIdx+1] = Xref[a][1]+disp[a][1];
    coords[xyzIdx+2] = Xref[a][2]+disp[a][2];
  }

  RegionMap::iterator
    iter_region     = domain.m_feElementManager.m_ElementRegions.begin(),
    end_region = domain.m_feElementManager.m_ElementRegions.end();

  for( ; iter_region != end_region ; ++iter_region)
  {
    ElementRegionT& region = iter_region->second;
    // Find out some basic information abou the element region mesh
    unsigned long int numElems = region.m_numElems;//domain.m_feElementManager.m_numElems;

    const unsigned numNodesPerElement = region.m_numNodesPerElem;
    lArray1d elemToNode_serialized(numNodesPerElement*numElems);

    // Figure out the element type for this element region
    int switchOnElement = numNodesPerElement;
    if (numNodesPerElement == 4 && region.m_ElementDimension == 3)
    {
      switchOnElement = 5;
    }
    switch (switchOnElement)
    {
    case 3:
      topo = ElementTypes::TRIS;
      break;
    case 4:
      topo = ElementTypes::QUADS;
      break;
    case 5:
      topo = ElementTypes::TETS;
      break;
    case 8:
      topo = ElementTypes::HEXES;
      break;
    default:
      cout << "There is no default case";
      throw GPException("There is no default case");
      break;
    }

    const FixedOneToManyRelation& elementToNodeMap = region.m_toNodesRelation;
    const array<integer> nodeOrdering = region.SiloNodeOrdering();

    // Interleave connectivity into array of integers such that for tri elements
    // you have:
    //     quads = [ elem1_node1, elem1_node2, elem1_node3, elem2_node1,
    // elem2_node2, elem2_node3, ...] and so on.
    nodeIdx = 0;
    for (localIndex a = 0 ; a < numElems ; ++a)
    {
      for (localIndex b = 0 ; b < numNodesPerElement ; ++b)
      {
        nodeIdx = a*numNodesPerElement + b;
        elemToNode_serialized[nodeIdx] = elementToNodeMap(a,nodeOrdering(b));
      }
    }

    MesquiteRelaxer::ConstructMesh(numNodes,
                                   coords.data(),
                                   isNodeExternal.data(),
                                   numElems,
                                   elemToNode_serialized.data(),
                                   topo);

    switch (solver)
    {
    case SolverTypes::SMOOTH_OFF:
      cout << "Why are you here?  No smoothing to be done...";
      return;
    case SolverTypes::LAPLACE:
      MesquiteRelaxer::RunLaplace(topo,num_iter);
      break;
    case SolverTypes::SMART_LAPLACE:
      MesquiteRelaxer::RunSmartLaplace(topo,num_iter);
      break;
    case SolverTypes::CONJUGATE_GRADIENT:
      MesquiteRelaxer::RunConjugateGradient(topo,num_iter,patch);
      break;
    case SolverTypes::FEASIBLE_NEWTON:
      MesquiteRelaxer::RunFeasibleNewton(topo,num_iter,patch);
      break;
    case SolverTypes::STEEPEST_DESCENT:
      MesquiteRelaxer::RunSteepestDescent(topo,num_iter,patch);
      break;
    case SolverTypes::BOUNDARY:
      cout << "Not Implemented yet";
      return;
    default:
      cout << "There is no default case";
      break;
    }
  }

  // Overwrite the reference nodal positions based on the movement from the
  // Mesquite library
  for (localIndex a = 0 ; a < numNodes ; ++a)
  {
    xyzIdx = 3*a;
    Xref[a][0] = coords[xyzIdx  ] - disp[a][0];
    Xref[a][1] = coords[xyzIdx+1] - disp[a][1];
    Xref[a][2] = coords[xyzIdx+2] - disp[a][2];
  }

  //Timer t;
  //double secs = t.since_birth();
  //std::cout << "Optimization completed in " << secs << " seconds" <<
  // std::endl;

  //std::cout << "New vertex location: ( " << endl;
  //for (int j=0; j<*num_nod; j++){
  //   cout << "                  ";
  //   for (int k=0; k<*dimen; k++){
  //          value[k] = coords[3*j+k];
  //          cout << value[k] << " ";
  //   }
  //   cout << endl;
  //}

// End of Main Routine
}

void MesquiteRelaxer::ConstructMesh( unsigned long num_nod,
                                     double* coords,
                                     const int* fixed,
                                     unsigned long num_elem,
                                     const unsigned long* quads,
                                     const ElementTypes topo)
{
  EntityTopology type;
  switch (topo)
  {
  case ElementTypes::TRIS:
    type=TRIANGLE;
    break;
  case ElementTypes::QUADS:
    type=QUADRILATERAL;
    break;
  case ElementTypes::TETS:
    type=TETRAHEDRON;
    break;
  case ElementTypes::HEXES:
    type=HEXAHEDRON;
    break;
  default:
    throw GPException("invalid topo");
  }
  ;

  mesh.set_mesh(3,num_nod,coords,fixed,num_elem,type,quads);
}

void MesquiteRelaxer::RunLaplace(const ElementTypes topo, const int num_iter)
{
  MsqError err;
  InstructionQueue queue1;
  PlanarDomain domain( PlanarDomain::XY );
  QualityAssessor assessor;
  TerminationCriterion outer, inner;

  // Metrics
  ConditionNumberQualityMetric shape_metric;
  EdgeLengthQualityMetric edge_metric;
  edge_metric.set_averaging_method(QualityMetric::RMS);

  //Objective Function

  //Solver
  LaplacianSmoother solver;
  if (err)
    cout << err << endl;
  //   Solver: set parameters

  //   Solver: set termination
  outer.add_iteration_limit(num_iter);
  solver.set_outer_termination_criterion(&outer);

  //Quality Assessment
  assessor.add_quality_assessment(&shape_metric);
  assessor.add_quality_assessment(&edge_metric);

  //Add quality and assessor schemes to queue
  queue1.add_quality_assessor(&assessor,err);
  if (err)
    cout << err << endl;
  queue1.set_master_quality_improver(&solver,err);
  if (err)
    cout << err << endl;

  switch (topo)
  {
  case ElementTypes::TRIS:
  case ElementTypes::QUADS:
  {
//     queue1.run_instructions( &mesh, &domain, err);
    //cout << "Im using tris";
    break;
  }
  case ElementTypes::TETS:
  case ElementTypes::HEXES:
  {
    queue1.run_instructions( &mesh, err);
    //cout << "Im using tets";
    break;
  }
  }
  if (err)
    cout << err << endl;

  queue1.clear();
  mesh.release();
}

void MesquiteRelaxer::RunSmartLaplace(const ElementTypes topo, const int num_iter)
{
  MsqError err;
  InstructionQueue queue1;
  PlanarDomain domain( PlanarDomain::XY );
  QualityAssessor assessor;
  TerminationCriterion outer, inner;

  // Metrics
  IdealWeightInverseMeanRatio inverse_metric;
  inverse_metric.set_averaging_method(QualityMetric::SUM, err);

  //Objective Function
  //LPtoPTemplate objective_function(inverse_metric, 1, err);
  LInfTemplate objective_function(&inverse_metric);

  //Solver
  SmartLaplacianSmoother solver( &objective_function );
  if (err)
    cout << err << endl;

  //   Solver: set parameters

  //   Solver: set termination
  outer.add_iteration_limit( num_iter );
  solver.set_outer_termination_criterion( &outer );

  //Quality Assessment
  assessor.add_quality_assessment( &inverse_metric );
  if (err)
    cout << err << endl;

  //Add quality and assessor schemes to queue
  queue1.set_master_quality_improver( &solver, err );
  if (err)
    cout << err << endl;
  queue1.add_quality_assessor( &assessor, err );
  if (err)
    cout << err << endl;

  switch (topo)
  {
  case ElementTypes::TRIS:
  case ElementTypes::QUADS:
  {
//     queue1.run_instructions( &mesh, &domain, err);
    //cout << "Im using tris";
    break;
  }
  case ElementTypes::TETS:
  case ElementTypes::HEXES:
  {
    queue1.run_instructions( &mesh, err);
    //cout << "Im using tets";
    break;
  }
  }
  if (err)
    cout << err << endl;

  queue1.clear();
  mesh.release();
}


void MesquiteRelaxer::RunConjugateGradient(const ElementTypes topo,
                                           const int num_iter,
                                           const PatchType patch)
{
  MsqError err;
  InstructionQueue queue1;
  PlanarDomain domain( PlanarDomain::XY );
  QualityAssessor assessor;
  TerminationCriterion outer, inner;

  // Metrics
  TShapeB1 target_metric;
  IdealShapeTarget W;
  TQualityMetric mu( &W, &target_metric );
  ConditionNumberQualityMetric shape_metric;

  //Objective Function
  PMeanPTemplate objective_function( 2.0, &mu);

  //Solver
  ConjugateGradient solver( &objective_function );
  if (err)
    cout << err << endl;

  //   Solver: set parameters
  switch (patch)
  {
  // Global Optimizer
  case PatchType::PATCH_GLOBAL:
  {
    solver.use_global_patch();
    outer.add_iteration_limit( 1 );
    inner.add_relative_successive_improvement(OF_value);
    break;
  }
  //Nash Game
  case PatchType::PATCH_NASH:
  {
    solver.use_element_on_vertex_patch();
    outer.add_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(num_iter);
    break;
  }
  // Block Coordinate Descent
  case PatchType::PATCH_BLOCK:
  {
    solver.use_element_on_vertex_patch();
    solver.do_block_coordinate_descent_optimization();
    inner.add_iteration_limit(num_iter);
    outer.add_relative_quality_improvement(OF_value);
    break;
  }

  // Culling
  case PatchType::PATCH_CULL:
  {
    solver.use_element_on_vertex_patch();
    inner.cull_on_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(2);
    break;
  }

  // Jacobi
  case PatchType::PATCH_JACOBI:
  {
    solver.use_element_on_vertex_patch();
    solver.do_jacobi_optimization();
    inner.add_iteration_limit( 2 );
    outer.add_absolute_vertex_movement(OF_value);
    break;
  }
  }

  //   Solver: set termination
  solver.set_inner_termination_criterion( &inner );
  solver.set_outer_termination_criterion( &outer );

  //Quality Assessment
  assessor.add_quality_assessment( &mu, 10 );
  assessor.add_quality_assessment( &shape_metric );

  //Add quality and assessor schemes to queue
  queue1.set_master_quality_improver( &solver, err );
  if (err)
    cout << err << endl;
  queue1.add_quality_assessor( &assessor, err );
  if (err)
    cout << err << endl;

  switch (topo)
  {
  case ElementTypes::TRIS:
  case ElementTypes::QUADS:
  {
    //    queue1.run_instructions( &mesh, &domain, err);
    //cout << "Im using tris";
    break;
  }
  case ElementTypes::TETS:
  case ElementTypes::HEXES:
  {
    queue1.run_instructions( &mesh, err);
    //cout << "Im using tets";
    break;
  }
  }
  if (err)
    cout << err << endl;

  queue1.clear();
  mesh.release();
}

void MesquiteRelaxer::RunFeasibleNewton(const ElementTypes topo,
                                        const int num_iter,
                                        const PatchType patch)
{
  MsqError err;
  InstructionQueue queue1;
  PlanarDomain domain( PlanarDomain::XY );
  QualityAssessor assessor;
  TerminationCriterion outer, inner;

  //Metric
  TShapeB1 target_metric;
  IdealShapeTarget W;
  TQualityMetric mu( &W, &target_metric );
  ConditionNumberQualityMetric shape_metric;

  //Objective Function
  PMeanPTemplate objective_function( 2.0, &mu);

  //Solver
  FeasibleNewton solver( &objective_function );
  if (err)
    cout << err << endl;

  //   Solver: set parameters
  switch (patch)
  {
  // Global Optimizer
  case PatchType::PATCH_GLOBAL:
  {
    solver.use_global_patch();
    outer.add_iteration_limit( 1 );
    inner.add_relative_successive_improvement(OF_value);
    break;
  }
  //Nash Game
  case PatchType::PATCH_NASH:
  {
    solver.use_element_on_vertex_patch();
    outer.add_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(num_iter);
    break;
  }
  // Block Coordinate Descent
  case PatchType::PATCH_BLOCK:
  {
    solver.use_element_on_vertex_patch();
    solver.do_block_coordinate_descent_optimization();
    inner.add_iteration_limit(num_iter);
    outer.add_relative_quality_improvement(OF_value);
    break;
  }

  // Culling
  case PatchType::PATCH_CULL:
  {
    solver.use_element_on_vertex_patch();
    inner.cull_on_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(2);
    break;
  }

  // Jacobi
  case PatchType::PATCH_JACOBI:
  {
    solver.use_element_on_vertex_patch();
    solver.do_jacobi_optimization();
    inner.add_iteration_limit( 2 );
    outer.add_absolute_vertex_movement(OF_value);
    break;
  }
  }

  //   Solver: set termination
  solver.set_inner_termination_criterion( &inner );
  solver.set_outer_termination_criterion( &outer );


  //Quality Assessment
  assessor.add_quality_assessment( &mu, 10 );
  assessor.add_quality_assessment( &shape_metric );

  //Add quality and assessor schemes to queue
  queue1.set_master_quality_improver( &solver, err );
  if (err)
    cout << err << endl;
  queue1.add_quality_assessor( &assessor, err );
  if (err)
    cout << err << endl;

  switch (topo)
  {
  case ElementTypes::TRIS:
  case ElementTypes::QUADS:
  {
    //    queue1.run_instructions( &mesh, &domain, err);
    //cout << "Im using tris";
    break;
  }
  case ElementTypes::TETS:
  case ElementTypes::HEXES:
  {
    queue1.run_instructions( &mesh, err);
    //cout << "Im using tets";
    break;
  }
  }
  if (err)
    cout << err << endl;

  queue1.clear();
  mesh.release();

}

void MesquiteRelaxer::RunSteepestDescent(const ElementTypes topo,
                                         const int num_iter,
                                         const PatchType patch)
{
  MsqError err;
  InstructionQueue queue1;
  PlanarDomain domain( PlanarDomain::XY );
  QualityAssessor assessor;
  TerminationCriterion outer, inner;

  //Metric
  TShapeB1 target_metric;
  IdealShapeTarget W;
  TQualityMetric mu( &W, &target_metric );
  ConditionNumberQualityMetric shape_metric;

  //Objective Function
  PMeanPTemplate objective_function( 2.0, &mu);

  //Solver
  SteepestDescent solver( &objective_function );
  if (err)
    cout << err << endl;

  //   Solver: set parameters
  switch (patch)
  {
  // Global Optimizer
  case PatchType::PATCH_GLOBAL:
  {
    solver.use_global_patch();
    outer.add_iteration_limit( 1 );
    inner.add_relative_successive_improvement(OF_value);
    break;
  }
  //Nash Game
  case PatchType::PATCH_NASH:
  {
    solver.use_element_on_vertex_patch();
    outer.add_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(num_iter);
    break;
  }
  // Block Coordinate Descent
  case PatchType::PATCH_BLOCK:
  {
    solver.use_element_on_vertex_patch();
    solver.do_block_coordinate_descent_optimization();
    inner.add_iteration_limit(num_iter);
    outer.add_relative_quality_improvement(OF_value);
    break;
  }

  // Culling
  case PatchType::PATCH_CULL:
  {
    solver.use_element_on_vertex_patch();
    inner.cull_on_absolute_vertex_movement(OF_value);
    inner.add_iteration_limit(2);
    break;
  }

  // Jacobi
  case PatchType::PATCH_JACOBI:
  {
    solver.use_element_on_vertex_patch();
    solver.do_jacobi_optimization();
    inner.add_iteration_limit( 2 );
    outer.add_absolute_vertex_movement(OF_value);
    break;
  }
  }

  //   Solver: set termination
  solver.set_inner_termination_criterion( &inner );
  solver.set_outer_termination_criterion( &outer );

  //Quality Assessment
  assessor.add_quality_assessment( &mu, 10 );
  assessor.add_quality_assessment( &shape_metric );

  //Add quality and assessor schemes to queue
  queue1.set_master_quality_improver( &solver, err );
  if (err)
    cout << err << endl;
  queue1.add_quality_assessor( &assessor, err );
  if (err)
    cout << err << endl;

  switch (topo)
  {
  case ElementTypes::TRIS:
  case ElementTypes::QUADS:
  {
//     queue1.run_instructions( &mesh, &domain, err);
    //cout << "Im using tris";
    break;
  }
  case ElementTypes::TETS:
  case ElementTypes::HEXES:
  {
    queue1.run_instructions( &mesh, err);
    //cout << "Im using tets";
    break;
  }
  }
  if (err)
    cout << err << endl;

  queue1.clear();
  mesh.release();
}
