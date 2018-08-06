/**
 * @file PoroelasticSolver.cpp
 *
 */

#include "PoroelasticSolver.hpp"
#include "../FiniteVolume/SinglePhaseFlow_TPFA.hpp"
#include "managers/DomainPartition.hpp"


namespace geosx
{

PoroelasticSolver::PoroelasticSolver( const std::string& name,
                                      ManagedGroup * const parent ):
  SolverBase(name,parent)
{
  this->RegisterViewWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0);
  this->RegisterViewWrapper(viewKeyStruct::fluidSolverNameString, &m_flowSolverName, 0);
}


void PoroelasticSolver::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example single phase flow solver");


  docNode->AllocateChildNode( viewKeyStruct::fluidSolverNameString,
                              viewKeyStruct::fluidSolverNameString,
                              -1,
                              "string",
                              "string",
                              "name of fluid solver",
                              "name of fluid solver",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::solidSolverNameString,
                              viewKeyStruct::solidSolverNameString,
                              -1,
                              "string",
                              "string",
                              "name of solid solver",
                              "name of solid solver",
                              "",
                              "",
                              0,
                              1,
                              0 );
}

void PoroelasticSolver::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  SolverBase::FillOtherDocumentationNodes( rootGroup );
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0);

    {


    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks( [&]( CellBlockSubRegion * const cellBlock ) -> void
      {
        cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

        docNode->AllocateChildNode( viewKeyStruct::effectiveStressString,
                                    viewKeyStruct::effectiveStressString,
                                    -1,
                                    "r2Sym_array",
                                    "r2Sym_array",
                                    "Effective Stress",
                                    "Effective Stress",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );
      });
    }
  }
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}


real64 PoroelasticSolver::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition * domain )
{
  return SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );

}

real64 PoroelasticSolver::SplitOperatorStep( real64 const& time_n,
                                             real64 const& dt,
                                             integer const cycleNumber,
                                             DomainPartition * const domain)
{

  SolverBase &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow_TPFA &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow_TPFA*>());

  SystemSolverParameters & solverParams = *(this->getSystemSolverParameters());

  fluidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  solidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  int iter = 0;
  bool balanced = false;
  while (iter < solverParams.maxIterNewton() )
  {
    fluidSolver.NonlinearImplicitStep( time_n,
                                       dt,
                                       cycleNumber,
                                       domain,
                                       getLinearSystemRepository() );
    balanced = fluidSolver.getSystemSolverParameters()->numNewtonIterations();
    if (balanced && iter > 0)
    {
      break;
    }

    solidSolver.NonlinearImplicitStep( time_n,
                                       dt,
                                       cycleNumber,
                                       domain,
                                       getLinearSystemRepository() );

    //fluidSolver.UpdateDeformation(domain);
    ++iter;
  }
  fluidSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
//  if (partition.m_rank == 0) std::cout << ". Finished in " << iter << " iterations." << std::endl;
  return dt;
}


REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolver, std::string const &, ManagedGroup * const )

} /* namespace geosx */
