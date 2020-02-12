/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior
 * University Copyright (c) 2018-2019 Total, S.A Copyright (c) 2019-     GEOSX
 * Contributors All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS
 * files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldDamageFEM.cpp
 *
 */

#include "PhaseFieldDamageFEM.hpp"

#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"

#include "managers/DomainPartition.hpp"

// this should be part of the input file

double myFunc2(double , double , double ) {
  return 0;
  // return pow(x, 2) + pow(y, 2) + pow(z, 2) + 6;
//  return x * (1 - x) * y * (1 - y) * z * (1 - z) -
//         2 * (x - 1) * x * (y - 1) * y - 2 * (x - 1) * x * (z - 1) * z -
//         2 * (y - 1) * y * (z - 1) * z;
}

///////////////////////////////////

namespace geosx {

namespace dataRepository {
namespace keys {}
} // namespace dataRepository

using namespace dataRepository;
using namespace constitutive;

PhaseFieldDamageFEM::PhaseFieldDamageFEM(const std::string &name,
                                           Group *const parent):
  SolverBase(name, parent),
  m_fieldName("primaryField"),
  m_solidModelName()
{

  registerWrapper<string>(PhaseFieldDamageFEMViewKeys.timeIntegrationOption.Key())->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("option for default time integration method");

  registerWrapper<string>(PhaseFieldDamageFEMViewKeys.fieldVarName.Key(), &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("name of field variable");

  registerWrapper(viewKeyStruct::localDissipationOption, &m_localDissipationOption, false)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Type of local dissipation function. Can be Linear or Quadratic");
//  registerWrapper(viewKeyStruct::coeffFieldName,
//                          &m_coeffFieldName, false)->
//    setInputFlag(InputFlags::REQUIRED)->
//    setDescription("name of field variable representing the diffusion coefficient");



}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM() {
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::RegisterDataOnMesh(Group *const MeshBodies) {
  for (auto &mesh : MeshBodies->GetSubGroups()) {

    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    NodeManager *const nodes = meshLevel->getNodeManager();

    nodes->registerWrapper<real64_array>(m_fieldName)
        ->setApplyDefaultValue(0.0)
        ->setPlotLevel(PlotLevel::LEVEL_0)
        ->setDescription("Primary field variable");

     ElementRegionManager * const elemManager = meshLevel->getElemManager();

     elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
     {
        subRegion->registerWrapper( viewKeyStruct::coeffName, &m_coeff, false )->
          setApplyDefaultValue(0.0)->
          setPlotLevel(PlotLevel::LEVEL_0)->
          setDescription("field variable representing the diffusion coefficient");
     });

  }
}

void PhaseFieldDamageFEM::PostProcessInput() {
  SolverBase::PostProcessInput();

  string tiOption = this->getReference<string>(
      PhaseFieldDamageFEMViewKeys.timeIntegrationOption);

  if (tiOption == "SteadyState") {
    this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
  } else if (tiOption == "ImplicitTransient") {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
  } else if (tiOption == "ExplicitTransient") {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
  } else {
    GEOS_ERROR("invalid time integration option");
  }

  if(m_localDissipationOption != "Linear" && m_localDissipationOption != "Quadratic"){
    GEOS_ERROR("invalid local dissipation option - must be Linear or Quadratic");
  }

  // Set basic parameters for solver
  m_linearSolverParameters.logLevel = 0;
  m_linearSolverParameters.solverType = "gmres";
  m_linearSolverParameters.krylov.tolerance = 1e-8;
  m_linearSolverParameters.krylov.maxIterations = 250;
  m_linearSolverParameters.krylov.maxRestart = 250;
  m_linearSolverParameters.preconditionerType = "amg";
  m_linearSolverParameters.amg.smootherType = "gaussSeidel";
  m_linearSolverParameters.amg.coarseType = "direct";
}

real64 PhaseFieldDamageFEM::SolverStep(real64 const &time_n, real64 const &dt,
                                        const int cycleNumber,
                                        DomainPartition *domain) {
  real64 dtReturn = dt;
  if (m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient) {
    dtReturn = ExplicitStep(time_n, dt, cycleNumber, domain);
  } else if (m_timeIntegrationOption ==
                 timeIntegrationOption::ImplicitTransient ||
             m_timeIntegrationOption == timeIntegrationOption::SteadyState) {
    dtReturn =
        this->LinearImplicitStep(time_n, dt, cycleNumber, domain, m_dofManager,
                                 m_matrix, m_rhs, m_solution);
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::ExplicitStep(
    real64 const &GEOSX_UNUSED_ARG(time_n), real64 const &dt,
    const int GEOSX_UNUSED_ARG(cycleNumber),
    DomainPartition *const GEOSX_UNUSED_ARG(domain)) {
  return dt;
}

void PhaseFieldDamageFEM::ImplicitStepSetup(
    real64 const &GEOSX_UNUSED_ARG(time_n), real64 const &GEOSX_UNUSED_ARG(dt),
    DomainPartition *const domain, DofManager &dofManager,
    ParallelMatrix &matrix, ParallelVector &rhs, ParallelVector &solution) {
  // Computation of the sparsity pattern
  SetupSystem(domain, dofManager, matrix, rhs, solution);
}

void PhaseFieldDamageFEM::ImplicitStepComplete(
    real64 const &GEOSX_UNUSED_ARG(time_n), real64 const &GEOSX_UNUSED_ARG(dt),
    DomainPartition *const GEOSX_UNUSED_ARG(domain)) {}

void PhaseFieldDamageFEM::SetupDofs(
    DomainPartition const *const GEOSX_UNUSED_ARG(domain),
    DofManager &dofManager) const {
  dofManager.addField(m_fieldName, DofManager::Location::Node,
                      DofManager::Connectivity::Elem);
}

void PhaseFieldDamageFEM::AssembleSystem( real64 const time_n,
                                          real64 const GEOSX_UNUSED_ARG(dt),
                                          DomainPartition *const domain,
                                          DofManager const &dofManager,
                                          ParallelMatrix &matrix,
                                          ParallelVector &rhs) {
  MeshLevel *const mesh =
      domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager *const nodeManager = mesh->getNodeManager();
  ElementRegionManager *const elemManager = mesh->getElemManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  NumericalMethodsManager const *numericalMethodManager =
      domain->getParent()->GetGroup<NumericalMethodsManager>(
          keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const *feDiscretizationManager =
      numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(
          keys::finiteElementDiscretizations);

  arrayView1d<globalIndex const> const &dofIndex =
      nodeManager->getReference<array1d<globalIndex>>(
          dofManager.getKey(m_fieldName));

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  //arrayView1d<R1Tensor> &X = nodeManager->referencePosition();

  // Initialize all entries to zero
  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  // begin region loop
  for (localIndex er = 0; er < elemManager->numRegions(); ++er) {
    ElementRegionBase *const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const *feDiscretization =
        feDiscretizationManager->GetGroup<FiniteElementDiscretization>(
            m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>(
        [&](localIndex const GEOSX_UNUSED_ARG(esr),
            CellElementSubRegion const *const elementSubRegion) {

//          localIndex m_solidMaterialFullIndex = 0;
//          SolidBase * solidModel =  constitutiveRelations[er][esr][m_solidMaterialFullIndex]->group_cast<SolidBase*>();

#if 0
          arrayView2d<real64 const> const & strainEnergy =
          elementSubRegion->GetConstitutiveModels()->GetGroup(m_solidModelName)->getReference<array2d<real64> >(SolidBase::viewKeyStruct::strainEnergyDensityString);
#else
#define strainEnergy( k, q ) m_coeff[k]
#endif
          array3d<R1Tensor> const &
          dNdX = elementSubRegion->getReference<array3d<R1Tensor>>(keys::dNdX);

          arrayView2d<real64> const &detJ =
              elementSubRegion->getReference<array2d<real64>>(keys::detJ);

          localIndex const numNodesPerElement =
              elementSubRegion->numNodesPerElement();
          arrayView2d<localIndex const,
                      CellBlock::NODE_MAP_UNIT_STRIDE_DIM> const &elemNodes =
              elementSubRegion->nodeList();

          // arrayView1d<real64 const> const &
          // coeff = elementSubRegion->getReference<array1d<real64> >(viewKeyStruct::coeffName);

          globalIndex_array elemDofIndex(numNodesPerElement);
          real64_array element_rhs(numNodesPerElement);
          real64_array2d element_matrix(numNodesPerElement, numNodesPerElement);

          integer_array const &elemGhostRank = elementSubRegion->m_ghostRank;
          localIndex const n_q_points =
              feDiscretization->m_finiteElement->n_quadrature_points();

          real64 ell = 0.05; //phase-field lenght scale
          real64 Gc = 1; //energy release rate
          double threshold = 3 * Gc / (16 * ell); //elastic energy threshold - use when LocalDissipation is Linear
          //real64 diffusion = 1.0;
          // begin element loop, skipping ghost elements
          for (localIndex k = 0; k < elementSubRegion->size(); ++k) {
            if (elemGhostRank[k] < 0) {
              element_rhs = 0.0;
              element_matrix = 0.0;
              for (localIndex q = 0; q < n_q_points; ++q) {
                double D = 0; //max between threshold and Elastic energy
                if(m_localDissipationOption == "Linear"){
                  D = max(threshold, strainEnergy(k,q));
                }
                /*real64 Xq = 0;
                real64 Yq = 0;
                real64 Zq = 0;*/
                //for (localIndex a = 0; a < numNodesPerElement; ++a) {
                  /*Xq = Xq + feDiscretization->m_finiteElement->value(a, q) *
                                X[elemNodes(k, a)][0];

                  Yq = Yq + feDiscretization->m_finiteElement->value(a, q) *
                                X[elemNodes(k, a)][1];

                  Zq = Zq + feDiscretization->m_finiteElement->value(a, q) *
                                X[elemNodes(k, a)][2];*/
                //}
                for (localIndex a = 0; a < numNodesPerElement; ++a) {
                  elemDofIndex[a] = dofIndex[elemNodes(k, a)];
                  //real64 diffusion = 1.0;
                  real64 Na = feDiscretization->m_finiteElement->value(a, q);
                  //element_rhs(a) += detJ[k][q] * Na * myFunc(Xq, Yq, Zq); //older reaction diffusion solver
                  if (m_localDissipationOption == "Linear"){
                    element_rhs(a) += -detJ[k][q] * Na * (- ell * D + 3 * Gc / 16 )/ Gc;
                  }
                  else{
                    element_rhs(a) += detJ[k][q] * Na * (2 * ell) * strainEnergy(k,q) / Gc;
                  }
                  for (localIndex b = 0; b < numNodesPerElement; ++b) {
                    real64 Nb = feDiscretization->m_finiteElement->value(b, q);
                    if(m_localDissipationOption == "Linear"){
                       element_matrix(a, b) += detJ[k][q] *
                       (0.375*pow(ell,2) * -Dot(dNdX[k][q][a], dNdX[k][q][b]) -
                        (ell * D/Gc) * Na * Nb);
                    }
                    else{
                        element_matrix(a, b) +=
                        detJ[k][q] *
                        (pow(ell,2) * -Dot(dNdX[k][q][a], dNdX[k][q][b]) -
                         Na * Nb * (1 + 2 * ell*strainEnergy(k,q)/Gc));
                    }
                }
              }
            }
              matrix.add(elemDofIndex, elemDofIndex, element_matrix);
              rhs.add(elemDofIndex, element_rhs);
            }
          }
        });
  }
  matrix.close();
  rhs.close();

  if (getLogLevel() == 2) {
    GEOS_LOG_RANK_0("After PhaseFieldDamageFEM::AssembleSystem");
    GEOS_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOS_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }

  if (getLogLevel() >= 3) {
    SystemSolverParameters *const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_" + std::to_string(time_n) + "_" +
                          std::to_string(newtonIter) + ".mtx";
    matrix.write(filename_mat, true);

    string filename_rhs = "rhs_" + std::to_string(time_n) + "_" +
                          std::to_string(newtonIter) + ".mtx";
    rhs.write(filename_rhs, true);

    GEOS_LOG_RANK_0("After PhaseFieldDamageFEM::AssembleSystem");
    GEOS_LOG_RANK_0("Jacobian: written to " << filename_mat);
    GEOS_LOG_RANK_0("Residual: written to " << filename_rhs);
  }
}

void PhaseFieldDamageFEM::ApplySystemSolution(DofManager const &dofManager,
                                               ParallelVector const &solution,
                                               real64 const scalingFactor,
                                               DomainPartition *const domain) {
  MeshLevel *const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  NodeManager *const nodeManager = mesh->getNodeManager();

  dofManager.copyVectorToField(solution, m_fieldName, scalingFactor,
                               nodeManager, m_fieldName);

  // Syncronize ghost nodes
  std::map<string, string_array> fieldNames;
  fieldNames["node"].push_back(m_fieldName);

  CommunicationTools::SynchronizeFields(
      fieldNames, mesh,
      domain->getReference<array1d<NeighborCommunicator>>(
          domain->viewKeys.neighbors));
}

void PhaseFieldDamageFEM::ApplyBoundaryConditions(
    real64 const time_n, real64 const dt, DomainPartition *const domain,
    DofManager const &dofManager, ParallelMatrix &matrix, ParallelVector &rhs)
{
  ApplyDirichletBC_implicit(time_n + dt, dofManager, *domain, m_matrix, m_rhs);



  if (getLogLevel() == 2) {
    GEOS_LOG_RANK_0("After PhaseFieldDamageFEM::ApplyBoundaryConditions");
    GEOS_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOS_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }

  if (getLogLevel() >= 3) {
    SystemSolverParameters *const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_bc_" + std::to_string(time_n) + "_" +
                          std::to_string(newtonIter) + ".mtx";
    matrix.write(filename_mat, true);

    string filename_rhs = "rhs_bc_" + std::to_string(time_n) + "_" +
                          std::to_string(newtonIter) + ".mtx";
    rhs.write(filename_rhs, true);

    GEOS_LOG_RANK_0("After PhaseFieldDamageFEM::ApplyBoundaryConditions");
    GEOS_LOG_RANK_0("Jacobian: written to " << filename_mat);
    GEOS_LOG_RANK_0("Residual: written to " << filename_rhs);
  }
}

void PhaseFieldDamageFEM::SolveSystem(DofManager const &dofManager,
                                       ParallelMatrix &matrix,
                                       ParallelVector &rhs,
                                       ParallelVector &solution) {
  rhs.scale(-1.0); // TODO decide if we want this here
  solution.zero();

  SolverBase::SolveSystem(dofManager, matrix, rhs, solution);

  if (getLogLevel() == 2) {
    GEOS_LOG_RANK_0("After PhaseFieldDamageFEM::SolveSystem");
    GEOS_LOG_RANK_0("\nSolution\n");
    std::cout << solution;
  }
}

void PhaseFieldDamageFEM::ApplyDirichletBC_implicit( real64 const time,
                                                     DofManager const &dofManager, DomainPartition &domain,
                                                     ParallelMatrix &matrix,
                                                     ParallelVector &rhs )

{
  FieldSpecificationManager const * const fsManager =
      FieldSpecificationManager::get();

  fsManager->Apply( time, &domain, "nodeManager", m_fieldName,
                    [&](FieldSpecificationBase const *const bc, string const &,
                        set<localIndex> const &targetSet,
                        Group *const targetGroup,
                        string const GEOSX_UNUSED_ARG(fieldName)) -> void
  {
                      bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual,
                      LAInterface>(
                          targetSet, false, time, targetGroup, m_fieldName,
                          dofManager.getKey(m_fieldName), 1, matrix, rhs);
  } );

  fsManager->ApplyFieldValue< serialPolicy >( time, &domain, "ElementRegions", viewKeyStruct::coeffName );

}

REGISTER_CATALOG_ENTRY(SolverBase, PhaseFieldDamageFEM, std::string const &,
                       Group *const)
} // namespace geosx
