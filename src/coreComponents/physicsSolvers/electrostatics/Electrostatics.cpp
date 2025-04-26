/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Electrostatics.cpp
 */

#include "Electrostatics.hpp"
#include "ElectrostaticsKernels.hpp"

#include "constitutive/electroChemistry/ElectroChemistryBase.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/TractionBoundaryCondition.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

Electrostatics::Electrostatics(const string& name, Group* const parent)
:
PhysicsSolverBase(name, parent),
m_fieldName("primaryField"),
m_timeIntegrationOption(TimeIntegrationOption::QuasiStatic)
{
  registerWrapper(viewKeyStruct::timeIntegrationOption(), &m_timeIntegrationOption).
    setApplyDefaultValue(m_timeIntegrationOption).
    setInputFlag(InputFlags:: OPTIONAL).
    setDescription("Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n*" ));

  registerWrapper(viewKeyStruct::fieldVarName(), &m_fieldName).
    setRTTypeName(rtTypes::CustomTypes::groupNameRef).
    setInputFlag(InputFlags::REQUIRED).
    setDescription("Name of field variable");
}

Electrostatics::~Electrostatics() {}

void Electrostatics::registerDataOnMesh(Group& meshBodies)
{
  forDiscretizationOnMeshTargets(meshBodies,
    [&](string const&, MeshLevel& mesh, string_array const& regionNames) {
      NodeManager& nodes = mesh.getNodeManager();

      nodes.registerWrapper<real64_array>(m_fieldName)
        .setApplyDefaultValue(0.0)
        .setPlotLevel(PlotLevel::LEVEL_0)
        .setDescription("Primary field variable");

      ElementRegionManager& elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< CellElementSubRegion >(regionNames,
        [&]( localIndex const, CellElementSubRegion& subRegion) {
          subRegion.registerWrapper<string>(viewKeyStruct::electroMaterialNamesString()).
            setPlotLevel(PlotLevel::NOPLOT).
            setRestartFlags(RestartFlags::NO_WRITE).
            setSizedFromParent(0);

          string& electroMaterialName = subRegion.getReference<string>(viewKeyStruct::electroMaterialNamesString());
          electroMaterialName = PhysicsSolverBase::getConstitutiveName<ElectroChemistryBase>(subRegion);
          GEOS_ERROR_IF( electroMaterialName.empty(), GEOS_FMT("{}: ElectroChemistryBase model not found on subregion {}",
                                                               getDataContext(), subRegion.getName()));
        });  
    });
}

real64 Electrostatics::solverStep(real64 const& time_n, real64 const& dt,
                                  const int cycleNumber, DomainPartition& domain)
{
  GEOS_MARK_FUNCTION;
  real64 dtReturn = dt;

  setupSystem(domain, m_dofManager, m_localMatrix, m_rhs, m_solution, false);

  dtReturn = linearImplicitStep(time_n, dt, cycleNumber, domain);

  return dtReturn;
}

void Electrostatics::implicitStepSetup(real64 const& GEOS_UNUSED_PARAM(time_n),
                                       real64 const& GEOS_UNUSED_PARAM(dt),
                                       DomainPartition& domain)
{
  Timestamp const meshModificationTimestamp = getMeshModificationTimestamp(domain);

  if (meshModificationTimestamp > getSystemSetupTimestamp())
  {
    setupSystem(domain, m_dofManager, m_localMatrix, m_rhs, m_solution);
    setSystemSetupTimestamp(meshModificationTimestamp);
  }
}

void Electrostatics::setupDofs(DomainPartition const& GEOS_UNUSED_PARAM(domain),
                               DofManager& dofManager) const
{
  GEOS_MARK_FUNCTION;
  dofManager.addField(m_fieldName, FieldLocation::Node, 1, getMeshTargets());
  dofManager.addCoupling(m_fieldName, m_fieldName, DofManager::Connector::Elem);
}

void Electrostatics::setupSystem(DomainPartition& domain, DofManager& dofManager,
                                 CRSMatrix< real64, globalIndex >& localMatrix,
                                 ParallelVector& rhs, ParallelVector& solution,
                                 bool const setSparsity)
{
  GEOS_LOG("Electrostatics::setupSystem");

  GEOS_MARK_FUNCTION;
  PhysicsSolverBase::setupSystem(domain, dofManager, localMatrix, rhs, solution, setSparsity);

  SparsityPattern<globalIndex> sparsityPattern(dofManager.numLocalDofs(), dofManager.numGlobalDofs(), 8*8*1.2);

  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
    [&](string const&, MeshLevel& mesh, string_array const& regionNames) {
      NodeManager const& nodeManager = mesh.getNodeManager();
      arrayView1d<globalIndex const> const dofNumber = nodeManager.getReference<globalIndex_array>(dofManager.getKey(m_fieldName));

      finiteElement::fillSparsity<CellElementSubRegion, ElectrostaticsKernel>(mesh, regionNames, getDiscretizationName(),
                                                                              dofNumber, dofManager.rankOffset(), sparsityPattern);
    }
  );

  sparsityPattern.compress();
  localMatrix.assimilate<parallelDevicePolicy<>>(std::move(sparsityPattern));
}

void Electrostatics::assembleSystem(real64 const GEOS_UNUSED_PARAM(time_n), real64 const dt,
                                    DomainPartition& domain, DofManager const& dofManager,
                                    CRSMatrixView<real64, globalIndex const> const& localMatrix,
                                    arrayView1d<real64> const& localRhs)
{
  GEOS_MARK_FUNCTION;

  localMatrix.zero();
  localRhs.zero();

  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
    [&](string const&, MeshLevel& mesh, string_array const& regionNames) {
      NodeManager& nodeManager = mesh.getNodeManager();
      string const dofKey = dofManager.getKey(m_fieldName);
      arrayView1d<globalIndex const> const& dofIndex = nodeManager.getReference<array1d<globalIndex>>(dofKey);

      ElectrostaticsKernelFactory kernelFactory(dofIndex, dofManager.rankOffset(), localMatrix, localRhs, dt, m_fieldName);

      finiteElement::regionBasedKernelApplication<parallelDevicePolicy<>, constitutive::ElectroChemistryBase, CellElementSubRegion>(
        mesh, regionNames, this->getDiscretizationName(), viewKeyStruct::electroMaterialNamesString(), kernelFactory);
    });
}

void Electrostatics::applyBoundaryConditions(real64 const time_n, real64 const dt,
                                             DomainPartition& domain, DofManager const& dofManager,
                                             CRSMatrixView<real64, globalIndex const> const& localMatrix,
                                             arrayView1d<real64> const& localRhs)
{
  applyPotentialBC(time_n + dt, dofManager, domain, localMatrix, localRhs);

  applyCurrentBC(time_n + dt, dofManager, domain, localRhs);
}

// real64 Electrostatics::calculateResidualNorm(real64 const& GEOS_UNUSED_PARAM(time_n), real64 const& GEOS_UNUSED_PARAM(dt),
//                                              DomainPartition const& domain, DofManager const& dofManager,
//                                              arrayView1d<real64 const> const& localRhs)
// {
//   GEOS_MARK_FUNCTION;
// 
//   real64 totalResidualNorm = 0.0;
// 
//   forDiscretizationOnMeshTargets(domain.getMeshBodies(),
//     [&](string const&, MeshLevel const& mesh, string_array const&){
//       NodeManager const& nodeManager = mesh.getNodeManager();
//       string const dofKey = dofManager.getKey(m_fieldName);
//       arrayView1d<globalIndex const> const dofNumber = nodeManager.getReference<array1d<globalIndex>>(dofKey);
// 
//       globalIndex const rankOffset = dofManager.rankOffset();
//       arrayView1d<integer const> const ghostRank = nodeManager.ghostRank();
// 
//       RAJA::ReduceSum<parallelDeviceReduce, real64> localSum(0.0);
// 
//       SortedArrayView<localIndex const> const& targetNodes = nodeManager.sets()
//     });
// }

void Electrostatics::applySystemSolution(DofManager const& dofManager, arrayView1d<real64 const> const& localSolution,
                                         real64 const scalingFactor, real64 const dt, DomainPartition& domain)
{
  GEOS_UNUSED_VAR(dt);
  dofManager.addVectorToField(localSolution, m_fieldName, m_fieldName, scalingFactor);

  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
    [&](string const&, MeshLevel& mesh, string_array const&) {
      FieldIdentifiers fieldsToBeSync;
      fieldsToBeSync.addFields(FieldLocation::Node, {m_fieldName});

      CommunicationTools::getInstance().synchronizeFields(fieldsToBeSync, mesh, domain.getNeighbors(), true);
    });
}

void Electrostatics::applyPotentialBC(real64 const time, DofManager const& dofManager, DomainPartition& domain,
                                      CRSMatrixView<real64, globalIndex const> const& localMatrix,
                                      arrayView1d<real64> const& localRhs)
{
  FieldSpecificationManager const& fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
    [&](string const&, MeshLevel& mesh, string_array const&)
    {
      fsManager.apply<NodeManager>(time, mesh, m_fieldName,
        [&](FieldSpecificationBase const& bc, string const&, SortedArrayView<localIndex const> const& targetSet,
            NodeManager& targetGroup, string const& GEOS_UNUSED_PARAM(fieldName))
        {
          bc.applyBoundaryConditionToSystem<FieldSpecificationEqual, parallelDevicePolicy<>>(
            targetSet, time, targetGroup, m_fieldName, dofManager.getKey(m_fieldName),
            dofManager.rankOffset(), localMatrix, localRhs);
        });
    });
}

void Electrostatics::applyCurrentBC(real64 const time, DofManager const& dofManager,
                                    DomainPartition& domain, arrayView1d<real64> const& localRhs)
{
  FieldSpecificationManager& fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets(domain.getMeshBodies(),
    [&](string const&, MeshLevel& mesh, string_array const&)
    {
      FaceManager const& faceManager = mesh.getFaceManager();
      NodeManager const& nodeManager = mesh.getNodeManager();

      string const dofKey = dofManager.getKey(m_fieldName);
      arrayView1d<globalIndex const> const blockLocalDofNumber = nodeManager.getReference<globalIndex_array>(dofKey);
      globalIndex const dofRankOffset = dofManager.rankOffset();

      fsManager.template apply<FaceManager, TractionBoundaryCondition>(time, mesh, TractionBoundaryCondition::catalogName(),
        [&](TractionBoundaryCondition const& bc, string const&, SortedArrayView<localIndex const> const& targetSet,
            Group&, string const&)
        {
          bc.launch(time, blockLocalDofNumber, dofRankOffset, faceManager, targetSet, localRhs);
        });
    });
}

void Electrostatics::updateState(DomainPartition& domain)
{
  GEOS_UNUSED_VAR(domain);
}

void Electrostatics::resetStateToBeginningOfStep(DomainPartition& GEOS_UNUSED_PARAM(domain)) {}

void Electrostatics::implicitStepComplete(real64 const& GEOS_UNUSED_PARAM(time_n),
                                          real64 const& GEOS_UNUSED_PARAM(dt),
                                          DomainPartition& GEOS_UNUSED_PARAM(domain))
{}

REGISTER_CATALOG_ENTRY(PhysicsSolverBase, Electrostatics, string const&, dataRepository::Group* const)
}