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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file LagrangeExplicitDynamicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef LAGRANGEEXPLICITDYNAMICSSOLVER_H_
#define LAGRANGEEXPLICITDYNAMICSSOLVER_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif

class LagrangeExplicitDynamicsSolver : public SolverBase
{
public:
  LagrangeExplicitDynamicsSolver( const std::string& name,
                                  ProblemManagerT* const pm );
  ~LagrangeExplicitDynamicsSolver();

  virtual double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const array<string>& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);

  virtual void SetMaxStableTimeStep( const realT& time,
                                     PhysicalDomainT& domain,
                                     const array<string>& namesOfSolverRegions,
                                     SpatialPartition& partition);

  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  virtual void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  //  void LagrangianNodeUpdatePart1of2( ObjectDataStructureBaseT& nodeManager,
  // const realT& time, const realT& dt );
  //  void LagrangianNodeUpdatePart2of2( ObjectDataStructureBaseT& nodeManager,
  // const realT& time, const realT& dt );


  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn);

  static const char*
  SolverName()
  {
    return "LagrangeExplicitDynamicsSolver";
  }
  ;


  void ApplyGapDamping( NodeManager& nodeManager,
                        const FaceManagerT& faceManager,
                        const realT dt );


  realT m_dampingM;
  realT m_gapdamping;
  array<lSet> m_KinematicConstraintNodes;
//  array<std::pair<localIndex,localIndex> > m_KinematicConstraintFaces;

  bool m_tiedNodesFlag;
  realT m_tiedNodeNormalRuptureStress;
  realT m_tiedNodeShearRuptureStress;
  realT m_tiedNodeTolerance;

  void ApplyForcesFromContact(PhysicalDomainT& domain,
                              StableTimeStep& timeStep,
                              const realT dt );


protected:


private:
  virtual void WriteSiloDerived( SiloFile& siloFile ) const;
  virtual void ReadSiloDerived( const SiloFile& siloFile );



};

namespace LagrangeExplicitDynamicsFunctions
{

void LinearPointUpdatePart1( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const bool clearForces = true);

void LinearPointUpdatePart2( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const realT& damping = 0.0 );

void RotationalPointUpdatePart1(DiscreteElementManagerBaseT& discreteElementManager,
                                const realT& time, const realT& dt);

void RotationalPointUpdatePart1b(DiscreteElementManagerT& discreteElementManager);


void RotationalPointUpdatePart2(ObjectDataStructureBaseT& objectManager,
                                const realT& time, const realT& dt);
}

#endif /* LAGRANGEEXPLICITDYNAMICSSOLVER_H_ */
