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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file ParallelPlateFlowSolverExplicit.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVEREXPLICIT_H_
#define PARALLELPLATEFLOWSOLVEREXPLICIT_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"

class ParallelPlateFlowSolverExplicit : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverExplicit( const std::string& name,
                           ProblemManagerT* const pm );
  virtual ~ParallelPlateFlowSolverExplicit();

  double TimeStep( const realT& time,
                 const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const array<string>& namesOfSolverRegions,
                 SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);
  virtual void RegisterFields( PhysicalDomainT& domain );

  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
          realT time,realT dt );
  
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;


  /// name of solver class
  /// NB: Name() function of parent class returns name of specific solver objects 
  static const char* SolverName(){return "ParallelPlateFlowSolverExplicit";};

  void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain );

  void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  void CalculateMatrixFlowLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain );

  void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
          BoundaryConditionBase* bc ,
          const lSet& set,
          realT time );

  void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain );

  void CalculateNodalPressure ( PhysicalDomainT& domain, SpatialPartition& partition);

private:
  realT m_bBarton;
  realT m_aBarton;
  realT m_wZeroStress;
  realT m_dT;
  realT m_farFieldPorePressure;
  int m_pressureDependentLeakoff;
  realT m_apertureMovingAverageCoeff;


};

#endif /* PARALLELPLATEFLOWSOLVEREXPLICIT_H_ */
