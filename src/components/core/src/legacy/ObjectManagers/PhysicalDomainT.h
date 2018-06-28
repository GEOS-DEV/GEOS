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
 * @file PhysicalDomainT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef PHYSICAL_DOMAIN_T_H_
#define PHYSICAL_DOMAIN_T_H_

#include "DataStructures/EncapsulatedFields/NodeManager.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "ElementManagerT.h"
#include "FaceManagerT.h"
#include "ExternalFaceManagerT.h"
#include "ContactManagerT.h"
#include "CartesianGridManagerT.h"
#include "EdgeManagerT.h"
#include "EllipsoidalContactManagerT.h"
#include "DiscreteElementManagerT.h"
#include "EllipsoidalDiscreteElementManagerT.h"

#ifdef SRC_INTERNAL2
#include "proprietary_internal/src/ObjectManagers/MicroseismicElementManagerT.h"
#include "proprietary_internal/src/SurfaceGeneration/JointPopulator2.h"
#include "proprietary_internal/src/SurfaceGeneration/XfemManager.h"
#endif
#include "EnergyT.h"

#include "CrackObject.h"
#include "CrackSurfaceVertex.h"

#ifdef SRC_EXTERNAL
#include "proprietary_Ok4Export/src_external/ObjectManagers/WellboreManagerT.h"
#include "ObjectManagers/FaultElementManagerT.h"
#endif

class PhysicalDomainT
{
public:
  PhysicalDomainT();

  ~PhysicalDomainT();


  int Initialize();

  void SetInterObjectRelations();

  int m_globalDomainNumber;

  //note: because discrete elements have nodes with much reduced
  // degrees-of-freedom,
  //as well as different solution requirements, we keep a separate collection
  // specifically
  //for discrete elements, which, in turn, necessitates a separate collection of
  // faces
  //the ExternalFaceManagerT ties the external faces of the finite elements
  //with those in the discrete elements to provide coupling between the methods

  EllipsoidalDiscreteElementManagerT m_ellipsoidalDiscreteElementManager;
  EllipsoidalContactManagerT m_ellipsoidalContactManager;

  FaceManagerT m_discreteElementSurfaceFaces;
  NodeManager m_discreteElementSurfaceNodes;
  DiscreteElementManagerT m_discreteElementManager;

  ContactManagerT m_contactManager;

  NodeManager m_feNodeManager;
  ElementManagerT m_feElementManager;
  FaceManagerT m_feFaceManager;
  EdgeManagerT m_feEdgeManager;

  CrackObjectManager* m_crackManager;
  CrackSurfaceVertex* m_crackSurfaceVertex;

  bool m_splitElemThisIter = false;
  bool m_WriteXFEM = false;

  NodeManager m_flowNodes;
  EdgeManagerT m_flowEdges;
  FaceManagerT m_flowFaces;

  ExternalFaceManagerT m_externalFaces;

  CartesianGridManagerT m_cartesianGridManager;

#ifdef SRC_EXTERNAL
  FaceManagerT m_faultPatchFaces;
  NodeManager m_faultPatchNodes;
  CartesianGridManagerT m_faultPorePressure;
  FaultElementManagerT m_faultElementManager;
  WellboreManagerT m_wellboreManager;
#endif
#ifdef SRC_INTERNAL2
  MicroseismicElementManagerT m_microseismicElementManager;
  JointPopulator2 m_jointSets2;
#endif

#ifdef SRC_INTERNAL2
  XfemManager* m_xfemManager;
#endif
  EnergyT m_energy;

  void WriteSilo( SiloFile& siloFile,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const bool writeXFEM,
                  const bool writeFEMMesh = true,
                  const bool writeFEMFaces = false,
                  const bool writeFEMEdges = false,
                  const bool writeDE = true,
                  const bool writeCP = true,
                  const bool writeCG = true,
                  const bool writeFaultElements = true,
                  const bool writeMicroseismicSources = true);



  void ReadSilo( const SiloFile& siloFile,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart,
                 const bool writeFEMMesh = true,
                 const bool writeFEMFaces = false,
                 const bool writeFEMEdges = false,
                 const bool writeDE = true,
                 const bool writeCP = true,
                 const bool writeCG = true,
                 const bool writeFaultElements = true,
                 const bool writeMicroseismicSources = true);

  void WriteFiniteElementMesh( SiloFile& siloFile,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart,
                               const bool writeFEMMesh,
                               const bool writeFEMFaces,
                               const bool writeFEMEdges);

  void ReadFiniteElementMesh( const SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart );

#ifdef SRC_EXTERNAL
  void WriteFaultElementMesh( SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const bool writeFaultElements);
  void ReadFaultElementMesh( const SiloFile& siloFile,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart );
#endif

  void WriteDiscreteElements( SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const bool writeDE );

  void WriteXFEMElements( SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart,
                          const bool writeXFEM );

  void ReadDiscreteElements( const SiloFile& siloFile,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart );

  void WriteEllipsoidalDiscreteElements( SiloFile& siloFile,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const bool writeEDE );

  void ReadEllipsoidalDiscreteElements( const SiloFile& siloFile,
                                        const int cycleNum,
                                        const realT problemTime,
                                        const bool isRestart );


  void WriteCommonPlanes( SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart,
                          const bool writeCP );

  void ReadCommonPlanes( const SiloFile& siloFile,
                         const int cycleNum,
                         const realT problemTime,
                         const bool isRestart );

  void WriteCartesianGrid( SiloFile& siloFile,
                           const int cycleNum,
                           const realT problemTime,
                           const bool isRestart,
                           const bool writeCG );

  void ReadCartesianGrid( const SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart );

  void RegisterBCFields();
  void RegisterBCFields(ObjectDataStructureBaseT& objectManager);

  void UpdateEnergy();

  //NOTE: the order of these is important; they must be in order of ascending
  // dependence
  //      e.g., faces are dependent on nodes, so faces must be of greater
  // ordinal value than nodes
  enum ObjectDataStructureKeys
  {
    EllipsoidalDiscreteElementManager  = 0,
    EllipsoidalContactManager = 1,
    DiscreteElementNodeManager = 2,
    DiscreteElementFaceManager = 3,
    DiscreteElementManager = 4,
    ContactManager = 5,
    FiniteElementNodeManager = 6,
    FiniteElementEdgeManager = 7,
    FiniteElementFaceManager = 8,
    FiniteElementElementManager = 9,
    FiniteElementElementRegion = 10,
    ExternalFaceManager = 11,
    CartesianGridManager = 12,
    FaultPatchNodeManager = 13,
    FaultPatchFaceManager = 14,
    FaultPatchElementManager = 15,
    VirtualEdgeManager = 16,
    VirtualFaceManager = 17,
    Last_ObjectDataStructureNames_Index = 18,
    numObjectDataStructureNames = 19,
    WellboreElementManager = 20
  };



  static const char* EllipsoidalDiscreteElementManagerStr() { return "EllipsoidalDiscreteElement"; }
  static const char* EllipsoidalContactManagerStr()         { return "EllipsoidalContactManager"; }
  static const char* EllipsoidalConditionStr()              { return "EllipsoidalDiscreteElement"; }

  static const char* DiscreteElementNodeManagerStr()        { return "DiscreteElement_NodeManager"; }
  static const char* DiscreteElementFaceManagerStr()        { return "DiscreteElement_FaceManager"; }
  static const char* DiscreteElementManagerStr()            { return "DiscreteElementManager"; }
  static const char* DiscreteElementConditionStr()          { return "DiscreteElement"; }

  static const char* ContactManagerStr()                    { return "ContactManager"; }
  static const char* CartesianGridManagerStr()              { return "CartesianGrid"; }

  static const char* FiniteElementNodeManagerStr()          { return "FiniteElement_NodeManager"; }
  static const char* FiniteElementNodeConditionStr()        { return "Node"; }
  static const char* FiniteElementEdgeManagerStr()          { return "FiniteElement_EdgeManager"; }
  static const char* FiniteElementEdgeConditionStr()        { return "Edge"; }
  static const char* FiniteElementFaceManagerStr()          { return "FiniteElement_FaceManager"; }
  static const char* FiniteElementFaceConditionStr()        { return "Face"; }
  static const char* FiniteElementElementManagerStr()       { return "FiniteElement_ElementManager"; }
  static const char* FiniteElementElementConditionStr()     { return "Element"; }
  static const char* FiniteElementElementRegionStr()        { return "FiniteElement_ElementRegion"; }

  static const char* ExternalFaceManagerStr()               { return "ExternalFaceManager"; }
  static const char* ExternalFaceConditionStr()             { return "ExternalFace"; }

  static const char* FaultNodeManagerStr()                  { return "Fault_NodeManager"; }
  static const char* FaultFaceManagerStr()                  { return "Fault_FaceManager"; }
  static const char* FaultElementManagerStr()               { return "Fault_ElementManager"; }
  static const char* FaultElementConditionStr()             { return "FaultElement"; }

  static const char* VirtualEdgeManagerStr()                { return "Virtual_EdgeManager"; }
  static const char* VirtualFaceManagerStr()                { return "Virtual_FaceManager"; }

  static const char* WellboreElementManagerStr()            { return "Wellbore"; }

  ObjectDataStructureBaseT& GetObjectDataStructure( ObjectDataStructureKeys key,
                                                    const std::string& regionName="" );

  static ObjectDataStructureKeys GetObjectDataStructureConditionKey(const std::string name);
  static ObjectDataStructureKeys GetObjectDataStructureKey(const std::string name);
  static std::string GetObjectDataStructureName(const ObjectDataStructureKeys key);

  template< typename T_indices >
  void Pack( const ObjectDataStructureKeys name,
             const T_indices& sendIndices,
             bufvector& buffer,
             const bool packConnectivityToGlobal,
             const bool packFields,
             const bool packMaps,
             const bool packSets );

  template< typename T_indices >
  void Unpack(const ObjectDataStructureKeys name,
              const char*& pbuffer,
              T_indices& indices,
              const bool unpackConnectivityToLocal,
              const bool unpackFields,
              const bool unpackMaps,
              const bool unpackSets  );



private:
  PhysicalDomainT( const PhysicalDomainT& );
  PhysicalDomainT& operator=( const PhysicalDomainT& rhs );

};

// the following facilitates looping over object data structure keys
// nb keys must be contiguous and 0 start.
inline
PhysicalDomainT::ObjectDataStructureKeys& operator++(PhysicalDomainT::ObjectDataStructureKeys& objDSK)
{
  // if ENUMS ARE CONTIGUOUS
  int i = static_cast<int>(objDSK);
  if( ++i > PhysicalDomainT::numObjectDataStructureNames )
    i = 0;
  objDSK = static_cast<PhysicalDomainT::ObjectDataStructureKeys>(i);
  return objDSK;
}


#endif /* PHYSICAL_DOMAIN_T_H_ */
