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

/**
 * @file ProblemManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef PROBLEMMANAGERT_H_
#define PROBLEMMANAGERT_H_

#include "../Common/GPException.h"
#include "DataStructures/Tables/Table.h"
#include "PhysicalDomainT.h"

#include "IO/FileManagerDataT.h"
#include "IO/FileManagerT.h"
#include "IO/silo/SiloFile.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "PhysicsSolvers/SolverApplicationSet.h"
#include "PhysicsSolvers/SolverBase.h"
#include "MPI_Communications/SpatialPartition.h"
#include "SurfaceGeneration/Fractunator.h"

#ifdef SRC_INTERNAL2
//@ annavarapusr1
#include "SurfaceGeneration/XfemManager.h"
#endif

#if 1
// This is for writing out node-based graph.  Pengcheng is working on it.
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#endif

#include "Constitutive/ConstitutivePropertiesTable.h"

#include "MeshUtilities/MeshGenerator.h"

#include "Epetra_ConfigDefs.h"
#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

using namespace TICPP;
class InitialConditionBase;

/**
 * @author Randolph Settgast
 */
class ProblemManagerT
{
public:
  ProblemManagerT();
  ~ProblemManagerT();

  void ParseCommandLineInput( const int& argc, char* const argv[]) throw(GPException);

  void ProblemSetup();

  void ReadGeometryInput( HierarchicalDataNode& hdn );

  void RegisterFields();

  void DisplayFields();

  void DisplayUnits();

  void SetInitialConditions();

  void SetInitialConstitutive();

  void CompleteObjectInitialization( );

  void OutputMeshInformation();

  void WriteSilo( const bool isRestart );

  void WriteFlowTxt (const int cycleNum, const realT time);

  void ReadSilo( const bool isRestart );

  void CalculateNextSiloTime( const std::string& intervalTableName,
                              realT& siloIntervals,
                              realT& nextSiloTime,
                              globalIndex& nextSiloIndex );



  void DisplayUsage( );

  void DisplayVersion( );

  void DisplaySplash( );

  void DisplaySolvers( );

  void VerifySolvers( );

  void UpdateOwnership( );

  void WriteNodeGraph( const localIndex nodeID,
                       PhysicalDomainT& domain);


  /// Run the problem
  void Run( double& outputTime, realT t_start );

  std::map<std::string,SolverBase*> m_solvers;
  std::vector<InitialConditionBase*> m_initialConditions;
  ConstitutivePropertiesTable m_initialConstitutive;

  GPAC_IO::FileManagerT m_FileManager;

  MeshGenerator* m_MeshGenerator;
  MeshGenerator myMeshGenerator;


  int m_numDomains;

  realT m_problemTime;
  realT m_dt;
  realT m_maxWallTime;
  bool m_forceOutputTime;

  int m_size;
  int m_rank;
  int m_cycleNumber;
  int m_cycleReportFreq;

  SiloFile m_siloFile;

  bool m_writeFEM;
  bool m_writeFEMFaces;
  bool m_writeFEMEdges;
  bool m_writeFlowText;
//  bool m_writeXFEM;
  array<string> m_nameFieldsToPlot;
  array<string> m_nameAdditionalFieldsToPlot;

  int m_fractureFlag;
  array<string> m_preFractureSets;

  bool m_writePlot;
  std::string m_visitFileGroupFile; // ".visit" file with list of time series
                                    // files (usually "geos.visit")
  realT m_plotIntervals;
  globalIndex m_nextPlotIndex;
  std::string m_plotIntervalTableName;
  realT m_nextPlotTime;

  bool m_writeRestart;
  realT m_restartIntervals;
  globalIndex m_nextRestartIndex;
  std::string m_restartIntervalTableName;
  realT m_nextRestartTime;


  //Fractunator3 m_surfaceSeparation;
  FractunatorBase* m_surfaceSeparation;

  //@annavarapusr1
  int m_xfem;
//  XfemManager* m_elementSplitting;

  std::vector<SolverApplicationSet> m_solverApplicationSets;

  SpatialPartition m_partition;

  std::map<std::string, std::string> m_simulationParameterMap;

  PhysicalDomainT m_Domains;

  // trilinos communicators
  #if GPAC_MPI
  Epetra_MpiComm m_epetraComm;
  #else
  Epetra_SerialComm m_epetraComm;
  #endif


private:
  bool m_beginFromRestart;
  std::string m_beginFromRestartFileName;

  bool m_doWriteXML;
  std::string m_xmlOutputFileName;

  bool m_hackInitialStress;

  bool m_useMetis = false;
  int m_trackEnergy;
  EnergyT m_energy;
  realT m_initialEnergy;

  bool m_displayFields;
  bool m_displaySplash;

  bool m_echoParameters;

  ProblemManagerT(const ProblemManagerT&);
  ProblemManagerT& operator=(const ProblemManagerT&);
  void RegisterFilesIncludedFromCommandLine(HierarchicalDataNode* hdn, array<string>& includedFiles);
  void ParseIncludedFiles(HierarchicalDataNode* hdn);
  void BuildSimulationParameterMap(HierarchicalDataNode* hdn);
  bool ReplaceMathematicalExpressions(std::string& valueStr);
  void ParseMetadata(HierarchicalDataNode* hdn,bool isRoot);
  bool EvaluateXMLIfThenElseStatement(HierarchicalDataNode* &hdn, HierarchicalDataNode* parentNode);
  void ReadXML(HierarchicalDataNode& hdn);
  void WriteXML(std::string& filename, HierarchicalDataNode& hdn);


  realT SetInitialTimeStep( const realT& time, std::vector<SolverApplicationSet>::iterator solverSet );

  //CompleteObjectInitialization helpers
  void SetDomainBoundaryObjects();
  void InitializeObjectManagers();
  void ResetGlobalToLocal();
  void SetExternal();
  void SetupNeighborListsPackUnpack(std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> >& syncedFields);
  void InitializeSurfaceSeparation();
  void InitializeElementSplitting();

  void UpdateEnergy( const bool set_initial = false );

  static realT SetLocalMaxTimeStep( const realT& hard_dt, const realT& solver_dt, const realT& current_max_dt );
  static realT SetGlobalMaxTimeStep( const realT& local_dt );

  void MarkFieldsToWrite();

};


#endif /* PROBLEMMANAGERT_H_ */
