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
 * @file ProblemManagerT.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "ProblemManagerT.h"

#if GPAC_MPI
#include <mpi.h>
#endif

#include <cstdlib>

#include <getopt.h>
#include <iostream>
#include <ios>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <sys/time.h>

#include "TableManager.h"
#include "FunctionManager.h"

#ifdef USE_CHEM
#include "ChemistryManager.h"
#endif

#include "UnitManager.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/SubstepSolver.h"
#include "InitialConditions/InitialConditions.h"
#include "IO/RestartFile.h"
#include "../../codingUtilities/Functions.hpp"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/IOUtilities.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#include "MeshUtilities/MeshUtilities.h"
#include "SurfaceGeneration/FractunatorFactory.h"
#include "Common/Version.h"

#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"

#ifdef SRC_INTERNAL
  #include "PhysicsSolvers/BackgroundAMR.h"

  #ifdef GPAC_GMM
    #include "ObjectManagers/GeodynMaterialModelManager.h"
  #endif
#endif

/**
 *
 */
ProblemManagerT::ProblemManagerT():
  m_solvers(),
  m_initialConditions(),
  m_initialConstitutive(),
  m_FileManager(),
  m_numDomains(0),
  m_problemTime(0.0),
  m_size(1),
  m_rank(0),
  m_cycleNumber(0),
  m_cycleReportFreq(1),
  m_siloFile(),
  m_writeFEM(true),
  m_writeFEMFaces(true),
  m_writeFEMEdges(false),
  m_writeFlowText(false),
  m_fractureFlag(0),
  m_preFractureSets(),
  m_writePlot(true),
  m_plotIntervals(0.0),
  m_nextPlotIndex(0),
  m_plotIntervalTableName(),
  m_nextPlotTime(0.0),
  m_writeRestart(true),
  m_restartIntervals(0.0),
  m_nextRestartIndex(0),
  m_restartIntervalTableName(),
  m_nextRestartTime(0.0),
  m_surfaceSeparation(NULL),
//m_elementSplitting(NULL),
  m_solverApplicationSets(),
  m_partition(),
  m_simulationParameterMap(),
  m_Domains(),
#if GPAC_MPI
  m_epetraComm(MPI_COMM_WORLD),
#endif
  m_beginFromRestart(false),
  m_beginFromRestartFileName(),
  m_doWriteXML(false),
  m_xmlOutputFileName(),
  m_hackInitialStress(false),
  m_trackEnergy(0),
  m_energy(),
  m_initialEnergy(0.0),
  m_displayFields(false),
  m_displaySplash(true),
  m_echoParameters(false)
{
#if GPAC_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &m_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
#endif

  FieldInfo::AllocateAttributes();

}

ProblemManagerT::~ProblemManagerT()
{
  FieldInfo::DeleteAttributes();

  for (std::map<std::string, SolverBase*>::iterator it_solver = m_solvers.begin() ; it_solver
       != m_solvers.end() ; ++it_solver)
  {
    delete it_solver->second;
  }

  for (std::vector<InitialConditionBase*>::iterator it_ic = m_initialConditions.begin() ; it_ic
       != m_initialConditions.end() ; ++it_ic)
  {
    delete *it_ic;
  }

#ifdef SRC_INTERNAL
  if (doBackgroundAMR)
  {
    int ret=p_exitbamr(0);
    if( ret )
      throw GPException( "Errors during BackgroundAMR exit!" );
  }
#endif
  if( m_surfaceSeparation!=NULL )
  {
    delete m_surfaceSeparation;
  }
}
/**
 * @brief Registers new include files added from the command line
 * @author walsh24
 * Note that files can be included from within the input file also via
 * <Include file="someFilePath.xml">
 */
void ProblemManagerT::RegisterFilesIncludedFromCommandLine(HierarchicalDataNode* hdn, array<string>& includedFiles)
{
  for( array<string>::size_type i =0 ; i < includedFiles.size() ; ++i )
  {
    HierarchicalDataNode* includeNode =  hdn->NewChild("Include");
    includeNode->AddAttributePair("file",includedFiles[i]);
  }
}

/**
 * @brief Convert the file to hierarchical data node format and traverses the
 * graph
 * @author walsh24
 * Performs a recursive search over the hierarchical data nodes for entries of
 * type
 * <Include file="someFilePath.xml">
 * Converts the file to hdn format, copies the result in the hdn, then continues
 * to traverse the (updated) graph.
 * Supports nested and multiple inclusions, but will loop indefinitely if two
 * files include each other.
 * i.e., WARNING: NO CHECK FOR CIRCULAR REFERENCES
 *
 * @param[in,out] hdn Hierarchical data node object
 */
void ProblemManagerT::ParseIncludedFiles(HierarchicalDataNode* hdn)
{

  int count = 0;
  for (HierarchicalDataNode* childNode = hdn->Next(true) ; childNode ; childNode = hdn->Next())
  {
    // load included file
    if (streq(childNode->Heading(), "Include"))
    {
      const std::string path = childNode->GetAttributeString("file");

      std::vector<HierarchicalDataNode> hdnv;
      bool success = TICPP::TinyXMLParser::Load(path.c_str(), hdnv, childNode->Level());

      if (success)
      {
        if (hdnv.size() > 0)
        {
          (*childNode) = hdnv[0];
          for (size_t i = 1 ; i < hdnv.size() ; ++i)
          {
            if (hdnv.size() > 1)
            {
              hdn->InsertAfter(childNode, hdnv.begin() + 1, hdnv.end());
              // recover pointer invalidated by InsertAfter
              int cc = 0;
              childNode = hdn->Next(true);
              while (cc < count)
              {
                childNode = hdn->Next();
                ++cc;
              }
            }
          }
        }
      }
      else
      {
        throw GPException("Error: Could not find included input file " + path);
      }

    }

    // depth first traversal
    ParseIncludedFiles(childNode);
    count++;
  }
}

/**
 * BuildSimulationParameterMap
 *
 * @author walsh24
 *
 * Searches the node for xml entries of the form
 * <Parameter a="0.1" b="some string" etc />
 * and inserts the values into m_simulationParameterMap
 *
 * Map insert will not overwrite existing keys, so the first value for each key
 * is stored in the map - this allows simulation parameters to be overwritten by
 * input from the command line (and previous defined params) simplifying
 * parameter space analysis.
 *
 **/
void ProblemManagerT::BuildSimulationParameterMap(HierarchicalDataNode* hdn)
{

  for (HierarchicalDataNode* childNode = hdn->Next(true) ; childNode ; childNode = hdn->Next())
  {
    if (streq(childNode->Heading(), "Parameters"))
    {
      //parametersFound = true;
      for (HierarchicalDataNode* grandChildNode = childNode->Next(true) ; grandChildNode ; grandChildNode = childNode->Next())
      {
        const std::map<std::string, std::string>& attributes = grandChildNode->GetAttributes();
        std::map<std::string, std::string>::const_iterator iter = attributes.find("name");
        if(iter == attributes.end())
          throw GPException("cannot find name for parameter");
        std::string name = iter->second;
        iter = attributes.find("value");
        if(iter == attributes.end())
          throw GPException("cannot find value for parameter");
        m_simulationParameterMap.insert( std::pair<std::string,std::string>(name,iter->second)); // nb
                                                                                                 // map
                                                                                                 // insert
                                                                                                 // will
                                                                                                 // not
                                                                                                 // overwrite
                                                                                                 // existing
                                                                                                 // keys.

        /** overwrite parameter in hdn **/
        // this allows command line parameter to be written to file in xml
        // output
        // Note that parameters are evaluated before xml output - so this change
        // is cosmetic
        // (the updated parameters should no longer impact the new xml input
        // file, but this will record what the parameters were set to)
        grandChildNode->AddAttributePair("value", m_simulationParameterMap[name]);

      }
    }
    //This next part is for Stuart's legacy format
    //Which looks nicer, but does not work with the schema
    else if (streq(childNode->Heading(), "Parameter"))
    {
      // const std::map<std::string, std::string>& attributes =
      // childNode->GetAttributes();
      std::map<std::string, std::string> attributes = childNode->GetAttributes(); // making
                                                                                  // a
                                                                                  // copy
                                                                                  // to
                                                                                  // avoid
                                                                                  // any
                                                                                  // funny
                                                                                  // business
      m_simulationParameterMap.insert(attributes.begin(), attributes.end());   // nb
                                                                               // map
                                                                               // insert
                                                                               // will
                                                                               // not
                                                                               // overwrite
                                                                               // existing
                                                                               // keys.

      /** overwrite parameters in hdn **/
      std::map<std::string, std::string>::const_iterator iter =attributes.begin();
      std::map<std::string, std::string>::const_iterator iend =attributes.end();
      for( ; iter != iend ; ++iter)
      {
        childNode->AddAttributePair(iter->first, m_simulationParameterMap[iter->first]);
      }
    }

  }

  if(m_echoParameters && m_rank == 0)
  {
    std::map<std::string, std::string>::iterator itr;
    std::cout << "Parameter values: " << std::endl;
    for( itr = m_simulationParameterMap.begin() ; itr != m_simulationParameterMap.end() ; ++itr)
      std::cout << "Key: " << (*itr).first << " Value: " << (*itr).second <<std::endl;
    // exit(0);
  }
}

/**
 * ParseMetadata
 *
 * @author walsh24
 *
 * Traverses nodes, replacing simulation parameters (in the form
 *$:parameter_name) in headings and attributes
 * with their assigned value.
 *
 * Simulation Parameters may contain other parameters but cannot be used to
 * change filenames in "include" statements.
 *
 **/
void ProblemManagerT::ParseMetadata(HierarchicalDataNode* hdn, bool isRoot)
{

  for (HierarchicalDataNode* childNode = hdn->Next(true) ; childNode ; childNode = hdn->Next())
  {

    // Evaluate if nodes
    bool needToReevaluate = true;
    while(needToReevaluate)
    {
      needToReevaluate = EvaluateXMLIfThenElseStatement(childNode, hdn); // needed
                                                                         // for
                                                                         // recursive
                                                                         // if
                                                                         // statements

      if(needToReevaluate && isRoot)
      {
        BuildSimulationParameterMap(hdn);    // need to check for new parameters
                                             // from within If then else
        hdn->ResetChildCursor();       // Build simulation param map moves child
                                       // node of hdn
        continue;      // will restart loop
        // FIXME - this will not catch
        // <if condition = "1" >
        // <then>
        // param = x
        // </then>
        // </if>
        // param = y
        //
        // i.e. param will = y as Build simulation parameter map not stopped by
        // if statement
        // could fix by forcing build sim param map to stop at if, but seems
        // clunky
        // ideally build sim param map and evaluate if then else should be in
        // same routine
      }

    }

    // Replace parameters outside of if conditions
    if(childNode) // node may have beem deleted after evaluating if statement
    {
      std::string heading = childNode->Heading();
      bool needsReplacement = ReplaceParameters(heading, m_simulationParameterMap);
      if (needsReplacement)
        childNode->SetHeading(heading.c_str());

      // attributes
      std::string value;
      std::string key;
      array<string> newKeys;
      array<string> newValues;
      childNode->ResetCursors();
      while (childNode->Next(key, value))
      {
        needsReplacement = ReplaceParameters(key, m_simulationParameterMap);
        needsReplacement += ReplaceParameters(value, m_simulationParameterMap);
        needsReplacement += ReplaceMathematicalExpressions(value);
        if (needsReplacement)
        { // nb value already updated if needed
          newKeys.push_back(key);
          newValues.push_back(value);
        }
      }
      for (array<string>::size_type i = 0 ; i < newKeys.size() ; ++i)
        childNode->AddAttributePair(newKeys[i], newValues[i]);

      // depth first traversal
      ParseMetadata(childNode, false);
    }
  }

}

/*
   bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
   }
 */

/**
 * EvaluateXMLIfThenElseStatement
 *
 * @author walsh24
 *
 *  <If condition="1" >
 *    <Then>
 *      <NodeA> fu fu fu </NodeA>
 *      <NodeB> bar bar bar </NodeB>
 *    </Then>
 *    <Else>
 *      <NodeC> fubar fubar fubar </NodeC>
 *    </Else>
 *  </If>
 *
 *  is equivalent to:
 *  <NodeA> fu fu fu </NodeA>
 *  <NodeB> bar bar bar </NodeB>
 *
 *  If statements can be nested and can contain parameters and expressions
 *
 **/
bool ProblemManagerT::EvaluateXMLIfThenElseStatement(HierarchicalDataNode* &hdn, HierarchicalDataNode* parentNode){

  bool isIfStatement =  ( hdn != 0) && streq(hdn->Heading(), "If");
  bool reevaluateNode = false;

  if(isIfStatement)
  {

    int ifLevel = hdn->Level();

    // replace parameters in condition string -> may not work for nested if
    // statements containing parameters
    std::string conditionString = hdn->GetAttributeString("condition");
    //bool dummy =
    ReplaceParameters(conditionString, m_simulationParameterMap);
    hdn->AddAttributePair("condition", conditionString);

    // evaluate condition
    bool conditionIsTrue =  hdn->GetAttributeOrDefault<bool>("condition",false);

    HierarchicalDataNode* thenNode = NULL;
    HierarchicalDataNode* elseNode = NULL;


    for (HierarchicalDataNode* childNode = hdn->Next(true) ; childNode ; childNode = hdn->Next())
    {
      if( streq(childNode->Heading(), "Then") )
      {
        thenNode = childNode;
      }
      else if( streq(childNode->Heading(), "Else") )
      {
        elseNode = childNode;
      }
      else
      {
        throw GPException( std::string("Unsupported node type ") + childNode->Heading() + " in If xml statment" );
      }
    }

    std::vector<HierarchicalDataNode> hdnv; // vector of nodes in then or else -
                                            // replaces if statement
    if(conditionIsTrue)
    {
      if(thenNode)
      {
        hdnv = thenNode->CopyChildren();
      }

    }
    else
    {
      if(elseNode)
      {
        hdnv = elseNode->CopyChildren();
      }
    }
    // reset level
    for (size_t i = 0 ; i < hdnv.size() ; ++i)
    {
      hdnv[i].SetLevel(ifLevel);
    }

    // replace if statement with vector of nodes (or empty node)
    if (hdnv.size() > 0)
    {
      (*hdn) = hdnv[0];
      for (size_t i = 1 ; i < hdnv.size() ; ++i)
      {
        if (hdnv.size() > 1)
        {
          int indx = parentNode->InsertAfter(hdn, hdnv.begin() + 1, hdnv.end());
          // recover pointer invalidated by InsertAfter
          indx -= 1; // index refers to index where insert starts - want the
                     // start of the hdnv instead
          int cc = 0;
          hdn = parentNode->Next(true);
          while (cc < indx)
          {
            hdn = parentNode->Next();
            ++cc;
          }
        }
      }
    }
    else     // hdnv.size() == 0
    { // (*hdn) = HierarchicalDataNode(); // empty node - ideally would want to
      // remove node in case surrounding environment does not li, but would need
      // to make sure that next child of parent was also evaluated.
      int indx = parentNode->RemoveChild(hdn);

      // recover pointer invalidated by InsertAfter and point to next node
      int cc = 0;
      hdn = parentNode->Next(true);
      while (cc < indx)
      {
        hdn = parentNode->Next();
        ++cc;
      }

    }

    reevaluateNode = true; // need to come back for recursive if statements.
  } // is if statement


  return reevaluateNode;
}

/**
 * ReplaceMathematicalExpressions
 *
 * @author walsh24
 *
 * @param[in/out] valueStr - a string containing real expressions enclosed in
 * grave accents (`) eg. `sin(3.14159265/4)`
 *
 * Parses mathematical expressions enclosed in grave accents and replaces the
 * substring with the result.
 * Returns true if a replacement is performed.
 */
bool ProblemManagerT::ReplaceMathematicalExpressions(std::string& valueStr)
{
  bool rv = false;

  size_t startIndx = valueStr.find('`');

  while (startIndx < valueStr.size())
  {
    rv = true;
    size_t endIndx = valueStr.find('`', startIndx + 1);

    std::string replaceStr = valueStr.substr(startIndx + 1, endIndx - startIndx - 1);
    replaceStr = toString(EvaluateStringFunction(replaceStr));

    valueStr.replace(startIndx, endIndx - startIndx + 1, replaceStr);

    startIndx = valueStr.find('`');
  }
  return rv;
}

void ProblemManagerT::ReadGeometryInput( HierarchicalDataNode& hdn )
{
  bool hasExternalMeshFile = true;
//  bool hasMetis = false;

  // Meshes
  ///////////
  {
    // may also be specified using -m or -f command-line options
    HierarchicalDataNode* meshNode = hdn.GetChild("Mesh");
    if (meshNode != NULL)
    {
      hasExternalMeshFile = m_FileManager.ReadMeshXML(meshNode);
      if( !hasExternalMeshFile )
      {
        m_MeshGenerator = &myMeshGenerator;
        m_MeshGenerator->ReadXML( *meshNode );
      }
      //else delete m_MeshGenerator;
      else
      {
        m_MeshGenerator = NULL;
      }
    }
  }

  {
    // Cartesian Grid
    HierarchicalDataNode* cgNode = hdn.GetChild("CartesianGrid");
    if (cgNode != NULL)
    {
      m_Domains.m_cartesianGridManager.ReadXML(cgNode);
    }
  }

  // Partitions
  //////////////
  {
    HierarchicalDataNode* PartitionNode = hdn.GetChild("Partition");
    if (PartitionNode != NULL)
    {
      HierarchicalDataNode* SpatialPartitionNode = PartitionNode->GetChild("SpatialPartition");
      if (SpatialPartitionNode != NULL)
      {
        m_useMetis = SpatialPartitionNode->HasAttribute("Metispar");
        m_partition.ReadXMLInput(*SpatialPartitionNode);
      }
    }

    // default to 1 partition unless already set
    if(m_partition.m_size  == 0)
    {
      m_partition.setPartitions( 1,1,1 );
      m_partition.setPeriodic( 0,0,0 );
    }

    // set radial partition boundary
    if( !hasExternalMeshFile )
    {
      if((m_MeshGenerator->isRadial()) && (m_partition.GetPartitions()[1]>2))
      {
        m_partition.setRadialPeriodic( 1 );
      }
    }
    if (m_useMetis)
    {
      m_partition.InitializeMetis();
    }
    else
    {
      m_partition.Initialize();
      m_partition.SetRankOfNeighborNeighbors();
      m_partition.SetDomain(m_Domains);
    }
  }

  if( !hasExternalMeshFile )
    m_MeshGenerator->GenerateElementRegions( m_Domains );


  bool has_fe = false;
  if( hasExternalMeshFile )  //!m_beginFromRestart && Fu note: this condition
                             // was removed because we need to read the mesh
                             // file to know whether this is a 2d or 3d problem
                             // in order to determine which fractunator to use.
  {
    has_fe = m_FileManager.ExtractElementRegionFromAbaqus( m_Domains, m_partition );
  }
  else
  {
    has_fe = true;
  }


  // Element regions
  //////////////////
  if(has_fe)
  {
    HierarchicalDataNode* regionsNode = hdn.GetChild("ElementRegions");
    if (regionsNode != NULL)
    {
      for (HierarchicalDataNode* erNode = regionsNode->Next(true) ; erNode ; erNode = regionsNode->Next())
      {
        //const std::string erName = erNode->Heading();
        const std::string erName = erNode->GetAttributeString("name");
        if(erName.length()==0)
        {
          throw GPException("Cannot have an ElementRegions child that does not have a name");
        }


        ElementRegionT* elemRegion = stlMapLookupPointer(m_Domains.m_feElementManager.m_ElementRegions,erName);

        if( elemRegion == NULL )
        {
          if( !m_beginFromRestart &&  erName != "FlowFaceRegion")
          {
            char msg[200];
            sprintf(msg, "ProblemManagerT::ReadGeometryInput(): region specified in XML does not exist in geometry: %s", erName.c_str());
            throw GPException(msg);
          }
          else if (!m_beginFromRestart && erName == "FlowFaceRegion")
          {
            m_Domains.m_feElementManager.InitializeFlowFaceRegion();
            elemRegion = stlMapLookupPointer(m_Domains.m_feElementManager.m_ElementRegions,erName);
          }
          else
          {
            m_Domains.m_feElementManager.m_ElementRegions[erName];
            elemRegion = stlMapLookupPointer(m_Domains.m_feElementManager.m_ElementRegions,erName);
          }
        }
        elemRegion->ReadXML( erNode, m_beginFromRestart );

        if (m_rank  == 0)
          std::cout << "    Element Region " << erName << " contains " << elemRegion->m_elementGeometryID
                    << " element types with " << elemRegion->m_elementType << " integration." << std::endl;
      }
    }
  }



  // Ellipsoidal DEM Boundary Contact
  /////////////////////////////////////
  //need to read the ellipsoidal DEM boundary contact model XML prior to any
  // resizing
//  {
//    HierarchicalDataNode* EDemNode = hdn.GetChild("EllipsoidalDEM");
//    if (EDemNode != NULL)
//    {
//      //set the boundary contact
//      HierarchicalDataNode* EDemNodeContact0 = EDemNode->GetChild("Contact");
//      if(EDemNodeContact0 == NULL)
//        throw GPException("Must define a boundary contact model for
// Elliposidal DEM");
//      HierarchicalDataNode* EDemNodeContact1 = EDemNodeContact0->Next(true);
//      if(EDemNodeContact1 == NULL)
//        throw GPException("Must define a boundary contact model for
// Elliposidal DEM (1)");
//      Interfaces::Allocate(m_Domains.m_ellipsoidalDiscreteElementManager.m_boundaryContact,
//                           Interfaces::StringToModelType(EDemNodeContact1->Heading()));
//      m_Domains.m_ellipsoidalDiscreteElementManager.m_boundaryContact->resize(0,1);
//      m_Domains.m_ellipsoidalDiscreteElementManager.m_boundaryContact->ReadXML(*EDemNodeContact1);
//    }
//  }



  if( !m_beginFromRestart )
  {
    if( hasExternalMeshFile )
    {
      //R.W. use Metis related functions if Metis is specified
      if (m_useMetis)
      {
        m_FileManager.ReadMeshforMetis(m_Domains, m_partition);
      }
      else
      {
        m_FileManager.ReadMesh(m_Domains, m_partition);
      }
    }
    else
    {
      m_MeshGenerator->GenerateMesh( m_partition,
                                     m_Domains );
    }

    HierarchicalDataNode* nodesetNode = hdn.GetChild("Nodesets");

    if( nodesetNode != NULL )
      MeshUtilities::GenerateNodesets( *nodesetNode, m_Domains.m_feNodeManager );


    HierarchicalDataNode* elementsetNode = hdn.GetChild("Elementsets");
    if( elementsetNode != NULL )
      MeshUtilities::GenerateElementsets( *elementsetNode, m_Domains.m_feNodeManager, m_Domains.m_feElementManager );
    //  std::cout<<"Initializing Domain...\n";

    // Fu: I am moving this to after faces are built.  We can only build the
    // flowFaceRegion after faces are build.  We need to call this function
    // then.
    //m_Domains.Initialize();

    //  std::cout<<"Done Initializing Domain.\n";
  }

}

/**
 * @brief Read XML parameters
 * @param[in,out] hdn Hierarchical data node object
 */
void ProblemManagerT::ReadXML(HierarchicalDataNode& hdn)
{

  // Units
  //////////////////////
  HierarchicalDataNode* unitsNode = hdn.GetChild("Units");
  UnitManager& unitManager = UnitManager::Instance();
  if (unitsNode)
  {
    if (m_rank == 0)
      std::cout << "Setting Default Units:" << std::endl;


    unitManager.ReadXML(unitsNode);

  }


  // Tables
  //////////
  HierarchicalDataNode* TablesNode = hdn.GetChild("Tables");
  TableManager& tableManager = TableManager::Instance();
  if (TablesNode != NULL)
  {
    if (m_rank  == 0)
      std::cout << "Reading Tables:" <<std::endl;
    tableManager.ReadXML(TablesNode);
  }

  // Functions
  /////////////
  HierarchicalDataNode* functionsNode = hdn.GetChild("Functions");
  FunctionManager& functionManager = FunctionManager::Instance();
  if (functionsNode)
  {
    if (m_rank  == 0)
      std::cout << "Reading Functions:" << std::endl;
    for (HierarchicalDataNode* fNode = functionsNode->Next(true) ; fNode ; fNode
           = functionsNode->Next())
    {
      std::string fType = fNode->Heading();
      std::string fName = fNode->GetAttributeString("name");
      //        std::cout << "    Function: " << fType << std::endl;
      //functionManager.Functions().insert(std::make_pair(fName,
      // newFunction(fType, fNode, this)));
      functionManager.AddFunction(fName, newFunction(fType, fNode, this));
    }
  }


#ifdef USE_CHEM
  // Chemistry
  //////////////
  HierarchicalDataNode* ChemNode = hdn.GetChild("ChemistryManager");
  if (ChemNode != NULL)
  {
    if (m_rank  == 0)
      std::cout << "Reading Chemistry Manager:" << std::endl;
    ChemistryManager& chemistryManager = ChemistryManager::Instance();
    chemistryManager.ReadXML(ChemNode);
  }
#endif


#ifdef SRC_INTERNAL
  #ifdef GPAC_GMM
  // Geodyn material model library
  //////////////
  HierarchicalDataNode* GMML_Node = hdn.GetChild("GeodynMaterialModelManager");
  if (GMML_Node != NULL)
  {
    std::cout << "Reading Geodyn Material Model Manager:" << std::endl;
    GeodynLibraryManager& theGMMLibManager = GeodynLibraryManager::Instance();
    theGMMLibManager.ReadXML(GMML_Node);

    theGMMLibManager.Initialize(); // must be initialized before materials are
                                   // requested
  }


  #endif
#endif


  // Meshes & Materials (inside element regions)
  ///////////////////////////////////////////////
  ReadGeometryInput(hdn);



  // Solvers
  ///////////
  HierarchicalDataNode* solversNode = hdn.GetChild("Solvers");
  {
    if (!solversNode)
      throw GPException("Problem Manager: Must have Solvers defined in the input file");

    if (m_rank  == 0)
      std::cout << "Initializing Solvers:" << std::endl;
    for (HierarchicalDataNode* solverNode = solversNode->Next(true) ; solverNode ; solverNode
           = solversNode->Next())
    {
      std::string solverType = solverNode->Heading();
      std::string solverName = solverNode->GetAttributeString("name");
      if(solverName.length() == 0)
        throw GPException("Solver must have a non-zero-length name");
      if (m_rank  == 0)
        std::cout << "    Solver: " << solverType << " Name: " << solverName << std::endl;
      m_solvers[solverName] = SolverFactory::NewSolver(solverType, solverNode, this);
      //if(!m_solvers[solverName])
      //  throw GPException("Solver not properly initialized");
      //m_solvers[solverName]->ReadXML(solverNode);
    }
  }

  // Solver applications
  ///////////////////////
  {
    HierarchicalDataNode* applicationsNode = hdn.GetChild("SolverApplications");
    if (!applicationsNode)
    {
      //throw GPException("Must have SolverApplications defined in the input
      // file");
      std::cout << "Warning: SolverApplications block not found. Solvers will be intialized but not run." << std::endl;
    }
    else
    {

      if (m_rank  == 0)
        std::cout << "Initializing Solver Applications:" << std::endl;

      realT endTimeOfLastApplication = 0.0;
      for (HierarchicalDataNode* applicationNode = applicationsNode->Next(true) ; applicationNode ; applicationNode
             = applicationsNode->Next())
      {

        //application node attributes
        const std::string name = applicationNode->GetAttributeString("name");
        const realT beginTime = applicationNode->GetAttributeValue<realT> ("begintime");
        const realT endTime = applicationNode->GetAttributeValue<realT> ("endtime");
        const realT dt = applicationNode->GetAttributeOrDefault<realT>("dt",std::numeric_limits<double>::max());

        if ( !isEqual(endTimeOfLastApplication,beginTime) )
          throw GPException("Solver application times are not contiguous\n");

        if (m_rank  == 0)
          std::cout << "   Application: " << name << " Time: " << beginTime << " to " << endTime
                    << std::endl;

        endTimeOfLastApplication = endTime;

        SolverApplicationSet currentSet;
        currentSet.m_beginTime = beginTime;
        currentSet.m_endTime = endTime;
        currentSet.m_deltaTime = dt;

        //get apply nodes and substeps
        HierarchicalDataNode* applyNode = 0;
        for (applyNode = applicationNode->Next(true) ; applyNode ; applyNode = applicationNode->Next())
        {

          const std::string& apType = applyNode->Heading();
          if( streq(apType,"Substep") )
          {
            // add new substep solver to the solver list
            std::string solverName = "SubstepSolver_" + toString(SubstepSolver::NumberOfInstances() );
            if(!applyNode->HasAttribute("name") )
            {
              applyNode->AddAttributePair("name",solverName);
            }
            m_solvers[solverName] = SolverFactory::NewSolver(SubstepSolver::SolverName(), applyNode,this);
            array<string> applyToRegionNames;

            currentSet.m_solverAppliedToRegion.push_back(SolverApplication(solverName,
                                                                           applyToRegionNames));

          }
          else
          {
            std::string solverName = applyNode->GetAttributeString("solver");
            std::string toRegionsNames = applyNode->GetAttributeString("toregions");

            if (m_rank  == 0)
              std::cout << "       Applying Solver: " << solverName << " to regions: "
                        << toRegionsNames << std::endl;

            std::istringstream iss(toRegionsNames);
            array<string> applyToRegionNames;

            copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                 std::back_inserter<array<string> >(applyToRegionNames));

            currentSet.m_solverAppliedToRegion.push_back(SolverApplication(solverName,
                                                                           applyToRegionNames));

          }

        }

        m_solverApplicationSets.push_back(currentSet);

      }
    }
  }



  // Boundary conditions
  //////////////////////
  HierarchicalDataNode* BoundaryConditionsNode = hdn.GetChild("BoundaryConditions");
  if (BoundaryConditionsNode != NULL)
  {
    //      std::cout << "Reading BoundaryConditions:" <<std::endl;
    for (HierarchicalDataNode* bcNode = BoundaryConditionsNode->Next(true) ; bcNode ; bcNode
           = BoundaryConditionsNode->Next())
    {

      //  std::string bcType = BCNode->Heading();
      //        std::cout << "    Boundary condition: " << bcType << std::endl;

      std::string bcType = bcNode->Heading();

      BoundaryConditionBase* bcPtr = newBoundaryCondition(bcType, bcNode, this);
      //ObjectDataStructureBaseT& objectManager =
      // m_Domains.GetObjectDataStructure(bcPtr->m_objectKey,"");
      ObjectDataStructureBaseT& objectManager = m_Domains.GetObjectDataStructure(bcPtr->m_objectKey,bcPtr->m_regionName);
      objectManager.m_bcData.push_back(bcPtr);

    }
  }

  // Initial conditions
  //////////////////////
  HierarchicalDataNode* initialConditionsNode = hdn.GetChild("InitialConditions");
  if (initialConditionsNode != NULL)
  {
    if (m_rank  == 0)
      std::cout << "Reading InitialConditions:" << std::endl;
    for (HierarchicalDataNode* icNode = initialConditionsNode->Next(true) ; icNode ; icNode
           = initialConditionsNode->Next())
    {
      std::string icType = icNode->Heading();
      //        std::cout << "    Initial Condition: " << icType << std::endl;
      m_initialConditions.push_back(newInitialCondition(icType, icNode, this));
    }
  }


  // Contact
  //////////////
  HierarchicalDataNode* ContactNode = hdn.GetChild("Contact");
  if (ContactNode != NULL)
  {
    m_Domains.m_externalFaces.ReadXML(ContactNode);
    m_Domains.m_contactManager.ReadXML(ContactNode);
  }


  // Ellipsoidal DEM
  //////////////
  HierarchicalDataNode* EDemNode = hdn.GetChild("EllipsoidalDEM");
  if (EDemNode != NULL)
  {
    m_Domains.m_ellipsoidalDiscreteElementManager.ReadXML(EDemNode);
    HierarchicalDataNode* EDemNodeContact = hdn.GetChild("Contact");
    if(EDemNodeContact == NULL)
      throw GPException("If you define an ellipsoidal DEM, you must also define a contact model (Contact)");
    m_Domains.m_ellipsoidalContactManager.ReadXML(EDemNodeContact);
  }


  // Joint sets
  //////////////
  #ifdef SRC_INTERNAL2
  HierarchicalDataNode* JointSetsNode = hdn.GetChild("JointSets");
  if (JointSetsNode != NULL)
  {
    m_Domains.m_jointSets2.ReadXML(JointSetsNode);
  }
  #endif

  // HierarchicalDataNode* OldJointSetsNode = hdn.GetChild("OldJointSets");
  // if (OldJointSetsNode != NULL)
  // {
  //   m_Domains.m_jointSets.ReadXML(OldJointSetsNode);
  //   m_Domains.m_jointSets.Populate(m_Domains.m_feElementManager.m_ElementRegions,
  //                                  m_Domains.m_feFaceManager,
  // m_Domains.m_feNodeManager);
  // }

#ifdef SRC_EXTERNAL
  // Wellbore mesh
  HierarchicalDataNode* WellboreNode = hdn.GetChild("WellboreMesh");
  if (WellboreNode != NULL)
  {
    m_Domains.m_wellboreManager.ReadXML(WellboreNode);
  }

  // Fault Elements
  ////////////
  HierarchicalDataNode* RiskNode = solversNode->GetChild("SeismicRiskSolver");
  {
    HierarchicalDataNode* FaultRuptureBEMNode = solversNode->GetChild("FaultRuptureBEMSolver");
    if (FaultRuptureBEMNode != NULL)
      m_Domains.m_faultElementManager.ReadXML(tableManager, FaultRuptureBEMNode);
    if (RiskNode != NULL)
      m_Domains.m_faultElementManager.ReadXML(tableManager, RiskNode);
  }
#endif
  // Fracture
  //////////////
  HierarchicalDataNode* FractureNode = hdn.GetChild("Fracture");
  if (FractureNode != NULL)
  {
    m_fractureFlag = FractureNode->GetAttributeOrDefault<int> ("fractureFlag", 0);

#ifdef SRC_INTERNAL2
    m_xfem = FractureNode->GetAttributeOrDefault<int> ("XFEM", 0);

    if(m_xfem != 0)
    {
      m_Domains.m_xfemManager = new XfemManager();
      m_Domains.m_crackManager = new CrackObjectManager();
      m_Domains.m_crackSurfaceVertex = new CrackSurfaceVertex();

      m_Domains.m_xfemManager->m_xfem = m_xfem;
      m_Domains.m_xfemManager->ReadXML(FractureNode, m_Domains);
//      m_elementSplitting = new XfemManager();
//      m_elementSplitting->ReadXML(FractureNode, m_Domains);
    }
#endif

    if( m_fractureFlag != 0 )
    {
      m_preFractureSets = FractureNode->GetStringVector("preFractureSetName");

      if( m_Domains.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension == 3 )
      {
        m_surfaceSeparation = FractunatorFactory::NewFractunator("Fractunator3", FractureNode);
      }
      else if( m_Domains.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension == 2 )
      {
        m_surfaceSeparation = FractunatorFactory::NewFractunator("Fractunator2D", FractureNode);
//        if(m_xfem != 0)
//        {
//          m_elementSplitting = new XfemManager();
//          m_elementSplitting->ReadXML(FractureNode, m_Domains);
//        }
      }
      else
      {
        throw GPException("Cannot find a valid fractunator!");
      }

      if(!FractureNode->GetAttributeString("type").empty())
        m_Domains.m_feFaceManager.m_cohesiveZone = CohesiveZoneFactory::NewCohesiveZone( FractureNode->GetAttributeString("type"), FractureNode );
    }
  }

  // HACKS
  //////////////
  HierarchicalDataNode* HackNode = hdn.GetChild("HACK");
  if (HackNode != NULL)
  {
    m_hackInitialStress = HackNode->GetAttributeOrDefault<int> ("hackInitStress", 0);
  }

  // Initial constitutive
  //////////////////////
  HierarchicalDataNode* initialConstitutiveNode = hdn.GetChild("InitialConstitutiveValues");
  if (initialConstitutiveNode != NULL)
  {
    if (m_rank  == 0)
      std::cout << "Reading InitialConstitutiveValues:" << std::endl;
    for (HierarchicalDataNode* icNode = initialConstitutiveNode->Next(true) ; icNode ; icNode
           = initialConstitutiveNode->Next())
    {
      const std::string tableName = icNode->GetAttributeString("tablename");
      realT value = 0.0;

      //check to make sure the table reference is valid
      if(!tableName.empty())
      {
        std::map<std::string,Table3D >::const_iterator it = tableManager.Tables<3>().find(tableName);
        if(it == tableManager.Tables<3>().end())
          throw GPException("Cannot find requested table in the table manager while attempting to set initial constitutive properties: " + tableName);
      }
      else
      {
        value = icNode->GetAttributeValue<realT>("value");
      }

      //get the ods
      PhysicalDomainT::ObjectDataStructureKeys key = m_Domains.GetObjectDataStructureKey(
        icNode->GetAttributeString("object"));

      //add a new table
      ObjectDataStructureBaseT& ods = m_Domains.GetObjectDataStructure(key, icNode->GetAttributeString("toregion"));

      const std::string fieldName = icNode->GetAttributeString("propertytype");
      if(fieldName.empty())
        throw GPException("Must specify a field name in the InitialConstitutiveValues node");

      NodeManager* nm = 0;
      if(key == PhysicalDomainT::FiniteElementElementRegion)
        nm = &m_Domains.m_feNodeManager;
      if(!tableName.empty())
        m_initialConstitutive.Add(fieldName, tableName, ods, nm);
      else
        m_initialConstitutive.Add(fieldName, value, ods, nm);
    }
#ifdef SRC_EXTERNAL
    if (RiskNode != NULL)
    {
      //because the RiskNode has to reset its states every aleatoric
      // realization,
      //we need to be able to cache the initial constitutive models in the fault
      //element manager
      m_Domains.m_faultElementManager.SetInitialConstitutive(m_initialConstitutive);
    }
#endif
  }


  // Output formatting
  ////////////////////
  {
    HierarchicalDataNode* output = hdn.GetChild("Output");
    if (!output)
      throw GPException("Must have Output defined in the input file");
    this->m_siloFile.m_numGroups = output->GetAttributeOrDefault<int> ("parallel_silo", 1);
    std::string temp = output->GetAttributeString("plotfile_root");
    if (!temp.empty())
      this->m_siloFile.m_fileRoot = temp;

    temp = output->GetAttributeString("restartfile_root");
    if (!temp.empty())
      this->m_siloFile.m_restartFileRoot = temp;

    temp = output->GetAttributeString("slave_directory");
    if (!temp.empty())
    {
      this->m_siloFile.m_slaveDirectory = temp;
      if (this->m_siloFile.m_numGroups > 1 && m_rank  == 0)
      {
        std::cout << "Creating the slave directory. System call returned: " << std::endl << system(("mkdir -p " + temp).c_str()) << std::endl;
        //This only works on Unix-like systems but should not harm a windows run
        // anyway.  It just won't be able to create the folder
      }
    }

    this->m_siloFile.m_markGhosts = output->GetAttributeOrDefault<int> ("markGhosts", 0);

    m_writeFEM = output->GetAttributeOrDefault<int> ("writeFEM", 1);
    m_writeFEMFaces = output->GetAttributeOrDefault<int> ("writeFEMFaces", 0);
    m_writeFEMEdges = output->GetAttributeOrDefault<int> ("writeFEMEdges", 0);
    m_writeFlowText = output->GetAttributeOrDefault<int> ("writeFlowText", 0);

    m_visitFileGroupFile = output->GetAttributeStringOrDefault( "visitGroupFile", "geos.visit");
    m_writePlot = output->GetAttributeOrDefault<bool>("writePlot",true);
    m_plotIntervals = output->GetAttributeOrDefault<realT> ("plot_interval",std::numeric_limits<realT>::max());
    m_plotIntervalTableName = output->GetAttributeString( "plotIntervalTable" );


    m_writeRestart = output->GetAttributeOrDefault<bool>("writeRestart",true);
    m_restartIntervals = output->GetAttributeOrDefault<realT> ("restart_interval",std::numeric_limits<realT>::max());
    m_restartIntervalTableName = output->GetAttributeString( "restartIntervalTable" );

    m_cycleReportFreq = output->GetAttributeOrDefault<int>("cycle_report_interval",1);


    m_nameFieldsToPlot = output->GetStringVector("fieldsToPlot");
    m_nameAdditionalFieldsToPlot = output->GetStringVector("additionalFieldsToPlot");

    m_maxWallTime = output->GetAttributeOrDefault<realT> ("maxWallTime",std::numeric_limits<realT>::max());
    m_forceOutputTime = output->GetAttributeOrDefault<bool>("forceOutputTime", false);

//    if( isEqual(m_restartIntervals,std::numeric_limits<realT>::max()) ){
//      m_nextRestartTime = m_restartIntervals;  // do this anyway - do we need
// to write a restart at time 0?
//    }

    for (HierarchicalDataNode* thNode = output->Next(true) ; thNode ; thNode = output->Next())
    {
      realT vval[nsdof];
      HierarchicalDataNode::StrVal<realT>(thNode->GetAttributeString("location"), vval, nsdof);// fixme
                                                                                               // -
                                                                                               // does
                                                                                               // nothing
      const std::string name = thNode->GetAttributeString("fieldname");// fixme
                                                                       // - does
                                                                       // nothing
    }

    m_trackEnergy = output->GetAttributeOrDefault("trackEnergy", m_trackEnergy);
  }



}

void ProblemManagerT::WriteXML(std::string& filename, HierarchicalDataNode& hdn){

  std::ofstream fid;
  fid.open( filename.c_str() );

  fid << hdn.ToXML();

  fid.close();


}

void ProblemManagerT::ParseCommandLineInput(const int& argc, char* const argv[])
throw (GPException)
{
  // Default values
  std::string fileRootString = "";
  std::string meshFileString, demeshFileString, edemeshFileString, inputFileString, fpmeshFileString;
  int displaySolvers_flag = 0;
  int displayFields_flag = 0;
  int displayUnits_flag = 0;
  int displayHistory_flag = 0;
  int reportParameter_flag = 0;

  int defaultVariableReportLevel = 0;

  bool setPartitions_flag = false;
#ifdef SRC_INTERNAL
  bool geodyn_pure =  true;
#endif
  unsigned int xPartitions = 1;
  unsigned int yPartitions = 1;
  unsigned int zPartitions = 1;

  array<string> commandLineIncludedFileList;

  // Get command line input
  while (true)
  {
    static struct option long_options[] =
    {
      { "help", no_argument, 0, 'h' },
      { "version", no_argument, 0, 'v' },
      { "solvers", no_argument, &displaySolvers_flag, 1 },
      { "fields", no_argument, &displayFields_flag, 1 },
      { "units", no_argument, &displayUnits_flag, 1 },
      { "history", no_argument, &displayHistory_flag, 1 },
      { "record_defaults", no_argument, &defaultVariableReportLevel, 1 },
      { "report_defaults", no_argument, &defaultVariableReportLevel, 2 },
      { "disable_defaults", no_argument, &defaultVariableReportLevel, 3 },
      { "report_parameters",no_argument,&reportParameter_flag,1},
      { "xpar", required_argument, 0, 0 },
      { "ypar", required_argument, 0, 0 },
      { "zpar", required_argument, 0, 0 },
      { "include", required_argument, 0, 0 },
      { "write_XML", required_argument, 0, 0 },
      { 0, 0, 0, 0 } };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long_only(argc, argv, "ahvf:p:m:d:i:e:r:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

#ifdef SRC_INTERNAL
    geodyn_pure = false; // a geos option has been detected
#endif

    switch (c)
    {
    case 0:
    {
      /* If option sets a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;

      /* long options without a short arg */
      if( streq( std::string("xpar"), long_options[option_index].name ) )
      {
        xPartitions = fromString<unsigned int>(optarg);
        setPartitions_flag = true;
      }
      else if( streq( std::string("ypar"), long_options[option_index].name ) )
      {
        yPartitions = fromString<int>(optarg);
        setPartitions_flag = true;
      }
      else if( streq( std::string("zpar"), long_options[option_index].name ) )
      {
        zPartitions = fromString<int>(optarg);
        setPartitions_flag = true;
      }
      else if( streq( std::string("include"), long_options[option_index].name ) )
      {
        commandLineIncludedFileList.push_back(optarg);
      }
      else if( streq( std::string("write_XML"), long_options[option_index].name ) )
      {
        m_doWriteXML = true;
        m_xmlOutputFileName = optarg;
      }

    }
    break;
    case 'a':   // Leave Empty: Included for totalview - does nothing
      break;
    case 'f':   // Record file root
    {
      fileRootString = optarg;
    }
    break;

    case 'i':   // Record input file
    {
      inputFileString = optarg;
    }
    break;

    case 'm':   // Record mesh file
    {
      meshFileString = optarg;
    }
    break;
    case 'd':   // Record discrete element mesh file
    {
      demeshFileString = optarg;
    }
    break;
    case 'e':   // Record ellipsoidal discrete element mesh file
    {
      edemeshFileString = optarg;
    }
    break;

    case 's':   // Record seismicity fault patch mesh file
    {
      fpmeshFileString = optarg;
    }
    break;

    case 'r':   // From restart
    {
      m_beginFromRestart = true;
      m_beginFromRestartFileName = optarg;
    }
    break;

    case 'p':   // Record model parameter key=value
    {
      std::string keyValStr = optarg;
      array<string> keyVal = Tokenize(keyValStr, "=");
      if (keyVal.size() == 2)
      {
        m_simulationParameterMap[keyVal[0]] = keyVal[1];
      }
      else
      {
        throw GPException("Error reading command line input: Parameter " + keyValStr);
      }
    }
    break;

    case 'h':   // help
      DisplayUsage();   // print help
      exit(0);
      break;

    case 'v':   // version
      DisplayVersion();
      exit(0);
      break;

    case '?':
      /* getopt_long has already printed an error message. */
      break;

    default:
      abort();
      break;
    }
  }

#ifdef SRC_INTERNAL
  if (p_initbamr && argc!=optind)
  {
    doBackgroundAMR=geodyn_pure ? -1 : 1;
    int argcg=argc-(optind-1); // passing wrong exe name (the last option)
    char** argvg=const_cast<char**>(argv+(optind-1));
    int ret=p_initbamr( argcg, argvg );
    if (ret)
      throw GPException("Errors during BackgroundAMR init");
  }
#endif

  if (displaySolvers_flag)
  {
    DisplaySolvers();
    exit(0);
  }


  if (displayUnits_flag)
  {
    DisplayUnits();
    exit(0);
  }


  if (displayHistory_flag)
  {
    DisplayVersionHistory(m_rank);
    exit(0);
  }

  if (displayFields_flag)
  {
    m_displayFields = true; // nb we want to initialize solvers, BC's initial
                            // conditions etc before reporting fields.
    m_displaySplash = false;
  }

  // this option concerns flags for the default values of variables
  switch(defaultVariableReportLevel)
  {
  case 0:
    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::silent);
    break;
  case 1:
    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::recordDefaults);
    break;
  case 2:
    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::reportDefaults);
    break;
  case 3:
    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::disableDefaults);
    break;
  }

  // this option sets the flag to report the current value of parameters
  if(reportParameter_flag)
  {
    m_echoParameters = true;
  }

  if(setPartitions_flag)
    m_partition.setPartitions(xPartitions,yPartitions,zPartitions );

  if (!fileRootString.empty())
    m_FileManager.SetRoot(fileRootString.c_str());
  if (!inputFileString.empty())
    m_FileManager.SetInputFilename(inputFileString.c_str());
  if (!meshFileString.empty())
    m_FileManager.SetGeometryFilename(meshFileString.c_str());
  if (!demeshFileString.empty())
    m_FileManager.SetDiscreteElementGeometryFilename(demeshFileString.c_str());
  if (!edemeshFileString.empty())
    m_FileManager.SetEllipsoidalDiscreteElementGeometryFilename(edemeshFileString.c_str());
#ifdef SRC_EXTERNAL
  if (!fpmeshFileString.empty())
    m_FileManager.SetFaultPatchElementGeometryFilename(fpmeshFileString.c_str());
#endif
  //////////////////////////////////////////////////////////////////////////////

#ifdef SRC_INTERNAL

  // Geodyn wrapper

  if(geodyn_pure && p_stepbamr)
  {
    // let GEOS act as a geodyn wrapper if regular input is not defined
    A2Gtype A2G;
    G2Atype G2A;
    for ( int done=0 ; !done ; done=p_stepbamr(A2G,G2A) )
      ;
    exit(p_exitbamr(1));
  }
#endif

  //////////////////////////////////////////////////////////////////////////////
  //this reads the given file into an intermediate
  //data structure (hdn), which can be interrogated for
  //values


  TICPP::HierarchicalDataNode hdn;
  RegisterFilesIncludedFromCommandLine(&hdn, commandLineIncludedFileList);
  bool success = m_FileManager.ReadXMLInput(hdn);
  if (success)
  {
    ParseIncludedFiles(&hdn);
    BuildSimulationParameterMap(&hdn);
    ParseMetadata(&hdn,true);
    ReadXML(hdn);

    if(m_doWriteXML)
    {
      WriteXML(m_xmlOutputFileName,hdn);
    }
  }
  else
  {
    std::cout << "Error reading xml file" << std::endl;
    exit(0);
  }


}

void ProblemManagerT::ProblemSetup()
{


#ifdef SRC_INTERNAL
  if (doBackgroundAMR==-1)
    return;
#endif

  RegisterFields();
  if( m_beginFromRestart )
  {
    ReadSilo( true );
    m_Domains.m_feNodeManager.Initialize( );
    m_Domains.m_discreteElementSurfaceNodes.Initialize();
#ifdef SRC_EXTERNAL
    m_Domains.m_faultPatchNodes.Initialize();
#endif
    MarkFieldsToWrite();
  }
  else
  {
    CompleteObjectInitialization();
  }

  if(m_displayFields)
  {
    DisplayFields();
    exit(0);
  }

  if(m_rank == 0)
    std::cout << "Completed problem setup." << std::endl;
}


void ProblemManagerT::RegisterFields(){

  if (m_rank  == 0)
    std::cout << "Registering Fields..." << std::endl;

  if (m_rank  == 0)
    std::cout << "    Initial Conditions" << std::endl;
  for (std::vector<InitialConditionBase*>::iterator it_ic = m_initialConditions.begin() ; it_ic
       != m_initialConditions.end() ; ++it_ic)
  {
    (*it_ic)->RegisterFields(m_Domains);
  }

  if (m_rank  == 0)
    std::cout << "    Boundary Conditions" << std::endl;
  m_Domains.RegisterBCFields();


  for (std::map<std::string, SolverBase*>::iterator it_solver = m_solvers.begin() ; it_solver
       != m_solvers.end() ; ++it_solver)
  {
    if (m_rank  == 0)
      std::cout << "    Solver: " << it_solver->first << std::endl;
    it_solver->second->RegisterFields(m_Domains);
  }

  if( m_fractureFlag!=0 )
  {
    if (m_rank  == 0)
      std::cout << "Registering fields and maps." << std::endl;
    m_surfaceSeparation->RegisterFieldsAndMaps(m_Domains.m_feNodeManager, m_Domains.m_feEdgeManager,
                                               m_Domains.m_feFaceManager);

#ifdef SRC_INTERNAL2
    if(m_Domains.m_xfemManager != NULL)
    {
//      m_elementSplitting->RegisterFieldsAndMaps(m_Domains.m_feNodeManager,
// m_Domains.m_feEdgeManager, m_Domains.m_feFaceManager,
// m_Domains.m_feElementManager, m_Domains.m_crackManager,
// m_elementSplitting->m_nonSeparableElementSet);
      m_Domains.m_xfemManager->RegisterFieldsAndMaps(m_Domains.m_feNodeManager, m_Domains.m_feEdgeManager, m_Domains.m_feFaceManager,
                                                     m_Domains.m_feElementManager, *(m_Domains.m_crackManager),
                                                     m_Domains.m_xfemManager->m_nonSeparableElementSet);
    }
#endif
  }

#ifdef SRC_EXTERNAL
  m_Domains.m_wellboreManager.RegisterFields();
#endif

#ifdef SRC_INTERNAL2
  m_Domains.m_microseismicElementManager.RegisterFields();
#endif


  if (m_rank  == 0)
    std::cout << "Done Registering Fields." << std::endl;

}



void ProblemManagerT::SetInitialConditions()
{
  if (m_rank  == 0)
    std::cout << "Setting Initial conditions." << std::endl;
  for (std::vector<InitialConditionBase*>::iterator it_ic = m_initialConditions.begin() ; it_ic
       != m_initialConditions.end() ; ++it_ic)
  {
    (*it_ic)->Apply(m_Domains);
  }
}

void ProblemManagerT::SetInitialConstitutive()
{
  // TODO: HACK
  // propagate initial conditions to material states

  //See if there is an initial condition related to stress
  bool initialStressSet(0);
  for (std::vector<InitialConditionBase*>::iterator it_ic = m_initialConditions.begin() ; it_ic
       != m_initialConditions.end() ; ++it_ic)
  {
    initialStressSet = initialStressSet||(*it_ic)->IsThisInitialStress();
  }

  if( m_hackInitialStress || initialStressSet)
    m_Domains.m_feElementManager.HACKInitialConditions(  );

  m_initialConstitutive.Apply();

  // If requested, deform the mesh after the materials have been applied
  if( m_MeshGenerator )
  {
    if (m_MeshGenerator->m_delayMeshDeformation == 1)
    {
      m_MeshGenerator->RemapMesh(m_Domains);
    }
  }
}


void ProblemManagerT::CompleteObjectInitialization()
{
  //-----------------------------
  //SETUP PARALLEL PARTITIONING
  //-----------------------------

  // initialize the object managers
  InitializeObjectManagers();

  // find what faces and edges are on the computational domain boundary
  SetDomainBoundaryObjects();

  // reset the global to local maps for FE
  ResetGlobalToLocal();

#ifdef SRC_EXTERNAL
  //NOTE: this must be done after resetting global indices, since the
  //BEM data-structure relies on immutable global indexing
  //because this initializes values, it should also be done before IC
  // application
  m_Domains.m_faultElementManager.Initialize();
#endif

  // set up the neighbor lists and pack/unpack ghosts
  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;
  SetupNeighborListsPackUnpack(syncedFields);

  if (m_useMetis)
  {
    // DEM has to use space-partitioning because there is not a graph.  The
    // following will mess up DEM neighbor list.
    // The check should really be based on FEM-DEM.  Will address it later.
    m_partition.DeleteExcessNeighbors();
    m_partition.GraphBasedColoring();
  }


  // set external faces
  SetExternal();

  // initialize communications for the solvers
  for (std::map<std::string, SolverBase*>::iterator i = m_solvers.begin() ; i != m_solvers.end() ; ++i)
  {
    i->second->InitializeCommunications(this->m_partition);
  }
  //-----------------------------
  //-----------------------------


  // set initial conditions
  this->SetInitialConditions();

  // set initial conditions
  this->SetInitialConstitutive();

  // initialize surface separation
  this->InitializeSurfaceSeparation();
#ifdef SRC_INTERNAL2
  if(m_Domains.m_xfemManager != NULL)
  {
    this->InitializeElementSplitting();
  }
#endif

#ifdef SRC_EXTERNAL
  // Build wellbore mesh
  m_Domains.m_wellboreManager.MeshWellbore( m_Domains.m_feFaceManager );
#endif

  // initialize solvers
  for (std::map<std::string, SolverBase*>::iterator i = m_solvers.begin() ; i != m_solvers.end() ; ++i)
  {
    i->second->Initialize(m_Domains, this->m_partition);
  }

//  for( auto&& region : m_Domains.m_feElementManager.m_ElementRegions )
//  {
//    region.second.CalculateNodalMasses( m_Domains.m_feNodeManager );
//  }

  // synchronize after solver initialization
  m_partition.SetBufferSizes(syncedFields);
  m_partition.SynchronizeFields(syncedFields);

  MarkFieldsToWrite();

  // WriteNodeGraph(13, m_Domains);
}

/// Check existence of solvers named in solver sets
void ProblemManagerT::VerifySolvers()
{
  if(m_rank == 0)
    std::cout << "Verifying Solvers." << std::endl;
  for (std::vector<SolverApplicationSet>::iterator solverSets = m_solverApplicationSets.begin() ; solverSets
       != m_solverApplicationSets.end() ; ++solverSets)
  {
    for (std::vector<SolverApplication>::iterator solverSet =
           solverSets->m_solverAppliedToRegion.begin() ; solverSet
         != solverSets->m_solverAppliedToRegion.end() ; ++solverSet)
    {
      const std::string& solverName = solverSet->first;
      if (!isMember(solverName, m_solvers))
      {
        std::cout << m_solvers.size() << std::endl;
        for( std::map<std::string,SolverBase*>::iterator itr = m_solvers.begin() ; itr != m_solvers.end() ; ++itr)
        {
          std::cout << itr->first << std::endl;
        }
        ;
        throw GPException("Error SolverApplicationSet: Solver " + solverName + " was not found.");
      }
    }
  }
  if(m_rank == 0)
    std::cout << "Solver verification complete." << std::endl;
}


void ProblemManagerT::InitializeObjectManagers()
{
  if (m_rank  == 0)
    std::cout << "Initializing managers." << std::endl;
  m_Domains.m_feNodeManager.Initialize( );
  m_Domains.m_discreteElementSurfaceNodes.Initialize();
#ifdef SRC_EXTERNAL
  m_Domains.m_faultPatchNodes.Initialize();
#endif

  m_Domains.m_ellipsoidalDiscreteElementManager.Initialize();
  m_Domains.m_discreteElementManager.Initialize();

  m_Domains.m_ellipsoidalDiscreteElementManager.RecalculatePhysicalProperties();

  // the elementToNode data is set when the mesh is read

  // build the faces, along with the faceToElement and faceToNode arrays
  if (m_rank  == 0)
    std::cout << "Building faces." << std::endl;
  m_Domains.m_feFaceManager.BuildFaces(m_Domains.m_feNodeManager, m_Domains.m_feElementManager);


  TICPP::HierarchicalDataNode hdn;
  bool success = m_FileManager.ReadXMLInput(hdn);
  if (success)
  {
    HierarchicalDataNode* facesetNode = hdn.GetChild("Facesets");
    if( facesetNode != NULL )
      MeshUtilities::GenerateFasesetsAndAssociatedNodesets( *facesetNode, m_Domains.m_feFaceManager, m_Domains.m_feNodeManager );
  }


  // build flow face regions for the porous media solver
  const std::string erName = "FlowFaceRegion";
  ElementRegionT* elemRegion = stlMapLookupPointer(m_Domains.m_feElementManager.m_ElementRegions, erName);
  if (elemRegion != NULL)
  {
    if (elemRegion->m_numElems == 0)
      m_Domains.m_feElementManager.GenerateFlowFaceRegion(m_Domains.m_feFaceManager);
  }
//
  m_Domains.Initialize();

  if (m_rank  == 0)
    std::cout << "Building edges." << std::endl;
  m_Domains.m_feEdgeManager.BuildEdges(m_Domains.m_feFaceManager, m_Domains.m_feNodeManager);

  // construct the nodeToElement map. This requires the elemToNode map to be
  // done already!
  if (m_rank  == 0)
    std::cout << "Constructing node maps." << std::endl;
  m_Domains.m_feNodeManager.ConstructNodeToElementMap(m_Domains.m_feElementManager);

  // construct the nodeToFace map.
  m_Domains.m_feNodeManager.ConstructNodeToFaceMap(m_Domains.m_feFaceManager);

/*
   if( m_fractureFlag!=0 )
   {
    if (m_rank  == 0)  std::cout << "Registering fields and maps." << std::endl;
    m_surfaceSeparation->RegisterFieldsAndMaps(m_Domains.m_feNodeManager,
       m_Domains.m_feEdgeManager,
                                              m_Domains.m_feFaceManager);
   }

 */
}

void ProblemManagerT::SetDomainBoundaryObjects()
{
  if (m_rank  == 0)
    std::cout << "Setting domain boundary objects." << std::endl;

  m_Domains.m_feFaceManager.SetDomainBoundaryObjects();
  m_Domains.m_feEdgeManager.SetDomainBoundaryObjects(&m_Domains.m_feFaceManager);
  m_Domains.m_feNodeManager.SetDomainBoundaryObjects(&m_Domains.m_feFaceManager);
  m_Domains.m_feElementManager.SetDomainBoundaryObjects(&m_Domains.m_feFaceManager);

  //TODO: be more discriminating about parallelization at some point
  //** from Scott:
  //Sample logic might be ...
  //  const array<R1Tensor>& x =
  // discreteElementSurfaceNodes.GetFieldData<FieldInfo::currentPosition>();
  //  IsCoordInContactGhostRange( x(a) )

  m_Domains.m_discreteElementManager.GetFieldData<FieldInfo::isDomainBoundary>() = 1;
  m_Domains.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::isDomainBoundary>() = 1;
  m_Domains.m_discreteElementSurfaceFaces.GetFieldData<FieldInfo::isDomainBoundary>() = 1;

  m_Domains.m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::isDomainBoundary>() = 1;

#ifdef SRC_EXTERNAL
  m_Domains.m_faultPatchNodes.GetFieldData<FieldInfo::isDomainBoundary>() = 1;
  m_Domains.m_faultPatchFaces.GetFieldData<FieldInfo::isDomainBoundary>() = 1;
#endif
  // set domain boundary objects for periodic boundary conditions
  m_partition.SetPeriodicDomainBoundaryObjects(m_Domains);

  // set the layers to be communicated based on a constant number of layers of
  // elements
  m_Domains.m_feNodeManager.SetLayersFromDomainBoundary( 2 );

}

void ProblemManagerT::ResetGlobalToLocal()
{
  // assign global numbers to faces and edges
  if (m_rank  == 0)
    std::cout << "Assigning global numbers." << std::endl;
  m_partition.AssignGlobalIndices(m_Domains);

  m_partition.ResetSinglePartitionGlobalToLocalMap(m_Domains);

  if (m_rank  == 0)
    std::cout << "Setting ghost objects for single periodic partitions." << std::endl;
  m_partition.CreateSinglePartitionGhostObjects(m_Domains,m_Domains.m_externalFaces.m_contactActive, 1);

  // set the globalToLocal map for each object
  if (m_rank  == 0)
    std::cout << "Setting global to local map." << std::endl;

  m_Domains.m_feFaceManager.ResetGlobalToLocalMap();
  m_Domains.m_feNodeManager.ResetGlobalToLocalMap();
  m_Domains.m_feElementManager.ResetGlobalToLocalMap();
  m_Domains.m_feEdgeManager.ResetGlobalToLocalMap();

  m_Domains.m_discreteElementManager.ResetGlobalToLocalMap();
  m_Domains.m_discreteElementSurfaceNodes.ResetGlobalToLocalMap();
  m_Domains.m_discreteElementSurfaceFaces.ResetGlobalToLocalMap();

  m_Domains.m_ellipsoidalDiscreteElementManager.ResetGlobalToLocalMap();

#ifdef SRC_EXTERNAL
  m_Domains.m_faultPatchNodes.ResetGlobalToLocalMap();
  m_Domains.m_faultPatchFaces.ResetGlobalToLocalMap();
#endif

  m_partition.ResetSinglePartitionGlobalToLocalMap(m_Domains);
}

void ProblemManagerT::SetupNeighborListsPackUnpack(std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> >& syncedFields)
{
  //***** SETTING UP NEIGHBOR LISTS, PACKING/UNPACKING GHOSTS
  // **********************************************************
  // now determine what the neighbor relations are for the elements

  if (m_rank  == 0)
    std::cout << "Setting neighbor lists." << std::endl;
  this->m_partition.SetUpNeighborLists(m_Domains, m_Domains.m_externalFaces.m_contactActive );

  if (m_rank  == 0)
    std::cout << "Setting ghost arrays." << std::endl;
  m_partition.SetGhostArrays(m_Domains);
  m_partition.SetSinglePartitionGhostArrays(m_Domains);

  if (m_rank  == 0)
    std::cout << "Enforcing periodic boundaries on node locations." << std::endl;
  m_partition.CorrectReferencePositionsForPeriodicBoundaries(m_Domains);

  if (m_rank  == 0)
    std::cout << "Setting external faces." << std::endl;
  m_Domains.m_feFaceManager.SetIsExternal();

  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::mass>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::isExternal>::Name());

  m_partition.SetBufferSizes(syncedFields);
  m_partition.SynchronizeFields(syncedFields);
}

void ProblemManagerT::SetExternal()
{
  m_Domains.m_feEdgeManager.SetIsExternal(&m_Domains.m_feFaceManager);
  m_Domains.m_feNodeManager.SetIsExternal(&m_Domains.m_feFaceManager);

  m_Domains.m_discreteElementSurfaceFaces.m_isExternal = 1;
  m_Domains.m_discreteElementSurfaceNodes.m_isExternal = 1;

#ifdef SRC_EXTERNAL
  m_Domains.m_faultPatchFaces.m_isExternal = 1;
  m_Domains.m_faultPatchNodes.m_isExternal = 1;
#endif

  if (m_rank  == 0)
    std::cout << "Building external faces." << std::endl;
  m_Domains.m_externalFaces.BuildExternalFaces();
}

void ProblemManagerT::InitializeSurfaceSeparation()
{
  if( m_fractureFlag!=0 )
  {
    m_surfaceSeparation->Initialize( m_Domains.m_feNodeManager,
                                     m_Domains.m_feEdgeManager,
                                     m_Domains.m_feFaceManager,
                                     m_Domains.m_feElementManager);
  }
  OutputMeshInformation();
  if( m_fractureFlag!=0 && m_preFractureSets.size()>0 )
  {
    for (array<string>::size_type i=0 ; i<m_preFractureSets.size() ; ++i)
    {
      m_Domains.m_feFaceManager.PreSeparateFaces(m_preFractureSets[i], 1);
      if (m_Domains.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension == 3 )
      {
        m_surfaceSeparation->m_virtualFaces.PreSeparateFaces(m_preFractureSets[i], 1);
      }
    }

    int numNodesSplit = 100;
    for (localIndex i = 0 ; i < 10 ; ++i)
    {
      numNodesSplit = this->m_surfaceSeparation->SeparationDriver( m_Domains.m_feNodeManager,
                                                                   m_Domains.m_feEdgeManager,
                                                                   m_Domains.m_feFaceManager,
                                                                   m_Domains.m_externalFaces,
                                                                   m_Domains.m_feElementManager,
                                                                   m_partition,
                                                                   true,
                                                                   realT(0.0));
      (void)numNodesSplit;
    }
  }

  if ( m_fractureFlag!=0 && m_Domains.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension == 2 )
  {
    this->m_surfaceSeparation->PreexistingFracture2D( m_Domains.m_feNodeManager,
                                                      m_Domains.m_feEdgeManager,
                                                      m_Domains.m_feFaceManager,
                                                      m_Domains.m_externalFaces,
                                                      m_Domains.m_feElementManager,
                                                      m_partition, true);
  }
}

void ProblemManagerT::InitializeElementSplitting()
{
//  this->m_elementSplitting->SeparationDriver( true,
//                                              realT(0.0),
//                                              m_rank,
//                                              m_Domains.m_feNodeManager,
//                                              m_Domains.m_feEdgeManager,
//                                              m_Domains.m_feFaceManager,
//                                              m_Domains.m_externalFaces,
//                                              m_Domains.m_feElementManager,
//                                              m_partition,
//                                              m_Domains.m_crackManager,
//                                              m_Domains.m_crackSurfaceVertex,
//                                              m_Domains.m_splitElemThisIter);

#ifdef SRC_INTERNAL2
  m_Domains.m_xfemManager->SeparationDriver( true,
                                             realT(0.0),
                                             m_rank,
                                             m_Domains.m_feNodeManager,
                                             m_Domains.m_feEdgeManager,
                                             m_Domains.m_feFaceManager,
                                             m_Domains.m_externalFaces,
                                             m_Domains.m_feElementManager,
                                             m_partition,
                                             *(m_Domains.m_crackManager),
                                             *(m_Domains.m_crackSurfaceVertex),
                                             m_Domains.m_splitElemThisIter);

  if(m_Domains.m_splitElemThisIter)
    m_Domains.m_WriteXFEM = true;
#endif
}

/**
 * @author everyone
 *
 */
void ProblemManagerT::Run( double& outputTime, realT t_start )
{

  if(m_rank == 0)
    std::cout << "Running problem." << std::endl;
  timeval timeVal;


  if( m_trackEnergy )
    UpdateEnergy( true );

  // reset (overwrite) preexisting visit group file
  if( !m_beginFromRestart )
  {
    if(m_rank == 0)
    {
      std::ofstream geosVisit;
      geosVisit.open(m_visitFileGroupFile.c_str(), std::ios::out| std::ios::trunc);
      geosVisit.close();
    }
  }

  // ***** loop over all application sets *****
  for (std::vector<SolverApplicationSet>::iterator solverSet = m_solverApplicationSets.begin() ; solverSet
       != m_solverApplicationSets.end() ; ++solverSet)
  {

    if( !m_beginFromRestart )
    {
      m_dt = SetInitialTimeStep( m_problemTime, solverSet );
    }


#ifdef SRC_INTERNAL
    if (doBackgroundAMR)
      m_dt=std::min(m_dt,BackgroundAMRdtmax);
#endif


    //***** now march forward in time range specified by the current solver set
    // *****
    while (m_problemTime < solverSet->m_endTime)
    {


      {
        gettimeofday(&timeVal, NULL);
        double time0, time1;
        time0=timeVal.tv_sec+(timeVal.tv_usec/1000000.0);

        // check walltime
        int myForceRestart(0), forceRestart;
        if (time0-t_start >= m_maxWallTime)
        {
          myForceRestart = 1;
        }

        MPI_Allreduce(&myForceRestart, &forceRestart, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        if (m_problemTime + m_dt > solverSet->m_endTime)
          m_dt = solverSet->m_endTime - m_problemTime;
        if (m_forceOutputTime && m_problemTime + m_dt > m_nextPlotTime &&  m_problemTime < m_nextPlotTime)
          m_dt = m_nextPlotTime - m_problemTime;
        if (m_forceOutputTime && m_problemTime + m_dt > m_nextRestartTime &&  m_problemTime < m_nextRestartTime)
          m_dt = m_nextRestartTime - m_problemTime;

        if( m_writePlot && m_problemTime > (m_nextPlotTime-m_dt*1.0e-2) )
        {

          for (std::vector<SolverApplication>::iterator solver = solverSet->m_solverAppliedToRegion.begin() ;
               solver != solverSet->m_solverAppliedToRegion.end() ; ++solver)
          {
            const std::string& solverName = solver->first;
            const array<string>& regionsSolved = solver->second;
            SolverBase* const solverPtr = m_solvers[solverName];
            solverPtr->PostProcess(m_Domains, m_partition, regionsSolved);
          }

          WriteSilo( false );
        }

        if( (m_writeRestart && m_problemTime >= (m_nextRestartTime-m_dt*1.0e-2)) || forceRestart)
        {
          for (std::vector<SolverApplication>::iterator solver = solverSet->m_solverAppliedToRegion.begin() ;
               solver != solverSet->m_solverAppliedToRegion.end() ; ++solver)
          {
            const std::string& solverName = solver->first;
            const array<string>& regionsSolved = solver->second;
            SolverBase* const solverPtr = m_solvers[solverName];
            solverPtr->PostProcess(m_Domains, m_partition, regionsSolved);
          }

          WriteSilo( true );

          if (forceRestart)
          {
            std::cout << "Time limit exceeded...  Exiting" << std::endl;
            return;
          }
        }

        gettimeofday(&timeVal, NULL);
        time1=timeVal.tv_sec+(timeVal.tv_usec/1000000.0);

        outputTime += time1-time0;

      }


      std::cout.precision(14);
      if (m_partition.m_rank == 0 && m_cycleNumber%m_cycleReportFreq == 0)
      {
        std::cout<<"Cycle Number=" << m_cycleNumber << ", Time=" << m_problemTime << ", dt=" << m_dt<< std::endl;
      }

#ifdef SRC_INTERNAL
      BackgroundAMRdt=m_dt; //shouldn't dt be a accessible by a function call?
                            // i.e. static data of ProblemManagerT
#endif
      // set the next dt as a really big number
      realT nextdt = std::numeric_limits<double>::max();

      // Track errors in solvers
      int mySolverErrors(0), solverErrors;

      // loop over all the solvers that are applied in this application set
      for (std::vector<SolverApplication>::iterator solver = solverSet->m_solverAppliedToRegion.begin() ;
           solver != solverSet->m_solverAppliedToRegion.end() ; ++solver)
      {

        const std::string& solverName = solver->first;
        const array<string>& regionsSolved = solver->second;
        SolverBase* const solverPtr = m_solvers[solverName];

//        if(m_elementSplitting != NULL)
//        {
//          m_dt = solverPtr->TimeStep( m_problemTime, m_dt, m_cycleNumber,
// m_Domains, regionsSolved, m_partition, m_surfaceSeparation,
// m_elementSplitting);
//        }
//        else
//        {
//          m_dt = solverPtr->TimeStep( m_problemTime, m_dt, m_cycleNumber,
// m_Domains, regionsSolved, m_partition, m_surfaceSeparation);
//        }

        m_dt = solverPtr->TimeStep( m_problemTime, m_dt, m_cycleNumber, m_Domains, regionsSolved, m_partition, m_surfaceSeparation);

        if (solverPtr->m_stabledt.m_maxdt > std::numeric_limits<double>::max()*0.9 &&  solverSet->m_deltaTime < std::numeric_limits<double>::max() * 0.99)
        {
          solverPtr->m_stabledt.m_maxdt = solverSet->m_deltaTime;
        }
        // If solver does not specify a maxdt, then we set to the hard dt.

        nextdt = SetLocalMaxTimeStep( solverSet->m_deltaTime, solverPtr->m_stabledt.m_maxdt, nextdt );
#ifdef SRC_INTERNAL
        if (doBackgroundAMR)
          nextdt=std::min(nextdt,BackgroundAMRdtmax);
#endif

        if (solverPtr->m_solverErrorFlag)
        {
          mySolverErrors++;
        }
      }

      // Write a restart and exit if solver errors were found
      MPI_Allreduce(&mySolverErrors, &solverErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (solverErrors > 0)
      {
        WriteSilo( true );
        throw GPException("Errors encountered in solvers!");
      }

      if( m_trackEnergy )
      {
        UpdateEnergy();
        if( m_cycleNumber%m_trackEnergy == 0 )
        {
          if( m_rank == 0 )
          {
            std::cout<<"    InitialEnergy   |     StressPower     | ElasticStrainEnergy |       NodeWork       |     KineticEnergy   |      SP+KE "<<std::endl;
            std::cout<<m_initialEnergy<<" | "<<
              this->m_energy.StressPower()<<" | "<<
              this->m_energy.StrainEnergy()<<" | "<<
              this->m_energy.WorkDoneOnNodes()<<" | "<<
              this->m_energy.KineticEnergy()<<" | "<<
              this->m_energy.StressPower()+this->m_energy.KineticEnergy()<<std::endl;
          }
        }
      }


      if(solverSet->m_solverAppliedToRegion.size() == 1)
      {
        m_dt = m_solvers[solverSet->m_solverAppliedToRegion.begin()->first]->UpdateTimeStepMid(m_dt);
      }

      m_problemTime += m_dt;

      //m_dt = nextdt;
      m_dt = SetGlobalMaxTimeStep( nextdt );

      m_cycleNumber++;

    }

    // The while loop above does not do this for the last step, so the
    // postprocessed info is not updated for the last output.
    // We have to do this update manually.
    if (m_problemTime >= solverSet->m_endTime)
    {
      for (std::vector<SolverApplication>::iterator solver = solverSet->m_solverAppliedToRegion.begin() ;
           solver != solverSet->m_solverAppliedToRegion.end() ; ++solver)
      {
        const std::string& solverName = solver->first;
        const array<string>& regionsSolved = solver->second;
        SolverBase* const solverPtr = m_solvers[solverName];
        solverPtr->PostProcess(m_Domains, m_partition, regionsSolved);
      }

    }
  }

  if( m_writePlot && m_problemTime > (m_nextPlotTime-m_dt*1.0e-2) )
  {
    WriteSilo( false );
  }

  if( m_writeRestart && m_problemTime >= (m_nextRestartTime-m_dt*1.0e-2) )
  {
    WriteSilo( true );
  }


  if( m_rank == 0 )
  {
    std::cout<<"    InitialEnergy   |     StressPower     | ElasticStrainEnergy |       NodeWork       |     KineticEnergy   |      SP+KE "<<std::endl;
    std::cout<<m_initialEnergy<<" | "<<
      this->m_energy.StressPower()<<" | "<<
      this->m_energy.StrainEnergy()<<" | "<<
      this->m_energy.WorkDoneOnNodes()<<" | "<<
      this->m_energy.KineticEnergy()<<" | "<<
      this->m_energy.StressPower()+this->m_energy.KineticEnergy()<<std::endl;
  }

//  std::cout<<"fracture time parameters"<<std::endl;
//  std::cout<<"  m_partition.m_tX = "<<m_partition.m_t1<<" "
//                                    <<m_partition.m_t2<<" "
//                                    <<m_partition.m_t3<<" "
//                                    <<m_partition.m_t4<<std::endl;
}



void ProblemManagerT::WriteSilo( const bool isRestart )
{


  if(m_rank==0)
  {
    if( isRestart )
      std::cout << "Writing silo restart file" << std::endl;
    else
      std::cout << "Writing silo plot file" << std::endl;

  }
  m_siloFile.Initialize(PMPIO_WRITE);
  m_siloFile.WaitForBaton(m_rank, m_cycleNumber, isRestart );


  if( isRestart )
  {
    m_partition.WriteSilo( m_siloFile );
  }

// This is for 2D hydrofrac postprocessing
  if (m_writeFlowText && m_Domains.m_feFaceManager.m_toNodesRelation[0].size() == 2)
  {
    const array<R1Tensor>* u0 = m_Domains.m_feNodeManager.GetFieldDataPointer<R1Tensor>("displacement0");
    if (u0 != NULL)
    {
      if (!isRestart)
        WriteFlowTxt(m_cycleNumber, m_problemTime);
    }
  }



  m_Domains.WriteSilo(m_siloFile, m_cycleNumber, m_problemTime, isRestart, m_Domains.m_WriteXFEM, m_writeFEM, m_writeFEMFaces,
                      m_writeFEMEdges);

  if( isRestart )
  {
    //Write solvers
    {
      std::string subDirectory =   "Solvers";
      DBMkDir( m_siloFile.m_dbFilePtr, subDirectory.c_str() );
      DBSetDir(m_siloFile.m_dbFilePtr, subDirectory.c_str());
      for(std::map<std::string,SolverBase*>::const_iterator iters = this->m_solvers.begin() ; iters != this->m_solvers.end() ; ++iters)
      {
        DBMkDir( m_siloFile.m_dbFilePtr, iters->second->m_name.c_str() );
        DBSetDir(m_siloFile.m_dbFilePtr, iters->second->m_name.c_str());
//        if(m_rank==0)
//        {
//            std::cout << "Writing solver restart: " << iters->second->m_name
// << std::endl;
//        }
        iters->second->WriteSilo(m_siloFile);
        DBSetDir(m_siloFile.m_dbFilePtr, "..");
      }
      DBSetDir(m_siloFile.m_dbFilePtr, "..");
    }

    m_siloFile.DBWriteWrapper("m_numDomains",m_numDomains);


    m_siloFile.DBWriteWrapper("m_problemTime",m_problemTime);
    m_siloFile.DBWriteWrapper("m_dt",m_dt);
    m_siloFile.DBWriteWrapper("m_size",m_size);
    m_siloFile.DBWriteWrapper("m_rank",m_rank);
    m_siloFile.DBWriteWrapper("m_cycleNumber",m_cycleNumber);
    m_siloFile.DBWriteWrapper("m_plotIntervals",m_plotIntervals);
    m_siloFile.DBWriteWrapper("m_nextPlotIndex",m_nextPlotIndex);
    m_siloFile.DBWriteWrapper("m_nextPlotTime",m_nextPlotTime);
    m_siloFile.DBWriteWrapper("m_restartIntervals",m_restartIntervals);
    m_siloFile.DBWriteWrapper("m_nextRestartTime",m_nextRestartTime);

    if( m_surfaceSeparation )
      m_surfaceSeparation->WriteSilo(m_siloFile, m_cycleNumber, m_problemTime, isRestart );

    m_siloFile.DBWriteWrapper("m_trackEnergy",m_trackEnergy);

    array<real64> energy( EnergyT::numVars );
    m_energy.Serialize( energy.data() );
    m_siloFile.DBWriteWrapper("energy",energy);

  }

  // Save version id
  if( !isRestart )  // saving version with restart would be nice but would mess
                    // up tests
  {
    DBSetDir(m_siloFile.m_dbFilePtr, "/");
    std::string versionStr = GetRepoVersionString();
    m_siloFile.DBWriteWrapper("version",versionStr);
  }

  m_siloFile.HandOffBaton();
  m_siloFile.ClearEmptiesFromMultiObjects(m_cycleNumber);

  m_siloFile.Finish();


  // Add silo file to visit file list (geos.visit or equivalent)"geos.visit"
  char fileName[256] = { 0 };
  if(m_rank == 0 && !isRestart)
  {
    std::ofstream geosVisit;
    geosVisit.open(m_visitFileGroupFile.c_str(), std::ios::out | std::ios::app);
    sprintf(fileName, "%s_%06d", m_siloFile.m_fileRoot.c_str(), m_cycleNumber);
    geosVisit << fileName << std::endl;
    geosVisit.close();
  }


  if( isRestart )
  {
    CalculateNextSiloTime( m_restartIntervalTableName, m_restartIntervals, m_nextRestartTime, m_nextRestartIndex );
  }
  else
  {
    CalculateNextSiloTime( m_plotIntervalTableName, m_plotIntervals, m_nextPlotTime, m_nextPlotIndex );
  }
}



void ProblemManagerT::WriteFlowTxt (const int cycleNum, const realT time)
{
  std::string filename = "Flow";
  filename += "_" + toString(cycleNum);

  //if(m_partition.m_size > 1) filename += "_" + toString(m_partition.m_rank);
  filename += ".txt";

  std::ofstream fStream;

  for (int iRank = 0 ; iRank < m_partition.m_size ; ++iRank)
  {
    if (m_partition.m_rank == iRank)
    {
      if (iRank == 0)
      {
        fStream.open(filename.c_str(),std::ios::out);
      }
      else
      {
        fStream.open(filename.c_str(),std::ios::app);
      }

      const array<R1Tensor>& u0 = m_Domains.m_feNodeManager.GetFieldData<R1Tensor>("displacement0");
      const array<integer>& isGhost = m_Domains.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
      const array<real64>& aperture = m_Domains.m_feFaceManager.GetFieldData<realT>( "Aperture" );
      const array<real64>& faceFluidPressure = m_Domains.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
      const array<real64>& faceFluidDensity = m_Domains.m_feFaceManager.GetFieldData<FieldInfo::density>();
      const array<real64>& flowRate = m_Domains.m_feFaceManager.GetFieldData<realT>("flowRate");
      const array<integer>& flowFaceType = m_Domains.m_feFaceManager.GetFieldData<int>("flowFaceType");

      if (iRank==0)
      {
        fStream << "t= " << time << std::endl;
        fStream << "CellID, xNodeA, yNodeA, xNodeB, yNodeB, fluidDensity, fluidPressure, apertureWidth, flowRate, rank" << std::endl;
      }
      for (localIndex i = 0 ; i != m_Domains.m_feFaceManager.DataLengths() ; ++i)
      {
        if (isGhost[i] < 0 && flowFaceType[i]>0)
        {
          fStream << i << ", ";
          //fStream << m_Domains.m_feFaceManager.m_localToGlobalMap[i] << ", ";

          for( lArray1d::iterator j = m_Domains.m_feFaceManager.m_toNodesRelation[i].begin() ;
               j!=m_Domains.m_feFaceManager.m_toNodesRelation[i].end() ; ++j )
          {
            fStream << (*m_Domains.m_feNodeManager.m_refposition)[*j][0] + u0[*j][0] << ", ";
            fStream << (*m_Domains.m_feNodeManager.m_refposition)[*j][1] + u0[*j][1] << ", ";
          }
          fStream << faceFluidDensity[i] << ", ";
          fStream << faceFluidPressure[i] << ", ";
          fStream << aperture[i] << ", ";
          fStream << flowRate[i] << ", ";
          fStream << m_partition.m_rank << std::endl;
        }
      }

      if (iRank == m_partition.m_size-1)
        fStream<<"EOF";

      fStream.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

  }
}

void ProblemManagerT::ReadSilo( const bool isRestart )
{

  if(m_rank==0)
  {
    if( isRestart )
      std::cout << "Reading silo restart file" << std::endl;
    else
      std::cout << "Reading silo plot file" << std::endl;

  }

  m_siloFile.Initialize(PMPIO_READ);
//  m_siloFile.WaitForBaton(m_rank, m_cycleNumber, isRestart );
  m_siloFile.WaitForBaton(m_rank, m_beginFromRestartFileName );

  if( isRestart )
  {
    m_partition.ReadSilo( m_siloFile );
    m_partition.SetDomain(m_Domains);
  }


  if( isRestart )
  {
    //Read solvers
    {
      std::string subDirectory =   "Solvers";
      DBSetDir(m_siloFile.m_dbFilePtr, subDirectory.c_str());
      DBtoc* silotoc = DBGetToc(m_siloFile.m_dbFilePtr);
      array<string> dirs(silotoc->ndir);
      for(int i = 0 ; i < silotoc->ndir ; i++)
        dirs[i] = silotoc->dir_names[i];
      for(localIndex i = 0 ; i < dirs.size() ; i++)
      {
        DBSetDir(m_siloFile.m_dbFilePtr, dirs[i].c_str());
        std::map<std::string,SolverBase*>::iterator iter = this->m_solvers.find(dirs[i]);
        if(iter == this->m_solvers.end())
          throw GPException("Attempt to read solver that does not exist");
        iter->second->ReadSilo(m_siloFile);
        DBSetDir(m_siloFile.m_dbFilePtr, "..");
      }
      DBSetDir(m_siloFile.m_dbFilePtr, "..");
    }

    m_siloFile.DBReadWrapper("m_numDomains",m_numDomains);

    m_siloFile.DBReadWrapper("m_problemTime",m_problemTime);
    m_siloFile.DBReadWrapper("m_dt",m_dt);
    m_siloFile.DBReadWrapper("m_size",m_size);
    m_siloFile.DBReadWrapper("m_rank",m_rank);
    m_siloFile.DBReadWrapper("m_cycleNumber",m_cycleNumber);
    // m_siloFile.DBReadWrapper("m_plotIntervals",m_plotIntervals);
    m_siloFile.DBReadWrapper("m_nextPlotIndex",m_nextPlotIndex);
    m_siloFile.DBReadWrapper("m_nextPlotTime",m_nextPlotTime);
    // m_siloFile.DBReadWrapper("m_restartIntervals",m_restartIntervals);
    m_siloFile.DBReadWrapper("m_nextRestartTime",m_nextRestartTime);

    m_Domains.ReadSilo( m_siloFile, m_cycleNumber, m_problemTime, isRestart,
                        m_writeFEM, m_writeFEMFaces, m_writeFEMEdges );

    if( m_surfaceSeparation )
    {
      m_surfaceSeparation->m_virtualFaces.m_elemManagerHACK = &(m_Domains.m_feElementManager);
      m_surfaceSeparation->ReadSilo(m_siloFile, m_cycleNumber, m_problemTime, isRestart);
    }

    m_siloFile.DBReadWrapper("m_trackEnergy",m_trackEnergy);

    array<real64> energy( EnergyT::numVars );
    m_siloFile.DBReadWrapper("energy",energy);
    m_energy.Deserialize( energy.data() );

  }



  m_siloFile.HandOffBaton();
  m_siloFile.ClearEmptiesFromMultiObjects(m_cycleNumber);

  m_siloFile.Finish();


  // Add silo file to geos.visit (or equivalent)
  char fileName[100] = { 0 };
  if(m_rank == 0)
  {

    std::ios_base::openmode omode;
    omode = std::ios::out | std::ios::app;
    std::ofstream geosVisit(m_visitFileGroupFile.c_str(),omode );
    sprintf(fileName, "%s_%06d", m_siloFile.m_fileRoot.c_str(), m_cycleNumber);
    geosVisit << fileName << std::endl;
    geosVisit.close();
  }


  if( isRestart )
  {
    CalculateNextSiloTime( m_restartIntervalTableName, m_restartIntervals, m_nextRestartTime, m_nextRestartIndex );
  }
  else
  {
    CalculateNextSiloTime( m_plotIntervalTableName, m_plotIntervals, m_nextPlotTime, m_nextPlotIndex );
  }
}


void ProblemManagerT::CalculateNextSiloTime( const std::string& intervalTableName,
                                             realT& siloIntervals,
                                             realT& nextSiloTime,
                                             globalIndex& nextSiloIndex )
{
  if( !(intervalTableName.empty()) )
  {
    const TableManager& tableManager = TableManager::Instance();
    siloIntervals = tableManager.LookupTable<1>(intervalTableName, &m_problemTime, TableInterpolation::zeroth);
    nextSiloTime += siloIntervals;
  }
  else if(siloIntervals > 0.0)
  {
    const globalIndex inext = (globalIndex)(this->m_problemTime / siloIntervals) + 1;
    nextSiloIndex = inext > nextSiloIndex ? inext : nextSiloIndex + 1;
    nextSiloTime = nextSiloIndex * siloIntervals;
  }

  /*
     const realT nextTimeish = this->m_problemTime + dt*0.5 ;
     if( m_nextPlotTime < nextTimeish )
     {
     m_nextPlotTime = this->m_problemTime + m_plotIntervals;
     }*/
}


void ProblemManagerT::DisplayUsage()
{
  std::cout << "Usage: GEOS.x [options]  \n" << "Options: \n"
            << "  -help, -h          : display this message.\n"
            << "  -f FILEROOT        : read problem input and mesh from files FILEROOT.xml and FILEROOT.geom.\n"
            << "  -i INPUTFILE       : read specified xml input file.\n"
            << "  -m MESHFILE        : read specified mesh file.\n"
            << "  -d DEMESHFILE      : read specified discrete element mesh file.\n"
            << "  -e EDEFILE         : read specified ellipsoidal discrete element mesh file.\n"
            << "  -p KEY=VALUE       : add model parameter 'KEY' with value 'VALUE'.\n"
            << "  -r RESTARTFILE     : restart using file 'RESTARTFILE' "
            << "  -include FILE      : include named xml file.\n"
            << "  -solvers           : report supported physics solvers.\n"
            << "  -units             : report supported units.\n"
            << "  -report_defaults   : report values of default fields when set.\n"
            << "  -disable_defaults  : disable default fields.\n"
            << "  -report_parameters : report all parameter values.\n"
            << "  -xpar N            : set number of x axis spatial partitions to N.\n"
            << "  -ypar N            : set number of y axis spatial partitions to N.\n"
            << "  -zpar N            : set number of z axis spatial partitions to N.\n"
            << "  -fields            : report active fields (typically used in conjunction with input file).\n"
            << "  -write_XML FILE    : write XML (after parsing includes, conditionals and parameters) to file.\n"
            << "  -version, -v       : display version.\n" << std::endl;
}


void ProblemManagerT::DisplayVersion()
{
  ::DisplayVersion(m_rank);
}

void ProblemManagerT::DisplaySplash()
{
  if( m_rank==0 && m_displaySplash )
  {
    ::DisplaySplash(m_rank);
  }
}


void ProblemManagerT::DisplaySolvers()
{
  std::vector<std::string> nameList;
  SolverFactory::GetSolverNames(nameList);
  for (size_t i = 0 ; i < nameList.size() ; ++i)
    std::cout << nameList[i] << std::endl;
}


void ProblemManagerT::DisplayUnits()
{

  HierarchicalDataNode emptyNode;
  UnitManager& unitManager = UnitManager::Instance();
  unitManager.ReadXML(&emptyNode);
  if (m_rank == 0)
  {
    unitManager.ReportAllUnits();
  }
}

/**
 * @brief Display all registered fields in a simulation
 * @author walsh24
 *  Display all registered fields in a simulation, then quit.
 *   Usage:
 *   GEOS.x -i SubstepTest.xml -fields
 *   N.B. a valid input file (and mesh etc.) is required.
 *   Empty objects (data length = 0) are not reported on.
 *   Prints the field name, object type, whether written to restart and plot.
 *
 **/
void ProblemManagerT::DisplayFields()
{

  std::cout << "                                    " << std::endl;
  std::cout << "####################################" << std::endl;
  std::cout << "##                                ##" << std::endl;
  std::cout << "##           Field Data           ##" << std::endl;
  std::cout << "##                                ##" << std::endl;
  std::cout << "####################################" << std::endl;

  // loop over data structures
  PhysicalDomainT::ObjectDataStructureKeys objKey;
  for(objKey=static_cast<PhysicalDomainT::ObjectDataStructureKeys>(0) ;
      objKey!=PhysicalDomainT::numObjectDataStructureNames ; ++objKey)
  {

    if(objKey != PhysicalDomainT::VirtualEdgeManager
       && objKey!= PhysicalDomainT::VirtualFaceManager  )     // virtual objects
                                                              // hold no data -
                                                              // can remove
                                                              // without risk if
                                                              // keys removed
    { // get object
      ObjectDataStructureBaseT& object  = m_Domains.GetObjectDataStructure( objKey);

      if(object.DataLengths() > 0)  // ignore empty objects
      { // get object name
        std::string dataStructureName = PhysicalDomainT::GetObjectDataStructureName(objKey);

        std::cout << std::endl;
        std::cout << "## Object: " << dataStructureName << std::endl
                  << "## Field    "<< "\t" << "Type" << "\t" << "Restart"<< "\t" << "Plot" << std::endl
                  << "####################################" << std::endl;
        std::cout << std::endl;


        // get field names
        array<string> fieldNames;
        object.GetAllFieldNames( fieldNames );

        // loop over fields and report status
        for(unsigned i =0 ; i < fieldNames.size() ; ++i )
        {

          std::string fieldName = fieldNames[i];

          // sanity check
          std::map<std::string, FieldBase*>::iterator itr = FieldInfo::AttributesByName.find(fieldName);
          if( itr == FieldInfo::AttributesByName.end() )
          {
            throw GPException("Field "+ fieldName + "not registered in FieldInfo::AttributesByName");
          }

          bool restart = FieldInfo::AttributesByName[fieldName]->m_WriteToRestart;
          bool plot = FieldInfo::AttributesByName[fieldName]->m_WriteToPlot;
          FieldInfo::FieldType type = object.GetFieldType(fieldName);

          std::string typeStr = "unrecognized";
          if(type != FieldInfo::numFieldTypes)
          {
            typeStr = FieldInfo::FieldTypeName(type);
          }

          // report
          std::cout << fieldName << "\t" << typeStr << "\t"<< restart << "\t" << plot <<"\n";

        } // loop over fields
      } // ignore empty objects
    } // ? is virtual object?
  } // loop over objects



  // Special treatment for element region fields
  std::cout << std::endl;
  std::cout << "## Object: " << "Element Region" << std::endl
            << "## Field    "<< "\t" << "Type" << "\t" << "Restart"<< "\t" << "Plot" << std::endl
            << "####################################" << std::endl;
  std::cout << std::endl;

  for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator elementRegionIter =
         m_Domains.m_feElementManager.m_ElementRegions.begin() ; elementRegionIter
       != m_Domains.m_feElementManager.m_ElementRegions.end() ; ++elementRegionIter)
  {
    ElementRegionT& elementRegion = elementRegionIter->second;
    // get field names
    array<string> fieldNames;
    elementRegion.GetAllFieldNames( fieldNames );

    // loop over fields and report status
    for(unsigned i =0 ; i < fieldNames.size() ; ++i )
    {

      std::string fieldName = fieldNames[i];

      // sanity check
      std::map<std::string, FieldBase*>::iterator itr = FieldInfo::AttributesByName.find(fieldName);
      if( itr == FieldInfo::AttributesByName.end() )
      {
        throw GPException("Field "+ fieldName + "not registered in FieldInfo::AttributesByName");
      }

      bool restart = FieldInfo::AttributesByName[fieldName]->m_WriteToRestart;
      bool plot = FieldInfo::AttributesByName[fieldName]->m_WriteToPlot;
      FieldInfo::FieldType type = elementRegion.GetFieldType(fieldName);

      std::string typeStr = "unrecognized";
      if(type != FieldInfo::numFieldTypes)
      {
        typeStr = FieldInfo::FieldTypeName(type);
      }

      // report
      std::cout << fieldName << "\t" << typeStr << "\t"<< restart << "\t" << plot <<"\n";

    } // loop over fields

  }

}

void ProblemManagerT::OutputMeshInformation()
{

  if (m_rank == 0)
  {
    printf(
      "   RANK    NUMELEMS(GHOST)     NUMNODES(GHOST)     NUMFACES(GHOST)     NUMEDGES(GHOST)   \n");
    printf(
      "  ------ ------------------- ------------------- ------------------- ------------------- \n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  localIndex numElems = 0;
  localIndex numGhostElems = 0;

  for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::const_iterator i =
         m_Domains.m_feElementManager.m_ElementRegions.begin() ; i
       != m_Domains.m_feElementManager.m_ElementRegions.end() ; ++i)
  {
    numElems += i->second.DataLengths();
    numGhostElems += i->second.GetNumGhosts();
  }

  const localIndex numNodes = m_Domains.m_feNodeManager.DataLengths();
  const localIndex numGhostNodes = m_Domains.m_feNodeManager.GetNumGhosts();
  const localIndex numFaces = m_Domains.m_feFaceManager.DataLengths();

  const localIndex numGhostFaces = m_Domains.m_feFaceManager.GetNumGhosts();
  const localIndex numEdges = m_Domains.m_feEdgeManager.DataLengths();

  const localIndex numGhostEdges = m_Domains.m_feEdgeManager.GetNumGhosts();

  printf(" %6i %8lu (%8lu) %8lu (%8lu) %8lu (%8lu) %8lu (%8lu)\n", m_rank, numElems, numGhostElems,
         numNodes, numGhostNodes, numFaces, numGhostFaces, numEdges, numGhostEdges);
}



/**
 *
 * @param solverSet
 * @return
 */
realT ProblemManagerT::SetInitialTimeStep( const realT& time, std::vector<SolverApplicationSet>::iterator solverSet )
{
  realT dt = solverSet->m_deltaTime;
  realT nextdt = dt;

  // if the timestep was not explicitly set, then we have to set it.
  if( solverSet->m_deltaTime > (0.99*std::numeric_limits<double>::max()) )
  {
    // loop over all solvers in this application set and get the timestep
    // requirement
    for (std::vector<SolverApplication>::iterator solver = solverSet->m_solverAppliedToRegion.begin() ;
         solver != solverSet->m_solverAppliedToRegion.end() ; ++solver)
    {
      const std::string& solverName = solver->first;
      const array<string>& regionsSolved = solver->second;
      SolverBase* const solverPtr = m_solvers[solverName];

      solverPtr->SetMaxStableTimeStep( time, m_Domains, regionsSolved, this->m_partition );

      if( solverPtr->m_stabledt.m_maxdt < dt )
        dt = solverPtr->m_stabledt.m_maxdt;
    }
  }

  nextdt = dt;
  MPI_Allreduce (&nextdt,&dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  return dt;
}

/**
 *
 * @param local_dt
 * @return
 */
realT ProblemManagerT::SetGlobalMaxTimeStep( const realT& local_dt )
{
  realT temp = local_dt;
  realT nextdt;
  MPI_Allreduce (&temp,&nextdt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  return nextdt;
}


// ***** STATIC MEMBER FUNCTIONS *****
realT ProblemManagerT::SetLocalMaxTimeStep( const realT& hard_dt, const realT& solver_dt, const realT& current_max_dt )
{
  realT max_dt = current_max_dt;
  // if a hard timestep was set, then just use it.
//  if( hard_dt<(0.99*std::numeric_limits<double>::max()) && solver_dt >=
// current_max_dt )
  if( hard_dt<(0.99*std::numeric_limits<double>::max()) )
  {
    max_dt = hard_dt;
  }
  else
  {
    // set the time step for the next step. If the solver returns a timestep,
    // and it is less than the value of nextdt, then
    // reset the value of nextdt
    if( ( solver_dt<(0.99*std::numeric_limits<double>::max()) ) )
    {
      if(solver_dt < max_dt)
        max_dt = solver_dt;
    }
    else
    {
      throw GPException( "timestep not properly set");
    }
  }

  if( max_dt > hard_dt )
    max_dt = hard_dt;

  return max_dt;
}

void ProblemManagerT::UpdateEnergy( const bool set_initial )
{
  m_Domains.UpdateEnergy();
  m_energy = m_Domains.m_energy;


#ifdef GPAC_MPI
  {

    realT ibuffer[EnergyT::numVars];
    realT obuffer[EnergyT::numVars];
    m_energy.Serialize(obuffer);
    MPI_Reduce(obuffer,ibuffer,EnergyT::numVars,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
//    MPI_Allreduce(obuffer,ibuffer,EnergyT::NumVars(),MPI_DOUBLE,
// MPI_SUM,MPI_COMM_WORLD);


    m_energy.Deserialize(ibuffer);


//    std::cout<<"Energy (" << elementRegionIter->second.m_regionName << ") " <<
// time << " " << tmp.TotalEnergy() << " " << tmp.StrainEnergy() << std::endl;
  }

#else
  EnergyT& tmp = this->elementEnergy[elementRegionIter->second.m_regionName];
  tmp.Increment(energyRegion);
  std::cout<<"Energy (" << elementRegionIter->second.m_regionName << ") " << time << " " << tmp.TotalEnergy() << " " << tmp.StrainEnergy() << std::endl;
#endif

  if( set_initial )
  {
    m_initialEnergy = this->m_energy.StressPower()+this->m_energy.KineticEnergy();
  }
/*
   if( m_rank == 0 )
   {
    std::cout<<"StressPower, ElasticStrainEnergy, KE, TE : "<<
               this->m_energy.StressPower()<<", "<<
               this->m_energy.StrainEnergy()<<", "<<
               this->m_energy.KineticEnergy()<<", "<<
               this->m_energy.WorkDoneOnNodes()<<", "<<
               this->m_energy.StressPower()+this->m_energy.KineticEnergy()<<",
                  "<<
               m_initialEnergy<<std::endl;
   }*/

}


void ProblemManagerT::UpdateOwnership()
{
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  m_Domains.m_feFaceManager.AddKeylessDataField<localIndex>("myRank", true, true);
  array<localIndex>& myRank = m_Domains.m_feFaceManager.GetFieldData<localIndex>( "myRank" );

  for (localIndex i = 0 ; i != m_Domains.m_feFaceManager.DataLengths() ; ++i)
  {
    myRank[i] = rank;
  }


}

void ProblemManagerT::MarkFieldsToWrite()
{

  if (m_nameFieldsToPlot.size() > 0)
  {
    {
      m_Domains.m_feNodeManager.UncheckAllFieldsForPlot();
      m_Domains.m_feFaceManager.UncheckAllFieldsForPlot();
      m_Domains.m_feEdgeManager.UncheckAllFieldsForPlot();
      m_Domains.m_feElementManager.UncheckAllFieldsForPlot();
      m_Domains.m_externalFaces.UncheckAllFieldsForPlot();
      m_Domains.m_discreteElementManager.UncheckAllFieldsForPlot();

      #ifdef SRC_EXTERNAL
      m_Domains.m_wellboreManager.UncheckAllFieldsForPlot();
      #endif

      #ifdef SRC_INTERNAL2
      m_Domains.m_microseismicElementManager.UncheckAllFieldsForPlot();
      #endif

      for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator elementRegionIter =
             m_Domains.m_feElementManager.m_ElementRegions.begin() ; elementRegionIter
           != m_Domains.m_feElementManager.m_ElementRegions.end() ; ++elementRegionIter)
      {
        ElementRegionT& elementRegion = elementRegionIter->second;
        elementRegion.UncheckAllFieldsForPlot();
      }
    }

    // Build list of field names using solver defined common fields
    for( localIndex i=0 ; i<m_nameFieldsToPlot.size() ; ++i )
    {
      if (m_nameFieldsToPlot[i] == "CommonFields")
      {
        m_nameFieldsToPlot.erase(m_nameFieldsToPlot.begin()+i);

        for (std::map<std::string, SolverBase*>::iterator it_solver = m_solvers.begin() ; it_solver!= m_solvers.end() ; ++it_solver)
        {
          for( array<string>::iterator cf=it_solver->second->m_commonFields.begin() ; cf!=it_solver->second->m_commonFields.end() ; ++cf )
          {
            m_nameFieldsToPlot.push_back(*cf);
          }
        }
      }
    }


    for (localIndex i=0 ; i<m_nameFieldsToPlot.size() ; ++i)
    {
      FieldBase*& fieldAttributesMap = stlMapLookup( FieldInfo::AttributesByName, m_nameFieldsToPlot[i], "" );
      if (fieldAttributesMap == NULL)
      {
        char msg[200];
        sprintf(msg, "Error! The following requested field does not exist: %s",  m_nameFieldsToPlot[i].c_str());
        throw GPException(msg);
      }
      fieldAttributesMap->m_WriteToPlot = true;

    }

    for (localIndex i=0 ; i<m_nameAdditionalFieldsToPlot.size() ; ++i)
    {
      FieldBase*& fieldAttributesMap = stlMapLookup( FieldInfo::AttributesByName, m_nameAdditionalFieldsToPlot[i], "" );
      if (fieldAttributesMap == NULL)
      {
        char msg[200];
        sprintf(msg, "Error! The following requested field does not exist: %s",  m_nameAdditionalFieldsToPlot[i].c_str());
        throw GPException(msg);
      }
      fieldAttributesMap->m_WriteToPlot = true;

    }

  }

}

void ProblemManagerT::WriteNodeGraph( const localIndex nodeID,
                                      PhysicalDomainT& domain)
{

//  int edgeIDToGraphNodeMap[domain.m_feEdgeManager.DataLengths()];
//  std::fill_n(edgeIDToGraphNodeMap, domain.m_feEdgeManager.DataLengths(), -1);
//  std::unique_ptr<int[]> edgeIDToGraphNodeMap(new
// int[domain.m_feEdgeManager.DataLengths()]);

  array<integer> edgeIDToGraphNodeMap(domain.m_feEdgeManager.DataLengths());
  edgeIDToGraphNodeMap = -1;

  localIndex nGNode = 0;
  R1Tensor xGNode;

  std::cout << "Writing nodal graph!" << std::endl;
  std::ofstream outfile, visitIn;
  outfile.open("NodeGraph.csv");
  visitIn.open("VisitIn.dat");

  outfile << "Id" << std::endl;
  visitIn << domain.m_feNodeManager.m_nodeToEdgeMap[nodeID].size() << std::endl;
  for( lSet::const_iterator iedge=domain.m_feNodeManager.m_nodeToEdgeMap[nodeID].begin() ;
       iedge!=domain.m_feNodeManager.m_nodeToEdgeMap[nodeID].end() ; ++iedge )
  {

    if (edgeIDToGraphNodeMap[*iedge] == -1)
    {
      ++nGNode;
      edgeIDToGraphNodeMap[*iedge] = nGNode;
    }

    xGNode = (*domain.m_feNodeManager.m_refposition)[domain.m_feEdgeManager.m_toNodesRelation(*iedge,0)];
    xGNode += (*domain.m_feNodeManager.m_refposition)[domain.m_feEdgeManager.m_toNodesRelation(*iedge,1)];
    xGNode /=2.0;
    xGNode -= (*domain.m_feNodeManager.m_refposition)[nodeID];

    visitIn << edgeIDToGraphNodeMap[*iedge] << " " << xGNode[0] << " " << xGNode[1]<< " " <<xGNode[2] << std::endl;
    outfile << *iedge << std::endl;

  }

  outfile.close();


  outfile.open("EdgeGraph.csv");

  outfile << "\"Source\";\"Target\";\"Label\"" << std::endl;
  visitIn << domain.m_feNodeManager.m_nodeToFaceMap[nodeID].size() << std::endl;
  for( lSet::const_iterator iface=domain.m_feNodeManager.m_nodeToFaceMap[nodeID].begin() ;
       iface!=domain.m_feNodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
  {
    for( lArray1d::iterator j = domain.m_feFaceManager.m_toEdgesRelation[*iface].begin() ;
         j!=domain.m_feFaceManager.m_toEdgesRelation[*iface].end() ; ++j )
    {
      if (domain.m_feEdgeManager.m_toNodesRelation(*j,0) == nodeID || domain.m_feEdgeManager.m_toNodesRelation(*j,1) == nodeID)
      {
        outfile << *j << ";";
        visitIn << edgeIDToGraphNodeMap[*j] << " ";
      }
    }
    outfile << "\"Face_" << *iface << "\""<<std::endl;
    visitIn << std::endl;


  }

  outfile.close();
}
