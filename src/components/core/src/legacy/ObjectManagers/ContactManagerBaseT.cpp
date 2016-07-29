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
 * @file ContactManagerBaseT.cpp
 * @author Randolph Settgast
 * @date created on July 14, 2011
 */

#include "ContactManagerBaseT.h"
#include <limits>
#include <cstdlib>

/**
 * @brief Constructor for the contacts / neighbor pairs
 * @author Scott Johnson
 * @date Jul 14, 2011
 * We are currently assuming that "contact" is the union of all "neighbor pairs" and "inactive + active physical contacts"
 * Note that neighbor pairs that are also contacts would be the intersection of these two sets
 * This data structure could be split into different structures later if it is
 * found that these are too cumbersome
 */
ContactManagerBaseT::ContactManagerBaseT() :
  ObjectDataStructureBaseT(ObjectDataStructureBaseT::ContactBaseManager),
  m_contact(NULL),
  m_sliding_law(0),
  m_yield_traction(1e16),
  m_coulomb_coefficient(1e16),
  m_glob_penalty_n(0.0),
  m_glob_penalty_t1(0.0),
  m_glob_penalty_t2(0.0),
  m_stab_magnifying_factor(1.0),
  m_nitsche_active(false),
  m_nitsche_symmetry_active(false),
  m_implicitContactActive(false),
  m_weighted_nitsche_active(false),
  m_nitsche_normal_active(false),
  m_nitsche_tau1_active(false),
  m_nitsche_tau2_active(false),
  m_traction_n_tol(1e-04)
{
  //---FACEPAIRS---
  this->AddKeylessDataField<localIndex>( "face1", true, true);
  this->AddKeylessDataField<localIndex>( "face2", true, true);
  this->AddKeylessDataField<int>( "active", true, true);
  this->AddKeylessDataField<R1Tensor>( "normal", true, true);
  this->AddKeylessDataField<realT>( "area", true, true);
  this->AddKeylessDataField<R1Tensor>( "applicationPoint", true, true);

  //normal
  this->AddKeylessDataField<realT>( "normalApproach", true, true);
  this->AddKeylessDataField<realT>( "normalApproachMax", true, true);

  //tangential
  this->AddKeylessDataField<R1Tensor>( "shearSlip", true, true);
  this->AddKeylessDataField<R1Tensor>( "velocity", true, true);
}

ContactManagerBaseT::ContactManagerBaseT( const ObjectType objectType ) :
  ObjectDataStructureBaseT(objectType),
  m_contact(NULL),
  m_sliding_law(0),
  m_yield_traction(1e16),
  m_coulomb_coefficient(1e16),
  m_glob_penalty_n(0.0),
  m_glob_penalty_t1(0.0),
  m_glob_penalty_t2(0.0),
  m_stab_magnifying_factor(1.0),
  m_nitsche_active(false),
  m_nitsche_symmetry_active(false),
  m_implicitContactActive(false),
  m_weighted_nitsche_active(false),
  m_nitsche_normal_active(false),
  m_nitsche_tau1_active(false),
  m_nitsche_tau2_active(false),
  m_traction_n_tol(1e-04)
{
  //---FACEPAIRS---
  this->AddKeylessDataField<localIndex>( "face1", true, true);
  this->AddKeylessDataField<localIndex>( "face2", true, true);
  this->AddKeylessDataField<int>( "active", true, true);
  this->AddKeylessDataField<R1Tensor>( "normal", true, true);

  this->AddKeylessDataField<realT>( "area", true, true);
  this->AddKeylessDataField<R1Tensor>( "applicationPoint", true, true);

  //normal
  this->AddKeylessDataField<realT>( "normalApproach", true, true);
  this->AddKeylessDataField<realT>( "normalApproachMax", true, true);

  //tangential
  this->AddKeylessDataField<R1Tensor>( "shearSlip", true, true);
  this->AddKeylessDataField<R1Tensor>( "velocity", true, true);
}

ContactManagerBaseT::~ContactManagerBaseT()
{
  if(m_contact)
    delete m_contact;
}

void ContactManagerBaseT::ReadXML(TICPP::HierarchicalDataNode* hdn)
{
  TICPP::HierarchicalDataNode* childNode = hdn->Next(true);
  if(childNode == NULL)
    throw GPException("Expected a sub-node to read in contact!\n");
  m_contact = InterfaceFactory::NewInterface(childNode->Heading(), childNode);

  //@author: Chandra Annavarapu
  //@brief: Penalty/stabilization parameter for interfacial contact
  m_implicitContactActive = hdn->GetAttributeOrDefault<bool>("ImplicitActive", false);
  m_use_contact_search = hdn->GetAttributeOrDefault<bool>("useContactSearch", false);
  bool threeD = hdn->GetAttributeOrDefault<bool>("ThreeDTrue", 0.0);
  if (m_implicitContactActive)
  {
    m_nitsche_active = hdn->GetAttributeOrDefault<realT> ("NitscheFlag", 0.0);

    m_sliding_law = hdn->GetAttributeOrDefault<int>("slidingLaw", 0);

    m_traction_n_tol = hdn->GetAttributeOrDefault<realT>("TractionCutoff", 1e-04);

    if(m_sliding_law == 1)
    {
      m_yield_traction = hdn->GetAttributeOrDefault<realT>("tracYield", 1e16);

      std::cout<<"Contact constitutive behavior: Tresca friction"<<std::endl;
      std::cout<<"Yield traction: "<<m_yield_traction<<std::endl;
    }
    if(m_sliding_law==2)
    {
      m_coulomb_coefficient = hdn->GetAttributeOrDefault<realT>("FrictionCoefficient", 1e16);
      std::cout<<"Contact constitutive behavior: Coulomb's friction"<<std::endl;
      std::cout<<"Coulomb's friction coefficient: "<<m_coulomb_coefficient<<std::endl;
    }

    if (m_nitsche_active == 0.0)
    {
      m_glob_penalty_n = hdn->GetAttributeOrDefault<realT>("PenaltyNormal", 0.0);
      m_glob_penalty_t1 = hdn->GetAttributeOrDefault<realT>("PenaltyTau1", 0.0);
      m_glob_penalty_t2 = hdn->GetAttributeOrDefault<realT>("PenaltyTau2", 0.0);

      if (m_glob_penalty_n == 0.0)
      {
        std::cout << "WARNING: Penalty parameter not set in n^1 direction."<<std::endl;
      }
      if (m_glob_penalty_t1 == 0.0)
      {
        std::cout << "WARNING: Penalty parameter not set in tau^1 direction."<<std::endl;
      }
      if(m_glob_penalty_t2 == 0.0 and threeD)
      {
        std::cout << "WARNING: Penalty parameter not set in tau^2 direction."<<std::endl;
      }
      std::cout<<"*****************************************************************************"<<std::endl;
      std::cout<<"          Penalty method for contact                                         "<<std::endl;
      std::cout<<"          Penalty stiffness in n^1 direction: "<<m_glob_penalty_n<<"             "<<std::endl;
      std::cout<<"          Penalty stiffness in tau^1 direction: "<<m_glob_penalty_t1<<"                   "<<std::endl;
      if (threeD)
        std::cout<<"          Penalty stiffness in tau^2 direction: "<<m_glob_penalty_t2<<"                   "<<std::endl;
      std::cout<<"*****************************************************************************"<<std::endl;

      // If Nitsche's method is not used all Nitsche-specific flags MUST be zero. This line overrides an oversight in the xml file.
      m_nitsche_symmetry_active = 0.0;
      m_weighted_nitsche_active = 0.0;
      m_nitsche_normal_active = 0.0;
      m_nitsche_tau1_active = 0.0;
      m_nitsche_tau2_active = 0.0;
    }
    else
    {
      m_weighted_nitsche_active = hdn->GetAttributeOrDefault<realT> ("NitscheWeightedFlag", 0.0);
      m_nitsche_normal_active = hdn->GetAttributeOrDefault<realT> ("NitscheNormalFlag", 1.0);
      m_nitsche_tau1_active = hdn->GetAttributeOrDefault<realT> ("NitscheTau1Flag", 1.0);
      m_nitsche_tau2_active = hdn->GetAttributeOrDefault<realT> ("NitscheTau2Flag", 1.0);

      if(m_weighted_nitsche_active)
      {
        std::cout<<"*****************************************************************************"<<std::endl;
        std::cout<<"                Weighted Nitsche's method for contact                        "<<std::endl;
        std::cout<<"*****************************************************************************"<<std::endl;
      }
      else
      {
        std::cout<<"*****************************************************************************"<<std::endl;
        std::cout<<"               Classical Nitsche's method for contact                        "<<std::endl;
        std::cout<<"*****************************************************************************"<<std::endl;
      }
      m_nitsche_symmetry_active = hdn->GetAttributeOrDefault<realT> ("NitscheSymmetryFlag", 0.0);
      if(m_nitsche_symmetry_active)
      {
        std::cout<<"Symmetric Nitsche's method"<<std::endl;
      }

      if(!m_nitsche_normal_active)
      {
        m_glob_penalty_n = hdn->GetAttributeOrDefault<realT>("PenaltyNormal", 0.0);

        std::cout<<"Penalty method used for normal contact"<<std::endl;
        if (m_glob_penalty_n == 0.0)
        {
          std::cout << "WARNING: Penalty parameter not set in n^1 direction."<<std::endl;
        }
      }
      if(!m_nitsche_tau1_active)
      {
        m_glob_penalty_t1 = hdn->GetAttributeOrDefault<realT>("PenaltyTau1", 0.0);

        std::cout<<"Penalty method used for tangential contact in tau^1 direction"<<std::endl;
        if (m_glob_penalty_t1 == 0.0)
        {
          std::cout << "WARNING: Penalty parameter not set in tau^1 direction."<<std::endl;
        }
      }
      if(!m_nitsche_tau2_active and threeD)
      {
        m_glob_penalty_t2 = hdn->GetAttributeOrDefault<realT>("PenaltyTau2", 0.0);

        std::cout<<"Penalty method used for tangential contact in tau^2 direction"<<std::endl;
        if (m_glob_penalty_t2 == 0.0)
        {
          std::cout << "WARNING: Penalty parameter not set in tau^2 direction."<<std::endl;
        }
      }

      m_stab_magnifying_factor = hdn->GetAttributeOrDefault<realT> ("MagnifyingFactor", 1.0);
    }
  }
  // @annavarapusr1
}

globalIndex ContactManagerBaseT::insert(const localIndex i, const bool assignGlobals )
{
  globalIndex gi = ObjectDataStructureBaseT::insert(i, assignGlobals);
  m_contact->insert(i);
  return gi;
}
void ContactManagerBaseT::erase( const localIndex i )
{
  ObjectDataStructureBaseT::erase(i);
  m_contact->erase(i);
}
globalIndex ContactManagerBaseT::resize( const localIndex size, const bool assignGlobals )
{
  globalIndex gi = ObjectDataStructureBaseT::resize( size, assignGlobals );
  m_contact->resize(size);
  return gi;
}

void ContactManagerBaseT::WriteVTKCellData(std::ofstream& out)
{
  sArray1d intVarNames;
  sArray1d realVarNames;
  sArray1d R1TensorVarNames;
  sArray1d R2TensorVarNames;
  sArray1d R2SymTensorVarNames;

  Array1dT<iArray1d*> intVars;
  Array1dT<rArray1d*> realVars;
  Array1dT<Array1dT<R1Tensor>*> R1Vars;
  Array1dT<Array1dT<R2Tensor>*> R2Vars;
  Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

  m_contact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_contact->Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  localIndex i = 0;
  for( sArray1d::const_iterator itn=realVarNames.begin() ; itn!=realVarNames.end() ; ++itn, ++i )
  {
    out << "SCALARS ln" << (*itn) << " double" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const rArray1d& scalar = *(realVars[i]);
    for(rArray1d::const_iterator it = scalar.begin(); it != scalar.end(); ++it)
    {
      out << (*it) << " ";
    }
    out << std::endl;
  }
  i = 0;
  for( sArray1d::const_iterator itn=R1TensorVarNames.begin() ; itn!=R1TensorVarNames.end() ; ++itn, ++i )
  {
    out << "SCALARS ln" << (*itn) << "_Magnitude double" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const Array1dT<R1Tensor>& scalar = *(R1Vars[i]);
    for(Array1dT<R1Tensor>::const_iterator it = scalar.begin(); it != scalar.end(); ++it)
    {
      out << it->L2_Norm() << " ";
    }
    out << std::endl;
  }
  i = 0;
  for( sArray1d::const_iterator itn=intVarNames.begin() ; itn!=intVarNames.end() ; ++itn, ++i )
  {
    out << "SCALARS ln" << (*itn) << " integer" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const iArray1d& scalar = *(intVars[i]);
    for(iArray1d::const_iterator it = scalar.begin(); it != scalar.end(); ++it)
    {
      out << (*it) << " ";
    }
    out << std::endl;
  }

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
}

void ContactManagerBaseT::WriteSilo( SiloFile& siloFile,
                                             const std::string& siloDirName,
                                             const std::string& meshname,
                                             const int centering,
                                             const int cycleNum,
                                             const realT problemTime,
                                             const bool isRestart,
                                             const std::string& regionName,
                                             const lArray1d& mask )
{
  std::string subDirectory = siloDirName;
  std::string rootDirectory = "/" + siloDirName;
  siloFile.MakeSubDirectory( subDirectory, rootDirectory );
  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());


  if(isRestart)
  {
    const int varParams = m_contact ? (m_contact->VariableParameters() ? 1 : 0) : 0;
    siloFile.DBWriteWrapper("m_contact_hasVariableParameters", varParams);
  }

  //-----------------------------------------------------
  //HANDLE THE VARIABLES ASSOCIATED WITH THE JOINT STATE
  sArray1d intVarNames;
  sArray1d realVarNames;
  sArray1d R1TensorVarNames;
  sArray1d R2TensorVarNames;
  sArray1d R2SymTensorVarNames;

  Array1dT<iArray1d*> intVars;
  Array1dT<rArray1d*> realVars;
  Array1dT<Array1dT<R1Tensor>*> R1Vars;
  Array1dT<Array1dT<R2Tensor>*> R2Vars;
  Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

  m_contact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_contact->Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  ObjectDataStructureBaseT::WriteSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
  DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
  DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

  WriteNonManagedDataMembersToSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DBSetDir(siloFile.m_dbFilePtr, "..");
}

void ContactManagerBaseT::ReadSilo( const SiloFile& siloFile,
                                            const std::string& siloDirName,
                                            const std::string& meshname,
                                            const int centering,
                                            const int cycleNum,
                                            const realT problemTime,
                                            const bool isRestart,
                                            const std::string& regionName,
                                            const lArray1d& mask )
{
  if( DBSetDir(siloFile.m_dbFilePtr, siloDirName.c_str()) != -1 )
  {

    if(isRestart)
    {
      int varParams;
      siloFile.DBReadWrapper("m_contact_hasVariableParameters", varParams);
      if(m_contact)
        m_contact->SetVariableParameters(varParams == 1);
    }

    sArray1d intVarNames;
    sArray1d realVarNames;
    sArray1d R1TensorVarNames;
    sArray1d R2TensorVarNames;
    sArray1d R2SymTensorVarNames;

    Array1dT<iArray1d*> intVars;
    Array1dT<rArray1d*> realVars;
    Array1dT<Array1dT<R1Tensor>*> R1Vars;
    Array1dT<Array1dT<R2Tensor>*> R2Vars;
    Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

    m_contact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

    AllocateTemporaryFields( intVarNames, intVars );
    AllocateTemporaryFields( realVarNames, realVars );
    AllocateTemporaryFields( R1TensorVarNames, R1Vars );
    AllocateTemporaryFields( R2TensorVarNames, R2Vars );
    AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

    ObjectDataStructureBaseT::ReadSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);

    ReadNonManagedDataMembersFromSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );

    m_contact->Deserialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

    DeallocateTemporaryFields<int>(intVarNames);
    DeallocateTemporaryFields<realT>(realVarNames);
    DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
    DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
    DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}



/** @brief Get the index in the list associated with the two body or face indices; -1 if not found
 * @param body1 Index of body 1 or face 1
 * @param body2 Index of body 2 or face 2
 * @return Index of the common plane
 */
int ContactManagerBaseT::GetIndex(const localIndex body1, const localIndex body2) const
{
  const lArray1d& s1 = this->GetFieldData<localIndex> ("face1");
  const lArray1d& s2 = this->GetFieldData<localIndex> ("face2");
  int imin = -1;
  for (localIndex i = 0; i < this->DataLengths(); ++i)
  {
    if (imin < 0)
    {
      if (s1[i] == body1)
      {
        imin = i;
        if (s2[i] == body2)
          return i;
      }
    }
    else
    {
      if (s1[i] != body1)
      {
        return -1;
      }
      else if (s2[i] == body2)
      {
        return i;
      }
    }
  }
  return -1;
}

/**
 * @brief Update the contact list using the given neighbor list
 * @author Scott Johnson
 * Sets contacts not in the neighbor list to "inactive" and adds in
 * any new entries from the given neighbor list
 * @param[in] neighborList List of external face indices associated with a particular external face
 */
size_t ContactManagerBaseT::Update(const Array1dT< lArray1d >& neighborList)
{
  lArray1d& face1 = this->GetFieldData<localIndex>("face1");
  lArray1d& face2 = this->GetFieldData<localIndex>("face2");
  iArray1d& activeC = this->GetFieldData<int>("active");
    
  //-------------------------------------------------
  //case A: neighbor list is empty, so remove all contacts
  if( neighborList.empty() )
  {
    // neighborList is empty, so clear out all common planes
    this->resize(0);
    return this->DataLengths();
  }

  //-------------------------------------------------
  //case B: common plane list is empty, so populate with neighbor list entries
  if( this->DataLengths() == 0 )
  {
    // contact manager contains no common planes, so copy the neighbor list completely

    // B.1: allocate space for the new common planes
    localIndex newsize = 0;
    for( Array1dT<lArray1d>::const_iterator i=neighborList.begin() ; i!=neighborList.end() ; ++i )
      newsize += i->size();
    this->resize( newsize );

    // B.2: now iterate over the neighborList and fill the common plane face entries
    localIndex icount = 0, ii = 0;
    for( Array1dT<lArray1d>::const_iterator i=neighborList.begin() ; i!=neighborList.end() ; ++i, ++icount)
    {
      // iterate over each entry in the neighborList for a given face1 value
      for( lArray1d::const_iterator j=i->begin() ; j!=i->end() ; ++j, ++ii)
      {
        face1[ii] = icount;
        face2[ii] = *j;
      }
    }

    //B.3: return result
    return this->DataLengths();
  }

  //-------------------------------------------------
  //case C: update contact list with neighbor list

  // define a counter to keep track of the index on the neighborList array

  // get the number of face pairs contained in the neighborList
  size_t numFacePairs = 0;
  for(size_t nl_indexf1 = 0; nl_indexf1<neighborList.size() ; ++nl_indexf1 )
    numFacePairs += neighborList[nl_indexf1].size();

  // check to see whether there are any face pairs in the neighbor list; return if there are none.
  if( numFacePairs == 0 )
  {
    // no face pairs..so no common planes
    this->resize(0);
    return this->DataLengths();
  }

  // there are face pairs defined in the neighbor list, so we have to check them against existing common planes.
  // common plane index (the index of the data in *this)
  {
    lArray1d::size_type cp_index = 0;
    for(Array1dT<lArray1d>::size_type nl_indexf1 = 0; nl_indexf1<neighborList.size() ; ++nl_indexf1 )
    {
      lArray1d::size_type i2=0;
      while(i2<neighborList[nl_indexf1].size())
      {
        lArray1d::size_type nl_indexf2 = neighborList[nl_indexf1][i2];

        //-------------------------------------------------
        //case C.0: **APPEND NEW** to end of current list's faces with remaining new entry's faces
        if( cp_index==this->DataLengths() )
        {
          // we got to the end of the list of existing common planes, so everything else in the neighborList needs to be
          // tacked on to the end....boy...a push_back() function would be nice.
          this->resize( cp_index+1 );
          face1[cp_index] = nl_indexf1;
          face2[cp_index] = nl_indexf2;
//          if(activeC[cp_index]!=0)
//          {
//            char buffer[100];
//            sprintf(buffer, "RESIZE Not 0 %d\n",activeC[cp_index]);
//            throw GPException(buffer);
//          }
          activeC[cp_index] = 0;
        }

        //-------------------------------------------------
        // case C.1: **UPDATE CURRENT** new face pair = current face pair, so just update
        else if( face1[cp_index]==nl_indexf1 && face2[cp_index]==nl_indexf2 )
        {
          //just go to the next ... nothing to do
        }

        //-------------------------------------------------
        // case C.2: **INSERT NEW** current lists's faces > new faces
        else if( face1[cp_index]>nl_indexf1 || ( face1[cp_index]==nl_indexf1 && face2[cp_index]>nl_indexf2 ) )
        {
          // the current list of common planes do not contain face pair specified in neighborList,
          // so we should INSERT the face pair from neighborList at the location of cp_index
          this->insert(cp_index);
          face1[cp_index] = nl_indexf1;
          face2[cp_index] = nl_indexf2;
//          if(activeC[cp_index]!=0)
//          {
//            char buffer[100];
//            sprintf(buffer, "INSERT Not 0 %d\n",activeC[cp_index]);
//            throw GPException(buffer);
//          }
          activeC[cp_index] = 0;
        }

        //-------------------------------------------------
        //case C.3: **REMOVE OLD** current list's faces < new entry's faces
        else if( face1[cp_index] < nl_indexf1 || (face1[cp_index]==nl_indexf1 && face2[cp_index] < nl_indexf2 ))
        {
          //the entry in the current list is not in the new list
          //so in the future ... (1) set the current list's entry to inactive
          //                 and (2) increment the current list to the next entry
          // ... for now ...
          //let's just take any "inactive" contact out ... this can be changed later
          this->erase( cp_index );
          //NOTE: we need to keep the neighbor list entry constant ... don't march on yet!!!
          continue;//necessary to prevent marching on in "i2"
        }

        //-------------------------------------------------
        // case C.default: something unanticipated happened
        else
        {
          throw GPException("ContactManagerBaseT::Update - case unaccounted for");
        }

        //if not REMOVE OLD or exception, march to next entries
        ++i2;
        ++cp_index;
      }
    }
  }

  //-------------------------------------------------
  //case D: error handling!

  if(this->DataLengths() < numFacePairs)
  {
    size_t neighbor_index = 0;
    std::cout << "CONTACTS: " << this->DataLengths() << " FACEPAIRS: " << numFacePairs << "\n";
    for( Array1dT<lArray1d>::size_type nl_indexf1=0 ; nl_indexf1<neighborList.size() ; ++nl_indexf1 )
    {
      for( lArray1d::size_type i2=0 ; i2<neighborList[nl_indexf1].size() ; ++i2, ++neighbor_index)
      {
        if(neighbor_index >= this->DataLengths())
          throw GPException("Failure in determining contact pairs: number in neighbor list does not equal number of contact!");
        std::cout << "index: " << neighbor_index << " face1: " << face1[neighbor_index] << " " << nl_indexf1 << " face2: " << face2[neighbor_index] << " " << neighborList[nl_indexf1][i2] << "\n";
        std::cout.flush();
      }
    }
  }

  return this->DataLengths();
}

