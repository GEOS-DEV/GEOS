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
/*
 * EllipsoidalDiscreteElementManagerT.cpp
 *
 *  Created on: July 14, 2011
 *      Author: Scott Johnson
 */
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "DiscreteElementManagerBaseT.h"
#include "EllipsoidalDiscreteElementManagerT.h"
#include "EllipsoidalContactManagerT.h"

#include "Constitutive/Interface/InterfaceFactory.h"

/**
 * @brief Constructor to set initial states
 * @author Scott Johnson
 */
EllipsoidalDiscreteElementManagerT::EllipsoidalDiscreteElementManagerT() :
  DiscreteElementManagerBaseT(ObjectDataStructureBaseT::EllipsoidalDiscreteElementManager),
  m_neighborList(m_VariableOneToManyMaps["neighborList"]),
  m_neighborListInverse(m_UnorderedVariableOneToManyMaps["neighborListInverse"]),
  sorted(false),
  m_sorter( NULL ),
  m_contactBufferOffset(0.0),
  m_contactBufferFactor(0.0),
  m_boundaryContact(NULL)
{
  m_boundary.m_mat = NULL;
  m_boundary.Reset();

  //add fields for rotation and translation
  this->AddBaseFields();

  //geometry
  this->AddKeylessDataField<R1Tensor>("principalRadii", true, true);

  //sorting
  this->AddKeylessDataField<realT>("boundingRadius", false, false);
  this->AddKeylessDataField<realT>("boundingRadiusLastSort", false, false);
  this->AddKeylessDataField<R1Tensor>("currentPositionLastSort", false, false);

  //boundary contact
  this->AddKeylessDataField<R1Tensor>("contactPointBoundary", true, true);
  this->AddKeylessDataField<R1Tensor>("shearSlipBoundary", true, true);
  this->AddKeylessDataField<R1Tensor>("normalBoundary", true, true);
  this->AddKeylessDataField<int>("activeBoundary", true, true);
  this->AddKeylessDataField<realT>("normalApproachBoundary", true, true);

  m_boundaryContact = InterfaceFactory::NewInterface("Hertzian");
  m_sorter = SpatialSorting::SpatialSorterFactory::NewSpatialSorter("N2");
}

EllipsoidalDiscreteElementManagerT::~EllipsoidalDiscreteElementManagerT()
{
  if(m_boundaryContact)
    delete m_boundaryContact;

#if USECPP11!=1
  if(m_boundary.m_mat)
    delete m_boundary.m_mat;
#endif

  if(m_sorter)
    delete m_sorter;
}

globalIndex EllipsoidalDiscreteElementManagerT::insert(const localIndex i, const bool assignGlobals )
{
  globalIndex gi = DiscreteElementManagerBaseT::insert(i, assignGlobals);
  if(m_boundary.Active())
    m_boundaryContact->insert(i);
  return gi;
}
void EllipsoidalDiscreteElementManagerT::erase( const localIndex i )
{
  DiscreteElementManagerBaseT::erase(i);
  if(m_boundary.Active())
    m_boundaryContact->erase(i);
}
globalIndex EllipsoidalDiscreteElementManagerT::resize( const localIndex size, const bool assignGlobals )
{
  globalIndex gi = DiscreteElementManagerBaseT::resize( size, assignGlobals );
  if(m_boundary.Active())
    m_boundaryContact->resize(size);
  return gi;
}

void EllipsoidalDiscreteElementManagerT::ReadXML(TICPP::HierarchicalDataNode* EDEMNode)
{
  //set boundary material
  m_boundary.ReadXML(EDEMNode);

  //set the boundary contact
  if(m_boundary.Active())
  {
    TICPP::HierarchicalDataNode* EDemNodeContact0 = EDEMNode->GetChild("Contact");
    if(EDemNodeContact0 == NULL)
      throw GPException("Must define a boundary contact model for Elliposidal DEM if boundary is active");
    TICPP::HierarchicalDataNode* EDemNodeContact1 = EDemNodeContact0->Next(true);
    if(EDemNodeContact1 == NULL)
      throw GPException("Must define a boundary contact model for Elliposidal DEM (1) if boundary is active");

    if(m_boundaryContact)
      delete m_boundaryContact;
    m_boundaryContact = InterfaceFactory::NewInterface(EDemNodeContact1->Heading(), EDemNodeContact1);
  }

  //set material
  {
    TICPP::HierarchicalDataNode* EDemNodeMaterial0 = EDEMNode->GetChild("Material");
    if(EDemNodeMaterial0 == NULL)
      throw GPException("Must define a material model for Elliposidal DEM");
    TICPP::HierarchicalDataNode* EDemNodeMaterial1 = EDemNodeMaterial0->Next(true);
    if(EDemNodeMaterial1 == NULL)
      throw GPException("Must define a material model for Elliposidal DEM (1)");

#if USECPP11!=1
    if(m_mat)
      delete m_mat;
#endif
    m_mat = MaterialFactory::NewMaterial(EDemNodeMaterial1->Heading(), EDemNodeMaterial1);
    m_mat->SetVariableParameters(true);
    m_mat->resize(DataLengths(), 1);
    for(localIndex i = 0; i < DataLengths(); i++)
    {
      m_mat->ParameterData(i)->ReadXML(*EDemNodeMaterial1);
    }
  }
}

void EllipsoidalDiscreteElementManagerT::WriteVTKPointData(std::ofstream& out)
{
  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  array<array<integer>*> intVars;
  array<array<real64>*> realVars;
  array<array<R1Tensor>*> R1Vars;
  array<array<R2Tensor>*> R2Vars;
  array<array<R2SymTensor>*> R2SymVars;

  m_boundaryContact->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_boundaryContact->Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  DiscreteElementManagerBaseT::WriteVTKPointData(out);

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
}

void EllipsoidalDiscreteElementManagerT::WriteSilo( SiloFile& siloFile,
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
    const int varParams = m_mat ? (m_mat->VariableParameters() ? 1 : 0) : 0;
    siloFile.DBWriteWrapper("m_mat_hasVariableParameters", varParams);
  }

  //-----------------------------------------------------
  //HANDLE THE VARIABLES ASSOCIATED WITH THE JOINT STATE
  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  array<array<integer>*> intVars;
  array<array<real64>*> realVars;
  array<array<R1Tensor>*> R1Vars;
  array<array<R2Tensor>*> R2Vars;
  array<array<R2SymTensor>*> R2SymVars;

  m_mat->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_mat->Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  //avoid writing R2SymVars to SILO ... silodiff doesn't understand them
  DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);
  ObjectDataStructureBaseT::WriteSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
  DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);

  WriteNonManagedDataMembersToSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DBSetDir(siloFile.m_dbFilePtr, "..");
}

void EllipsoidalDiscreteElementManagerT::ReadSilo( const SiloFile& siloFile,
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
      siloFile.DBReadWrapper("m_mat_hasVariableParameters", varParams);
      if(m_mat)
        m_mat->SetVariableParameters(varParams == 1);
    }

    array<string> intVarNames;
    array<string> realVarNames;
    array<string> R1TensorVarNames;
    array<string> R2TensorVarNames;
    array<string> R2SymTensorVarNames;

    array<array<integer>*> intVars;
    array<array<real64>*> realVars;
    array<array<R1Tensor>*> R1Vars;
    array<array<R2Tensor>*> R2Vars;
    array<array<R2SymTensor>*> R2SymVars;

    m_mat->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

    AllocateTemporaryFields( intVarNames, intVars );
    AllocateTemporaryFields( realVarNames, realVars );
    AllocateTemporaryFields( R1TensorVarNames, R1Vars );
    AllocateTemporaryFields( R2TensorVarNames, R2Vars );

    ObjectDataStructureBaseT::ReadSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);

    ReadNonManagedDataMembersFromSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );

    //avoid reading R2SymVars from SILO ... silodiff doesn't understand them
    AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );
    m_mat->Deserialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

    DeallocateTemporaryFields<int>(intVarNames);
    DeallocateTemporaryFields<realT>(realVarNames);
    DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
    DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
    DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}


/**
 * @brief Recalculate the moment of inertia for each ellipsoid
 * @author Scott Johnson
 */
void EllipsoidalDiscreteElementManagerT::RecalculatePhysicalProperties()
{
  array<realT>& mass = GetFieldData<FieldInfo::mass>();
  const array<R1Tensor>& radii = GetFieldData<R1Tensor>("principalRadii");
  array<R1Tensor>& I = GetFieldData<FieldInfo::rotationalInertia>();

  for (localIndex i = 0; i < this->DataLengths(); i++)
  {
    mass[i] = m_mat->ParameterData(i)->init_density * Volume(radii[i]);
    MomentOfInertia(radii[i], mass[i], I[i]);
  }
}

/**
 * @brief Recalculate the list of closest neighbors if necessary
 * @author Scott Johnson
 * @param[in] dt Current timestep
 * @return Whether a resort was performed
 */
bool EllipsoidalDiscreteElementManagerT::RecalculateNeighborList(const realT dt)
{
  bool sort = false;
  localIndex num = this->DataLengths();
  if(num == 0)
    return sort;

  const array<R1Tensor>& pradii = this->GetFieldData<R1Tensor>("principalRadii");

  const array<R1Tensor>& disp = this->GetFieldData<FieldInfo::displacement>();
  const array<R1Tensor>& refpos = this->GetFieldData<FieldInfo::referencePosition>();
  array<R1Tensor>& centers = this->GetFieldData<FieldInfo::currentPosition>();

  const array<R1Tensor>& vel = this->GetFieldData<FieldInfo::velocity>();
  array<real64>& radii = this->GetFieldData<realT>("boundingRadius");

  if(this->sorted)
  {
    const array<real64>& lradii = this->GetFieldData<realT>("boundingRadiusLastSort");
    const array<R1Tensor>& lcenters = this->GetFieldData<R1Tensor>("currentPositionLastSort");
    R1Tensor dx;

    //iterate through
    for (localIndex kf0 = 0; kf0 < this->DataLengths(); ++kf0)
    {
      centers[kf0] = refpos[kf0];
      centers[kf0] += disp[kf0];

      //get new radius of influence
      radii[kf0] = pradii[kf0].MaxVal();
      radii[kf0] += dt * vel[kf0].L2_Norm();
      realT buffer = this->m_contactBufferFactor * lradii[kf0] + this->m_contactBufferOffset;
      radii[kf0] += buffer;

      //get new distance between potentials
      dx = centers[kf0];
      dx -= lcenters[kf0];

      //get the square of the centroid-to-centroid distance
      realT dd = Dot(dx, dx);

      //get the allowable travel distance - squared
      realT dd0 = buffer - (radii[kf0] - lradii[kf0]);
      dd0 *= dd0;
      if(dd > dd0)
        sort = true;
    }
  }
  else
  {
    sort = true;
    sorted = true;

    //iterate through
    for (localIndex kf0 = 0; kf0 < this->DataLengths(); ++kf0)
    {
      centers[kf0] = refpos[kf0];
      centers[kf0] += disp[kf0];

      //get new radius of influence
      radii[kf0] = pradii[kf0].MaxVal();
      radii[kf0] += dt * vel[kf0].L2_Norm();
      radii[kf0] *= 1.0 + this->m_contactBufferFactor;
      radii[kf0] += this->m_contactBufferOffset;
    }
  }
  if(!sort)
    return false;

  //reset the values of the last sort
  {
    array<real64>& lradii = this->GetFieldData<realT>("boundingRadiusLastSort");
    array<R1Tensor>& lcenters = this->GetFieldData<R1Tensor>("currentPositionLastSort");
    std::copy(centers.begin(), centers.end(), lcenters.begin());
    std::copy(radii.begin(), radii.end(), lradii.begin());
  }

//  //
//  if(num < 100)
//    Sort = &SpatialSorting::NSquared::Sort;
//  else
//    Sort = &SpatialSorting::CellVerlet::Sort;

  m_sorter->Sort(radii, centers, this->m_neighborList, this->m_neighborListInverse);

  return true;
}

/**
 * @brief Update the contact properties and then the forces/stresses from the contact
 * @author Scott Johnson
 * @param[in] dt Current timestep
 * @param[in] ecm Ellipsoidal contact manager object
 */
void EllipsoidalDiscreteElementManagerT::UpdateAndApplyContactStresses( StableTimeStep& maxdt, const realT dt,
                                                                        EllipsoidalContactManagerT& ecm)
{
  //TODO make not so hard-wired for the boundary ... need to more gracefully handle boundaries in general
  //this->m_boundaryContactParameters = ecm.m_contactParameters;
  if(this->m_DataLengths == 0)
    return;

  std::map<std::string, int>          pi0,   pi1;
  std::map<std::string, realT>        pr0,   pr1;
  std::map<std::string, R1Tensor>     pr10,  pr11;
  std::map<std::string, R2Tensor>     pr20,  pr21;
  std::map<std::string, R2SymTensor>  pr2s0, pr2s1;

  const array<R1Tensor>& pradii = this->GetFieldData<R1Tensor>("principalRadii");
  const array<R1Tensor>& centers = this->GetFieldData<FieldInfo::currentPosition>();
  const array<R1Tensor>& velocity = this->GetFieldData<FieldInfo::velocity>();
  array<R1Tensor>& frc = this->GetFieldData<FieldInfo::force>();
  array<R1Tensor>& mom = this->GetFieldData<FieldInfo::moment>();
  const array<real64>& mass = this->GetFieldData<FieldInfo::mass>();

  //contact points
  array<R1Tensor>& contactPointsC = ecm.GetFieldData<R1Tensor>("applicationPoint");
  array<R1Tensor>& normalC = ecm.GetFieldData<R1Tensor>("normal");
  array<integer>& activeC = ecm.GetFieldData<int>("active");
  array<real64>& normalApproachC = ecm.GetFieldData<realT>("normalApproach");
  array<R1Tensor>& velocityC = ecm.GetFieldData<R1Tensor>("velocity");
  array<R1Tensor>& shearSlipC = ecm.GetFieldData<R1Tensor>("shearSlip");

  //-----ITERATE NEIGHBORLIST AND UPDATE WHETHER IT IS ACTIVE, ETC-----
  localIndex index = 0;
  {
    bool isActive;
    R1Tensor dx, cpLocal;
    R2Tensor rotation0, rotation1;
    for (array<lArray1d>::size_type kf0 = 0; kf0 < this->m_neighborList.size(); ++kf0)
    {
      this->RotationTensor(kf0, rotation0);
      for (lArray1d::size_type it = 0; it < this->m_neighborList[kf0].size(); ++it, ++index)
      {
        const localIndex kf1 = this->m_neighborList[kf0][it];
        this->RotationTensor(kf1, rotation1);

        InterfaceBaseStateData& state = *ecm.m_contact->StateData(index,0);

        R1Tensor& normal = normalC[index];
        R1Tensor& contactPoint = contactPointsC[index];
        R1Tensor& vc = velocityC[index];

        isActive = InContact(pradii[kf0], centers[kf0], rotation0, pradii[kf1], centers[kf1],
                             rotation1, state.normalApproach, normal, contactPoint);
        activeC[index] = isActive ? 1 : 0;
        normalApproachC[index] = state.normalApproach;

        vc = velocity[kf1];
        vc -= velocity[kf0];

        if (isActive)
        {
          this->GlobalToLocal(kf0, contactPoint, cpLocal);
          const realT r0 = GaussianCurvature(cpLocal, pradii[kf0]);
          this->GlobalToLocal(kf1, contactPoint, cpLocal);
          const realT r1 = GaussianCurvature(cpLocal, pradii[kf1]);

          {
            m_mat->ParameterData(kf0)->GetVariableValues(pi0, pr0, pr10, pr20, pr2s0);
            m_mat->ParameterData(kf1)->GetVariableValues(pi1, pr1, pr11, pr21, pr2s1);
            pr0["radius"] = r0;
            pr1["radius"] = r1;
            pr0["mass"] = mass[kf0];
            pr1["mass"] = mass[kf1];
            ecm.m_contact->UpdateProperties(index, pr0, pr1);
          }

          R1Tensor dxs;
          state.UpdateOrientation(normalApproachC[index], dt, normalC[index], velocityC[index], dxs,
                                   shearSlipC[index]);

          //-----------------------------------------------------------
          //CONTACT 4: ONLY FOR FE-DE OR FE-FE CONTACT
          //           resolve the stresses
          ecm.m_contact->StrainDrivenUpdate(index);

          R1Tensor force(normal);
          force *= state.stress;
          force += state.stressShearVector;

          DiscreteElementManagerBaseT::ForceToCouple(contactPoint, force, centers[kf0], rotation0,
                                                     frc[kf0], mom[kf0]);
          force *= -1;
          DiscreteElementManagerBaseT::ForceToCouple(contactPoint, force, centers[kf1], rotation1,
                                                     frc[kf1], mom[kf1]);

          //----GET THE STABLE TIMESTEP ESTIMATE ----
          maxdt.SetIfSmaller(
              StableTimestep(
                  mass[kf0],
                  mass[kf1],
                  ecm.m_contact->StiffnessProjected(index)));
          //-----------------------------------------

          //UPDATE THE SHEAR SLIP
          shearSlipC[index] += dxs;
        }
      }
    }
  }
}

void EllipsoidalDiscreteElementManagerT::UpdateCylindricalBoundary(StableTimeStep& maxdt, const realT dt, const realT time)
{
  if (this->DataLengths() == 0 || !m_boundary.Active())
    return;

  std::map<std::string, int>          pi0,   pi1;
  std::map<std::string, realT>        pr0,   pr1;
  std::map<std::string, R1Tensor>     pr10,  pr11;
  std::map<std::string, R2Tensor>     pr20,  pr21;
  std::map<std::string, R2SymTensor>  pr2s0, pr2s1;

  m_boundary.Update(time);

  const array<R1Tensor>& pradii = this->GetFieldData<R1Tensor>("principalRadii");
  const array<R1Tensor>& centers = this->GetFieldData<FieldInfo::currentPosition>();
  const array<R1Tensor>& velocity = this->GetFieldData<FieldInfo::velocity>();
  array<R1Tensor>& frc = this->GetFieldData<FieldInfo::force>();
  array<R1Tensor>& mom = this->GetFieldData<FieldInfo::moment>();
  const array<real64>& mass = this->GetFieldData<FieldInfo::mass>();

  //Boundary contact
  array<R1Tensor>& shearSlips = this->GetFieldData<R1Tensor>("shearSlipBoundary");
  array<R1Tensor>& contactPoints = this->GetFieldData<R1Tensor>("contactPointBoundary");
  array<R1Tensor>& normals = this->GetFieldData<R1Tensor>("normalBoundary");
  array<real64>& normalApproaches = this->GetFieldData<realT>("normalApproachBoundary");
  array<integer>& actives = this->GetFieldData<int>("activeBoundary");

  //-----ITERATE SHAPE LIST FOR BOUNDARY CONTACT AND UPDATE WHETHER IT IS ACTIVE, ETC-----
  {
    m_boundary.m_mat->ParameterData(0)->GetVariableValues(pi0, pr0, pr10, pr20, pr2s0);
    pr0["radius"] = m_boundary.radiusCylinder;

    bool isActive;
    R1Tensor dx, cpLocal;
    R2Tensor rotation;
    for (localIndex a = 0; a < this->m_DataLengths; ++a)
    {
      this->RotationTensor(a, rotation);

      InterfaceBaseStateData& state = *m_boundaryContact->StateData(a, 0);

      R1Tensor& contactPoint = contactPoints[a];
      R1Tensor& normal = normals[a];

      isActive = InContact(pradii[a], centers[a], rotation, state.normalApproach, normal,
                           contactPoint);

      normalApproaches[a] = state.normalApproach;
      actives[a] = isActive ? 1 : 0;

      if (isActive)
      {
        this->GlobalToLocal(a, contactPoint, cpLocal);
        const realT radius = GaussianCurvature(cpLocal, pradii[a]);

        {
          //get material properties for this pair
          m_mat->ParameterData(a)->GetVariableValues(pi1, pr1, pr11, pr21, pr2s1);
          pr1["radius"] = radius;
          pr0["mass"] = pr1["mass"] = mass[a];
          m_boundaryContact->UpdateProperties(a, pr0, pr1);
        }

        R1Tensor dxs;
        state.UpdateOrientation(normalApproaches[a], dt, normals[a], velocity[a], dxs,
                                 shearSlips[a]);

        //-----------------------------------------------------------
        //CONTACT 4: ONLY FOR FE-DE OR FE-FE CONTACT
        //           resolve the stresses
        m_boundaryContact->StrainDrivenUpdate(a);

        R1Tensor force(normal);
        force *= state.stress;
        force += state.stressShearVector;

        DiscreteElementManagerBaseT::ForceToCouple(contactPoint, force, centers[a], rotation,
                                                   frc[a], mom[a]);

        //----GET THE STABLE TIMESTEP ESTIMATE ----
        maxdt.SetIfSmaller(
            StableTimestep(
                mass[a],
                mass[a],
                m_boundaryContact->StiffnessProjected(a)));

        //-----------------------------------------

        //UPDATE THE SHEAR SLIP
        shearSlips[a] += dxs;
      }
    } //for each ellipsoid against the cylindrical boundary
  }

//  //TODO: add ability to iterate to find a particular applied stress
//  //this would involve Newton solver on contact with a sub-set of particles
//  //determined during the above loop (probably just grab those particles within a certain distance of
//  //the current axial bound)
//  if (!isZero(this->stressRateCylinder))
//  {
//    throw GPException("Cannot yet handle stress rate application in uniaxial cylinder");
//    const realT err = 1.0e-12;
//    realT tol = 0.1 * err; //TODO: change this to an error metric
//    while (tol > err)
//    {
//      //do Newton iteration here
//      //need to have ordered list of particles according to fabs(dx1) descending
//      //to complete the procedure
//    }
//  }
}

bool EllipsoidalDiscreteElementManagerT::InContact(const R1Tensor& radii,
                                                   const R1Tensor& center,
                                                   const R2Tensor&,
                                                   realT& overlap,
                                                   R1Tensor& normal1,
                                                   R1Tensor& contactPoint)
{
  if(!isEqual(radii(0),radii(1)) ||
      !isEqual(radii(1),radii(2)))
    throw GPException("Cannot current handle non-spherical bodies in contact with cylindrical boundaries");

  const realT radius = radii.MaxVal();
  return m_boundary.InContact(radius, center, overlap, normal1, contactPoint);
}

bool EllipsoidalDiscreteElementManagerT::InContact(const R1Tensor& radii0,
                                                    const R1Tensor& center0,
                                                    const R2Tensor& rotation0,
                                                    const R1Tensor& radii1,
                                                    const R1Tensor& center1,
                                                    const R2Tensor& rotation1,
                                                    realT& overlap,
                                                    R1Tensor& normal1,
                                                    R1Tensor& contactPoint)
{
  //exclude obviously out of contact case
  R1Tensor dist(center0);
  dist -= center1;
  const realT rsum = radii0.MaxVal() + radii1.MaxVal();
  {
    if(rsum*rsum <= Dot(dist, dist))
      return false;
  }

  //GET CONTACT POINT
  contactPoint = dist;
  contactPoint *= radii1(0) / rsum;
  contactPoint += center1;

  //handle spheres
  if(isEqual(radii0(0),radii1(0)) &&
      isEqual(radii0(1),radii1(1)) &&
      isEqual(radii0(2),radii1(2)) &&
      isEqual(radii0(0),radii1(1)))
  {
    //GET NORMAL (1 to 0)
    normal1 = dist;
    normal1.Normalize();

    //GET OVERLAP
    overlap = rsum - dist.L2_Norm();
    return true;
  }

  //ITERATE SOLUTION FOR ELLIPSOID CONTACT
  {
    //get initial estimates
    R1Tensor contactPoint0, contactPoint1;
    PointAtSurface(radii0, center0, rotation0, contactPoint, contactPoint0);
    PointAtSurface(radii1, center1, rotation1, contactPoint, contactPoint1);

    //iterate solution
    {
      const realT tol = 1e-5;
      realT error = std::numeric_limits<realT>::max();
      realT lastPotential = PotentialValue(radii0, center0, rotation0,contactPoint);

      while(error>tol)
      {
        //update normal to face 1 and contactPoint
        const realT currentPotential = Iterate(radii0, center0, rotation0,
                                               radii1, center1, rotation1,
                                               contactPoint0, contactPoint1);
        error = currentPotential - lastPotential;
        lastPotential = currentPotential;
      }
    }

    //GET CONTACT POINT
    contactPoint = contactPoint0;
    contactPoint += contactPoint1;
    contactPoint *= 0.5;

    //GET NORMAL
    {
      R1Tensor tmp(contactPoint);
      tmp -= center1;
      R1Tensor local;
      local.AijBj(rotation1, tmp);
      NormalCartesian(radii1, local, tmp);
      normal1.AijBi(rotation1, tmp);
    }

    //GET OVERLAP
    dist = contactPoint0;
    dist -= contactPoint1;
    overlap = dist.L2_Norm() * (Dot(normal1, dist) > 0 ? -1.0 : 1.0);
  }
  return overlap > 0;
}

/**
 * @brief Iterate the ellipsoid contact solution
 * @author Scott Johnson
 *
 * @param[in] radii0 Principal radii of the first ellipsoid
 * @param[in] center0 Centroid of the first ellipsoid
 * @param[in] rotation0 Rotation tensor of the first ellipsoid
 * @param[in] radii1 Principal radii of the second ellipsoid
 * @param[in] center1 Centroid of the second ellipsoid
 * @param[in] rotation1 Rotation tensor of the second ellipsoid
 * @param[in,out] normal Normal to the surface of the second ellipsoid
 * @param[in,out] contactPoint Current contact point estimate
 * @return New potential evaluation
 */

realT EllipsoidalDiscreteElementManagerT::Iterate(const R1Tensor& radii0,
                                                  const R1Tensor& center0,
                                                  const R2Tensor& rotation0,
                                                  const R1Tensor& radii1,
                                                  const R1Tensor& center1,
                                                  const R2Tensor& rotation1,
                                                  R1Tensor& contactPoint0,
                                                  R1Tensor& contactPoint1)
{
  //use the points to get the new estimates of the normals to
  //surfaces 0 and 1
  R1Tensor normal0, normal1;
  NormalCartesian(radii0, center0, rotation0, contactPoint0, normal0);
  NormalCartesian(radii1, center1, rotation1, contactPoint1, normal1);
  normal1 -= normal0;
  normal1.Normalize();
  PointAtNormalCartesian(radii1, center1, rotation1, normal1, contactPoint1);
  normal1 *= -1.0;
  PointAtNormalCartesian(radii0, center0, rotation0, normal1, contactPoint0);
  R1Tensor contactPoint(contactPoint1);
  contactPoint += contactPoint0;
  contactPoint *= 0.5;
  const realT potentialValue = PotentialValue(radii1, center1, rotation1, contactPoint);
  PointAtSurface(radii0, center0, rotation0, contactPoint, contactPoint0);
  PointAtSurface(radii1, center1, rotation1, contactPoint, contactPoint1);
  return potentialValue;
}




/**
 * @brief Zero out the discrete element forces and accelerations
 * @author Scott Johnson
 */
void EllipsoidalDiscreteElementManagerT::UpdateNodalStatesZeroForcesAndAccelerations()
{
  //zero the DE accelerations and forces
  array<R1Tensor>& acceleration = this->GetFieldData<FieldInfo::acceleration> ();
  array<R1Tensor>& force = this->GetFieldData<FieldInfo::force> ();

  //zero the DE rotation accelerations and moments
  array<R1Tensor>& racc = this->GetFieldData<FieldInfo::rotationalAcceleration> ();
  array<R1Tensor>& moment = this->GetFieldData<FieldInfo::moment> ();

  acceleration = 0.0;
  force = 0.0;
  racc = 0.0;
  moment = 0.0;
}

/**
 * @brief Moment of inertia of a triaxial ellipsoid in local frame
 * @author Scott Johnson
 * @param[in] radii Principal radii of the ellipsoid
 * @param[in] mass Mass of the ellipsoid
 * @param[out] moi Moment of inertia of the ellipsoid in local frame
 */
void EllipsoidalDiscreteElementManagerT::MomentOfInertia(const R1Tensor& radii, const realT mass, R1Tensor& moi)
{
  const realT fct = 0.2 * mass;
  moi(0) = fct * (radii(1) * radii(1) + radii(2) * radii(2));
  moi(1) = fct * (radii(0) * radii(0) + radii(2) * radii(2));
  moi(2) = fct * (radii(0) * radii(0) + radii(1) * radii(1));
}







/**
 * @brief Get the point projected to the ellipsoid surface in the global frame
 * @author Scott Johnson
 * Description
 * @param[in] radii Principal radii
 * @param[in] center Center of the ellipsoid
 * @param[in] rotation Rotation tensor for the ellipsoid
 * @param[in] point Point in global frame
 * @param[out] pointAtSurface Point on the potential surface in the global frame
 */
void EllipsoidalDiscreteElementManagerT::PointAtSurface( const R1Tensor& radii,
                                                         const R1Tensor& center,
                                                         const R2Tensor& rotation,
                                                         const R1Tensor& point,
                                                         R1Tensor& pointAtSurface)
{
  R1Tensor tmp(point);
  tmp -= center;
  R1Tensor local;
  local.AijBj(rotation, tmp);
  PointAtSurface(radii, local, tmp);
  pointAtSurface.AijBi(rotation, tmp);
  pointAtSurface += center;
}

/**
 * @brief Get the point projected to the ellipsoid surface in the local frame
 * @author Scott Johnson
 *
 * @param[in] radii Principal radii
 * @param[in] point Point in local frame
 * @param[out] pointAtSurface Point on the potential surface in the local frame
 */
void EllipsoidalDiscreteElementManagerT::PointAtSurface( const R1Tensor& radii,
                                                         const R1Tensor& point,
                                                         R1Tensor& local)
{
  const realT alpha = 1.0/PotentialValue(radii, point);
  local = point;
  local *= alpha;
}


/**
 * @brief Get the normal to the point in the local frame on the scaled potential surface of the ellipsoid
 * @author Scott Johnson
 * \f[\frac{x}{a}^2 + \frac{y}{b}^2 + \frac{z}{c}^2 = 1\f]
 * \f[\hat{n} = \left[\begin{array}{ccc} \frac{x}{a^2} & \frac{y}{b^2} & \frac{z}{c^2} \end{array}\right] / | \left[\begin{array}{ccc} \frac{x}{a^2} & \frac{y}{b^2} & \frac{z}{c^2} \end{array}\right] | \f]
 * @param[in] local Point in space in the local frame
 * @param[in] radii Principal radii of the ellipsoid
 * @param[out] normal Normal to the potential surface in local frame
 */
void EllipsoidalDiscreteElementManagerT::NormalCartesian(const R1Tensor& radii,
                                                         const R1Tensor& local,
                                                         R1Tensor& normal)
 {
  realT mag = 0.;
  normal = local;
  for (unsigned int i = 0; i < nsdof; i++)
  {
    normal(i) /= (radii(i) * radii(i));
    mag += normal(i) * normal(i);
  }
  if (mag > 0)
  {
    mag = 1. / sqrt(mag);
    normal *= mag;
  }
}

/**
 * @brief Get the normal to the point in the glboal frame on the scaled potential surface of the ellipsoid
 * @author Scott Johnson
 * @param[in] radii Principal radii of the ellipsoid
 * @param[in] center Center of the ellipsoid
 * @param[in] rotation Rotation tensor
 * @param[in] point Point in space in the global frame
 * @param[out] normal Normal to the potential surface in local frame
 */
void EllipsoidalDiscreteElementManagerT::NormalCartesian(const R1Tensor& radii,
                                                         const R1Tensor& center,
                                                         const R2Tensor& rotation,
                                                         const R1Tensor& point,
                                                         R1Tensor& normal)
{
  R1Tensor tmp(point);
  tmp -= center;
  R1Tensor local;
  local.AijBj(rotation, tmp);
  NormalCartesian(radii, local, tmp);
  normal.AijBi(rotation, tmp);
}


void EllipsoidalDiscreteElementManagerT::PointAtNormalCartesian(const R1Tensor& radii,
                                                                const R1Tensor& center,
                                                                const R2Tensor& rotation,
                                                                const R1Tensor& normal,
                                                                R1Tensor& point)
{
  R1Tensor tmp;
  tmp.AijBj(rotation, normal);
  R1Tensor local;
  PointAtNormalCartesian(radii, tmp, local);
  point.AijBi(rotation, local);
  point += center;
}

void EllipsoidalDiscreteElementManagerT::PointAtNormalCartesian(const R1Tensor& radii,
                                                                const R1Tensor& normalLocal,
                                                                R1Tensor& local)
{
  local = normalLocal;
  //realT mag = 0.;
  for (unsigned int i = 0; i < nsdof; i++)
  {
    local(i) *= radii(i);
  }
}

/**
 * @brief Get the potential evaluated at the given point in the global frame
 * @author Scott Johnson
 *
 * @param[in] radii Principal radii of the ellipsoid
 * @param[in] center Center of the ellipsoid
 * @param[in] rotation Rotation tensor transform for the ellipsoid
 * @param[in] point Point in global frame at which to evaluate the potential
 * @return Value of the potential
 */
realT EllipsoidalDiscreteElementManagerT::PotentialValue(const R1Tensor& radii,
                                                         const R1Tensor& center,
                                                         const R2Tensor& rotation,
                                                         const R1Tensor& point)
{
  R1Tensor tmp(point);
  tmp -= center;
  R1Tensor local;
  local.AijBj(rotation, tmp);
  return PotentialValue(radii, local);
}

/**
 * @brief Get the potential evaluated at the given point in the local frame
 * @author Scott Johnson
 *
 * @param[in] radii Principal radii of the ellipsoid
 * @param[in] local Point in local frame at which to evaluate the potential
 * @return Value of the potential
 */
realT EllipsoidalDiscreteElementManagerT::PotentialValue(const R1Tensor& radii,
                                                         const R1Tensor& local)
{
  realT ret = 0;
  for(unsigned int i = 0; i < 3; i++)
  {
    realT ff = local(i) / radii(i);
    ff *= ff;
    ret += ff;
  }
  return ret;
}



/**
 * @brief Get the angular coordinates
 * @author Scott Johnson
 * @param[in] point Point in space in the local frame
 * @param[in] r Principal radii of the ellipsoid
 * @param[out] cu U-angle cosine
 * @param[out] cv V-angle cosine
 * @param[out] su U-angle sine
 * @param[out] sv V-angle sine
 */
void EllipsoidalDiscreteElementManagerT::AngularCoordinates(const R1Tensor& point, const R1Tensor& r, realT& cu, realT& su, realT& cv, realT& sv)
{
  R1Tensor spt = point;
  for(localIndex i = 0; i < nsdof; i++)
    spt(i) /= r(i);
  spt.Normalize();

  sv = spt(2);
  cv = cos(asin(sv));
  if(isZero( cv ))
  {
    su = 0;
    cu = 1;
  }
  else
  {
    su = spt(1) / cv;
    cu = spt(0) / cv;
  }
}

/**
 * @brief Return the Gaussian curvature of the surface
 * @author Scott Johnson
 * @param[in] r Principal radii
 * @param[in] cu U-angle cosine
 * @param[in] cv V-angle cosine
 * @param[in] su U-angle sine
 * @param[in] sv V-angle sine
 * @return Gaussian curvature of the surface at the given local angular coordinate
 */
realT EllipsoidalDiscreteElementManagerT::GaussianCurvature(const R1Tensor& r, const realT cu, const realT su, const realT cv, const realT sv)
{
  realT r2 = r.ProductOfSquares();
  realT fct0 = r(2) * cv;
  {
    if( !isZero( fct0 ) )
    {
      fct0 *= fct0;
      {
        realT fct1 = r(1) * cu;
        fct1 *= fct1;
        {
          realT fct2 = r(0) * su;
          fct2 *= fct2;
          fct1 += fct2;
        }
        fct0 *= fct1;
      }
    }
    {
      realT fct1b = r(0) * r(1) * sv;
      fct1b *= fct1b;
      fct0 += fct1b;
    }
  }
  fct0 *= fct0;
  r2 /= fct0;
  return r2;
}

/**
 * @brief Get the curvature at the given point in the local frame
 * @author Scott Johnson
 * @param[in] point Point in the local frame
 * @param[in] r Principal radii of the potential
 * @return Curvature at the point
 */
realT EllipsoidalDiscreteElementManagerT::GaussianCurvature(const R1Tensor& point, const R1Tensor& r)
{
  //check whether this is a sphere
  if(  isEqual(r(0),r(1))&&isEqual(r(1),r(2)) )
    return 1.0/r(0);
  realT cu, cv, su, sv;
  AngularCoordinates(point, r, cu, su, cv, sv);
  return GaussianCurvature(r, cu, su, cv, sv);
}

realT EllipsoidalDiscreteElementManagerT::Volume(const R1Tensor& radii)
{
  return (4.0/3.0) * 3.1415926535897932385 * radii(0) * radii(1) * radii(2);
}
