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
 * @file EllipsoidalDiscreteElementManagerT.h
 * @author Scott Johnson
 * @date created on July 14, 2011
 */

#ifndef ELLIPSOIDALDISCRETELEMENTMANAGERT_H_
#define ELLIPSOIDALDISCRETELEMENTMANAGERT_H_

#include "DiscreteElementManagerBaseT.h"
#include "EllipsoidalContactManagerT.h"

#include "DataStructures/VectorFields/StableTimeStep.h"

#include "Contact/SpatialSorterFactory.h"

#include "Constitutive/Interface/InterfaceFactory.h"
#include "Constitutive/Material/MaterialFactory.h"

/**
 * @author Scott Johnson
 * @brief Struct to manage the cylindrical boundary
 */
struct DiscreteElementCylindricalBoundary
{
  //boundary geometry
  realT radiusCylinder, heightCylinder, heightCylinderNow;
  R1Tensor centerCylinder;
  realT strainRateCylinder, stressRateCylinder, stopRateTime;
#if USECPP11==1
  std::unique_ptr<MaterialBase> m_mat;
#else
  MaterialBase* m_mat;
#endif

  bool Active() const
  {
    return heightCylinder > 0 && radiusCylinder > 0;
  }

  void Reset()
  {
    radiusCylinder = -1.0;
    heightCylinder = -1.0;
    heightCylinderNow = -1.0;
    centerCylinder = 0.0;
    strainRateCylinder = 0.0;
    stressRateCylinder = 0.0;
    stopRateTime=std::numeric_limits<realT>::max();

#if USECPP11!=1
    if(m_mat)
      delete m_mat;
#endif
    m_mat = MaterialFactory::NewMaterial("LinearElasticDEMMaterial");
    m_mat->resize(0,1);
  }

  void Update(const realT time)
  {
    if (time > stopRateTime)
      return;
    heightCylinderNow = 0.5 * heightCylinder * (1 - (strainRateCylinder * time));
  }

  bool InContact(const realT radius,
                 const R1Tensor& center,
                 realT& overlap,
                 R1Tensor& normal1,
                 R1Tensor& contactPoint)
  {
    R1Tensor dist(center);
    dist -= centerCylinder;
    const bool axialContact = fabs(dist(1)) > (heightCylinderNow - radius);

    const realT rsum = radiusCylinder - radius;
    const realT ddRadial = Dot(dist, dist) - dist(1)*dist(1);
    const bool radialContact = rsum*rsum < ddRadial;

    if((!radialContact) && (!axialContact))
      return false;

    realT overlapRadial = 0;
    R1Tensor normalRadial(dist);
    R1Tensor contactPointRadial(0.0);
    if(radialContact)
    {
      normalRadial(1) = 0.;
      normalRadial.Normalize();
      normalRadial *= -1.0;

      overlapRadial = sqrt(ddRadial) - rsum;

      contactPointRadial = normalRadial;
      contactPointRadial *= -radiusCylinder;
    }

    realT overlapAxial = 0;
    R1Tensor normalAxial(0.0);
    R1Tensor contactPointAxial(0.0);
    if(axialContact)
    {
      normalAxial(1) = 1.0;
      normalAxial *= dist(1) > 0 ? -1.0 : 1.0;

      overlapAxial = radius + fabs(dist(1)) - heightCylinderNow;

      contactPointAxial = normalAxial;
      contactPointAxial *= -heightCylinderNow;
    }

    //here we will just do a simple procedure to homogenize the contact with a concave body
    normalRadial *= overlapRadial;
    normalAxial *= overlapAxial;
    normal1 = normalRadial;
    normal1 += normalAxial;
    normal1.Normalize();

    overlap = overlapAxial + overlapRadial;

    contactPointRadial *= overlapRadial;
    contactPointAxial *= overlapAxial;
    contactPoint = contactPointRadial;
    contactPoint += contactPointAxial;
    contactPoint *= 0.5;
    contactPoint += centerCylinder;
    return true;
  }

  void ReadXML(TICPP::HierarchicalDataNode* EDEMNode)
  {
    radiusCylinder = EDEMNode->GetAttributeOrDefault<realT> ("radius_cylinder", -1.0);
    heightCylinder = EDEMNode->GetAttributeOrDefault<realT> ("height_cylinder", -1.0);
    strainRateCylinder = EDEMNode->GetAttributeOrDefault<realT> ("strain_rate_cylinder", 0.0);
    stressRateCylinder = EDEMNode->GetAttributeOrDefault<realT> ("stress_rate_cylinder", 0.0);
    stopRateTime = EDEMNode->GetAttributeOrDefault<realT> ("stop_rate_at_time",1.01e20);

    if(!Active())
      return;

    TICPP::HierarchicalDataNode* EDemNodeMaterial0 = EDEMNode->GetChild("BoundaryMaterial");
    if(EDemNodeMaterial0 == NULL)
      throw GPException("Must define a boundary material model for Elliposidal DEM if boundary active");
    TICPP::HierarchicalDataNode* EDemNodeMaterial1 = EDemNodeMaterial0->Next(true);
    if(EDemNodeMaterial1 == NULL)
      throw GPException("Must define a boundary material model for Elliposidal DEM (1) if boundary active");

#if USECPP11!=1
    if(m_mat)
      delete m_mat;
#endif
    m_mat = MaterialFactory::NewMaterial(EDemNodeMaterial1->Heading(), EDemNodeMaterial1);
    m_mat->resize(0,1);
  }
};

/**
 * @author Scott Johnson
 * @brief Class to manager the collection of ellipsoidal DEM elements
 */
class EllipsoidalDiscreteElementManagerT: public DiscreteElementManagerBaseT
{
public:
  EllipsoidalDiscreteElementManagerT();

  virtual ~EllipsoidalDiscreteElementManagerT();

  void erase( const localIndex i );
  globalIndex resize( const localIndex size, const bool assignGlobals = false );

  void RecalculatePhysicalProperties();

  virtual void WriteSilo( SiloFile& siloFile,
                          const std::string& siloDirName,
                          const std::string& meshname,
                          const int centering,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart,
                          const std::string& regionName = "none",
                          const lArray1d& mask = lArray1d());

  virtual void ReadSilo( const SiloFile& siloFile,
                         const std::string& siloDirName,
                         const std::string& meshname,
                         const int centering,
                         const int cycleNum,
                         const realT problemTime,
                         const bool isRestart,
                         const std::string& regionName = "none",
                         const lArray1d& mask  = lArray1d());

  void ReadXML(TICPP::HierarchicalDataNode* EDEMNode);

  bool RecalculateNeighborList(const realT dt);

  void UpdateAndApplyContactStresses( StableTimeStep& maxdt, const realT dt,
                                      EllipsoidalContactManagerT& ecm);

  void UpdateCylindricalBoundary(StableTimeStep& maxdt, const realT dt, const realT time);

  void UpdateNodalStatesZeroForcesAndAccelerations();

  static void MomentOfInertia(const R1Tensor& radii, const realT mass, R1Tensor& moi);

  static void AngularCoordinates(const R1Tensor& point, const R1Tensor& r, realT& cu, realT& su, realT& cv, realT& sv);

  static realT GaussianCurvature(const R1Tensor& r, const realT cu, const realT su, const realT cv, const realT sv);

  static realT GaussianCurvature(const R1Tensor& point, const R1Tensor& r);

  static realT Volume(const R1Tensor& radii);

protected:
  void WriteVTKPointData(std::ofstream& out);

  virtual globalIndex insert( const localIndex i, const bool assignGlobals = false );

  static realT StableTimestep(const realT mass1, const realT mass2, const realT stiffness)
  {
    const realT mass = mass2 > mass1 ? mass1 : mass2;
    const realT dttmp = isZero(stiffness, 10*std::numeric_limits<realT>::min()) ?
        std::numeric_limits<realT>::max() :
        0.3 * sqrt(mass/stiffness);
    return dttmp;
  }

public:
  //search attributes
  OrderedVariableOneToManyRelation& m_neighborList;
  UnorderedVariableOneToManyRelation& m_neighborListInverse;

  bool sorted;
  SpatialSorting::SpatialSorterBase* m_sorter;

  realT m_contactBufferOffset;
  realT m_contactBufferFactor;

  //boundary contact states - need to be public to ReadXML
  DiscreteElementCylindricalBoundary m_boundary;
  InterfaceBase* m_boundaryContact;

private:

  bool InContact(const R1Tensor& radii,
                 const R1Tensor& center,
                 const R2Tensor& rotation,
                 realT& overlap,
                 R1Tensor& normal1,
                 R1Tensor& contactPoint);

  static bool InContact(const R1Tensor& radii0,
                         const R1Tensor& center0,
                         const R2Tensor& rotation0,
                         const R1Tensor& radii1,
                         const R1Tensor& center1,
                         const R2Tensor& rotation1,
                         realT& overlap,
                         R1Tensor& normal1,
                         R1Tensor& contactPoint);

  static realT Iterate(const R1Tensor& radii0,
                       const R1Tensor& center0,
                       const R2Tensor& rotation0,
                       const R1Tensor& radii1,
                       const R1Tensor& center1,
                       const R2Tensor& rotation1,
                       R1Tensor& contactPoint0,
                       R1Tensor& contactPoint1);

  static void PointAtSurface(const R1Tensor& radii,
                             const R1Tensor& center,
                             const R2Tensor& rotation,
                             const R1Tensor& point,
                             R1Tensor& pointAtSurface);

  static void PointAtSurface(const R1Tensor& radii,
                             const R1Tensor& point,
                             R1Tensor& local);

  static realT PotentialValue(const R1Tensor& radii,
                              const R1Tensor& center,
                              const R2Tensor& rotation,
                              const R1Tensor& point);

  static realT PotentialValue(const R1Tensor& radii,
                              const R1Tensor& local);

  static void NormalCartesian(const R1Tensor& radii,
                              const R1Tensor& center,
                              const R2Tensor& rotation,
                              const R1Tensor& point,
                              R1Tensor& normal);

  static void NormalCartesian(const R1Tensor& radii,
                              const R1Tensor& local,
                              R1Tensor& normal);

  static void PointAtNormalCartesian(const R1Tensor& radii,
                              const R1Tensor& center,
                              const R2Tensor& rotation,
                              const R1Tensor& normal,
                              R1Tensor& point);

  static void PointAtNormalCartesian(const R1Tensor& radii,
                                     const R1Tensor& normalLocal,
                                     R1Tensor& local);

};

#endif /* ELLIPSOIDALDISCRETELEMENTMANAGERT_H_ */
