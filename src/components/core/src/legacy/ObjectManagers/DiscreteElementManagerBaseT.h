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
 * @file DiscreteElementManagerBaseT.h
 * @author Scott Johnson
 * @date created on July 12, 2011
 */

#ifndef DISCRETELEMENTMANAGERBASET_H_
#define DISCRETELEMENTMANAGERBASET_H_

#include "../../dataRepository/Group.hpp"
#include "../IO/ticpp/HierarchicalDataNode.h.old"
#include "Common/Common.h"
#include "ContactManagerBaseT.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "Constitutive/Material/MaterialFactory.h"

/**
 * @author Scott Johnson
 * @brief Class to manager the collection of basic DEM elements
 */
class DiscreteElementManagerBaseT : public ObjectDataStructureBaseT
{
public:
  DiscreteElementManagerBaseT();
  DiscreteElementManagerBaseT( const ObjectType objectType );
  virtual ~DiscreteElementManagerBaseT();

  void erase( const localIndex i );
  globalIndex resize( const localIndex size, const bool assignGlobals = false );

  void AddBaseFields();

  virtual void Initialize(){AddBaseFields();}

  virtual unsigned int Unpack( const char*& buffer, lArray1d& elementReceiveLocalIndices );

  virtual unsigned int Pack( const lArray1d& sendElements, bufvector& buffer ) const;

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject  = NULL) { (void)referenceObject; }
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) { (void)referenceObject; }
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("DiscreteElementManager::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  virtual void RecalculatePhysicalProperties() {}

  virtual void ReadXML(TICPP::HierarchicalDataNode*);

  virtual void WriteVTK(const int cycleNum, ContactManagerBaseT& contacts);

protected:
  virtual void WriteVTKPointData(std::ofstream& out);

  virtual globalIndex insert( const localIndex i, const bool assignGlobals = false );

public:

  bool writeVTK;

  //boundary contact states - need to be public to ReadXML

#if USECPP11==1
  std::unique_ptr<MaterialBase> m_mat;
#else
  MaterialBase* m_mat;
#endif



  static void ForceToCouple(const R1Tensor& position,
                            const R1Tensor& cforce,
                            const R1Tensor& currentPosition,
                            const R2Tensor& rotation,
                            R1Tensor& force,
                            R1Tensor& moment);

  /**
   * @brief Transform the given local point to the global frame
   * @author Scott Johnson
   * @param[in] a Index of the discrete element
   * @param[in] local Local position vector to transform
   * @param[out] global Global position vector result
   */
  inline void LocalToGlobal(const localIndex a, const R1Tensor& local, R1Tensor& global) const
  {
    LocalToGlobalDirection(a, local, global);
    global += this->GetFieldData<FieldInfo::currentPosition> () [a];
  }

  /**
   * @brief Transform the given global point to the local frame
   * @author Scott Johnson
   * @param[in] a Index of the discrete element
   * @param[in] global Global direction vector to transform
   * @param[out] local Local direction vector result
   */
  inline void GlobalToLocal(const localIndex a, const R1Tensor& global, R1Tensor& local) const
  {
    R1Tensor tmp = global;
    tmp -= this->GetFieldData<FieldInfo::currentPosition> () [a];
    GlobalToLocalDirection(a, tmp, local);
  }

  /**
   * @brief Transform the given local vector to the global frame
   * @author Scott Johnson
   * @param[in] a Index of the discrete element
   * @param[in] local Local position vector to transform
   * @param[out] global Global position vector result
   */
  inline void LocalToGlobalDirection(const localIndex a, const R1Tensor& local, R1Tensor& global) const
  {
    R2Tensor rotation;
    RotationTensor(a, rotation);
    global.AijBi(rotation, local);
  }

  /**
   * @brief Transform the given global vector to the local frame
   * @author Scott Johnson
   * @param[in] a Index of the discrete element
   * @param[in] global Global direction vector to transform
   * @param[out] local Local direction vector result
   */
  inline void GlobalToLocalDirection(const localIndex a, const R1Tensor& global, R1Tensor& local) const
  {
    R2Tensor rotation;
    RotationTensor(a, rotation);
    R1Tensor tmp = global;
    local.AijBj(rotation, tmp);
  }

  /**
   * @brief Convert quaternion to rotation tensor
   * @author Scott Johnson
   * Converts the quaternion to an equivalent rotation tensor
   * @param[in] w Magnitude component of quaternion
   * @param[in] xyz Directional component of quaternion
   * @param[out] rotation Rotation tensor
   */
  static inline void QuaternionToRotation(const realT& w, const R1Tensor& xyz, R2Tensor& rotation)
  {
    const int nsdofp1 = nsdof + 1;
    R2SymTensorT<nsdofp1> x2;
    {
      R1TensorT<nsdofp1> q;
      q[0] = xyz[0];
      q[1] = xyz[1];
      q[2] = xyz[2];
      q[3] = w;
      x2.dyadic_aa(q);
    }

    rotation(0,0) = 1.0 - 2.0 * (x2(1,1) + x2(2,2));
    rotation(1,0) = 2.0 * (x2(0,1) - x2(2,3));
    rotation(2,0) = 2.0 * (x2(0,2) + x2(1,3));
    rotation(0,1) = 2.0 * (x2(0,1) + x2(2,3));
    rotation(1,1) = 1.0 - 2.0 * (x2(0,0) + x2(2,2));
    rotation(2,1) = 2.0 * (x2(1,2) - x2(0,3));
    rotation(0,2) = 2.0 * (x2(0,2) - x2(1,3));
    rotation(1,2) = 2.0 * (x2(1,2) + x2(0,3));
    rotation(2,2) = 1.0 - 2.0 * (x2(0,0) + x2(1,1));
  }

  /**
   * @brief Convert rotation tensor to quaternion
   * @author Scott Johnson
   * Converts the rotation tensor to an equivalent quaternion
   * @param[in] rotation Rotation tensor
   * @param[out] w Magnitude component of quaternion
   * @param[out] xyz Directional component of quaternion
   */
  static inline void RotationToQuaternion(const R2Tensor& rotation, realT& w, R1Tensor& xyz)
  {
    realT trace = 1.0 + rotation.Trace();
    realT S = 0;
    if (trace > 1.0e-10)
    {
      S = sqrt(trace) * 2.0;
      realT tmp = 1.0/S;
      xyz[0] = (rotation(1,2) - rotation(2,1)) * tmp;
      xyz[1] = (rotation(2,0) - rotation(0,2)) * tmp;
      xyz[2] = (rotation(0,1) - rotation(1,0)) * tmp;
      w = 0.25 * S;
    }
    else
    {
      if (rotation(0,0) > rotation(1,1) && rotation(0,0) > rotation(2,2))
      {
        // Column 0:
        S = sqrt(1.0 + rotation(0,0) - rotation(1,1) - rotation(2,2)) * 2.0;
        realT tmp = 1.0/S;
        xyz[0] = 0.25 * S;
        xyz[1] = (rotation(0,1) + rotation(1,0)) * tmp;
        xyz[2] = (rotation(2,0) + rotation(0,2)) * tmp;
        w = (rotation(1,2) - rotation(2,1)) * tmp;
      }
      else if (rotation(1,1) > rotation(2,2))
      {
        // Column 1:
        S = sqrt(1.0 + rotation(1,1) - rotation(0,0) - rotation(2,2)) * 2.0;
        realT tmp = 1.0/S;
        xyz[0] = (rotation(0,1) + rotation(1,0)) * tmp;
        xyz[1] = 0.25 * S;
        xyz[2] = (rotation(1,2) + rotation(2,1)) * tmp;
        w = (rotation(2,0) - rotation(0,2)) * tmp;
      }
      else
      {
        // Column 2:
        S = sqrt(1.0 + rotation(2,2) - rotation(0,0) - rotation(1,1)) * 2.0;
        realT tmp = 1.0/S;
        xyz[0] = (rotation(2,0) + rotation(0,2)) * tmp;
        xyz[1] = (rotation(1,2) + rotation(2,1)) * tmp;
        xyz[2] = 0.25 * S;
        w = (rotation(0,1) - rotation(1,0)) * tmp;
      }
    }

    //normalize
    {
      realT tmp = w;
      tmp *= w;
      tmp += Dot(xyz,xyz);
      if(!isZero( tmp ) && !isEqual( tmp, 1.0 ) )
      {
        tmp = 1./sqrt(tmp);
        w *= tmp;
        xyz *= tmp;
      }
    }
  }

  /**
   * @brief Convert the given rotation tensor to the equivalent quaternion
   * @author Scott Johnson
   * @description see also
   * http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
   * @param[in] a Index of the discrete element
   * @param[in] rotation Rotation tensor
   */
  inline void SetRotation(const localIndex a, const R2Tensor& rotation)
  {
    //get array references
    R1Tensor& rotationAxis = GetFieldData<FieldInfo::rotationAxis> ()[a];//qx,qy,qz
    realT& rotationMagnitude = GetFieldData<FieldInfo::rotationMagnitude> ()[a];//qw
    RotationToQuaternion(rotation, rotationMagnitude, rotationAxis);
    //rotationAxis *= -1.0;
  }

  /**
   * @brief Get the rotation tensor
   * @author Scott Johnson
   * The global to local rotational transform tensor
   * @param[in] a  Discrete element index
   * @param[out] rotation Rotation tensor
   */
  inline void RotationTensor(const localIndex a, R2Tensor& rotation) const
  {
    //note: we are storing and operating on the quaternion in the following
    // manner
    //R1Tensor raxis = GetFieldData<FieldInfo::rotationAxis> ()[a];
    //raxis *= -1.0;
    QuaternionToRotation(
      GetFieldData<FieldInfo::rotationMagnitude> ()[a],
      GetFieldData<FieldInfo::rotationAxis> ()[a],
      rotation);
  }

  /**
   * @brief Change in quaternion with respect to time
   * @author Scott Johnson
   * @param[in] a Discrete element index
   * @param[in] qh Quaternion representation of the rotation
   * @param[out] dqdt Change in quaternion with respect to time
   */
  inline void Calculate_dqdt(const localIndex a, const realT* const qh, realT* const dqdt) const
  {
    const array<R1Tensor>& rotationalVelocity   = GetFieldData<FieldInfo::rotationalVelocity> ();

    //dq[0] = dot([-q1, -q2, -q3], [w0, w1, w2])
    dqdt[0] = -1. * (qh[1] * rotationalVelocity[a][0] + qh[2]
                     * rotationalVelocity[a][1] + qh[3] * rotationalVelocity[a][2]);
    //dq[1] = dot([q0, -q3, q2], [w0, w1, w2])
    dqdt[1] = (qh[0] * rotationalVelocity[a][0] - qh[3]
               * rotationalVelocity[a][1] + qh[2] * rotationalVelocity[a][2]);
    //dq[2] = dot([q3, q0, -q1], [w0, w1, w2])
    dqdt[2] = ( qh[3] * rotationalVelocity[a][0]
                + qh[0] * rotationalVelocity[a][1]
                - qh[1] * rotationalVelocity[a][2] );
    //dq[3] = dot([-q2, q1, q0], [w0, w1, w2])
    dqdt[3] = (-qh[2] * rotationalVelocity[a][0] + qh[1]
               * rotationalVelocity[a][1] + qh[0] * rotationalVelocity[a][2]);
  }
private:

};

#endif /* DISCRETELEMENTMANAGERBASET_H_ */
