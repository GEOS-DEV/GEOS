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
 * @file DiscreteElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 20, 2010
 */

#ifndef DISCRETELEMENTMANAGERT_H_
#define DISCRETELEMENTMANAGERT_H_

#include "Common/Common.h"
#include "DiscreteElementManagerBaseT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "FaceManagerT.h"
#include "ContactManagerT.h"
#include "IO/silo/SiloFile.h"

struct DiscreteElementContact
{
  realT m_area;
  realT m_r0, m_r1;
  realT m_overlap;
  R1Tensor m_contactPoint;
  R1Tensor m_normal;

  void Reset()
  {
    m_area = m_r0 = m_r1 = m_overlap = 0.0;
    m_contactPoint = 0.0;
  }


  void AddContribution(const realT contactArea, const realT dn, const realT r0, const realT r1,
                       const R1Tensor& cp, const R1Tensor& normal)
  {
    //weight the contact point and normal
    {
      R1Tensor tmp(cp);
      tmp *= contactArea;
      this->m_contactPoint += tmp;

      tmp = normal;
      tmp *= contactArea;
      this->m_normal += tmp;
    }

    //weight the area
    this->m_area += contactArea;
    this->m_r0 += contactArea * r0;
    this->m_r1 += contactArea * r1;
    this->m_overlap += contactArea * dn;
  }

  bool Values(realT& dn, realT& reff, R1Tensor& cp, R1Tensor& normal) const
  {
    if( isZero( this->m_area, 0.0 ) )
      return false;

    const realT ia = 1.0 / this->m_area;

    //normal approach
    dn = this->m_overlap * ia;

    //contact point
    cp = this->m_contactPoint;
    cp *= ia;

    //effective radius
    reff = this->m_area * (1.0/this->m_r0 + 1.0/this->m_r1);

    //normal
    normal = this->m_normal;
    normal.Normalize();

    return dn > 0;
  }

  void ReadSilo(const SiloFile& siloFile)
  {
    siloFile.DBReadWrapper("m_area",m_area);
    siloFile.DBReadWrapper("m_r0",m_r0);
    siloFile.DBReadWrapper("m_r1",m_r1);
    siloFile.DBReadWrapper("m_overlap",m_overlap);
    siloFile.DBReadWrapper("m_contactPoint_x",m_contactPoint[0]);
    siloFile.DBReadWrapper("m_contactPoint_y",m_contactPoint[1]);
    siloFile.DBReadWrapper("m_contactPoint_z",m_contactPoint[2]);
    siloFile.DBReadWrapper("m_normal_x",m_normal[0]);
    siloFile.DBReadWrapper("m_normal_y",m_normal[1]);
    siloFile.DBReadWrapper("m_normal_z",m_normal[2]);
  }

  void WriteSilo(SiloFile& siloFile) const
  {
    siloFile.DBWriteWrapper("m_area",m_area);
    siloFile.DBWriteWrapper("m_r0",m_r0);
    siloFile.DBWriteWrapper("m_r1",m_r1);
    siloFile.DBWriteWrapper("m_overlap",m_overlap);
    siloFile.DBWriteWrapper("m_contactPoint_x",m_contactPoint[0]);
    siloFile.DBWriteWrapper("m_contactPoint_y",m_contactPoint[1]);
    siloFile.DBWriteWrapper("m_contactPoint_z",m_contactPoint[2]);
    siloFile.DBWriteWrapper("m_normal_x",m_normal[0]);
    siloFile.DBWriteWrapper("m_normal_y",m_normal[1]);
    siloFile.DBWriteWrapper("m_normal_z",m_normal[2]);
  }

  void Serialize(std::map<std::string, array<real64> >& realFields,
                 std::map<std::string, array<R1Tensor> >& R1Fields) const
  {
    realFields["CDEM_Area"].push_back(m_area);
    realFields["CDEM_R1"].push_back(m_r0);
    realFields["CDEM_R2"].push_back(m_r1);
    realFields["CDEM_Overlap"].push_back(m_overlap);

    R1Fields["CDEM_ContactPoint"].push_back(m_contactPoint);
    R1Fields["CDEM_Normal"].push_back(m_normal);
  }

  static void WriteSilo(SiloFile& siloFile,
                        std::map<std::string, array<real64> >& realFields,
                        std::map<std::string, array<R1Tensor> >& R1Fields,
                        const std::string& siloDirName,
                        const std::string& meshname,
                        const int centering,
                        const int cycleNum,
                        const realT problemTime,
                        const bool )
  {
    std::string subDirectory = siloDirName;
    std::string rootDirectory = "/" + siloDirName;
    siloFile.MakeSubDirectory( subDirectory, rootDirectory );
    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

    for(std::map<std::string, array<real64> >::iterator it = realFields.begin() ; it != realFields.end() ; ++it)
      siloFile.WriteDataField<realT>(meshname.c_str(), it->first, it->second, centering, cycleNum, problemTime, rootDirectory, "none" );

    for(std::map<std::string, array<R1Tensor> >::iterator it = R1Fields.begin() ; it != R1Fields.end() ; ++it)
      siloFile.WriteDataField<realT>(meshname.c_str(), it->first, it->second, centering, cycleNum, problemTime, rootDirectory, "none" );

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
};

/**
 * @author Scott Johnson
 * @brief Class to manager the collection of DEM elements
 */
class DiscreteElementManagerT : public DiscreteElementManagerBaseT
{
public:
  DiscreteElementManagerT();
  DiscreteElementManagerT(NodeManager* nm, FaceManagerT* fm);
  virtual ~DiscreteElementManagerT();

  unsigned int Unpack( const char*& buffer, lArray1d& elementReceiveLocalIndices );

  unsigned int Pack( const lArray1d& sendElements, bufvector& buffer ) const;

  void UpdateNodalStatesZeroForcesAndAccelerations();

  void ApplyNodalForce(const localIndex nodeIndex, const R1Tensor& f);

  void ApplyFaceForce( const localIndex faceIndex,
                       const R1Tensor& x,
                       const R1Tensor& f);

  void ApplyDiscreteElementContactForces(const realT youngs, const realT poissons);

  void CalculateEnergy();

  void ApplyNodalForces();

protected:
  void WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                         const std::string& siloDirName,
                                         const std::string& meshname,
                                         const int centering,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const std::string& multiRoot,
                                         const std::string& regionName = "none",
                                         const lArray1d& mask = lArray1d());

  void ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                          const std::string& siloDirName,
                                          const std::string& meshname,
                                          const int centering,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const std::string& regionName = "none",
                                          const lArray1d& mask = lArray1d());
public:
  NodeManager* m_nodeManager;
  FaceManagerT* m_faceManager;

  ///maps the discrete element to the indices in the _external_ face manager
  OrderedVariableOneToManyRelation& m_discreteElementToExternalFacesMap;

  ///maps the discrete element to the indices in the _external_ node manager
  OrderedVariableOneToManyRelation& m_discreteElementToExternalNodesMap;

  //the following needs to be written and read for restart
  std::map<localIndex, std::map<localIndex, DiscreteElementContact> > m_discreteElementToDiscreteElementContactsMap;



  /**
   * @brief Global nodal velocity
   * @author Scott Johnson
   * @param[in] idem Discrete element index
   * @param[in] a Local nodal index in the discrete element's list
   * @param[out] v Velocity of the node in the global frame
   */
  inline void VelocityAtPoint(const localIndex idem, const R1Tensor& point, R1Tensor& v) const
  {
    v = 0.0;

    //get the local frame velocity
    R1Tensor vl;
    {
      const array<R1Tensor>& rotationalVelocity = this->GetFieldData< FieldInfo::rotationalVelocity> ();
      R1Tensor lpt;
      this->GlobalToLocal(idem, point, lpt);
      vl.Cross(lpt, rotationalVelocity[idem]);
    }

    //get the global frame velocity
    this->LocalToGlobalDirection(idem, vl, v);

    //add in the translational component
    v += this->GetFieldData<FieldInfo::velocity>()[idem];
  }

  /**
   * @brief Global nodal velocity
   * @author Scott Johnson
   * @param[in] idem Discrete element index
   * @param[in] a Local nodal index in the discrete element's list
   * @param[out] v Velocity of the node in the global frame
   */
  inline void NodalVelocity(const localIndex idem, const localIndex a, R1Tensor& v) const
  {
    v = 0.0;

    //get the local frame velocity
    R1Tensor vl;
    {
      const localIndex inode = m_discreteElementToExternalNodesMap[idem][a];
      const array<R1Tensor>& rotationalVelocity = this->GetFieldData< FieldInfo::rotationalVelocity> ();
      const array<R1Tensor>& relativePosition = m_nodeManager->GetFieldData< FieldInfo::relativePosition> ();
      vl.Cross(relativePosition[inode], rotationalVelocity[idem]);
    }

    //get the global frame velocity
    this->LocalToGlobalDirection(idem, vl, v);

    //add in the translational component
    v += this->GetFieldData<FieldInfo::velocity>()[idem];
  }

  /**
   * @brief Global nodal velocity
   * @author Scott Johnson
   * @param[in] i Discrete element index
   * @param[in] j Local nodal index in the discrete element's list
   * @param[in] r Rotation matrix for discrete element i
   * @param[out] v Velocity of the node in the global frame
   */
  inline void NodalVelocity(const int i, const int j, const R2Tensor r, R1Tensor& v) const
  {
    v = 0.0;

    //get the local frame velocity
    R1Tensor vl;
    {
      int inode = m_discreteElementToExternalNodesMap[i][j];
      const array<R1Tensor>& rotationalVelocity = this->GetFieldData< FieldInfo::rotationalVelocity> ();
      const array<R1Tensor>& relativePosition = m_nodeManager->GetFieldData< FieldInfo::relativePosition> ();
      vl.Cross(relativePosition[inode], rotationalVelocity[i]);
    }

    //get the global frame velocity
    v.AijBi(r, vl);

    //add in the translational component
    v += this->GetFieldData<FieldInfo::velocity>()[i];
  }

  void RecalculatePhysicalProperties();
private:

  /**
   * @brief Get the geometric properties of the tetrahedron described by the
   * four points
   * @author Scott Johnson
   * @param[in] x0 First point in global frame (about which rotational inertia
   * is taken)
   * @param[in] x1 Second point in global frame
   * @param[in] x2 Third point in global frame
   * @param[in] x3 Fourth point in global frame
   * @param[out] volume Volume of the tetrahedron
   * @param[out] momentOfInertia Moment of inertia of the tetrahedron about x0
   * @param[out] centroid Centroid of the tetrahedron
   * @return return
   */
  void TetrahedronProperties( const R1Tensor& x0,
                              const R1Tensor& x1,
                              const R1Tensor& x2,
                              const R1Tensor& x3,
                              realT& volume,
                              R2Tensor& momentOfInertia,
                              R1Tensor& centroid)
  {
    R1Tensor v1 = x1;
    v1 -= x0;
    R1Tensor v2 = x2;
    v1 -= x0;
    R1Tensor v12;
    v12.Cross(v2, v1);
    R1Tensor v3 = x3;
    v1 -= x0;
    //volume is 1/6 of the abs(triple product)
    volume = (1./6.)*fabs(Dot(v12, v3));
    //centroid is x0 + 0.25*sum(x1..3)
    centroid = x1;
    centroid += x2;
    centroid += x3;
    centroid *= 0.25;
    centroid += x0;
    //moment of inertia is
    //I = [1]*tr(C) - C
  }

  /**
   * @brief Helper function for
   * @author Scott Johnson
   * From David Eberly's "Geometric Tools"
   * see
   * http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
   * @param[in] w0
   * @param[in] w1
   * @param[in] w2
   * @param[out] f1
   * @param[out] f2
   * @param[out] f3
   * @param[out] g0
   * @param[out] g1
   * @param[out] g2
   */
  static void Subexpressions(const realT w0, const realT w1, const realT w2,
                             realT& f1,realT& f2,realT& f3,
                             realT& g0, realT& g1, realT& g2)
  {
    realT temp0 = w0 + w1;
    f1 = temp0 + w2;
    realT temp1 = w0 * w0;
    realT temp2 = temp1 + w1 * temp0;
    f2 = temp2 + w2 * f1;
    f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
    g0 = f2 + w0 * (f1 + w0);
    g1 = f2 + w1 * (f1 + w1);
    g2 = f2 + w2 * (f1 + w2);
  }

  /**
   * @brief Compute the moment of inertia of the polyhedron
   * @author Scott Johnson
   * From David Eberly's "Geometric Tools"
   * see
   * http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
   * @param[in] p Points on the polyhedron
   * @param[in] tmax Number of triangular facets
   * @param[in] index Indices of the nodes of the triangles in an ordered array
   * of 3-element tuples
   * @param[out] mass Mass of the polyhedron
   * @param[out] cm Centroid of the polyhedron
   * @param[out] inertia Moment of inertia of the polyhedron
   */
  static void CalculatePhysicalProperties(
    const array<R1Tensor>& p,
    const unsigned int tmax,
    const lArray1d& index,
    realT& mass,
    R1Tensor& cm,
    R2SymTensor& inertia)
  {
    const realT mult[10] = { 1./6,1./24,1./24,1./24,1./60,1./60,1./60,1./120,1./120,1./120};
    realT intg[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx

    for (unsigned int t = 0 ; t < tmax ; ++t)
    {
      //get the coordinates of the vertices
      R1Tensor x[3];
      for(unsigned int ii = 0 ; ii < 3 ; ii++)
        x[ii] = p[index[3*t + ii]];

      //get the cross product of the triangular face vectors
      R1Tensor d;
      {
        R1Tensor a1 = x[1];
        a1 -= x[0];
        R1Tensor a2 = x[2];
        a2 -= x[0];
        d.Cross(a1, a2);
      }

      // compute integral terms
      R1Tensor f[3];
      R1Tensor g[3];
      for(unsigned int ii = 0 ; ii < 3 ; ii++)
        Subexpressions(x[0](ii), x[1](ii), x[2](ii), f[0](ii), f[1](ii), f[2](ii), g[0](ii), g[1](ii), g[2](ii));

      // update integrals
      intg[0] += d(0) * f[0](0);

      intg[1] += d(0) * f[1](0);
      intg[2] += d(1) * f[1](1);
      intg[3] += d(2) * f[1](2);

      intg[4] += d(0) * f[2](0);
      intg[5] += d(1) * f[2](1);
      intg[6] += d(2) * f[2](2);

      intg[7] += d(0) * (x[0](1) * g[0](0) + x[1](1) * g[1](0) + x[2](1) * g[2](0));
      intg[8] += d(1) * (x[0](2) * g[0](1) + x[1](2) * g[1](1) + x[2](2) * g[2](1));
      intg[9] += d(2) * (x[0](0) * g[0](2) + x[1](0) * g[1](2) + x[2](0) * g[2](2));
    }
    for (unsigned int i = 0 ; i < 10 ; i++)
      intg[i] *= mult[i];
    mass = intg[0];

    // center of mass
    for (unsigned int i = 0 ; i < 3 ; i++)
      cm(i) = intg[i+1] / mass;

    // inertia tensor relative to center of mass
    inertia(0,0) = intg[5] + intg[6] - mass * (cm(1) * cm(1) + cm(2) * cm(2));
    inertia(1,1) = intg[4] + intg[6] - mass * (cm(2) * cm(2) + cm(0) * cm(0));
    inertia(2,2) = intg[4] + intg[5] - mass * (cm(0) * cm(0) + cm(1) * cm(1));
    inertia(0,1) = -(intg[7] - mass * cm(0) * cm(1));
    inertia(1,2) = -(intg[8] - mass * cm(1) * cm(2));
    inertia(0,2) = -(intg[9] - mass * cm(2) * cm(0));
  }
};



/**
 * @brief Compute various integrations over projection of face
 * @author Scott Johnson
 *
 * Adapted from volInt.c at
 * http://www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar
 * Disclaimer with volInt.c:
 *
 *  This code computes volume integrals needed for
 *  determining mass properties of polyhedral bodies.
 *
 *  For more information, see the accompanying READM
 *  file, and the paper
 *
 *  Brian Mirtich, "Fast and Accurate Computation of
 *  Polyhedral Mass Properties," journal of graphics
 *  tools, volume 1, number 1, 1996.
 *
 *  This source code is public domain, and may be used
 *  in any way, shape or form, free of charge.
 *
 *  Copyright 1995 by Brian Mirtich
 *  mirtich@cs.berkeley.edu
 *  http://www.cs.berkeley.edu/~mirtich
 *
 * @param[in] p Points on the polyhedron
 * @param[in] tmax Number of triangular facets
 * @param[in] index Indices of the nodes of the triangles in an ordered array of
 * 3-element tuples
 * @param[in] faceIndex Index of the triangular facet
 * @param[in] A Alpha
 * @param[in] B Beta
 * @param[in] C Gamma
 * @param[out] P1 Projection integral
 * @param[out] Pa Projection integral
 * @param[out] Pb Projection integral
 * @param[out] Paa Projection integral
 * @param[out] Pab Projection integral
 * @param[out] Pbb Projection integral
 * @param[out] Paaa Projection integral
 * @param[out] Paab Projection integral
 * @param[out] Pabb Projection integral
 * @param[out] Pbbb Projection integral
 */
/*
   static void ComputeProjectionIntegrals(const array<R1Tensor>& p,
                                    const unsigned int tmax,
                                    const lArray1d& index,
                                    const unsigned int faceIndex,
                                    const int A,
                                    const int B,
                                    const int C,
                                    realT& P1,
                                    realT& Pa,
                                    realT& Pb,
                                    realT& Paa,
                                    realT& Pab,
                                    realT& Pbb,
                                    realT& Paaa,
                                    realT& Paab,
                                    realT& Pabb,
                                    realT& Pbbb)
   {
   realT a0, a1, da;
   realT b0, b1, db;
   realT a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
   realT a1_2, a1_3, b1_2, b1_3;
   realT C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
   realT Cab, Kab, Caab, Kaab, Cabb, Kabb;

   P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

   //note: here we assume TRIANGULAR ELEMENTS
   for (int i = 0; i < 3; i++)
   {
    a0 = p[index[3*faceIndex+i]](A);
    b0 = p[index[3*faceIndex+i]](B);

    a1 = p[index[3*faceIndex+((i + 1) % 3)]](A);
    b1 = p[index[3*faceIndex+((i + 1) % 3)]](B);

    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0;
    a0_3 = a0_2 * a0;
    a0_4 = a0_3 * a0;
    b0_2 = b0 * b0;
    b0_3 = b0_2 * b0;
    b0_4 = b0_3 * b0;
    a1_2 = a1 * a1;
    a1_3 = a1_2 * a1;
    b1_2 = b1 * b1;
    b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1 * C1 + a0_2;
    Caa = a1 * Ca + a0_3;
    Caaa = a1 * Caa + a0_4;
    Cb = b1 * (b1 + b0) + b0_2;
    Cbb = b1 * Cb + b0_3;
    Cbbb = b1 * Cbb + b0_4;
    Cab = 3 * a1_2 + 2 * a1 * a0 + a0_2;
    Kab = a1_2 + 2 * a1 * a0 + 3 * a0_2;
    Caab = a0 * Cab + 4 * a1_3;
    Kaab = a1 * Kab + 4 * a0_3;
    Cabb = 4 * b1_3 + 3 * b1_2 * b0 + 2 * b1 * b0_2 + b0_3;
    Kabb = b1_3 + 2 * b1_2 * b0 + 3 * b1 * b0_2 + 4 * b0_3;

    P1 += db * C1;
    Pa += db * Ca;
    Paa += db * Caa;
    Paaa += db * Caaa;
    Pb += da * Cb;
    Pbb += da * Cbb;
    Pbbb += da * Cbbb;
    Pab += db * (b1 * Cab + b0 * Kab);
    Paab += db * (b1 * Caab + b0 * Kaab);
    Pabb += da * (a1 * Cabb + a0 * Kabb);
   }

   P1 /= 2.0;
   Pa /= 6.0;
   Paa /= 12.0;
   Paaa /= 20.0;
   Pb /= -6.0;
   Pbb /= -12.0;
   Pbbb /= -20.0;
   Pab /= 24.0;
   Paab /= 60.0;
   Pabb /= -60.0;
   };
 */
/**
 * @brief Compute various integrations over projection of face
 * @author Scott Johnson
 *
 * Adapted from volInt.c at
 * http://www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar
 * Disclaimer with volInt.c:
 *
 *  This code computes volume integrals needed for
 *  determining mass properties of polyhedral bodies.
 *
 *  For more information, see the accompanying READM
 *  file, and the paper
 *
 *  Brian Mirtich, "Fast and Accurate Computation of
 *  Polyhedral Mass Properties," journal of graphics
 *  tools, volume 1, number 1, 1996.
 *
 *  This source code is public domain, and may be used
 *  in any way, shape or form, free of charge.
 *
 *  Copyright 1995 by Brian Mirtich
 *  mirtich@cs.berkeley.edu
 *  http://www.cs.berkeley.edu/~mirtich
 *
 * @param[in] p Points on the polyhedron
 * @param[in] tmax Number of triangular facets
 * @param[in] index Indices of the nodes of the triangles in an ordered array of
 * 3-element tuples
 * @param[in] faceIndex Index of the triangular facet
 * @param[in] A Alpha
 * @param[in] B Beta
 * @param[in] C Gamma
 * @param[in] w Component of v0 in the direction opposite the normal
 * @param[in] n Normal to the face
 * @param[out] Fa Face integral
 * @param[out] Fb Face integral
 * @param[out] Fc Face integral
 * @param[out] Faa Face integral
 * @param[out] Fbb Face integral
 * @param[out] Fcc Face integral
 * @param[out] Faaa Face integral
 * @param[out] Fbbb Face integral
 * @param[out] Fccc Face integral
 * @param[out] Faab Face integral
 * @param[out] Fbbc Face integral
 * @param[out] Fcca Face integral
 */
/*
   static void ComputeFaceIntegrals(const array<R1Tensor>& p,
                          const unsigned int tmax,
                          const lArray1d& index,
                          const unsigned int faceIndex,
                          const realT w,
                          const R1Tensor& n,
                          const int A,
                          const int B,
                          const int C,
                          realT& Fa,
                          realT& Fb,
                          realT& Fc,
                          realT& Faa,
                          realT& Fbb,
                          realT& Fcc,
                          realT& Faaa,
                          realT& Fbbb,
                          realT& Fccc,
                          realT& Faab,
                          realT& Fbbc,
                          realT& Fcca)
   {
   realT k1, k2, k3, k4;
   realT P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

   ComputeProjectionIntegrals(p, tmax, index, faceIndex,
                             A, B, C,
                             P1, Pa, Pb,
                             Paa, Pab, Pbb,
                             Paaa, Paab, Pabb, Pbbb);

   k1 = 1 / n(C);
   k2 = k1 * k1;
   k3 = k2 * k1;
   k4 = k3 * k1;

   Fa = k1 * Pa;
   Fb = k1 * Pb;
   Fc = -k2 * (n(A)*Pa + n(B)*Pb + w*P1);

   Faa = k1 * Paa;
   Fbb = k1 * Pbb;
   Fcc = k3 * ((n(A)*n(A))*Paa + 2*n(A)*n(B)*Pab + (n(B)*n(B))*Pbb
 + w*(2*(n(A)*Pa + n(B)*Pb) + w*P1));

   Faaa = k1 * Paaa;
   Fbbb = k1 * Pbbb;
   Fccc = -k4 * ((n(A)*n(A)*n(A))*Paaa + 3*(n(A)*n(A))*n(B)*Paab
 + 3*n(A)*(n(B)*n(B))*Pabb + (n(B)*n(B)*n(B))*Pbbb
 + 3*w*((n(A)*n(A))*Paa + 2*n(A)*n(B)*Pab + (n(B)*n(B))*Pbb)
 + w*w*(3*(n(A)*Pa + n(B)*Pb) + w*P1));

   Faab = k1 * Paab;
   Fbbc = -k2 * (n(A)*Pabb + n(B)*Pbbb + w*Pbb);
   Fcca = k3 * ((n(A)*n(A))*Paaa + 2*n(A)*n(B)*Paab + (n(B)*n(B))*Pabb
 + w*(2*(n(A)*Paa + n(B)*Pab) + w*Pa));
   };
 */
/**
 * @brief Compute various integrations over projection of face
 * @author Scott Johnson
 *
 * Adapted from volInt.c at
 * http://www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar
 * Disclaimer with volInt.c:
 *
 *  This code computes volume integrals needed for
 *  determining mass properties of polyhedral bodies.
 *
 *  For more information, see the accompanying READM
 *  file, and the paper
 *
 *  Brian Mirtich, "Fast and Accurate Computation of
 *  Polyhedral Mass Properties," journal of graphics
 *  tools, volume 1, number 1, 1996.
 *
 *  This source code is public domain, and may be used
 *  in any way, shape or form, free of charge.
 *
 *  Copyright 1995 by Brian Mirtich
 *  mirtich@cs.berkeley.edu
 *  http://www.cs.berkeley.edu/~mirtich
 *
 * @param[in] p Points on the polyhedron
 * @param[in] tmax Number of triangular facets
 * @param[in] index Indices of the nodes of the triangles in an ordered array of
 * 3-element tuples
 * @param[out] T0 Volume integral
 * @param[out] T Volume integral
 *//*
   static void ComputeVolumeIntegrals( const array<R1Tensor>& p,
                             const unsigned int tmax,
                             const lArray1d& index,
                             realT& T0,
                             R2Tensor& T)
   {
   T0  = 0;
   T = 0.;

   for (unsigned int i = 0; i < tmax; i++)
   {
    // compute face normal and offset w from first 3 vertices
    R1Tensor n;
    {
      R1Tensor dx = p[index[3*i+1]];
      dx -= p[index[3*i]];
      R1Tensor dx2 = p[index[3*i+2]];
      dx2 -= p[index[3*i+1]];
      n.Cross(dx, dx2);
      n.Normalize();
    }

   //    //display
   //    std::cout << "FACE " << i << " of " << tmax << "\n";
   //    for(int j = 0; j < 3; j++)
   //      std::cout << "  v" << j << ": " << p[index[3*i+j]](0) << " " <<
 * p[index[3*i+j]](1) << " " << p[index[3*i+j]](2) << "\n";
   //    std::cout << "  n: " << n(0) << " " << n(1) << " " << n(2) << "\n";

    //get the component of v0 in the opposite normal direction
    realT w = -1. * Dot(n, p[index[3*i]]);

    //get the absolute value of the norm
    R1Tensor nabs = n;
    for(int j = 0; j < 3; j++)
      nabs(j) = fabs(nabs(j));

    int C = 0;
    if (nabs(0) < nabs(1) || nabs(0) < nabs(2))
      C = (nabs(1) > nabs(2)) ? 1 : 2;

    int A = (C + 1) % 3;
    int B = (A + 1) % 3;

    realT Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
    ComputeFaceIntegrals(p, tmax, index, i, w, n,
                         A, B, C,
                         Fa, Fb, Fc,
                         Faa, Fbb, Fcc,
                         Faaa, Fbbb, Fccc,
                         Faab, Fbbc, Fcca);

    T0 += n(0) * ((A == 0) ? Fa : ((B == 0) ? Fb : Fc));

    T(0,A) += n(A) * Faa;
    T(0,B) += n(B) * Fbb;
    T(0,C) += n(C) * Fcc;
    T(1,A) += n(A) * Faaa;
    T(1,B) += n(B) * Fbbb;
    T(1,C) += n(C) * Fccc;
    T(2,A) += n(A) * Faab;
    T(2,B) += n(B) * Fbbc;
    T(2,C) += n(C) * Fcca;
   }

   T(0,0) /= 2;
   T(0,1) /= 2;
   T(0,2) /= 2;
   T(1,0) /= 3;
   T(1,1) /= 3;
   T(1,2) /= 3;
   T(2,0) /= 2;
   T(2,1) /= 2;
   T(2,2) /= 2;
   };
 */
/**
 * @brief Compute the moment of inertia of the polyhedron
 * @author Scott Johnson
 * Adapted from Mirtrich's volInt.c at
 * http://www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar
 * @param[in] p Points on the polyhedron
 * @param[in] tmax Number of triangular facets
 * @param[in] index Indices of the nodes of the triangles in an ordered array of
 * 3-element tuples
 * @param[out] mass Mass of the polyhedron
 * @param[out] cm Centroid of the polyhedron
 * @param[out] inertia Moment of inertia of the polyhedron
 */
/*
   static void CalculatePhysicalProperties2(
    const array<R1Tensor>& p,
    const unsigned int tmax,
    const lArray1d& index,
    realT& mass,
    R1Tensor& cm,
    R2SymTensor& inertia)
   {
   realT density = 1.;

   R2Tensor T;
   ComputeVolumeIntegrals(p, tmax, index, mass, T);

   // compute center of mass
   for(int i = 0; i < 3; i++)
    cm(i) = T(0,i) / mass;
   mass *= density;

   // compute inertia tensor
   inertia(0,0) = density * (T(1,1) + T(1,2));
   inertia(1,1) = density * (T(1,2) + T(1,0));
   inertia(2,2) = density * (T(1,0) + T(1,1));
   inertia(0,1) = inertia(1,0) = - density * T(2,0);
   inertia(1,2) = inertia(2,1) = - density * T(2,1);
   inertia(2,0) = inertia(0,2) = - density * T(2,2);

   // translate inertia tensor to center of mass
   inertia(0,0) -= mass * (cm(1)*cm(1) + cm(2)*cm(2));
   inertia(1,1) -= mass * (cm(2)*cm(2) + cm(0)*cm(0));
   inertia(2,2) -= mass * (cm(0)*cm(0) + cm(1)*cm(1));
   inertia(0,1) = inertia(1,0) += mass * cm(0) * cm(1);
   inertia(1,2) = inertia(2,1) += mass * cm(1) * cm(2);
   inertia(2,0) = inertia(0,2) += mass * cm(2) * cm(0);
   };
 */
#endif /* DISCRETELEMENTMANAGERT_H_ */
