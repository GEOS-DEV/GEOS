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


/*
 * InterfaceBase.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "InterfaceBase.h"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

InterfaceBase::InterfaceBase( const int paramSize, const int stateSize ):
  ConstitutiveBase(),
  m_paramSize( paramSize ),
  m_stateSize( stateSize )
{}

InterfaceBase::~InterfaceBase()
{}

void
InterfaceBaseStateData::UpdateOrientation(const realT normalApproachIn,
                                          const realT dtIn,
                                          const R1Tensor& normal,
                                          const R1Tensor& velocity,
                                          R1Tensor& dShearSlip,
                                          R1Tensor& shearSlip)
{
  dt = dtIn;
  normalApproach = normalApproachIn;

  dxndt = Dot(velocity, normal);

  //transform the tangential stress and slip
  const realT shearSlipNormal = Dot(shearSlip, normal);
  if (!isZero(shearSlipNormal))
  {
    //project shear slip and shear stress to current coordinates
    const realT shearSlipMagnitude = shearSlip.L2_Norm();
    R1Tensor slipNormal(normal);
    slipNormal *= shearSlipNormal;
    shearSlip -= slipNormal;
    shearSlip.Normalize();
    stressShearVector = shearSlip;
    shearSlip *= shearSlipMagnitude;
  }
  else
  {
    stressShearVector = shearSlip;
    stressShearVector.Normalize();
  }
  stressShearVector *= -stressShear;

  //update the contact velocity
  GeometryUtilities::OrthogonalVectorComponent(velocity, normal, dShearSlip);

  //next statement assumes ever-increasing slip ... should handle reversal of
  // slip, too!
  dxsdt = dShearSlip.L2_Norm();
  dShearSlip *= dt;

  //update the slip magnitude
  xs = shearSlip.L2_Norm();
}
void
InterfaceBase::Initialize( const localIndex index,
                           const realT stressNormal,
                           const realT stressShear)
{
  throw GPException("Cannot call InterfaceBase::Initialize\n");
}

void
InterfaceBase::UpdateProperties(const localIndex, std::map<std::string, realT>&, std::map<std::string, realT>& )
{}

realT
InterfaceBase::SetPermeabilityTerm(const localIndex index)
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  if(matState.normalApproach <= 0)
  {
    matState.kappa = matParams.kappa0 - (matState.normalApproach * matState.normalApproach * matState.normalApproach);
  }
  else
  {
    matState.kappa = matParams.kappa0 +
                     matParams.dkappadnFct *
                     (isEqual(matParams.dkappadnExp, 1.0) ? matState.normalApproach : pow(matState.normalApproach, matParams.dkappadnExp)) +
                     matParams.dkappadsFct * (isEqual(matParams.dkappadsExp, 1.0) ? matState.xs : pow(matState.xs, matParams.dkappadsExp));
  }
  return matState.kappa;
}

realT
InterfaceBase::StiffnessProjected( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  const realT normalApproachEstimate = matState.normalApproach +
                                       (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);

  return NormalStiffness(matParams, matState, normalApproachEstimate, false);
}

void
InterfaceBase::StrainDrivenUpdate( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  //evolve the normal force ("stress" ... bad terminology, but I like the
  // inheritance)
  const realT stiffness = NormalStiffness(matParams, matState, matState.normalApproach, true);

  //evolve the friction
  UpdateFriction(matParams, matState);

  //evolve the shear stress
  {
    matState.stressShear += -stiffness * matParams.ston * matState.dxsdt * matState.dt;
    matState.xs += matState.dxsdt * matState.dt;
  }

  //determine whether strength exceeded and return to failure surface if
  // necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  //update vector representation of shear stress
  //NOTE: JointBaseStateData::UpdateOrientation sets the direction as -(shear
  // slip direction)
  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}
realT InterfaceBase::ShearStrength(const InterfaceBaseParameterData& matParams,
                                   InterfaceBaseStateData& matState) const
{
  return matState.stress * matParams.mu0;
}

void InterfaceBase::UpdateFriction(const InterfaceBaseParameterData&,
                                   InterfaceBaseStateData&) const
{}

realT
InterfaceBase::NormalStiffness(const InterfaceBaseParameterData&,
                               InterfaceBaseStateData&,
                               const realT,
                               const bool ) const
{
  return 0.0;
}

void
InterfaceBase::ThresholdToFailureSurface(const InterfaceBaseParameterData& matParams,
                                         InterfaceBaseStateData& matState,
                                         const realT dxs)
{
  //determine whether strength exceeded and return to failure surface if
  // necessary
  const realT strength = ShearStrength( matParams, matState);
  if(strength < fabs(matState.stressShear))
  {
    matState.DissipatedEnergy += fabs(dxs) * strength;
    matState.stressShear = matState.stressShear < 0 ? -strength : strength;
  }
}
