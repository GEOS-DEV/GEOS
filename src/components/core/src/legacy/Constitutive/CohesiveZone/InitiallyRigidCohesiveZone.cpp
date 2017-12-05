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
 * InitiallyRigidCohesiveZone.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "DataStructures/VectorFields/ElementRegionT.h"
#include "Constitutive/Material/MaterialBase.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "InitiallyRigidCohesiveZone.h"
#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

InitiallyRigidCohesiveZone::InitiallyRigidCohesiveZone( ):
  CohesiveZoneBase( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

InitiallyRigidCohesiveZone::~InitiallyRigidCohesiveZone()
{}

int
InitiallyRigidCohesiveZone::UpdateCohesiveZone( const localIndex index,
                                                const R1Tensor& gap,
                                                const R1Tensor& N,
                                                const std::pair< ElementRegionT*, localIndex >& elem0,
                                                const std::pair< ElementRegionT*, localIndex >& elem1,
                                                R1Tensor& traction,
                                                R2Tensor& stiffness )
{
  const InitiallyRigidCohesiveZoneParameterData& params = *ParameterData(index);
  InitiallyRigidCohesiveZoneStateData& state      = *StateData(index,0);

  int rval = 2;


  realT gapMag = gap.L2_Norm();

  realT k;

  const realT s = gapMag / params.failGap < 1.0 ?  gapMag / params.failGap : 1.0;
  bool mainline = false;

  state.maxGap.SetMax( gap );

  if( state.separationCoeff<s )
  {
    state.separationCoeff = s;
    mainline = true;
  }


  R1Tensor direction;

  traction = gap;
  traction.Normalize();

  if( isZero(gapMag,1e-14) )
  {
    R2SymTensor stress0;
    R2SymTensor stress1;


    ElementRegionT& er0 = *(elem0.first);
    ElementRegionT& er1 = *(elem1.first);
    const localIndex elemIndex0 = elem0.second;
    const localIndex elemIndex1 = elem1.second;

    realT pressure;
    er0.m_mat->MeanPressureDevStress(elemIndex0,pressure, stress0);
    stress0.PlusIdentity(pressure);

    er1.m_mat->MeanPressureDevStress(elemIndex1,pressure, stress1);
    stress1.PlusIdentity(pressure);

    stress0 += stress1;
    stress0 *= 0.5;

    direction.AijBj( stress0, N );
    direction.Normalize();
    gapMag = 1.0e-14;
  }
  else
  {
    direction = gap;
    direction /= gapMag;
  }

  direction[0] = 0;
  direction[2] = 0;
  //todo: why is this constrained this way??
  traction = direction;

  realT scale = 0.0;
  if( state.separationCoeff >=1.0 )
  {
    // fully separated
    state.maxTraction = 0.0;
    k = 0.0;
    traction = 0.0;
    rval = 3;
  }
  else if( isZero(state.separationCoeff) )
  {
    // initial
    traction *= params.failStress;
    k = -params.failStress/params.failGap;
  }
  else if( mainline || params.unloadFlag==0 )
  {
    // going down the curve
    scale = 1.0 - state.separationCoeff;
    state.maxTraction = params.failStress * scale;
    traction *= state.maxTraction;
    k = -params.failStress/params.failGap;
  }
  else
  {
    // unload/reload line
    scale = s/state.separationCoeff;
    traction *= scale * state.maxTraction;

    k = state.maxTraction/(state.separationCoeff*params.failGap);

  }


  realT dterm = traction.L2_Norm() / gapMag;

  R2Tensor kI;
  kI.PlusIdentity( dterm );

  stiffness.dyadic_ab( direction, direction );
  stiffness *= k - dterm;

  stiffness += kI;

  state.traction = traction;
  state.stiffness = stiffness;

  return rval;

}

/// Register class in the class factory
REGISTER_COHESIVEZONE( InitiallyRigidCohesiveZone )
