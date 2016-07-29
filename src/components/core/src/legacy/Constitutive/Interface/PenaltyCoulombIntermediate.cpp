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
 * PenaltyCoulombIntermediate.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "PenaltyCoulombIntermediate.h"
#include <typeinfo>
#include <assert.h>

PenaltyCoulombIntermediate::PenaltyCoulombIntermediate( const int paramSize, const int stateSize ):
InterfaceBase( paramSize, stateSize )
{
  // TODO Auto-generated constructor stub

}

PenaltyCoulombIntermediate::~PenaltyCoulombIntermediate()
{

}

realT
PenaltyCoulombIntermediateParameterData::Stiffness(const realT normalApproach) const
{
  realT arealStiffness = 0.0;

  if(normalApproach <= 0)
    arealStiffness = 0.0;
  else if(normalApproach < normalApproachYield)
    arealStiffness = ktildeAperture / (aperture - normalApproach);
  else if(normalApproach < normalApproachSoften)
    arealStiffness = kyield;
  else
    arealStiffness = ksoften;

  return arealStiffness;
}

realT
PenaltyCoulombIntermediateParameterData::Stress(const realT normalApproach) const
{
  realT normalStress = 0.0;

  if(normalApproach <= 0)
    normalStress = 0;
  else if(normalApproach < normalApproachYield)
    normalStress = ktildeAperture * log(aperture / (aperture - normalApproach));
  else if(normalApproach < normalApproachSoften)
    normalStress = stressYield + kyield * (normalApproach - normalApproachYield);
  else
    normalStress = stressSoften + ksoften * (normalApproach - normalApproachSoften);

  return normalStress;
}
void
PenaltyCoulombIntermediateParameterData::Initialize()
{
  //check user settings
  const realT tol = 1e-6;
  if(normalApproachYield > aperture)
    normalApproachYield = (1.0 - tol) * aperture;
  if(stressYield > stressSoften)
    stressYield = stressSoften;

  //set derived values
  const realT factor = 1.0/(aperture - normalApproachYield);
  ktildeAperture = stressYield / log(aperture * factor);
  kyield = ktildeAperture * factor;
  if(isEqual(stressSoften, std::numeric_limits<realT>::max()))
    normalApproachSoften = std::numeric_limits<realT>::max();
  else
    normalApproachSoften = normalApproachYield + (stressSoften - stressYield) / kyield;
}

void
PenaltyCoulombIntermediateParameterData::PostReadXML( TICPP::HierarchicalDataNode& )
{
  if(aperture <= 0)
    aperture = 1e-2;

  //normal approach at the onset of yielding
  if(normalApproachYield <= 0)
    throw GPException("You must set a normal approach for the onset of yielding > 0");
  if(normalApproachYield >= aperture)
    throw GPException("You must set a normal approach for the onset of yielding < aperture");
  if(stressYield <= 0)
    throw GPException("You must set a stress for the onset of yielding > 0");
  if(stressSoften <= 0)
    stressSoften = std::numeric_limits<realT>::max();
  if(stressSoften < stressYield)
    stressSoften = stressYield;

  //set all derived values
  Initialize();

  //set the rest ... arealStiffnessSoften is not calculated or used in SimpleInitialize
  if(ksoften <= 0)
    ksoften = (1e-5) * kyield;
  if(kshear <= 0)
    kshear = 0.7 * kyield;
}

void
PenaltyCoulombIntermediateParameterData::PostSetValue()
{
  TICPP::HierarchicalDataNode node;
  PostReadXML(node);
}
void
PenaltyCoulombIntermediate::StrainDrivenUpdate( const localIndex index )
{
  //get temporary references
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);

  //(1) evolve normal stress
  matState.stress = matParams.Stress(matState.normalApproach);

  //(2) evolve friction
  UpdateFriction(matParams, matState);

  //(3) evolve shear stress
  const realT dssdt = matParams.kshear * matState.dxsdt;
  matState.stressShear += dssdt * matState.dt;

  //(4) return shear stress to the failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}

realT
PenaltyCoulombIntermediate::StiffnessProjected(const localIndex index)
{
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  const PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);
  const realT normalApproachEstimate = matState.normalApproach +
      (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);
  return matParams.Stiffness(normalApproachEstimate);
}
