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
 * MaterialBase.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "MaterialBase.h"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

MaterialBase::MaterialBase( const int paramSize, const int stateSize ):
ConstitutiveBase(),
m_paramSize( paramSize ),
m_stateSize( stateSize )
{

}

MaterialBase::~MaterialBase()
{

}

void
MaterialBaseStateData::TotalStress(R2SymTensor& totalStress) const
{
  totalStress = devStress;
  totalStress.PlusIdentity(pressure);
}

void
MaterialBaseStateData::RotateState( const R2Tensor& Rot )
{
  R2SymTensor temp;

  devStress.PlusIdentity( pressure );
  temp.QijAjkQlk(devStress,Rot);
  pressure = temp.Trace() / 3.0;
  devStress = temp;
  devStress.PlusIdentity(-pressure);
}
void
MaterialBase::StrainDrivenUpdateMember( const localIndex index0,
                            const localIndex index1,
                            const R2SymTensorT < 3 >& Ddt,
                            const R2TensorT < 3 >& L,
                            const R2Tensor& Rot,
                            const realT dt )
{
  throw GPException("Cannot call MaterialBase::StrainDrivenUpdateMember; must have derived method\n");
}

void
MaterialBase::StrainDrivenUpdateMember( const localIndex index0,
                            const localIndex index1,
                            const R2SymTensorT < 3 >& Ddt,
                            const R2TensorT < 3 >& L,
                            const R2Tensor& Rot,
                            const realT& volume_n,
                            const realT& volume_np1,
                            const realT dt)
{
  throw GPException("Cannot call MaterialBase::StrainDrivenUpdateMember; must have derived method\n");
}

void
MaterialBase::MeanPressureDevStress( const localIndex index,
                                     realT& pressure, R2SymTensor& devStress) const
{
  throw GPException("Cannot call MaterialBase::MeanPressureDevStress; must have derived method\n");
}

