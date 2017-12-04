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
 * @file ParallelPlateFluidModels.cpp
 * @author walsh24
 * @date March 10, 2014
 */

#include "ParallelPlateFluidModels.h"

///////////////////////////////////////////////////////////////////////////////

/**
 *  Newtonian fluid model
 */

NewtonianFluidModel::NewtonianFluidModel():
  ParallelPlateFluidModelBase(), m_mu(){ /** empty **/};

void NewtonianFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
}

realT NewtonianFluidModel::CalculatePermeability(const realT la, const realT lb,
                                                 const realT apa, const realT apb,
                                                 const realT w, const realT qMag,
                                                 const realT SHP_FCT){
  return PPFS::CalculatePermeability(la,lb,apa,apb, w, m_mu, SHP_FCT);
}

// one sided permeability
realT NewtonianFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                 const realT w, const realT qMag,
                                                 const realT SHP_FCT){
  return PPFS::CalculatePermeability(l,ap, w,  m_mu, SHP_FCT);
}

/// Register Fluid model in the solver factory
REGISTER_PARALLEL_PLATE_FLUID_MODEL( NewtonianFluidModel )

///////////////////////////////////////////////////////////////////////////////

/**
 *  Power law fluid model
 */


PowerlawFluidModel::PowerlawFluidModel():
  ParallelPlateFluidModelBase(), m_n(),m_M(),m_phiM(){ /** empty **/};

void PowerlawFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  m_M = hdn->GetAttributeOrDefault("ConsistencyIndex","1.0e-3 N.s/m^2");
  m_n = hdn->GetAttributeOrDefault("FluidBehaviorIndex","1");  // newtonian by
                                                               // default
  m_phiM = PPFS::CalculateModifiedConsistencyIndex(m_M,m_n);
}

realT PowerlawFluidModel::CalculatePermeability(const realT la, const realT lb,
                                                const realT apa, const realT apb,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
  return PPFS::CalculatePermeability_PowerLawFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, SHP_FCT);
}

// one sided permeability
realT PowerlawFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
  return PPFS::CalculatePermeability_PowerLawFluid(l,ap, w, qMag, m_phiM, m_n, SHP_FCT);
}

/// Register Fluid model in the solver factory
REGISTER_PARALLEL_PLATE_FLUID_MODEL( PowerlawFluidModel )


///////////////////////////////////////////////////////////////////////////////

/**
 *  Herschel Bulkley Parallel Plate Fluid Model
 */
HerschelBulkleyParallelPlateFluidModel::HerschelBulkleyParallelPlateFluidModel( ):
  PowerlawFluidModel(),
  m_tau_y(0)
{}

void HerschelBulkleyParallelPlateFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
  PowerlawFluidModel::ReadXML(hdn);
  m_tau_y = hdn->GetAttributeOrDefault("YieldStress","0.0");
  m_phi = m_phiM/m_M;

  /* build lookup table */
  array<real64> xs;
  array<real64> values;
  xs.resize(16);
  values.resize(16);
  for(int i =0 ; i < 16 ; i++)
  {
    realT x = -3 + 0.2*i;
    xs[i] = x;
    realT kappa = pow(10,x);
    values[i] = CalculateZp(kappa,m_n);
  }
  m_zp_kappaLookup.SetAxisValues(0,xs);
  m_zp_kappaLookup.SetValues(values);

  for(int i =0 ; i < 16 ; i++)
  {
    std::cout <<  "k: " <<  pow(10,xs[i])
              << " v: " <<  values[i]
              << "\n";
  }


}

// analytical solution relating dimensionless parameter V to zp (dimensionless
// plug thickness)
realT HerschelBulkleyParallelPlateFluidModel::analyticalVFunc(realT zp, realT n){
  realT V = pow(1-zp,n+1)*pow(n/(n+1)*zp + 1,n);
  return V;
}
// derivative of analyticalVFunc wrt zp
realT HerschelBulkleyParallelPlateFluidModel::dVdz(realT zp, realT n){

  realT deriv = -(n+1)*pow(1-zp,n) *pow(n/(n+1)*zp + 1,n) +
                pow(1-zp,n+1)*(n*n/(n+1))*pow(n/(n+1)*zp + 1,n-1);
  return deriv;
}

// calculate Zp as a function of kappa = Kq^n
realT HerschelBulkleyParallelPlateFluidModel::CalculateZp(realT Kqn, realT n){

  realT zp = 0.5;
  realT dzp = 1;
  int numIters = 0;
  realT tol = 1e-8;
  // newton's method
  while(fabs(dzp) > tol*zp && numIters < 500)
  {
    realT f = Kqn*zp - analyticalVFunc(zp,n);
    realT dfdz = Kqn - dVdz(zp,n);
    dzp = -f/dfdz;
    zp += dzp;
    numIters += 1;
  }
  if(numIters == 500)
    throw GPException("HerschelBulkleyParallelPlateFluidModel: Zp calculation failed to converge");
  return zp;
};

realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT la, const realT lb,
                                                                    const realT apa, const realT apb,
                                                                    const realT w, const realT qMag,
                                                                    const realT SHP_FCT){
  realT kappa = m_M*pow(qMag,m_n);
  realT zp(0.0);
  if(kappa > 1000.0)
  {
    zp = 1.0/kappa;
  }
  else if(kappa < 1e-3)
  {
    zp = 1 - pow(kappa/(m_n/(m_n+1) -1), 1.0/(m_n+1));
  }
  else
  {
    realT xx[1];
    xx[0] = log10(kappa);
    zp = m_zp_kappaLookup.Lookup(xx);
  }
  return PPFS::CalculatePermeability_HerschelBulkleyFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, zp, SHP_FCT);
}

// one sided permeability
realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                                    const realT w, const realT qMag,
                                                                    const realT SHP_FCT){
  realT kappa = m_M*pow(qMag,m_n);
  realT zp(0.0);
  if(kappa > 1000.0)
  {
    zp = 1.0/kappa;
  }
  else if(kappa < 1e-3)
  {
    zp = 1 - pow(kappa/(m_n/(m_n+1) -1), 1.0/(m_n+1));
  }
  else
  {
    realT xx[1];
    xx[0] = log10(kappa);
    zp = m_zp_kappaLookup.Lookup(xx);
  }
  return PPFS::CalculatePermeability_HerschelBulkleyFluid(l,ap, w, qMag, m_phiM, m_n,zp, SHP_FCT);
}

REGISTER_PARALLEL_PLATE_FLUID_MODEL( HerschelBulkleyParallelPlateFluidModel )
