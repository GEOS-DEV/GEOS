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
 * @file ParallelPlateFluidModels.h
 * @author walsh24
 * @date March 10, 2014
 */

#ifndef PARALLELPLATEFLUIDMODELS_H_
#define PARALLELPLATEFLUIDMODELS_H_

#include "ParallelPlateFluidModelBase.h"

/*
 * Newtonian Fluid Model
 *
 */

class NewtonianFluidModel : public ParallelPlateFluidModelBase
{

public:
  NewtonianFluidModel();
  ~NewtonianFluidModel();

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
  virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                      const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
  virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                      const realT qMag, const realT SHP_FCT); // one
                                                                              // sided
  static const char* FluidModelName(){return "Newtonian";};

  realT m_mu; // viscosity
};


/*
 * Powerlaw Fluid Model
 *
 */

class PowerlawFluidModel : public ParallelPlateFluidModelBase
{

public:
  PowerlawFluidModel();
  ~PowerlawFluidModel();

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
  virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                      const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
  virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                      const realT qMag, const realT SHP_FCT); // one
                                                                              // sided
  static const char* FluidModelName(){return "PowerlawFluid";};

  realT m_n; // fluid behavior index
  realT m_M; // consistency index
  realT m_phiM;// modified consistency index (phi accounts for geometry of flow
               // between parallel plates)
};

/*
 * The Herschel Bulkley model describes a non-Newtonian fluid with a Bingham
 * law-like yield stress
 * accompanied by a power-law fluid like growth in stress as a function of shear
 * strain rate:
 *
 * \tau = \tau_y + k |du/dr|^n   for |\tau| > \tau_y
 * du/dr = 0                     for |\tau| < \tau_y
 *
 * The parallel plate model implemented here follows the approximate solution
 * provided in
 * Wang and Gordaninejad 1999, Flow analysis of field controllable, electro- and
 * magneto-rheological fluids
 * using Herschel-Bulkley model.
 *
 */
class HerschelBulkleyParallelPlateFluidModel : public PowerlawFluidModel
{

public:
  HerschelBulkleyParallelPlateFluidModel();
  ~HerschelBulkleyParallelPlateFluidModel(){ /** empty **/};

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
  virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                      const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
  virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                      const realT qMag, const realT SHP_FCT); // one
                                                                              // sided
  static const char* FluidModelName(){return "HerschelBulkleyFluid";};


  realT m_tau_y; // fluid yield stress
  realT m_phi; // modified consistency index factor

private:
  Table<1, realT> m_zp_kappaLookup;

  realT CalculateZp(realT kappa, realT n);
  /*
     lookup table relating the dimensionless plug thickness (zp) to the
        dimensionless variable
     kappa = K q^n where
     K = ((2n+1)/n)^n 2^n/h^(n+1) * M/\tau_y = phi/2* 1/h^(n+1) *M/\tau_y
   */
  // functions used in Newton's method
  realT analyticalVFunc(realT zp, realT n);
  realT dVdz(realT zp, realT n);

};

#endif
