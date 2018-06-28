/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
