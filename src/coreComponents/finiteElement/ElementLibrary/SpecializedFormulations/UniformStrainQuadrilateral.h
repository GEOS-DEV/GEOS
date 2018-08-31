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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file UniformStrainQuadrilateral.h
 * @author settgast1
 * @date Jun 6, 2011
 */

#include "legacy/ElementLibrary/FiniteElement.h"

#ifndef UniformStrainQuadrilateral_H_
#define UniformStrainQuadrilateral_H_

class UniformStrainQuadrilateral : public FiniteElement<3>
{
public:
  UniformStrainQuadrilateral();
  virtual ~UniformStrainQuadrilateral();

  void reinit(const std::vector<R1TensorT<3> > &mapped_support_points);

  void zero_energy_mode_control( const array1d<R1Tensor>& dNdx,
                                 const realT& volume,
                                 const array1d<R1Tensor>& x,
                                 const array1d<R1Tensor>& vel,
                                 const realT& dampcoef,
                                 const realT& stiffcoef,
                                 const realT& rho,
                                 const realT& modulus,
                                 const realT& dt,
                                 array1d<R1Tensor>& Qstiffness,
                                 array1d<R1Tensor>& force );

};

#endif /* UniformStrainQuadrilateral_H_ */
