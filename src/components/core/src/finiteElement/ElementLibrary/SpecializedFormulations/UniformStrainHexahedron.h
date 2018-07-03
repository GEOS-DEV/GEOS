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
 * @file UniformStrainHexahedron.h
 * @author settgast1
 * @date Jun 6, 2011
 */

#include "legacy/ElementLibrary/FiniteElement.h"

#ifndef UNIFORMSTRAINHEXAHEDRON_H_
#define UNIFORMSTRAINHEXAHEDRON_H_

class UniformStrainHexahedron : public FiniteElement<3>
{
public:
  UniformStrainHexahedron();
  virtual ~UniformStrainHexahedron();

  void reinit(const std::vector<R1TensorT<3> > &mapped_support_points);

  void zero_energy_mode_control( const array<R1Tensor>& dNdx,
                                 const realT& volume,
                                 const array<R1Tensor>& x,
                                 const array<R1Tensor>& vel,
                                 const realT& dampcoef,
                                 const realT& stiffcoef,
                                 const realT& rho,
                                 const realT& modulus,
                                 const realT& dt,
                                 array<R1Tensor>& Qstiffness,
                                 array<R1Tensor>& force );

};

#endif /* UNIFORMSTRAINHEXAHEDRON_H_ */
