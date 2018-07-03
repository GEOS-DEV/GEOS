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

//FUNCTION_BEGIN_PARSE
void
MaterialBaseStateData::TotalStress(R2SymTensor& totalStress) const
{
  totalStress = devStress;
  totalStress.PlusIdentity(pressure);
}

//FUNCTION_BEGIN_PARSE
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
