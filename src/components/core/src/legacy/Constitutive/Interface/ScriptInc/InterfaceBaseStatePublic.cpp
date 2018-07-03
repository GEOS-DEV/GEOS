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
