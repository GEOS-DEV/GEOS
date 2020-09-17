/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AuxiliaryFunctionsSpectral.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_
#include <algorithm>
#include <cmath> 
#include "LvArray/src/tensorOps.hpp"
namespace geosx
{

//get the positive part only of tensor T using spectral split
GEOSX_HOST_DEVICE inline
void PositivePartOfTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&positivePart)[6] )
{
  real64 positiveEigs[6];
  for (int i=0; i < 3; i++)
    {
      positiveEigs[i] = std::max(0.0, eigs[i]); 
    }
  LvArray::tensorOps::AikSymBklAjl<3>(positivePart, eigvecs, positiveEigs); 
}

GEOSX_HOST_DEVICE inline
void NegativePartOfTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&negativePart)[6] )
{
  real64 negativeEigs[6];
  for (int i=0; i < 3; i++)
    {
      negativeEigs[i] = std::min(0.0, eigs[i]); 
    }
  LvArray::tensorOps::AikSymBklAjl<3>(negativePart, eigvecs, negativeEigs); 
}

GEOSX_HOST_DEVICE inline
void recoverStrainFromStress(arraySlice1d<real64> const stress, real64 (&strain)[6], real64 const K, real64 const mu)
{
  real64 E = 9*K*mu / (3*K + mu);
  real64 nu = (3*K - 2*mu) / (6*K + 2*mu);
  strain[0] = (stress[0] - nu*(stress[1] + stress[2]))/E;
  strain[1] = (stress[1] - nu*(stress[0] + stress[2]))/E;
  strain[2] = (stress[2] - nu*(stress[0] + stress[1]))/E;
  strain[3] = (1 + nu)*stress[3]/E;
  strain[4] = (1 + nu)*stress[4]/E;
  strain[5] = (1 + nu)*stress[5]/E;
}

GEOSX_HOST_DEVICE inline
int voigt(int voigtIndex, int pairIndex)
{
  if (pairIndex != 1 && pairIndex != 2){
    GEOSX_ERROR("Index of the pair must be 1 or 2");
  }
  switch (voigtIndex) {
    case 0:
      if (pairIndex == 1){
	return 1;
      }
      else {
	return 1;
      }
      break;
    case 1:
      if (pairIndex == 1){
	return 2;
      }
      else {
	return 2;
      }
      break;
    case 2:
      if (pairIndex == 1){
	return 3;
      }
      else {
	return 3;
      }
      break;
    case 3:
      if (pairIndex == 1){
	return 2;
      }
      else {
	return 3;
      }
      break;
    case 4:
      if (pairIndex == 1){
	return 1;
      }
      else {
	return 3;
      }
      break;
    case 5:
      if (pairIndex == 1){
	return 1;
      }
      else {
	return 2;
      }
      break;
    default:  
      GEOSX_ERROR("Voigt index must be from 0 to 5");
      return 1000000;
    }    
}

GEOSX_HOST_DEVICE inline
real64 heaviside(real64 x)
{
  real64 tol = 1e-12;
  if (std::abs(x) < tol) {
    return 0.5;
  }
  if (x > 0) {
    return 1;
  }
  if (x < 0) {
    return 0;
  }
  GEOSX_ERROR("This function was not supposed to reach this line");
  return 1000000;
}

GEOSX_HOST_DEVICE inline
void QTensor( real64 (&eigvector)[3], real64 (&Q)[6][6] )
{
  real64 M[3][3];
  LvArray::tensorOps::AiBj<3,3>(M, eigvector, eigvector);
  for (int i = 0; i<6; i++)
  {
    for (int j = 0; j<6; j++)
      {
	Q[i][j] = M[voigt(i,1)][voigt(i,2)]*M[voigt(j,1)][voigt(j,2)];
      }
  }  	  
}

GEOSX_HOST_DEVICE inline
void GTensor( real64 (&eigvec1)[3], real64 (&eigvec2)[3], real64 (&G)[6][6] )
{
  real64 M1[3][3];
  real64 M2[3][3];
  LvArray::tensorOps::AiBj<3,3>(M1, eigvec1, eigvec1);
  LvArray::tensorOps::AiBj<3,3>(M2, eigvec2, eigvec2);
  for (int i = 0; i<6; i++)
  {
    for (int j = 0; j<6; j++)
      {
	G[i][j] = M1[voigt(i,1)][voigt(j,1)]*M2[voigt(i,2)][voigt(j,2)] + M1[voigt(i,1)][voigt(j,2)]*M2[voigt(i,2)][voigt(j,1)];
      }
  }  	  
}

GEOSX_HOST_DEVICE inline
void PositiveProjectorTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&PositiveProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-8;
  if ( std::abs(eigs[0] - eigs[1]) < tol || std::abs(eigs[0]-eigs[2]) < tol || std::abs(eigs[1]-eigs[2]) < tol ) {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6];
  //init GVoigt
  real64 Gsym[6][6];
  real64 Gji[6][6];

  //compute projector
  for (int i = 0; i < 3; i++){
    real64 ithEigenVector[3];
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor(ithEigenVector, Qi);
    //do update
    LvArray::tensorOps::scale<6,6>(Qi, heaviside(eigs[i]));
    LvArray::tensorOps::add<6,6>(PositiveProjector, Qi);
    if (!repeatedEigenvalues) { 

      for(int j = 0; j < 3; j++){
	real64 jthEigenVector[3];
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	LvArray::tensorOps::add<6,6>(Gsym,Gji);
	//Do update
	LvArray::tensorOps::scale<6,6>(Gsym, 0.5 * (std::max(eigs[i],0.0) + std::max(eigs[j],0.0))/(2*(eigs[i]+eigs[j])));
	LvArray::tensorOps::add<6,6>(PositiveProjector, Gsym); 
      }
    }
    else {
      for(int j = 0; j < 3; j++){
	real64 jthEigenVector[3];
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	LvArray::tensorOps::add<6,6>(Gsym,Gji);
	LvArray::tensorOps::scale<6,6>(Gsym, 0.5 * (heaviside(eigs[i]) + heaviside(eigs[j]))/4);
	//do update
	LvArray::tensorOps::add<6,6>(PositiveProjector, Gsym); 
      }
    }
  }
  								      
}

GEOSX_HOST_DEVICE inline
void NegativeProjectorTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&NegativeProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-8;
  if ( std::abs(eigs[0] - eigs[1]) < tol || std::abs(eigs[0]-eigs[2]) < tol || std::abs(eigs[1]-eigs[2]) < tol ) {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6];
  //init GVoigt
  real64 Gsym[6][6];
  real64 Gji[6][6];

  //compute projector
  for (int i = 0; i < 3; i++){
    real64 ithEigenVector[3];
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor(ithEigenVector, Qi);
    //do update
    LvArray::tensorOps::scale<6,6>(Qi, heaviside(-eigs[i]));
    LvArray::tensorOps::add<6,6>(NegativeProjector, Qi);
    if (!repeatedEigenvalues) { 

      for(int j = 0; j < 3; j++){
	real64 jthEigenVector[3];
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	LvArray::tensorOps::add<6,6>(Gsym,Gji);
	//Do update
	LvArray::tensorOps::scale<6,6>(Gsym, 0.5 * (std::min(eigs[i],0.0) + std::min(eigs[j],0.0))/(2*(eigs[i]+eigs[j])));
	LvArray::tensorOps::add<6,6>(NegativeProjector, Gsym); 
      }
    }
    else {
      for(int j = 0; j < 3; j++){
	real64 jthEigenVector[3];
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	LvArray::tensorOps::add<6,6>(Gsym,Gji);
	LvArray::tensorOps::scale<6,6>(Gsym, 0.5 * (heaviside(-eigs[i]) + heaviside(-eigs[j]))/4);
	//do update
	LvArray::tensorOps::add<6,6>(NegativeProjector, Gsym); 
      }
    }
  }
  								      
}

} //namespacegeosx 
  
#endif /* GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_ */
