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
//get the positive part only of tensor T using spectral split
GEOSX_HOST_DEVICE
void PositivePartOfTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&positivePart)[6] )
{
  real64 positiveEigs[6];
  for (int i=0; i < 3; i++)
    {
      positiveEigs[i] = std::max(0, eigs[i]); 
    }
  AikSymBklAjl<3>(positivePart, eigvecs, positiveEigs); 
}

GEOSX_HOST_DEVICE
void NegativePartOfTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&negativePart)[6] )
{
  real64 negativeEigs[6];
  for (int i=0; i < 3; i++)
    {
      negativeEigs[i] = std::min(0, eigs[i]); 
    }
  AikSymBklAjl<3>(negativePart, eigvecs, negativeEigs); 
}

GEOSX_HOST_DEVICE
void recoverStrainFromStress(array1d const stress, real64 (&strain)[6], real64 const K, real64 const mu)
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

GEOSX_HOST_DEVICE
int voigt(int voigtIndex, int pairIndex)
{
  if (pairIndex != 1 && pairIndex != 2){
    GEOSX_ERROR("Index of the pair must be 1 or 2");
  }
  switch voigtIndex {
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
      return NaN;
    }    
}

GEOSX_HOST_DEVICE
real64 heaviside(real64 x)
{
  real64 tol = 1e-12;
  if (abs(x) < tol) {
    return 0.5;
  }
  if (x > 0) {
    return 1;
  }
  if (x < 0) {
    return 0;
  }
  GEOSX_ERROR("This function was not supposed to reach this line")
  return NaN;
}

GEOSX_HOST_DEVICE
void QTensor( real64 (&eigvector)[3], real64 (&Q)[6][6] )
{
  real64 M[3][3];
  AiBj<3,3>(M, eigvector, eigvector);
  for (int i = 0; i<6; i++)
  {
    for (int j = 0; j<6; j++)
      {
	Q[i][j] = M[voigt(i,1)][voigt(i,2)]*M[voigt(j,1)][voigt(j,2)];
      }
  }  	  
}

GEOSX_HOST_DEVICE
void GTensor( real64 (&eigvec1)[3], real64 (&eigvec2)[3], real64 (&G)[6][6] )
{
  real64 M1[3][3];
  real64 M2[3][3];
  AiBj<3,3>(eigvec1, eigvec1);
  AiBj<3,3>(eigvec2, eigvec2);
  for (int i = 0; i<6; i++)
  {
    for (int j = 0; j<6; j++)
      {
	G[i][j] = M1[voigt(i,1)][voigt(j,1)]*M2[voigt(i,2)][voigt(j,2)] + M1[voigt(i,1)][voigt(j,2)]*M2[voigt(i,2)][voigt(j,1)];
      }
  }  	  
}

GEOSX_HOST_DEVICE
void PositiveProjectorTensor( real64 (&eigs)[3], real64 (&eigvecs)[3][3], real64 (&PositiveProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  tol = 1e-8
  if ( abs(eigs[0] - eigs[1]) < tol || abs(eigs[0]-eigs[2]) < tol || abs(eigs[1]-eigs[2]) < tol ) {
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
    scaledAdd(PositiveProjector, Qi, heaviside(eigs[i]));
    if (!repeatedEigenvalues) { 

      for(int j = 0; j < 3; j++){
	real64 jthEigenVector;
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	add(Gsym,Gji);
	//Do update
	scaledAdd(PositiveProjector, Gsym, 0.5 * (max(eigs[i],0) + max(eigs[j],0))/(2*(eigs[i]+eigs[j]))); 
      }
    }
    else {
      for(int j = 0; j < 3; j++){
	real64 jthEigenVector;
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	add(Gsym,Gji);
	//do update
	scaledAdd(PositiveProjector, Gsym, 0.5 * (heaviside(eigs[i]) + heaviside(eigs[j]))/4 ); 
      }
    }
  }
  								      
}

GEOSX_HOST_DEVICE
void NegativeProjectorTensor( real64[3] eigs, real64[3][3] eigvecs, voigtTensor P )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  tol = 1e-8
  if ( abs(eigs[0] - eigs[1]) < tol || abs(eigs[0]-eigs[2]) < tol || abs(eigs[1]-eigs[2]) < tol ) {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6];
  //init GVoigt
  real64 Gij[6][6];
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
    scaledAdd(NegativeProjector, Qi, heaviside(-eigs[i]) );
    if (!repeatedEigenvalues) { 

      for(int j = 0; j < 3; j++){
	real64 jthEigenVector;
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	add(Gsym,Gji);
	//Do update
	scaledAdd(NegativeProjector, Gsym, 0.5 * (min(eigs[i],0) + min(eigs[j],0))/(2*(eigs[i]+eigs[j]))); 
      }
    }
    else {
      for(int j = 0; j < 3; j++){
	real64 jthEigenVector;
	jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
	if (i == j) {
	  continue;
	}
	//compute Gij and Gji
	GTensor(ithEigenVector, jthEigenVector, Gsym);
	GTensor(jthEigenVector, ithEigenVector, Gji);
	add(Gsym,Gji);
	//Do update
	scaledAdd(NegativeProjector, Gsym, 0.5 * (heaviside(-eigs[i]) + heaviside(-eigs[j]))/4);  
      }
    }
  }
}

#endif /* GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_ */
