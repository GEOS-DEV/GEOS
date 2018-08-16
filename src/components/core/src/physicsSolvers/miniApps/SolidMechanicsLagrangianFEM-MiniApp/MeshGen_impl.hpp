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

#ifndef __MESH_GEN_HPP__
#define __MESH_GEN_HPP__

#include "MatrixMath_impl.hpp"
#include "RAJA/RAJA.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

#define VX(id,i) VX[i + 3*id]
#define elemToNodes(k, i) elemToNodes[i + 8*k]

int myrandom (int i) { return std::rand()%i;}

using RAJA::Index_type;

//Creates a domain on the bi-unit cube, i.e.  [-1, 1]^3

//Input: 
// Kx - number of elements in a cartesian dimension

//Populates:
// VX with physical vertices
// elemToNodes with the element to nodes list
void meshGen(real64 * VX, Index_type * elemToNodes,
             Index_type * RAJA_RESTRICT constitutiveMap, Index_type Kx){             

  std::srand ( unsigned ( std::time(0) ) );
  size_t nx = Kx+1;
  size_t K = Kx*Kx*Kx; //total no of elements
  
  //Create list of nodes
  for(unsigned int z=0; z<nx; ++z){
    for(unsigned int y=0; y<nx; ++y){
      for(unsigned int x=0; x<nx; ++x){
        
        unsigned int id = x + nx*(y + z*nx);
        VX(id,0) = x + (2.0*x)/Kx;
        VX(id,1) = y + (2.0*y)/Kx;
        VX(id,2) = z + (2.0*z)/Kx;
        
      }
    }
  }
  
  //shuffle array
  std::vector<int> rangeSeg(K);
  for(size_t i=0; i<K; ++i) rangeSeg[i] = i;
#if !defined(STRUCTURED_GRID)           
 std::random_shuffle(rangeSeg.begin(),rangeSeg.end(),myrandom);
#endif

  //-----------
  //Generate element to nodes array
  size_t iter=0; 
  for(unsigned int z=0; z<Kx; ++z){
    for(unsigned int y=0; y<Kx; ++y){
      for(unsigned int x=0; x<Kx; ++x){

        size_t tid = rangeSeg[iter];
        elemToNodes(tid, 0) = x + nx*(y + z*nx);
        elemToNodes(tid, 1) = (x+1) + nx*(y + z*nx);
        elemToNodes(tid, 2) = (x+1) + nx*((y+1) + z*nx);
        elemToNodes(tid, 3) = x + nx*((y+1) + z*nx);

        elemToNodes(tid, 4) = x + nx*(y + (z+1)*nx);
        elemToNodes(tid, 5) = (x+1) + nx*(y + (z+1)*nx);
        elemToNodes(tid, 6) = (x+1) + nx*((y+1) + (z+1)*nx);
        elemToNodes(tid, 7) = x + nx*((y+1) + (z+1)*nx);

        for(int q=0; q<8; ++q){
          Index_type id = q + 8*tid; 
          constitutiveMap[id] = id;
        }
        
        iter++;
      }
    }
  }
  
  
}

#endif
