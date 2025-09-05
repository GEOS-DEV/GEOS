/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Pk_TriangleFace_BCDBB.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKTRIANGLEFACEBCDBB_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKTRIANGLEFACEBCDBB_HPP_

#include "FiniteElementBase.hpp"
#include "BB_Tetrahedron.hpp"
#include "Pk_Pyramid_BCD.hpp"
#include <utility>



namespace geos
{
namespace finiteElement
{


/**
 * This class is the basis class for the coupling of pyramid finite element cells with tetrahedron finite element cells
 * It will contains all the matrices on the triangle faces and use both function to compute shape functions and 
 * their derivatives from Pk_TriangleFace_BCDBB and BB_Tetrahedron.
 */
template< int ORDER >
class Pk_TriangleFace_BCDBB final : public FiniteElementBase
{
public:

  /// The order of the finite element.
  static constexpr int order = ORDER;
  static constexpr int numNodesTet = BB_Tetrahedron< ORDER >::numNodes;
  static constexpr int numNodesPyr = Pk_Pyramid_BCD< ORDER >::numNodes;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~Pk_TriangleFace_BCDBB() = default;


  /*
    * @brief Compute the quadrature points and weights on a triangle defined by its three vertices
    * @param[in] X The coordinates of the triangle vertices
    * @param[out] coords The coordinates of the quadrature points
    * @param[out] weights The weights associated to the quadrature points
  */
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computeQuadraturePoints( real64 const (&X)[3][3], real64 (&coords)[9][3], real64 weights[9] )
  {


    //Gauss-Legendre points and weights on [0,1]
    real64 GLPts[3] = {0.1127016654, 0.5, 0.8872983346};
    real64 GlWeights[3] = {0.2777777778, 0.4444444444, 0.2777777778};

    // Compute the area of the triangle defined by the three vertices
    real64 v1[3] = {X[1][0]-X[0][0], X[1][1]-X[0][1], X[1][2]-X[0][2]};
    real64 v2[3] = {X[2][0]-X[0][0], X[2][1]-X[0][1], X[2][2]-X[0][2]};
    real64 cross[3] = {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]};
    real64 area = 0.5*LvArray::math::sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);  

    // Compute the quadrature points and weights on a triangle in the physical space.
    // To do so, we first put the 1D point on the reference triangle thanks to the Duffy transformation
    // then we put them on the real triangle using barycentrique coordinates

    localIndex count = 0;
    for(localIndex i = 0 ; i<3 ; ++i)
    {
      for(localIndex j = 0 ; j<3 ; ++j)
      {        
        real64 r = GLPts[i];
        real64 s = GLPts[j];
  
        //Duffy transformation + add (1-z) weight
        real64 xi = r;
        real64 eta = s*(1.0-r);

        weights[count] = GlWeights[i]*GlWeights[j]*(1.0-r)*(2.0*area);

        //Barycentric coordinates
        real64 l1 = 1.0 - xi - eta;
        real64 l2 = xi;
        real64 l3 = eta;

        //Quadratures points in the physical triangle

        coords[count][0] = l1*X[0][0] + l2*X[1][0] + l3*X[2][0];
        coords[count][1] = l1*X[0][1] + l2*X[1][1] + l3*X[2][1];
        coords[count][2] = l1*X[0][2] + l2*X[1][2] + l3*X[2][2];


        ++count;
      }

    }

  }

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void transformBarycentricToPhysical( real64 const (&X)[4][3], real64 const (&gradBary)[numNodesTet][4], real64 (&gradPhys)[numNodesTet][3] )
  {
    real64 J[3][3];
    for(localIndex k=0;k<3;++k)
    {
      J[k][0] = X[1][k] - X[0][k];
      J[k][1] = X[2][k] - X[0][k];
      J[k][2] = X[3][k] - X[0][k];
    }

    real64 const detJ = LvArray::tensorOps::determinant< 3 >( J );
    // Inverse of the transpose of J
    real64 invT[3][3];
    invT[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
    invT[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1])/detJ;
    invT[0][2] =  (J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;

    invT[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0])/detJ;
    invT[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
    invT[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0])/detJ;

    invT[2][0] =  (J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
    invT[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0])/detJ;
    invT[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;

    // Barycentric to physical
    // gradNtet_xyz = invT * gradNtet_bary
    //Only three component are independant 
    for(localIndex i=0; i<numNodesTet; ++i)
    {
      // Seules 3 composantes sont indépendantes
      real64 grad_lambda[3];
      grad_lambda[0] = gradBary[i][1];
      grad_lambda[1] = gradBary[i][2];
      grad_lambda[2] = gradBary[i][3];
      // Produit invT * grad_lambda pour obtenir le gradient physique
      for(localIndex k=0; k<3; ++k)
      {
        gradPhys[i][k] = 0.0;
        for(localIndex m=0; m<3; ++m)
        {
          gradPhys[i][k] += invT[k][m] * grad_lambda[m];
        }
      }
    }


  }


  /*
    * @brief Compute the penalization terms for the coupling between a tetrahedron and a pyramid on a triangle face
    * @param[in] xTet The coordinates of the tetrahedron vertices
    * @param[in] X The coordinates of the triangle face vertices
    * @param[in] VDM_inv The inverse of the Vandermonde matrix for the pyramid element
    * @param[in] func The function to call with the computed penalization terms
  */
  template <typename FUNC >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computePenalizationTermsTetPyr( real64 const (&xTet)[4][3], real64 const (&X)[3][3], arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    //Compute the quadrature points and weights on the triangle face
    real64 coords[9][3] = {{0.0}};
    real64 weights[9] = {0.0};
    computeQuadraturePoints(X, coords, weights);

    real64 Ntet[numNodesTet] = {0.0};
    real64 Npyr[numNodesPyr] = {0.0};

    for (localIndex i = 0; i < numNodesTet; ++i)
    {
      for (localIndex j = 0; j < numNodesPyr; ++j)
      {
        real64 val = 0.0;
        for (localIndex q = 0; q < 9; ++q)
        {
          //Compute the shape functions of the tetrahedron at the quadrature point
          BB_Tetrahedron<ORDER>::calcN( xTet, coords[q], Ntet );

          //Normalize Ntet to ensure partition of unity
          real64 sumNtet = 0.0;
          for(localIndex k=0; k<numNodesTet; k++)
          { 
            sumNtet += Ntet[k];
          }
          for(localIndex k=0; k<numNodesTet; k++) 
          {
            Ntet[k] /= sumNtet;
          }

          //Compute the shape functions of the pyramid at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcN( coords[q], VDM_inv, Npyr );

          val += weights[q] * Ntet[i] * Npyr[j];
        }
        func( i, j, val );
      }
    }
  

  }

  /*
    * @brief Compute the penalization terms for the coupling between two pyramids on a triangle face
    * @param[in] X The coordinates of the triangle face vertices
    * @param[in] VDM_inv The inverse of the Vandermonde matrix for the pyramid element
    * @param[in] func The function to call with the computed penalization terms
  */
  template <typename FUNC >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computePenalizationTermsPyrPyr( real64 const (&X)[3][3], arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    //Compute the quadrature points and weights on the triangle face
    real64 coords[9][3] = {{0.0}};
    real64 weights[9] = {0.0};
    computeQuadraturePoints(X, coords, weights);

    real64 Npyr[numNodesPyr] = {0.0};

    for (localIndex i = 0; i < numNodesPyr; ++i)
    {
      for (localIndex j = 0; j < numNodesPyr; ++j)
      {
        real64 val = 0.0;
        for (localIndex q = 0; q < 9; ++q)
        {
          //Compute the shape functions of the pyramid at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcN( coords[q], VDM_inv, Npyr );

          val += weights[q] * Npyr[i] * Npyr[j];
        }
        func( i, j, val );
      }
    }
  

  }

  /*
    * @brief Compute the flux terms for the coupling between a tetrahedron and a pyramid on a triangle face (derivatives on tetrahedron side)
    * @param[in] xTet The coordinates of the tetrahedron vertices
    * @param[in] X The coordinates of the triangle face vertices
    * @param[in] VDM_inv The inverse of the Vandermonde matrix for the pyramid element
    * @param[in] func The function to call with the computed flux terms
  */
  template <typename FUNC >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computeFluxTermsTetPyr( real64 const (&xTet)[4][3], real64 const (&X)[3][3], arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    //Compute the quadrature points and weights on the triangle face
    real64 coords[9][3] = {{0.0}};
    real64 weights[9] = {0.0};
    computeQuadraturePoints(X, coords, weights);

    real64 Ntet[numNodesTet] = {0.0};
    real64 Npyr[numNodesPyr] = {0.0};
    real64 gradNtet[numNodesTet][4] = {{0.0}};
    real64 gradPhys[numNodesTet][3] = {{0.0}};

    for (localIndex i = 0; i < numNodesTet; ++i)
    {
      for (localIndex j = 0; j < numNodesPyr; ++j)
      {
        real64 val[3] = {0.0};
        for (localIndex q = 0; q < 9; ++q)
        {
          //Compute the shape functions of the tetrahedron at the quadrature point
          BB_Tetrahedron<ORDER>::calcNandGradN( xTet, coords[q], Ntet, gradNtet );

          transformBarycentricToPhysical( xTet, gradNtet, gradPhys );

          //Compute the shape functions of the pyramid at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcN( coords[q], VDM_inv, Npyr );

          for(localIndex k=0; k<3; ++k)
          {
            val[k] += weights[q] * gradPhys[i][k] * Npyr[j];
          }
        }
        func( i, j, val );
      }
    }
  

  }

  /*
   * @brief Compute the flux terms for the coupling between a tetrahedron and a pyramid on a triangle face (derivatives on pyramid side)
   * @param[in] xTet The coordinates of the tetrahedron vertices
   * @param[in] X The coordinates of the triangle face vertices
   * @param[in] VDM_inv The inverse of the Vandermonde matrix for the pyramid element
   * @param[in] func The function to call with the computed flux terms
  */
  template <typename FUNC >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computeFluxTermsPyrTet( real64 const (&xTet)[4][3], real64 const (&X)[3][3], arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    //Compute the quadrature points and weights on the triangle face
    real64 coords[9][3] = {{0.0}};
    real64 weights[9] = {0.0};
    computeQuadraturePoints(X, coords, weights);

    real64 Ntet[numNodesTet] = {0.0};
    real64 Npyr[numNodesPyr] = {0.0};
    real64 gradNPyr[numNodesPyr][3] = {{0.0}};

    for (localIndex i = 0; i < numNodesTet; ++i)
    {
      for (localIndex j = 0; j < numNodesPyr; ++j)
      {
        real64 val[3] = {0.0};
        for (localIndex q = 0; q < 9; ++q)
        {
          //Compute the shape functions of the tetrahedron at the quadrature point
          BB_Tetrahedron<ORDER>::calcN( xTet, coords[q], Ntet );


          //Compute the shape functions of the pyramid at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcgradN( coords[q], VDM_inv, gradNPyr );

          for(localIndex k=0; k<3; ++k)
          {
            val[k] += weights[q] * Ntet[i] * gradNPyr[j][k]; 
          }
        }
        func( i, j, val );
      }
    }
  

  }

  /*
   * @brief Compute the flux terms for the coupling between two pyramids on a triangle face
   * @param[in] X The coordinates of the triangle face vertices
   * @param[in] VDM_inv The inverse of the Vandermonde matrix for the pyramid element
   * @param[in] func The function to call with the computed flux terms
  */ 
  template <typename FUNC >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static constexpr void computeFluxTermsPyrPyr( real64 const (&X)[3][3], arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    //Compute the quadrature points and weights on the triangle face
    real64 coords[9][3] = {{0.0}};
    real64 weights[9] = {0.0};
    computeQuadraturePoints(X, coords, weights);

    real64 Npyr[numNodesPyr] = {0.0};
    real64 gradNPyr[numNodesPyr][3] = {{0.0}};

    for (localIndex i = 0; i < numNodesPyr; ++i)
    {
      for (localIndex j = 0; j < numNodesPyr; ++j)
      {
        real64 val[3] = {0.0};
        for (localIndex q = 0; q < 9; ++q)
        {
          //Compute the shape functions of the tetrahedron at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcN( coords[q], VDM_inv, Npyr );


          //Compute the shape functions of the pyramid at the quadrature point
          Pk_Pyramid_BCD<ORDER>::calcgradN( coords[q], VDM_inv, gradNPyr );

          for(localIndex k=0; k<3; ++k)
          {
            val[k] += weights[q] * Npyr[i] * gradNPyr[j][k]; 
          }
        }
        func( i, j, val );
      }
    }
  

  }


};

/**
 *  Pyramid element with BCD basis functions of order 1.
 */
using P1_TriangleFace_BCDBB = Pk_TriangleFace_BCDBB< 1 >;


}

}

#endif // GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKPYRAMIDBCD_HPP_
