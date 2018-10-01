#ifndef __SHAPE_FUN_HPP__
#define __SHAPE_FUN_HPP__
#include "RAJA/RAJA.hpp"
#include "Layout.hpp"
#include "MatrixMath_impl.hpp"
//
//GEOSX shape function derivatives are defined on [0, 1]
// shape functions { 1-x, x}
// a in {0,1}, b in {+1,-1}
// x in [0,1]
//
RAJA_HOST_DEVICE
RAJA_INLINE real64 shapefun_prime(localIndex ix, real64 ay, real64 by, real64 y,
                          real64 az, real64 bz, real64 z)
{  
  return ix*(ay-by*y)*(az-bz*z);
}


//
//P is a tensor which holds the values of the shape function derivatives
//at the quadrature points w.r.t to the parent coordinate system
//
RAJA_HOST_DEVICE
RAJA_INLINE void generateP(P_Wrapper &P, localIndex nQpts, localIndex nDofs)
{

  localIndex sign[2] = {-1,1};
  real64 a_val[2] = {1.0,0.0};
  real64 b_val[2] = {1,-1};
  real64 intPoints[2]  = {0.211324865405187, 0.788675134594813};
  
  for(localIndex a=0; a < nQpts; ++a)
    {
      //Index for quadrature point (8 different combinations)
      localIndex ix = a % 2;
      localIndex iy = (a / 2) % 2;
      localIndex iz = a / (2 * 2);
      //dNi/dxi ..... 
      //dNi/nu  .....
      //dNi/mu  .....
      //Evaluate point across all basis functions
      for(localIndex q=0; q < 8; ++q)
        {
          localIndex qx = q % 2;
          localIndex qy = (q / 2) % 2;
          localIndex qz = q / (2 * 2);
                    
          P(a,0,q) = shapefun_prime(sign[qx],a_val[qy],b_val[qy],intPoints[iy],
                                    a_val[qz],b_val[qz],intPoints[iz]);
          P(a,1,q) = shapefun_prime(sign[qy],a_val[qx],b_val[qx],intPoints[ix],
                                    a_val[qz],b_val[qz],intPoints[iz]);
          P(a,2,q) = shapefun_prime(sign[qz],a_val[qx],b_val[qx],intPoints[ix],
                                    a_val[qy],b_val[qy],intPoints[iy]);
        }      
    }//loop over quadrature points
      
}

//
//Computes shape function derivatives for a specific element
//
template<typename T, typename U>
RAJA_HOST_DEVICE //x,y,z, format
RAJA_INLINE void make_dNdX(T nodeList, const real64 * X,
                           U dNdX[QUAD_PTS][NUM_NODES][LOCAL_DIM],
                           P_Wrapper P, 
                           localIndex nQpts, localIndex nDofs)
{
    real64 xyz_loc[24];
  for(localIndex i=0; i<nDofs; ++i)
    {
      xyz_loc[0 + LOCAL_DIM*i] = X[0 + LOCAL_DIM * nodeList[i]];
      xyz_loc[1 + LOCAL_DIM*i] = X[1 + LOCAL_DIM * nodeList[i]];
      xyz_loc[2 + LOCAL_DIM*i] = X[2 + LOCAL_DIM * nodeList[i]];            
    }
    
  //loop over quadrature points
  for(localIndex a=0; a < nQpts; ++a)
    {
      
      //dNi/dxi ..... 
      //dNi/nu  .....
      //dNi/mu  .....
      //Evaluate point across all basis functions
            
      //--------[Form the Jacobian]---------------------
      real64 J[LOCAL_DIM][LOCAL_DIM];
      real64 Jinv[LOCAL_DIM][LOCAL_DIM];
      
      //Form Jacobian Matrix
      for(int row=0; row<3; ++row){
        real64 dot[3] = {0.0,0.0,0.0};
        for(int col=0; col<8; ++col){
          dot[0] += P(a,row,col)*xyz_loc[0 + LOCAL_DIM*col];
          dot[1] += P(a,row,col)*xyz_loc[1 + LOCAL_DIM*col];
          dot[2] += P(a,row,col)*xyz_loc[2 + LOCAL_DIM*col];
        }
        J[row][0] = dot[0];
        J[row][1] = dot[1];
        J[row][2] = dot[2];
      }
      
      //Compute Jacobian inverse
      Finverse(J,Jinv);
        
      //Matrix - vector multiplication across the evaluation of
      //the points at the different basis functions 
      for(int q=0; q < 8; ++q){          
        real64 dX[3] = {0.0,0.0,0.0};                  
        for(int row=0; row<3; ++row){                    
          for(int col=0; col<3; ++col){
            dX[row] += Jinv[row][col]*P(a,col,q);
          }
        }

        //write out to data structure
        dNdX[a][q][0] = dX[0];
        dNdX[a][q][1] = dX[1];
        dNdX[a][q][2] = dX[2];
      }
      
    }//loop over quadrature points
  
}

//
//Computes shape function derivatives for a specific element
//
template<typename T, typename U>
RAJA_HOST_DEVICE //x,y,z, format
RAJA_INLINE void make_dNdX(T nodeList, const real64 * X,
                           U dNdX[QUAD_PTS][NUM_NODES][LOCAL_DIM], localIndex nQpts, localIndex nDofs)
{
  
  localIndex sign[2] = {-1,1};
  real64 a_val[2] = {1.0,0.0};
  real64 b_val[2] = {1,-1};
  real64 intPoints[2]  = {0.211324865405187, 0.788675134594813};

  real64 xyz_loc[24];
  for(localIndex i=0; i<nDofs; ++i)
    {
      xyz_loc[0 + LOCAL_DIM*i] = X[0 + LOCAL_DIM * nodeList[i]];
      xyz_loc[1 + LOCAL_DIM*i] = X[1 + LOCAL_DIM * nodeList[i]];
      xyz_loc[2 + LOCAL_DIM*i] = X[2 + LOCAL_DIM * nodeList[i]];            
    }
    
  //loop over quadrature points
  for(localIndex a=0; a < nQpts; ++a)
    {
      //Index for quadrature point (8 different combinations)
      localIndex ix = a % 2;
      localIndex iy = (a / 2) % 2;
      localIndex iz = a / (2 * 2);
      
      real64 P[3][8];   //store evaluation of points at all basis functions
      //dNi/dxi ..... 
      //dNi/nu  .....
      //dNi/mu  .....
      //Evaluate point across all basis functions
      
      //Evaluate point across all basis functions
      for(localIndex q=0; q < 8; ++q)
        {
          localIndex qx = q % 2;
          localIndex qy = (q / 2) % 2;
          localIndex qz = q / (2 * 2);
                    
          P[0][q] = shapefun_prime(sign[qx],a_val[qy],b_val[qy],intPoints[iy],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
          P[1][q] = shapefun_prime(sign[qy],a_val[qx],b_val[qx],intPoints[ix],
                           a_val[qz],b_val[qz],intPoints[iz]);                    
          P[2][q] = shapefun_prime(sign[qz],a_val[qx],b_val[qx],intPoints[ix],
                           a_val[qy],b_val[qy],intPoints[iy]);                    
        }
      
      //--------[Form the Jacobian]---------------------
      real64 J[LOCAL_DIM][LOCAL_DIM];
      real64 Jinv[LOCAL_DIM][LOCAL_DIM];
      
      //Form Jacobian Matrix
      for(int row=0; row<3; ++row){
        real64 dot[3] = {0.0,0.0,0.0};
        for(int col=0; col<8; ++col){
          dot[0] += P[row][col]*xyz_loc[0 + LOCAL_DIM*col];
          dot[1] += P[row][col]*xyz_loc[1 + LOCAL_DIM*col];
          dot[2] += P[row][col]*xyz_loc[2 + LOCAL_DIM*col];
        }
        J[row][0] = dot[0];
        J[row][1] = dot[1];
        J[row][2] = dot[2];
      }
      
      //Compute Jacobian inverse
      Finverse(J,Jinv);

        
      //Matrix - vector multiplication across the evaluation of
      //the points at the different basis functions 
      for(int q=0; q < 8; ++q){          
        real64 dX[3] = {0.0,0.0,0.0};                  
        for(int row=0; row<3; ++row){                    
          for(int col=0; col<3; ++col){
            dX[row] += Jinv[row][col]*P[col][q];
          }
        }

        //write out to data structure
        dNdX[a][q][0] = dX[0];
        dNdX[a][q][1] = dX[1];
        dNdX[a][q][2] = dX[2];
      }
      
    }//loop over quadrature points
  
}

//
//Computes shape function derivatives for a specific element
//
template<typename T, typename U>
RAJA_HOST_DEVICE
RAJA_INLINE void make_dNdX(T nodeList, const real64 * X, 
                           U dNdX_x[QUAD_PTS][NUM_NODES],
                           U dNdX_y[QUAD_PTS][NUM_NODES],
                           U dNdX_z[QUAD_PTS][NUM_NODES],
                           localIndex nQpts, localIndex nDofs)
{

  localIndex sign[2] = {-1,1};
  real64 a_val[2] = {1.0,0.0};
  real64 b_val[2] = {1,-1};
  real64 intPoints[2]  = {0.211324865405187, 0.788675134594813};

  real64 x_loc[8];
  real64 y_loc[8];
  real64 z_loc[8];

  for(localIndex i=0; i<nDofs; ++i)
    {
      x_loc[i]  = X[0 + LOCAL_DIM * nodeList[i]];
      y_loc[i]  = X[1 + LOCAL_DIM * nodeList[i]];
      z_loc[i]  = X[2 + LOCAL_DIM * nodeList[i]];      
    }
  
  //loop over quadrature points
  for(localIndex a=0; a < nQpts; ++a)
    {
      //Index for quadrature point (8 different combinations)
      localIndex ix = a % 2;
      localIndex iy = (a / 2) % 2;
      localIndex iz = a / (2 * 2);
      
      real64 P[3][8];   //store evaluation of points at all basis functions
      //dNi/dxi ..... 
      //dNi/nu  .....
      //dNi/mu  .....
      //Evaluate point across all basis functions
      
      //Evaluate point across all basis functions
      for(localIndex q=0; q < 8; ++q)
        {
          localIndex qx = q % 2;
          localIndex qy = (q / 2) % 2;
          localIndex qz = q / (2 * 2);
          
          P[0][q] = shapefun_prime(sign[qx],a_val[qy],b_val[qy],intPoints[iy],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
          P[1][q] = shapefun_prime(sign[qy],a_val[qx],b_val[qx],intPoints[ix],
                           a_val[qz],b_val[qz],intPoints[iz]);                    
          P[2][q] = shapefun_prime(sign[qz],a_val[qx],b_val[qx],intPoints[ix],
                           a_val[qy],b_val[qy],intPoints[iy]);                    
        }
      
      //--------[Form the inverse of the Jacobian]---------
      real64 J[LOCAL_DIM][LOCAL_DIM];
      real64 Jinv[LOCAL_DIM][LOCAL_DIM];
      
      //Form Jacobian Matrix
      for(int row=0; row<3; ++row){
        real64 dot[3] = {0.0,0.0,0.0};
        for(int col=0; col<8; ++col){
          dot[0] += P[row][col]*x_loc[col];
          dot[1] += P[row][col]*y_loc[col];
          dot[2] += P[row][col]*z_loc[col];
        }
        J[row][0] = dot[0];
        J[row][1] = dot[1];
        J[row][2] = dot[2];
      }
      
      //Compute Jacobian inverse
      Finverse(J,Jinv);
      //---------[Done forming the inverse of the Jacobian]----------
        
      //Matrix - vector multiplication across the evaluation of
      //the points at the different basis functions 
      for(int q=0; q < 8; ++q){          
        real64 dX[3] = {0.0,0.0,0.0};                  
        for(int row=0; row<3; ++row){                    
          for(int col=0; col<3; ++col){
            dX[row] += Jinv[row][col]*P[col][q];
          }
        }

        //write out to data structure
        dNdX_x[a][q] = dX[0];
        dNdX_y[a][q] = dX[1];
        dNdX_z[a][q] = dX[2];
      }
      
    }//loop over quadrature points
  
}

template<typename T, typename U>
RAJA_HOST_DEVICE
RAJA_INLINE void make_dNdX(T nodeList, const real64 * X, 
                           U dNdX_x[QUAD_PTS][NUM_NODES],
                           U dNdX_y[QUAD_PTS][NUM_NODES],
                           U dNdX_z[QUAD_PTS][NUM_NODES],
                           P_Wrapper P,
                           localIndex nQpts, localIndex nDofs)
{

  real64 x_loc[8];
  real64 y_loc[8];
  real64 z_loc[8];

  for(localIndex i=0; i<nDofs; ++i)
    {
      x_loc[i]  = X[0 + LOCAL_DIM * nodeList[i]];
      y_loc[i]  = X[1 + LOCAL_DIM * nodeList[i]];
      z_loc[i]  = X[2 + LOCAL_DIM * nodeList[i]];      
    }
  
  //loop over quadrature points
  for(localIndex a=0; a < nQpts; ++a)
    {
      //dNi/dxi ..... 
      //dNi/nu  .....
      //dNi/mu  .....
      //Evaluate point across all basis functions
      
      //--------[Form the inverse of the Jacobian]---------
      real64 J[LOCAL_DIM][LOCAL_DIM];
      real64 Jinv[LOCAL_DIM][LOCAL_DIM];
      
      //Form Jacobian Matrix
      for(int row=0; row<3; ++row){
        real64 dot[3] = {0.0,0.0,0.0};
        for(int col=0; col<8; ++col){
          dot[0] += P(a,row,col)*x_loc[col];
          dot[1] += P(a,row,col)*y_loc[col];
          dot[2] += P(a,row,col)*z_loc[col];
        }
        J[row][0] = dot[0];
        J[row][1] = dot[1];
        J[row][2] = dot[2];
      }
      
      //Compute Jacobian inverse
      Finverse(J,Jinv);
      //---------[Done forming the inverse of the Jacobian]----------
        
      //Matrix - vector multiplication across the evaluation of
      //the points at the different basis functions 
      for(int q=0; q < 8; ++q){          
        real64 dX[3] = {0.0,0.0,0.0};                  
        for(int row=0; row<3; ++row){                    
          for(int col=0; col<3; ++col){
            dX[row] += Jinv[row][col]*P(a,col,q);
          }
        }

        //write out to data structure
        dNdX_x[a][q] = dX[0];
        dNdX_y[a][q]= dX[1];
        dNdX_z[a][q] = dX[2];
      }
      
    }//loop over quadrature points
  
}

RAJA_HOST_DEVICE
RAJA_INLINE void make_dNdX(geosxData dNdX, const real64 * X,  const localIndex * elemsToNodes,
                           localIndex NoElem, localIndex nQpts, localIndex nDofs)
{

  localIndex sign[2] = {-1,1};
  real64 a_val[2] = {1.0,0.0};
  real64 b_val[2] = {1,-1};
  real64 intPoints[2]  = {0.211324865405187, 0.788675134594813};

  for(localIndex k=0; k < NoElem; ++k){

    real64 x_loc[8];
    real64 y_loc[8];
    real64 z_loc[8];
    
    //Copy local dofs
    for(localIndex a=0; a<nDofs; ++a)
      {
        localIndex id = elemsToNodes[a + NUM_NODES*k];
        x_loc[a] = X[0 + LOCAL_DIM * id];
        y_loc[a] = X[1 + LOCAL_DIM * id];
        z_loc[a] = X[2 + LOCAL_DIM * id]; 
      }

    //loop over quadrature points
    for(localIndex a=0; a < nQpts; ++a)
      {
        //Index for quadrature point (8 different combinations)
        localIndex ix = a % 2;
        localIndex iy = (a / 2) % 2;
        localIndex iz = a / (2 * 2);

        real64 P[3][8];   //store evaluation of points at all basis functions
        //dNi/dxi ..... 
        //dNi/nu  .....
        //dNi/mu  .....
        //Evaluate point across all basis functions
        
        //Evaluate point across all basis functions
        for(localIndex q=0; q < 8; ++q)
          {
            localIndex qx = q % 2;
            localIndex qy = (q / 2) % 2;
            localIndex qz = q / (2 * 2);
            
            P[0][q] = shapefun_prime(sign[qx],a_val[qy],b_val[qy],intPoints[iy],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
            P[1][q] = shapefun_prime(sign[qy],a_val[qx],b_val[qx],intPoints[ix],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
            P[2][q] = shapefun_prime(sign[qz],a_val[qx],b_val[qx],intPoints[ix],
                             a_val[qy],b_val[qy],intPoints[iy]);                    
          }

        //--------[Form the inverse of the Jacobian]---------
        real64 J[LOCAL_DIM][LOCAL_DIM];
        real64 Jinv[LOCAL_DIM][LOCAL_DIM];

        //Form Jacobian Matrix
        for(int row=0; row<3; ++row){
          real64 dot[3] = {0.0,0.0,0.0};
          for(int col=0; col<8; ++col){
            dot[0] += P[row][col]*x_loc[col];
            dot[1] += P[row][col]*y_loc[col];
            dot[2] += P[row][col]*z_loc[col];
          }
          J[row][0] = dot[0];
          J[row][1] = dot[1];
          J[row][2] = dot[2];
        }
        
        //Compute Jacobian inverse
        Finverse(J,Jinv);
        //---------[Done forming the inverse of the Jacobian]----------
        
        //Matrix - vector multiplication across the evaluation of
        //the points at the different basis functions 
        for(int q=0; q < 8; ++q){          
          real64 dX[3] = {0.0,0.0,0.0};                  
          for(int row=0; row<3; ++row){                    
            for(int col=0; col<3; ++col){
              dX[row] += Jinv[row][col]*P[col][q];
            }
          }
          
          localIndex id = q + QUAD_PTS*(a + NUM_NODES*k);
          dNdX[0 + LOCAL_DIM*id] = dX[0];
          dNdX[1 + LOCAL_DIM*id] = dX[1];
          dNdX[2 + LOCAL_DIM*id] = dX[2];
        }
        
      }//loop over quadrature points
        
  } //loop over element list 
  

}

RAJA_HOST_DEVICE
RAJA_INLINE void make_dNdX(geosxData dNdX_x, geosxData dNdX_y,  geosxData dNdX_z,
                           const real64 * X, const localIndex * elemsToNodes,
                           localIndex NoElem, localIndex nQpts, localIndex nDofs)
{

  localIndex sign[2] = {-1,1};
  real64 a_val[2] = {1.0,0.0};
  real64 b_val[2] = {1,-1};
  real64 intPoints[2]  = {0.211324865405187, 0.788675134594813};

  for(localIndex k=0; k < NoElem; ++k){

    real64 x_loc[8];
    real64 y_loc[8];
    real64 z_loc[8];
    
    //Copy local dofs
    for(localIndex a=0; a<nDofs; ++a)
      {
        localIndex id = elemsToNodes[a + NUM_NODES*k];
        x_loc[a] = X[0 + LOCAL_DIM * id];
        y_loc[a] = X[1 + LOCAL_DIM * id];
        z_loc[a] = X[2 + LOCAL_DIM * id]; 
      }

    //loop over quadrature points
    for(localIndex a=0; a < nQpts; ++a)
      {
        //Index for quadrature point (8 different combinations)
        localIndex ix = a % 2;
        localIndex iy = (a / 2) % 2;
        localIndex iz = a / (2 * 2);

        real64 P[3][8];   //store evaluation of points at all basis functions
        //dNi/dxi ..... 
        //dNi/nu  .....
        //dNi/mu  .....
        //Evaluate point across all basis functions       
        for(localIndex q=0; q < 8; ++q)
          {
            localIndex qx = q % 2;
            localIndex qy = (q / 2) % 2;
            localIndex qz = q / (2 * 2);
            
            P[0][q] = shapefun_prime(sign[qx],a_val[qy],b_val[qy],intPoints[iy],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
            P[1][q] = shapefun_prime(sign[qy],a_val[qx],b_val[qx],intPoints[ix],
                             a_val[qz],b_val[qz],intPoints[iz]);                    
            P[2][q] = shapefun_prime(sign[qz],a_val[qx],b_val[qx],intPoints[ix],
                             a_val[qy],b_val[qy],intPoints[iy]);                    
          }

        //--------[Form the inverse of the Jacobian]---------
        real64 J[LOCAL_DIM][LOCAL_DIM];
        real64 Jinv[LOCAL_DIM][LOCAL_DIM];

        //Form Jacobian Matrix
        for(int row=0; row<3; ++row){
          real64 dot[3] = {0.0,0.0,0.0};
          for(int col=0; col<8; ++col){
            dot[0] += P[row][col]*x_loc[col];
            dot[1] += P[row][col]*y_loc[col];
            dot[2] += P[row][col]*z_loc[col];
          }
          J[row][0] = dot[0];
          J[row][1] = dot[1];
          J[row][2] = dot[2];
        }
        
        //Compute Jacobian inverse
        Finverse(J,Jinv);
        //---------[Done forming the inverse of the Jacobian]----------
        
        //Matrix - vector multiplication across the evaluation of
        //the points at the different basis functions 
        for(int q=0; q < 8; ++q){          
          real64 dX[3] = {0.0,0.0,0.0};                  
          for(int row=0; row<3; ++row){                    
            for(int col=0; col<3; ++col){
              dX[row] += Jinv[row][col]*P[col][q];
            }
          }
          
          localIndex id = q + QUAD_PTS*(a + NUM_NODES*k);
          dNdX_x[id] = dX[0];
          dNdX_y[id] = dX[1];
          dNdX_z[id] = dX[2];
        }
        
      }//loop over quadrature points
        
  } //loop over element list 

  
}
#endif
