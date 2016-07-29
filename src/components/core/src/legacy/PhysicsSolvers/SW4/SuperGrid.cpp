#include "SuperGrid.h"
#include <iostream>

using namespace std;

SuperGrid::SuperGrid()
{
  m_left = false;
  m_right = false;
  m_x0=0.;
  m_x1=1.;
  m_width=0.1;
  m_const_width=0.;
  m_epsL = 1e-4;
}

void SuperGrid::print_parameters() const
{
   cout << "SuperGrid parameters left= " << m_left << " right= " << m_right
	<< " x0= " << m_x0 << " x1= " << m_x1 << " width= " << m_width
	<< " transition= " << m_trans_width << " epsL= " << m_epsL << endl;
}

void SuperGrid::define_taper(bool left, realT leftStart, bool right, realT rightEnd, realT width)
{
  m_left = left;
  m_x0 = leftStart;
  m_right = right;
  m_x1 = rightEnd;
  m_width = width;
// always use the full width for the transition, making m_const_width=0
  m_trans_width = width;
//  m_trans_width = transWidth;
  m_const_width = m_width - m_trans_width;
  
// sanity checks
  if (m_left || m_right)
  {
     realT dlen = m_x1-m_x0;
     if( m_width <= 0 )
	cout << "The supergrid taper width must be positive, not = " << m_width << endl;
     if( m_width >= dlen )
	cout << "The supergrid taper width must be smaller than the domain, not = " << m_width << endl;
     if( m_trans_width <= 0 )
        cout <<  "The supergrid taper transition width must be positive, not = " << m_trans_width << endl;
     if( m_const_width < 0 )
	cout << "The supergrid const_width = width - trans_width must be non-negative, not = " << m_const_width << endl;
  }
  
  if (m_left && m_right)
  {
    if (m_x0+m_width > m_x1-m_width)
    {
      print_parameters();
      cout << "The supergrid taper functions at the left and right must be separated. Here x0+width = " << m_x0+m_width << 
	 " and x1-width = " << m_x1-m_width << endl;
    }
    
  }
  else if( m_left )
  {
    if (m_x0+m_width > m_x1 )
    {
      print_parameters();
      cout << "The supergrid taper functions at the left must be smaller than the domain. Here x0+width = " << m_x0+m_width << 
	 " and x1 = " << m_x1 << endl;
    }
  }    
  else if( m_right )
  {
    if (m_x0 > m_x1-m_width )
    {
      print_parameters();
      cout << "The supergrid taper functions at the right must be smaller than the domain. Here x0 = " << m_x0 << 
	 " and x1-width = " << m_x1-m_width << endl;
    }
  }
}

realT SuperGrid::dampingCoeff(realT x) const
{
  realT phi = stretching(x);
// should be equivalent to sigmaScale/(1-(1-epsL)*sigmaScale)
  realT f=(1-phi)/phi/(1-m_epsL);
  return f;
}

realT SuperGrid::sigmaScale(realT x) const
{ 
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  realT f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 < const_width <= x <= const_width+trans_width = m_width
// constant damping in 0 <= x <= const_width
    f=sigma( (m_x0 + m_width - x)/m_trans_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 - const_width < m_x1
// constant damping in m_x1 - const_width <= x <= m_x1
    f=sigma( (x - (m_x1-m_width) )/m_trans_width);
  return f;
}

realT SuperGrid::linTaper(realT x) const
{ 
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  realT f=0.;
  if (m_left && x < m_x0+m_width)
//  linear taper from 0 to 1
    f= (m_x0 + m_width - x)/m_width; 
  else if (m_right && x > m_x1-m_width)
// linear taper from 0 to 1
    f= (x - (m_x1-m_width) )/m_width;
  return f;
}


// used for damping coefficient
realT SuperGrid::sigma(realT xi) const
{
   realT f;
   if (xi<=0.)
      f = 0;
   else if (xi>=1.)
      f = 1.0;
   else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
// C4 function
//    f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
//      126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
// C5 function
      f =  xi*xi*xi*xi*xi*xi*(
    462-1980*xi+3465*xi*xi-3080*xi*xi*xi+1386*xi*xi*xi*xi-252*xi*xi*xi*xi*xi);
   return f;
}

realT SuperGrid::stretching( realT x ) const
{ // this function satisfies 0 < epsL <= f <= 1
   return 1-(1-m_epsL)*sigmaScale(x);
}

realT SuperGrid::cornerTaper( realT x ) const
{ // this function is 1 in the interior and tapers linearly to 1/2 in the SG layers
  const realT cmin=0.33;
  return 1.0 - (1.0-cmin)*linTaper(x);
}
