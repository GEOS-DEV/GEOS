#ifndef EW_QSPLINE_H
#define EW_QSPLINE_H

class Qspline
{
   int m_npts;
   double* m_polcof;
   double m_tmin, m_dt, m_dti;
public:
   Qspline( int npts, double* fun, double tmin, double dt, int bclow=1, int bchigh=1,
	    double s1=0, double t1=0, double sn=0, double tn=0 );
   void Qsplineold( int npts, double* fun, double tmin, double dt );
   ~Qspline();
   double* get_polycof_ptr() { return m_polcof;}
   void evalf( double t, double& f );
   //   void evald( double t, double& f, double& f1, double& f2 );
   void evaldd( double t, double& f, double& f1, double& f2, double& f3, double& f4 );
};

#endif
