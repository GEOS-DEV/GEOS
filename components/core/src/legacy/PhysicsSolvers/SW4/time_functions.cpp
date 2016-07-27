#include <cmath>
//#include "Require.h"
//#include <iostream>

double VerySmoothBump(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = - 1024*pow(t*freq,10) + 5120*pow(t*freq,9) - 10240*pow(t*freq,8) + 10240*pow(t*freq,7) - 5120*pow(t*freq,6) + 1024*pow(t*freq,5);
  return tmp;
}

double VerySmoothBump_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

double VerySmoothBump_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

double VerySmoothBump_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_tom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 5120*pow(t*freq,4)*(-20*pow(t*freq,5) + 81*pow(t*freq,4) -
			       128*pow(t*freq,3) + 98*pow(t*freq,2) - 36*t*freq + 5 );
  return tmp;
}

double VerySmoothBump_omom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*t*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) );
  return tmp;
}

double VerySmoothBump_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*t*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) )
+2*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_tttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = pow(freq,4)*122880*
	(-42*pow(t*freq,6)+126*pow(t*freq,5)-140*pow(t*freq,4)+70*pow(t*freq,3)-15*pow(t*freq,2)+t*freq);
  return tmp;
}

double VerySmoothBump_tttom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*61440*
	(-120*pow(t*freq,7)+378*pow(t*freq,6)-448*pow(t*freq,5)+245*pow(t*freq,4)-60*pow(t*freq,3)+5*t*freq*t*freq);
  return tmp;
}

double VerySmoothBump_ttomom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*t*t*122880*
	(-42*pow(t*freq,6)+126*pow(t*freq,5)-140*pow(t*freq,4)+70*pow(t*freq,3)-15*pow(t*freq,2)+t*freq) +
20480*freq*freq*freq*t*t*t*(-153*pow(t*freq,5)+540*pow(t*freq,4)-728*pow(t*freq,3)+462*pow(t*freq,2)-135*t*freq+14);
  return tmp;
}


double RickerWavelet(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return (2*factor - 1)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,2)*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

double RickerWavelet_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*t*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

double RickerWavelet_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*( 6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,4)*t*( -60+80*factor-16*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*(12-108*factor+96*factor*factor-16*factor*factor*factor)*exp(-factor);
 
  else
    return 0;
}

double RickerInt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return -t*exp(-factor);
  else
    return 0;
}

double RickerInt_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return (2*factor-1)*exp(-factor);
  else
    return 0;
}

double RickerInt_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return 2*t*t*t*freq*M_PI*M_PI*exp(-factor);
  else
    return 0;
}

double RickerInt_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*t*(6-4*factor)*exp(-factor);
  else
    return 0;
}

double RickerInt_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*(6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerInt_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return t*M_PI*M_PI*freq*(12-28*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq*freq*freq*t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return (1-2*factor)/ sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

double Gaussian_tom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_omom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*t*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq*t / sqrt(2*M_PI)*(3-2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*(12*factor-3-4*factor*factor)/sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_tttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq / sqrt(2*M_PI)*(3-12*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_tttom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*t / sqrt(2*M_PI)*(15-20*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_ttomom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)*(-6+54*factor-48*factor*factor+8*factor*factor*factor)*exp(-factor);
  else
    return 0;
}

double Erf( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  return 0.5*(1+erf( freq*t/sqrt(2.0)) );
}

double Erf_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Erf_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Erf_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq / sqrt(2*M_PI)* freq*freq*t*exp(-factor);
  else
    return 0;
}

double Erf_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

double Erf_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return t / sqrt(2*M_PI)* freq*freq*(2*factor-3)*exp(-factor);
  else
    return 0;
}

double Ramp(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 1.0;
  else
    tmp = 0.5*(1 - cos(M_PI*t*freq));
  
  return tmp;
}

double Ramp_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*freq*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*t*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*M_PI*freq*freq*cos(M_PI*t*freq);
  
  return tmp;
}

double Ramp_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = -0.5*M_PI*M_PI*M_PI*freq*freq*freq*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = M_PI*M_PI*freq*(cos(M_PI*t*freq)-0.5*M_PI*t*freq*sin(M_PI*t*freq));
  
  return tmp;
}

double Triangle(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*(sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49);

  return tmp; 
}

double Triangle_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*M_PI*freq*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7);

  return tmp; 
}

double Triangle_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*8./pow(M_PI,2)*(M_PI*freq*t*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7) + (sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49) );

  return tmp; 
}

double Triangle_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*8./pow(M_PI,2)*(-M_PI*M_PI*freq*freq)*
	( sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq)) + 
          sin(5*M_PI*(t*freq)) - sin(7*M_PI*(t*freq)) );

  return tmp; 
}

double Triangle_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*8./pow(M_PI,2)*(-M_PI*M_PI*M_PI*freq*freq*freq)*
	( cos(M_PI*(t*freq)) - 3*cos(3*M_PI*(t*freq)) + 
          5*cos(5*M_PI*(t*freq)) - 7*cos(7*M_PI*(t*freq)) );

  return tmp; 
}

double Triangle_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*freq*M_PI*M_PI*8/pow(M_PI,2)*( 
	(-3)*( sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq)) + 
	       sin(5*M_PI*(t*freq)) - sin(7*M_PI*(t*freq)) ) 
	-freq*t*M_PI*(cos(M_PI*(t*freq)) - 3*cos(3*M_PI*(t*freq)) + 
		 5*cos(5*M_PI*(t*freq)) - 7*cos(7*M_PI*(t*freq)) ));
  return tmp; 
}

double Sawtooth(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 8./pow(M_PI,2)*(sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq))/9 + sin(5*M_PI*(2*t*freq))/25 - sin(7*M_PI*(2*t*freq))/49);

  return tmp; 
}

double Sawtooth_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*freq)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

double Sawtooth_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*t)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

double Sawtooth_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(-M_PI*M_PI*2*2*freq*freq)*
              (sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq)) +
	       sin(5*M_PI*(2*t*freq)) - sin(7*M_PI*(2*t*freq)));

  return tmp; 
}

double Sawtooth_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = -64.*(M_PI*freq*freq*freq)*
              (cos(M_PI*(2*t*freq)) - 3*cos(3*M_PI*(2*t*freq)) +
	       5*cos(5*M_PI*(2*t*freq)) - 7*cos(7*M_PI*(2*t*freq)));

  return tmp; 
}

double Sawtooth_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = -64*freq*(sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq)) +
		     sin(5*M_PI*(2*t*freq)) - sin(7*M_PI*(2*t*freq))) 
          -64*M_PI*freq*freq*t*
               (cos(M_PI*(2*t*freq)) - 3*cos(3*M_PI*(2*t*freq)) +
		5*cos(5*M_PI*(2*t*freq)) - 7*cos(7*M_PI*(2*t*freq)));

  return tmp; 
}

double SmoothWave(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = (c0*pow(t*freq,3)+c1*pow(t*freq,4)+c2*pow(t*freq,5)+c3*pow(t*freq,6)+c4*pow(t*freq,7));
  
  return tmp;
}

double SmoothWave_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = freq*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

double SmoothWave_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = t*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

double SmoothWave_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5));
  
  return tmp;
}

double SmoothWave_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

double SmoothWave_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5)) +
         freq*freq*t*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

double Brune( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	return 1-exp(-tf)*(1+tf);
      else
	return 1;
    }
}

double Brune_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*freq*exp(-tf);
      else
	return 0;
    }
}

double Brune_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*t*exp(-tf);
      else
	return 0;
    }
}

double Brune_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

double Brune_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

double Brune_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return freq*(2-4*tf+tf*tf)*exp(-tf);
      else
	return 0;
    }
}

double DBrune( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	return tf*freq*exp(-tf);
      else
	return 0;
    }
}

double DBrune_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

double DBrune_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*(2-tf)*exp(-tf);
      else
	return 0;
    }
}

double DBrune_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

double DBrune_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (3-tf)*freq*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

double DBrune_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (6*tf-6-tf*tf)*freq*freq*exp(-tf);
      else
	return 0;
    }
}

double BruneSmoothed( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
    return 1-exp(-tf)*(1 + tf + 0.5*tf*tf - 1.5*hi*tf*tf*tf + 
		       1.5*hi*hi*tf*tf*tf*tf -0.5*hi*hi*hi*tf*tf*tf*tf*tf);
  else
    {
      if( -tf > par[0] )
	return 1-exp(-tf)*(1+tf);
      else
	return 1;
    }
}

double BruneSmoothed_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*((0.5-3*c3)*tf*tf+(c3-4*c4)*tf*tf*tf+(c4-5*c5)*tf*tf*tf*tf+c5*tf*tf*tf*tf*tf);
  }
  else
    {
      if( -tf > par[0] )
	 return tf*freq*exp(-tf);
      else
	return 0;
    }
}

double BruneSmoothed_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
     return exp(-tf)*t*((0.5-3*c3)*tf*tf+(c3-4*c4)*tf*tf*tf+(c4-5*c5)*tf*tf*tf*tf+c5*tf*tf*tf*tf*tf);
  }
  else
    {
      if( -tf > par[0] )
	 return tf*t*exp(-tf);
      else
	return 0;
    }
}

double BruneSmoothed_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
     return exp(-tf)*( freq*freq*( (1-6*c3)*tf+(-0.5+6*c3-12*c4)*tf*tf+(-c3+8*c4-20*c5)*tf*tf*tf+
				   (-c4+10*c5)*tf*tf*tf*tf -c5*tf*tf*tf*tf*tf));
				
  }
  else
    {
      if( -tf > par[0] )
	return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

double BruneSmoothed_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*freq*freq*( (1-6*c3) + (18*c3-2-24*c4)*tf+(0.5-9*c3+36*c4-60*c5)*tf*tf+(c3-12*c4+60*c5)*tf*tf*tf+
					(c4-15*c5)*tf*tf*tf*tf +c5*tf*tf*tf*tf*tf);
				
  }
  else
    {
      if( -tf > par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

double BruneSmoothed_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*( tf*( (1-6*c3) + (12*c3-1-24*c4)*tf+(-3*c3+24*c4-60*c5)*tf*tf+
				 (-4*c4+40*c5)*tf*tf*tf -5*c5*tf*tf*tf*tf ) + 
			    (2-tf)*((1-6*c3)*tf+(-0.5+6*c3-12*c4)*tf*tf+(-c3+8*c4-20*c5)*tf*tf*tf+
				    (-c4+10*c5)*tf*tf*tf*tf -c5*tf*tf*tf*tf*tf) );
  }
  else
    {
      if( -tf > par[0] )
	 return freq*(2-4*tf+tf*tf)*exp(-tf);
      else
	return 0;
    }
}


double GaussianWindow( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
    return sin(tf)*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (freq*cos(tf)-freq*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

double GaussianWindow_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (t*cos(tf)-t*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

double GaussianWindow_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( (-freq*freq-freq*freq*incyc2+freq*freq*tf*tf*incyc2*incyc2)*sin(tf)-
	      tf*2*freq*freq*incyc2*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*freq*freq*(3*tf*incyc2*(1+incyc2)-tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      freq*freq*freq*( 3*tf*tf*incyc2*incyc2-3*incyc2-1)*cos(tf))*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*(-2-2*incyc2 + 3*incyc2*tf*tf +5*tf*tf*incyc2*incyc2 -tf*tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      t*freq*freq*( 3*tf*tf*incyc2*incyc2-7*incyc2-1)*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double Liu( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 1;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*t-0.7*tau1*ipi*sin(M_PI*t/tau1)-1.2*tau1*ipi*(cos(0.5*M_PI*t/tau1)-1));
      else if( t <= 2*tau1 )
	 return cn*(1.0*t-0.3*tau1+1.2*tau1*ipi - 0.7*tau1*ipi*sin(M_PI*t/tau1)+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3*t+1.1*tau1+1.2*tau1*ipi+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = t*1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)+0.3*M_PI*cos(0.5*M_PI*t/tau1))/tau1;
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)/tau1-0.3*M_PI*sin(M_PI*(t-tau1)/tau2)/tau2);
      else if( t <= tau )
	 return cn*(-0.3*M_PI*sin(M_PI*(t-tau1)/tau2))/tau2;
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI/tau1*t)-0.15*M_PI*M_PI*sin(0.5*M_PI/tau1*t))/(tau1*tau1);
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI*t/tau1)/(tau1*tau1)-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2));
      else if( t <= tau )
	 return cn*(-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2))/(tau2*tau2);
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
      if( t <= tau1 )
	 return cn*(2*(0.7*M_PI*sin(M_PI/tau1*t)+0.3*M_PI*cos(0.5*M_PI/tau1*t)) + 
        (0.7*M_PI*M_PI*t/tau1*cos(M_PI/tau1*t)-0.15*M_PI*M_PI*t/tau1*sin(0.5*M_PI/tau1*t)) )/(tau1);
      else if( t <= 2*tau1 )
	 return cn*(2*(0.7*M_PI*sin(M_PI*t/tau1)/tau1-0.3*M_PI*sin(M_PI*(t-tau1)/tau2)/tau2) + 
		    t*(0.7*M_PI*M_PI*cos(M_PI*t/tau1)/(tau1*tau1)-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2)));
      else if( t <= tau )
	 return -0.3*M_PI*cn*( 2*sin(M_PI*(t-tau1)/tau2)/tau2 + M_PI*t*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2) );
   }
   return 0.; // should never get here, but keeps compiler happy
}

double NullFunc( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
// this function should never be called
//  CHECK_INPUT(false,"The NullFunc time function was called!");
  return 0.;  
}

double Dirac( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   const double c1=2.59765625;
   const double c2=10.0625;
   const double c3=22.875;
   const double c4=26.6666666666667;
   const double c5=11.6666666666667;
   const double o12=0.0833333333333333;
   const double o6=0.166666666666667;
   const double a1=0.41015625;
   const double a2=2.140625;
   const double a3=3.4609375;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt

   double kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   //   std::cout << "t="<< t << " kc=" << kc << " k0= " << k0 << std::endl;
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      double alpha =(-t*freq-k0);
      double alpha2=alpha*alpha;
      double pol=alpha2*alpha2*( c1 - c2*alpha2 + c3*alpha2*alpha2 -
	    c4*alpha2*alpha2*alpha2+c5*alpha2*alpha2*alpha2*alpha2);
      double wgh;
      if( k0 == 2 )
         wgh = o12*alpha*(1-alpha2)-a1*alpha2+pol;
      else if( k0 == 1 )
         wgh = o6*alpha*(-4+alpha2)+a2*alpha2-4*pol;
      else if( k0 == 0 )
         wgh = 1-a3*alpha2 + 6*pol;
      else if( k0 == -1 )
         wgh = o6*alpha*(4-alpha2)+a2*alpha2-4*pol;
      else if( k0 == -2 )
         wgh = o12*alpha*(-1+alpha2)-a1*alpha2+pol;
      //      std::cout << "wgh = " << wgh << std::endl;
      return freq*wgh;
   }
}

double Dirac_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   const double c1=10.390625;
   const double c2= 60.375;
   const double c3=183.0;
   const double c4=266.666666666667;
   const double c5=140;
   const double o12=0.0833333333333333;
   const double o6=0.166666666666667;
   const double a1=0.8203125;
   const double a2=4.281250;
   const double a3=6.921875;
   //   delta'(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   double kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      double wgh;
      double alpha =(-t*freq-k0);
      double alpha2=alpha*alpha;
      double polp=alpha2*alpha*( c1 - c2*alpha2 + c3*alpha2*alpha2 -
	    c4*alpha2*alpha2*alpha2+c5*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = o12*(1-3*alpha2)-a1*alpha + polp;
      else if( k0 == 1 )
         wgh = o6*(-4+3*alpha2)+a2*alpha-4*polp;
      else if( k0 == 0 )
         wgh = -a3*alpha + 6*polp;
      else if( k0 == -1 )
         wgh = o6*(4-3*alpha2)+a2*alpha-4*polp;
      else if( k0 == -2 )
         wgh = o12*(-1+3*alpha2)-a1*alpha+polp;
      return freq*freq*wgh;
   }
}

double Dirac_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   const double c1=31.171875;
   const double c2=301.875;
   const double a1=0.8203125;
   const double a2=4.281250;
   const double a3=6.921875;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   double kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      double wgh;
      double alpha =(-t*freq-k0);
      double alpha2=alpha*alpha;
      double polpp=alpha2*( c1 - c2*alpha2 + 1281.0*alpha2*alpha2 -
	    2400.0*alpha2*alpha2*alpha2+1540*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = -0.5*alpha-a1 + polpp;
      else if( k0 == 1 )
         wgh = alpha + a2 - 4*polpp;
      else if( k0 == 0 )
         wgh = -a3 + 6*polpp;
      else if( k0 == -1 )
         wgh = -alpha+a2-4*polpp;
      else if( k0 == -2 )
         wgh = 0.5*alpha-a1+polpp;
      return freq*freq*freq*wgh;
   }
}
double Dirac_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   const double c1=62.34375;
   const double c2=1207.5;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   double kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      double wgh;
      double alpha =(-t*freq-k0);
      double alpha2=alpha*alpha;
      double polppp=alpha*( c1 - c2*alpha2 + 7686.0*alpha2*alpha2 -
	    19200.0*alpha2*alpha2*alpha2+15400*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = -0.5 + polppp;
      else if( k0 == 1 )
         wgh = 1.0 - 4*polppp;
      else if( k0 == 0 )
         wgh =  6*polppp;
      else if( k0 == -1 )
         wgh = -1.0-4*polppp;
      else if( k0 == -2 )
         wgh = 0.5+polppp;
      return freq*freq*freq*freq*wgh;
   }
}
double Dirac_tttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   const double c1=62.34375;
   const double c2=3622.5;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   double kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      double wgh;
      double alpha =(-t*freq-k0);
      double alpha2=alpha*alpha;
      double polpppp=( c1 - c2*alpha2 + 38430.0*alpha2*alpha2 -
	    134400.0*alpha2*alpha2*alpha2+138600.0*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = polpppp;
      else if( k0 == 1 )
         wgh = -4*polpppp;
      else if( k0 == 0 )
         wgh =  6*polpppp;
      else if( k0 == -1 )
         wgh = -4*polpppp;
      else if( k0 == -2 )
         wgh = polpppp;
      return freq*freq*freq*freq*freq*wgh;
   }
}
double Dirac_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Dirac_tom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Dirac_omom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Dirac_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Dirac_tttom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Dirac_ttomom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}

double Discrete( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
// freq holds 1/dt
   double tstart = par[0];
   int npts = ipar[0];

   int k = static_cast<int>(floor((t-tstart)*freq));

   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }

   double arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
//std::cout <<  "t= " << t << " npts " << npts << " k= " << k << "arg = " << arg <<  std::endl;
   return par[6*k+1] + par[2+6*k]*arg + par[3+6*k]*arg*arg + par[4+6*k]*arg*arg*arg +
       par[5+6*k]*arg*arg*arg*arg + par[6+6*k]*arg*arg*arg*arg*arg; 
}

double Discrete_t( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
// freq holds 1/dt
   double tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   double arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (par[2+6*k] + 2*par[3+6*k]*arg + 3*par[4+6*k]*arg*arg + 4*par[5+6*k]*arg*arg*arg+
	5*par[6+6*k]*arg*arg*arg*arg)*freq;
}

double Discrete_tt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   double arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (2*par[3+6*k] + 6*par[4+6*k]*arg + 12*par[5+6*k]*arg*arg + 20*par[6+6*k]*arg*arg*arg)*freq*freq;
}

double Discrete_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   double arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (6*par[4+6*k] + 24*par[5+6*k]*arg + 60*par[6+6*k]*arg*arg)*freq*freq*freq;
}

double Discrete_tttt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   double tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      t = tstart;
      k = 0;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   double arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (24*par[5+6*k] + 120*par[6+6*k]*arg)*freq*freq*freq*freq;
}

double Discrete_om( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Discrete_tom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Discrete_omom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Discrete_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Discrete_tttom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
double Discrete_ttomom( double freq, double t, double* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}

double C6SmoothBump(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 51480*pow(t*freq*(1-t*freq),7);
  //    tmp = 16384*pow(t*freq*(1-t*freq),7);
  return tmp;
}

double C6SmoothBump_t(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
     tmp = 51480*freq*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
  return tmp;
}

double C6SmoothBump_om(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*t*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
     tmp = 51480*t*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
  return tmp;
}

double C6SmoothBump_tt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
     tmp = 51480*freq*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
  return tmp;
}

double C6SmoothBump_tom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*(7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6)+
     //	          t*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6) ));
     tmp = 51480*(7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6)+
	          t*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6) ));
  return tmp;
}

double C6SmoothBump_omom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*t*t*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
     tmp = 51480*t*t*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
  return tmp;
}

double C6SmoothBump_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*freq*42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //					    5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4));
     tmp = 51480*freq*freq*freq*42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
					    5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4));
  return tmp;
}

double C6SmoothBump_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*(2*freq*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)) ) 
     //            + freq*freq*t*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //					     5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4)) ) );
     tmp = 51480*(2*freq*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)) ) 
            + freq*freq*t*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
					     5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4)) ) );
  return tmp;
}

double C6SmoothBump_tttt(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*freq*freq*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) );
     tmp = 51480*freq*freq*freq*freq*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) );
  return tmp;
}

double C6SmoothBump_tttom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*(3*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //		 5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
     //	    freq*t*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+
     //                    840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) ) );
     tmp = 51480*freq*freq*(3*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
		 5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
	    freq*t*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+
                    840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) ) );
  return tmp;
}

double C6SmoothBump_ttomom(double freq, double t, double* par, int npar, int* ipar, int nipar )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*( 2*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)))+
     //		   4*freq*t*( 42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //				5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
     //	   freq*freq*t*t*( (12*42*(pow(t*freq*(1-t*freq),5)-
     //	   5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) )));
     tmp = 51480*( 2*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)))+
		   4*freq*t*( 42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
				5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
	   freq*freq*t*t*( (12*42*(pow(t*freq*(1-t*freq),5)-
	   5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) )));
  return tmp;
}
