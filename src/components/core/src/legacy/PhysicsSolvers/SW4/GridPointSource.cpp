#include "GridPointSource.h"

#include <cmath>

using namespace std;

#include "time_functions.h"

//-----------------------------------------------------------------------
GridPointSource::GridPointSource( double frequency, double t0,
                                  int N, int M, int L, int G,
                                  double Fx, double Fy, double Fz,
                                  timeDep tDep,
                                  int ncyc, double* pars, int npar, int* ipars, int nipar ):
  m_i0(N), m_j0(M), m_k0(L),m_grid(G),
  mFreq(frequency),
  mT0(t0),
  mTimeDependence(tDep),
  mNcyc(ncyc)
{
  // Copy only pointers, not memory.
  // mPar, mIpar will no longer be changed by this class, they should
  // be set correctly in Source.
  mPar   = pars;
  mNpar  = npar;
  mIpar  = ipars;
  mNipar = nipar;

  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;

  initializeTimeFunction();
}

//-----------------------------------------------------------------------
GridPointSource::~GridPointSource()
{
  //   if( mPar != NULL )
  //  delete[] mPar;
}

//-----------------------------------------------------------------------
void GridPointSource::initializeTimeFunction()
{
  //   if( mTimeDependence != iDiscrete )
  //      mPar[0] = m_min_exponent;
  switch(mTimeDependence)
  {
  case iRicker:
    mTimeFunc = RickerWavelet;
    mTimeFunc_t = RickerWavelet_t;
    mTimeFunc_tt = RickerWavelet_tt;
    mTimeFunc_ttt = RickerWavelet_ttt;
    mTimeFunc_om = RickerWavelet_om;
    mTimeFunc_omtt = RickerWavelet_omtt;
    break;
  case iGaussian:
    mTimeFunc   = Gaussian;
    mTimeFunc_t = Gaussian_t;
    mTimeFunc_tt = Gaussian_tt;
    mTimeFunc_ttt = Gaussian_ttt;
    mTimeFunc_om = Gaussian_om;
    mTimeFunc_omtt = Gaussian_omtt;
    break;
  case iRamp:
    mTimeFunc = Ramp;
    mTimeFunc_t = Ramp_t;
    mTimeFunc_tt = Ramp_tt;
    mTimeFunc_ttt = Ramp_ttt;
    mTimeFunc_om = Ramp_om;
    mTimeFunc_omtt = Ramp_omtt;
    break;
  case iTriangle:
    mTimeFunc = Triangle;
    mTimeFunc_t = Triangle_t;
    mTimeFunc_tt = Triangle_tt;
    mTimeFunc_ttt = Triangle_ttt;
    mTimeFunc_om = Triangle_om;
    mTimeFunc_omtt = Triangle_omtt;
    break;
  case iSawtooth:
    mTimeFunc = Sawtooth;
    mTimeFunc_t = Sawtooth_t;
    mTimeFunc_tt = Sawtooth_tt;
    mTimeFunc_ttt = Sawtooth_ttt;
    mTimeFunc_om = Sawtooth_om;
    mTimeFunc_omtt = Sawtooth_omtt;
    break;
  case iSmoothWave:
    mTimeFunc = SmoothWave;
    mTimeFunc_t = SmoothWave_t;
    mTimeFunc_tt = SmoothWave_tt;
    mTimeFunc_ttt = SmoothWave_ttt;
    mTimeFunc_om = SmoothWave_om;
    mTimeFunc_omtt = SmoothWave_omtt;
    break;
  case iErf:
    mTimeFunc = Erf;
    mTimeFunc_t = Erf_t;
    mTimeFunc_tt = Erf_tt;
    mTimeFunc_ttt = Erf_ttt;
    mTimeFunc_om = Erf_om;
    mTimeFunc_omtt = Erf_omtt;
    break;
  case iVerySmoothBump:
    mTimeFunc = VerySmoothBump;
    mTimeFunc_t = VerySmoothBump_t;
    mTimeFunc_tt = VerySmoothBump_tt;
    mTimeFunc_ttt = VerySmoothBump_ttt;
    mTimeFunc_om = VerySmoothBump_om;
    mTimeFunc_omtt = VerySmoothBump_omtt;
    break;
  case iRickerInt:
    mTimeFunc = RickerInt;
    mTimeFunc_t = RickerInt_t;
    mTimeFunc_tt = RickerInt_tt;
    mTimeFunc_ttt = RickerInt_ttt;
    mTimeFunc_om = RickerInt_om;
    mTimeFunc_omtt = RickerInt_omtt;
    break;
  case iBrune:
    mTimeFunc = Brune;
    mTimeFunc_t = Brune_t;
    mTimeFunc_tt = Brune_tt;
    mTimeFunc_ttt = Brune_ttt;
    mTimeFunc_om = Brune_om;
    mTimeFunc_omtt = Brune_omtt;
    break;
  case iBruneSmoothed:
    mTimeFunc = BruneSmoothed;
    mTimeFunc_t = BruneSmoothed_t;
    mTimeFunc_tt = BruneSmoothed_tt;
    mTimeFunc_ttt = BruneSmoothed_ttt;
    mTimeFunc_om = BruneSmoothed_om;
    mTimeFunc_omtt = BruneSmoothed_omtt;
    break;
  case iDBrune:
    mTimeFunc = DBrune;
    mTimeFunc_t = DBrune_t;
    mTimeFunc_tt = DBrune_tt;
    mTimeFunc_ttt = DBrune_ttt;
    mTimeFunc_om = DBrune_om;
    mTimeFunc_omtt = DBrune_omtt;
    break;
  case iGaussianWindow:
    //      mPar[1] = mNcyc;
    mTimeFunc = GaussianWindow;
    mTimeFunc_t = GaussianWindow_t;
    mTimeFunc_tt = GaussianWindow_tt;
    mTimeFunc_ttt = GaussianWindow_ttt;
    mTimeFunc_om = GaussianWindow_om;
    mTimeFunc_omtt = GaussianWindow_omtt;
    break;
  case iLiu:
    mTimeFunc = Liu;
    mTimeFunc_t = Liu_t;
    mTimeFunc_tt = Liu_tt;
    mTimeFunc_ttt = Liu_ttt;
    mTimeFunc_om = Liu_om;
    mTimeFunc_omtt = Liu_omtt;
    break;
  case iDirac:
    mTimeFunc = Dirac;
    mTimeFunc_t = Dirac_t;
    mTimeFunc_tt = Dirac_tt;
    mTimeFunc_ttt = Dirac_ttt;
    mTimeFunc_om = Dirac_om;
    mTimeFunc_omtt = Dirac_omtt;
    break;
  case iDiscrete:
    mTimeFunc = Discrete;
    mTimeFunc_t = Discrete_t;
    mTimeFunc_tt = Discrete_tt;
    mTimeFunc_ttt = Discrete_ttt;
    mTimeFunc_om = Discrete_om;
    mTimeFunc_omtt = Discrete_omtt;
    break;
  case iC6SmoothBump:
    mTimeFunc = C6SmoothBump;
    mTimeFunc_t = C6SmoothBump_t;
    mTimeFunc_tt = C6SmoothBump_tt;
    mTimeFunc_ttt = C6SmoothBump_ttt;
    mTimeFunc_om = C6SmoothBump_om;
    mTimeFunc_omtt = C6SmoothBump_omtt;
    break;
  default:
    std::cout << "incorrect argument to GridPointSource constructor : default RickerWavelet used " << std::endl;
    mTimeFunc = RickerWavelet;
    mTimeFunc_t = RickerWavelet_t;
    mTimeFunc_tt = RickerWavelet_tt;
    mTimeFunc_ttt = RickerWavelet_ttt;
    mTimeFunc_om = RickerWavelet_om;
    mTimeFunc_omtt = RickerWavelet_omtt;
  }
  // Treat fourth derivatives in special 'switch', because not (yet?)
  // implemented for all time functions
  switch( mTimeDependence )
  {
  case iVerySmoothBump:
    mTimeFunc_tttt = VerySmoothBump_tttt;
    mTimeFunc_tttom = VerySmoothBump_tttom;
    mTimeFunc_ttomom = VerySmoothBump_ttomom;
    mTimeFunc_tom = VerySmoothBump_tom;
    mTimeFunc_omom = VerySmoothBump_omom;
    break;
  case iGaussian:
    mTimeFunc_tttt = Gaussian_tttt;
    mTimeFunc_tttom = Gaussian_tttom;
    mTimeFunc_ttomom = Gaussian_ttomom;
    mTimeFunc_tom = Gaussian_tom;
    mTimeFunc_omom = Gaussian_omom;
    break;
  case iDirac:
    mTimeFunc_tttt = Dirac_tttt;
    mTimeFunc_tttom = Dirac_tttom;
    mTimeFunc_ttomom = Dirac_ttomom;
    mTimeFunc_tom = Dirac_tom;
    mTimeFunc_omom = Dirac_omom;
    break;
  case iDiscrete:
    mTimeFunc_tttt = Discrete_tttt;
    mTimeFunc_tttom = Discrete_tttom;
    mTimeFunc_ttomom = Discrete_ttomom;
    mTimeFunc_tom = Discrete_tom;
    mTimeFunc_omom = Discrete_omom;
    break;
  default:
// tmp
// std::cout << "High derivatives not implemented for time fuction:" <<
// mTimeDependence <<
//   " default Gaussian used for tttt, ttt-omega derivatives, etc " <<
// std::endl;
    mTimeFunc_tttt = Gaussian_tttt;
    mTimeFunc_tttom = Gaussian_tttom;
    mTimeFunc_ttomom = Gaussian_ttomom;
    mTimeFunc_tom = Gaussian_tom;
    mTimeFunc_omom = Gaussian_omom;
  }
}

//-----------------------------------------------------------------------
void GridPointSource::getFxyz( double t, double* fxyz ) const
{
  double afun = mTimeFunc(mFreq,t-mT0,mPar, mNpar, mIpar, mNipar );
  fxyz[0] = mForces[0]*afun;
  fxyz[1] = mForces[1]*afun;
  fxyz[2] = mForces[2]*afun;
}

//-----------------------------------------------------------------------
void GridPointSource::getFxyz_notime( double* fxyz ) const
{
// For source spatial discretization testing
  fxyz[0] = mForces[0];
  fxyz[1] = mForces[1];
  fxyz[2] = mForces[2];
}

//-----------------------------------------------------------------------
void GridPointSource::getFxyztt( double t, double* fxyz ) const
{
  double afun = mTimeFunc_tt(mFreq,t-mT0,mPar, mNpar, mIpar, mNipar);
  fxyz[0] = mForces[0]*afun;
  fxyz[1] = mForces[1]*afun;
  fxyz[2] = mForces[2]*afun;
}

//-----------------------------------------------------------------------
void GridPointSource::limitFrequency(double max_freq)
{
  if (mFreq > max_freq)
    mFreq=max_freq;
}

//-----------------------------------------------------------------------
double GridPointSource::getTimeFunc(double t) const
{
  return mTimeFunc(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
double GridPointSource::evalTimeFunc_t(double t) const
{
  return mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
double GridPointSource::evalTimeFunc_tt(double t) const
{
  return mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
double GridPointSource::evalTimeFunc_ttt(double t) const
{
  return mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
double GridPointSource::evalTimeFunc_tttt(double t) const
{
  return mTimeFunc_tttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const GridPointSource& s )
{
  output << "GridPointSource at (i,j,k) = " << s.m_i0 << "," << s.m_j0 << "," << s.m_k0 <<
    " in grid no " << s.m_grid << endl;
  //   output << "   Strength " << s.mAmp;
  output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[2] << endl;
  output << " freq = " << s.mFreq << " t0 = " <<  s.mT0 << endl;
  output << " npar = " <<  s.mNpar << " nipar= " <<  s.mNipar  << endl;
  if(  s.mNpar > 0 )
    output << " mpar[0] = " <<  s.mPar[0];
  if(  s.mNipar > 0 )
    output << " mipar[0] = " <<  s.mIpar[0];
  if(  s.mNpar >0 ||  s.mNipar > 0 )
    cout << endl;


  return output;
}

//-----------------------------------------------------------------------
void GridPointSource::print_info() const
{
  cout << "-----------------------------------------------------------------------"<<endl;
  cout << " position " << m_i0 << " " << m_j0 << " " << m_k0 << endl;
  cout << "Forces = " << mForces[0] << " " << mForces[1] << " " << mForces[2] << endl;
  cout << "Time dep " << mTimeDependence << endl;
  cout << "-----------------------------------------------------------------------"<<endl;

}
