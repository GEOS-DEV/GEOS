#ifndef EW_TIMEFUNCTIONS_H
#define EW_TIMEFUNCTIONS_H

double VerySmoothBump(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_tom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_omom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_tttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_tttom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double VerySmoothBump_ttomom(double freq, double t, double* par, int npar, int* ipar, int nipar );

double RickerWavelet(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerWavelet_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerWavelet_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerWavelet_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerWavelet_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerWavelet_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double RickerInt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerInt_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerInt_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerInt_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerInt_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double RickerInt_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Gaussian(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_tom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_omom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_tttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_tttom(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Gaussian_ttomom(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Erf( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Erf_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Erf_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Erf_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Erf_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Erf_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Ramp(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Ramp_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Ramp_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Ramp_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Ramp_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Ramp_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Triangle(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Triangle_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Triangle_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Triangle_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Triangle_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Triangle_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Sawtooth(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Sawtooth_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Sawtooth_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Sawtooth_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Sawtooth_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double Sawtooth_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double SmoothWave(double freq, double t, double* par, int npar, int* ipar, int nipar );
double SmoothWave_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double SmoothWave_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double SmoothWave_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double SmoothWave_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double SmoothWave_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

double Brune( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Brune_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Brune_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Brune_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Brune_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Brune_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );

double DBrune( double freq, double t, double* par, int npar, int* ipar, int nipar );
double DBrune_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double DBrune_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double DBrune_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double DBrune_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double DBrune_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );

double BruneSmoothed( double freq, double t, double* par, int npar, int* ipar, int nipar );
double BruneSmoothed_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double BruneSmoothed_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double BruneSmoothed_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double BruneSmoothed_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double BruneSmoothed_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );

double GaussianWindow( double freq, double t, double* par, int npar, int* ipar, int nipar );
double GaussianWindow_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double GaussianWindow_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double GaussianWindow_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double GaussianWindow_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double GaussianWindow_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );

double Liu( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Liu_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Liu_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Liu_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Liu_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Liu_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );

double NullFunc( double freq, double t, double* par, int npar, int* ipar, int nipar );

double Dirac( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_tttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_tom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_omom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_tttom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Dirac_ttomom( double freq, double t, double* par, int npar, int* ipar, int nipar );

double Discrete( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_t( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_tt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_ttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_tttt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_om( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_tom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_omom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_omtt( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_tttom( double freq, double t, double* par, int npar, int* ipar, int nipar );
double Discrete_ttomom( double freq, double t, double* par, int npar, int* ipar, int nipar );

double C6SmoothBump(double freq, double t, double* par, int npar, int* ipar, int nipar );
double C6SmoothBump_t(double freq, double t, double* par, int npar, int* ipar, int nipar );
double C6SmoothBump_om(double freq, double t, double* par, int npar, int* ipar, int nipar );
double C6SmoothBump_tt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double C6SmoothBump_ttt(double freq, double t, double* par, int npar, int* ipar, int nipar );
double C6SmoothBump_omtt(double freq, double t, double* par, int npar, int* ipar, int nipar );

#endif
