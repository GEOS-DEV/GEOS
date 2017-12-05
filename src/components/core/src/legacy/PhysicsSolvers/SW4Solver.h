#ifndef SW4SOLVER_H_
#define SW4SOLVER_H_

#include "SolverBase.h"
#include "SW4/SuperGrid.h"

class Source;
class Filter;
class GridPointSource;
class MaterialData;
class TimeSeries;

class SW4Solver : public SolverBase
{
public:

  SW4Solver(const std::string& name,
            ProblemManagerT* const pm );

  ~SW4Solver();

  double TimeStep( const realT& time,
                   const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions,
                   SpatialPartition& partition,
                   FractunatorBase* const fractunator );

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  //  void PostProcess (PhysicalDomainT& domain,
  //                    SpatialPartition& partition,
  //                    const array<string>& namesOfSolverRegions);

  void RegisterFields( PhysicalDomainT& domain );

  static const char* SolverName() {return "SW4Solver";}

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  void SetMaxStableTimeStep( const realT& time,
                             PhysicalDomainT& domain,
                             const array<string>& namesOfSolverRegions,
                             SpatialPartition& partition  );

  void GetStencilCoefficients( realT* _acof, realT* _ghcof, realT* _bope, realT* _sbop );

  void sw4discretization( realT* displacement, realT* lame_mu, realT* lame_lambda,
                          realT* metric, realT* jacobian, realT* lhs,
                          char op, realT* _strx, realT* _stry, realT* _strz,
                          realT* _acof, realT* _ghcof, realT* _bope );

private:
  void Forcing( realT time, array<R1Tensor>& F );

  void Forcing_tt( realT time, array<R1Tensor>& F );

  void evalPredictor( array<R1Tensor>& Up, array<R1Tensor>& U, array<R1Tensor>& Um,
                      array<realT>& rho, array<R1Tensor>& Lu, array<R1Tensor>& F, realT dt );
  void evalDpDmInTime( array<R1Tensor>& Up, array<R1Tensor>& U,  array<R1Tensor>& Um,
                       array<R1Tensor>& Uacc, realT dt );

  void evalCorrector( array<R1Tensor>& Up,  array<realT>& rho, array<R1Tensor>& Lu,
                      array<R1Tensor>& F, realT dt );

  void enforceBC( array<R1Tensor>& Up, array<R2SymTensor>& metric,
                  array<realT>& lame_mu, array<realT>& lame_lambda );

  void addSuperGridDamping( realT* Up, realT* U, realT* Um, realT* _rho,
                            realT* _dcx, realT* _dcy, realT* _dcz, realT* _strx, realT* _stry, realT* _strz,
                            realT* _jac, realT* _cox, realT* _coy, realT* _coz, realT coeff );

  void cycleSolutionArrays( array<R1Tensor>& Up, array<R1Tensor>& U, array<R1Tensor>& Um );

  void setupSupergrid( realT* strx, realT* stry, realT* strz, realT* dcx, realT* dcy, realT* dcz,
                       realT* cox, realT* coy, realT* coz );

  void setupCartesianMetric( array<R2SymTensor>& metric, array<realT>& jacobian );

  void preProcessSources( realT dt, int numberOfTimeSteps, realT tstart );

  bool startswith(const char begin[], char *line );

  void parseInputFile( std::string fname );
  void processSource(char* buffer );
  void processMaterialBlock(char* buffer);
  void processMaterialPfile(char* buffer);
  void processSupergrid(char * buffer);
  void processPrefilter(char* buffer);
  void processTimeSeries( char* buffer );
  void convert_to_mulambda( array<realT>& rho, array<realT>& mu,
                            array<realT>& lambda );
  void Forcing( realT time, realT* F );
  void Forcing_tt( realT time, realT* F );
  void evalPredictor( realT* Up, realT* U, realT* Um, realT* rho,
                      realT* Lu, realT* F, realT dt );
  void evalDpDmInTime( realT* Up, realT* U,  realT* Um,
                       realT* Uacc, realT dt );
  void evalCorrector( realT* Up, realT* rho, realT* Lu,
                      realT* F, realT dt );

  void create_output_directory( );
  int mkdirs(const std::string& path);

  // tw-testing:
  void init_smooth( realT* u, realT* um, realT* coords, realT dt );
  void debug_output( realT* Lu );

  // Stencil coefficient
  realT *m_acof, *m_ghcof, *m_bope, *m_sbop;

  // Variables used for supergrid sponge layers at far-field boundaries
  realT *m_strx, *m_stry, *m_strz, *m_dcx, *m_dcy, *m_dcz, *m_cox, *m_coy, *m_coz;
  SuperGrid m_supergrid_taper_x, m_supergrid_taper_y, m_supergrid_taper_z;
  int   m_sg_ptsinlayer;
  realT m_sg_width, m_sg_damping_coefficient;

  Filter* m_source_filter;
  std::vector<Source*> m_sources;
  std::vector<GridPointSource*> m_point_sources;
  std::vector<MaterialData*> m_material;
  std::vector<TimeSeries*> m_stations;

// Index limits, local to my processor.
  int m_ifirst, m_ilast, m_jfirst, m_jlast, m_kfirst, m_klast;
  // Total number of points
  size_t m_npts;

// Index limits without ghost points, local to my processor.
  int m_ifirst_int, m_ilast_int, m_jfirst_int, m_jlast_int, m_kfirst_int, m_klast_int;

// k-index such that the grid can be considered Cartesian for k >= m_kcart
  int m_kcart;

  int m_nghost; // Number of ghost points
  int m_npadding; // Number of parallel overlap points

  // Global size is 1..m_ni_global etc. for j and k.
  int m_ni_global, m_nj_global, m_nk_global;

  //Grid spacing, same in all directions.
  realT m_dx;

  // Grid is given by x_i= (i-1)*dx + m_xmin_global, i=-nghost+1,..,N+nghost,
  // where
  // the domain is i=1,..,N, and x_1=m_xmin_global, x_N = m_xmax_global.
  // Similarly for y_j and z_k.
  realT m_xmin_global, m_xmax_global, m_ymin_global, m_ymax_global, m_zmin_global, m_zmax_global;

  // Boundary conditions 0-no boundary (parallel overlap), 1-Far field supergrid
  // boundary, 2-free surface.
  int m_bctype[6];
  bool m_onesided; // One sided difference formula at the upper boundary?

  // Supergrid far-field boundary condition

  std::string m_input_file;
  std::string m_output_path;
  bool m_output_timing;
};
#endif
