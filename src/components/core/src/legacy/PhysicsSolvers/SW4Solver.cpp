#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "SW4Solver.h"
#include "SolverFactory.h"
#include "SW4/Source.h"
#include "SW4/GridPointSource.h"
#include "SW4/Filter.h"
#include "SW4/MaterialBlock.h"
#include "SW4/MaterialPfile.h"
#include "SW4/TimeSeries.h"

//-----------------------------------------------------------------------
SW4Solver::SW4Solver( const std::string& name,
                      ProblemManagerT* const pm ) :
    SolverBase( name, pm )
{
  m_acof = new realT[6 * 8 * 8];
  m_ghcof = new realT[6];
  m_bope = new realT[6 * 8];
  m_sbop = new realT[5];
  m_source_filter = static_cast<Filter*>( 0 );
  GetStencilCoefficients( m_acof, m_ghcof, m_bope, m_sbop );
  m_nghost = 2;
  m_npadding = 2;
  m_sg_damping_coefficient = 0.02;
  m_sg_ptsinlayer = 20;
  m_output_timing = false;
  // Assume Cartesian grid
  m_kcart = 5;
}

//-----------------------------------------------------------------------
SW4Solver::~SW4Solver()
{
  for( unsigned int s = 0 ; s < m_stations.size() ; s++ )
    m_stations[s]->writeFile();

  delete[] m_acof;
  delete[] m_ghcof;
  delete[] m_bope;
  delete[] m_sbop;
}

//-----------------------------------------------------------------------
double SW4Solver::TimeStep( const realT& time,
                            const realT& dt,
                            const int cycleNumber,
                            PhysicalDomainT& domain,
                            const sArray1d& namesOfSolverRegions,
                            SpatialPartition& partition,
                            FractunatorBase* const fractunator )
{
  realT dt_taken = dt;
  //   std::cout << "SW4 Step no. " <<  cycleNumber << std::endl;
  //   std::cout << " stable dt = " <<   m_stabledt.m_maxdt << std::endl;
  //   std::cout << "        dt = " <<   dt << std::endl;
  //   std::cout << "       cfl = " << m_courant << std::endl;

  ObjectDataStructureBaseT& objectManager = domain.m_feNodeManager;

  Array1dT<realT>& rho = objectManager.GetFieldData<realT>( "rho" );
  Array1dT<realT>& mu = objectManager.GetFieldData<realT>( "mu" );
  Array1dT<realT>& lambda = objectManager.GetFieldData<realT>( "lambda" );
  Array1dT<R1Tensor>& Lu = objectManager.GetFieldData<R1Tensor>( "Lu" );
  Array1dT<R1Tensor>& F = objectManager.GetFieldData<R1Tensor>( "Src" );
  Array1dT<R1Tensor>& U = objectManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>& Uacc = objectManager.GetFieldData<R1Tensor>( "Uacc" );
  Array1dT<R1Tensor>& Um = objectManager.GetFieldData<R1Tensor>( "Um" );
  Array1dT<R1Tensor>& Up = objectManager.GetFieldData<R1Tensor>( "Up" );
  Array1dT<R2SymTensor>& metric = objectManager.GetFieldData<R2SymTensor>( "metric" );
  Array1dT<realT>& jacobian = objectManager.GetFieldData<realT>( "jacobian" );

  std::vector<realT> t( 13 );

  t[0] = MPI_Wtime();
  Forcing( time, &F[0][0] );

  t[1] = MPI_Wtime();
  sw4discretization( &U[0][0], &mu[0], &lambda[0], &metric[0].Data()[0], &jacobian[0], &Lu[0][0], '=', m_strx,
                     m_stry,
                     m_strz, m_acof, m_ghcof, m_bope );

  //   debug_output( &Lu[0][0]);

  t[2] = MPI_Wtime();
  evalPredictor( &Up[0][0], &U[0][0], &Um[0][0], &rho[0], &Lu[0][0], &F[0][0], dt );

  t[3] = MPI_Wtime();
// Communicate Up here

  enforceBC( Up, metric, mu, lambda );

  t[4] = MPI_Wtime();
  Forcing_tt( time, &F[0][0] );

  t[5] = MPI_Wtime();
  evalDpDmInTime( &Up[0][0], &U[0][0], &Um[0][0], &Uacc[0][0], dt ); // store result in Uacc

  t[6] = MPI_Wtime();
  sw4discretization( &Uacc[0][0], &mu[0], &lambda[0], &metric[0].Data()[0], &jacobian[0], &Lu[0][0], '=', m_strx,
                     m_stry,
                     m_strz, m_acof, m_ghcof, m_bope );

  t[7] = MPI_Wtime();
  evalCorrector( &Up[0][0], &rho[0], &Lu[0][0], &F[0][0], dt );

  t[8] = MPI_Wtime();
  addSuperGridDamping( &Up[0][0], &U[0][0], &Um[0][0], &rho[0], m_dcx, m_dcy, m_dcz,
                       m_strx,
                       m_stry, m_strz, &jacobian[0],
                       m_cox,
                       m_coy, m_coz, m_sg_damping_coefficient );

  t[9] = MPI_Wtime();
// Communicate Up here

  enforceBC( Up, metric, mu, lambda );
  t[10] = MPI_Wtime();

  localIndex nj = ( m_jlast - m_jfirst + 1 );
  localIndex nk = ( m_klast - m_kfirst + 1 );
  for( unsigned int s = 0 ; s < m_stations.size() ; s++ )
  {
    localIndex ind = nk * nj * ( m_stations[s]->m_i0 - m_ifirst ) + nk * ( m_stations[s]->m_j0 - m_jfirst ) + ( m_stations[s]
        ->m_k0 - m_kfirst );
    vector<realT> urec( 3 );
    urec[0] = U[ind][0];
    urec[1] = U[ind][1];
    urec[2] = U[ind][2];
    m_stations[s]->recordData( urec );
  }

  t[11] = MPI_Wtime();
  cycleSolutionArrays( Up, U, Um );

  t[12] = MPI_Wtime();

  if( m_output_timing )
  {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Timing for one step " << std::endl;
    std::cout << " Forcing    " << t[1] - t[0] << " sec." << std::endl;
    std::cout << " Disc.1     " << t[2] - t[1] << " sec." << std::endl;
    std::cout << " Predictor  " << t[3] - t[2] << " sec." << std::endl;
    std::cout << " BC         " << t[4] - t[3] << " sec." << std::endl;
    std::cout << " Forcing_tt " << t[5] - t[4] << " sec." << std::endl;
    std::cout << " DpDm intime" << t[6] - t[5] << " sec." << std::endl;
    std::cout << " Disc.2     " << t[7] - t[6] << " sec." << std::endl;
    std::cout << " Corrector  " << t[8] - t[7] << " sec." << std::endl;
    std::cout << " Sg.damping " << t[9] - t[8] << " sec." << std::endl;
    std::cout << " BC         " << t[10] - t[9] << " sec." << std::endl;
    std::cout << " save SAC   " << t[11] - t[10] << " sec." << std::endl;
    std::cout << " cycle array" << t[12] - t[11] << " sec." << std::endl;
  }

  return dt_taken;
}

//-----------------------------------------------------------------------
void SW4Solver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  std::cout << "SW4 initializing  " << std::endl;

  ObjectDataStructureBaseT& objectManager = domain.m_feNodeManager;

  Array1dT<realT>& jac = objectManager.GetFieldData<realT>( "jacobian" );
  Array1dT<R2SymTensor>& metric = objectManager.GetFieldData<R2SymTensor>( "metric" );
  Array1dT<realT>& rho = objectManager.GetFieldData<realT>( "rho" );
  Array1dT<realT>& mu = objectManager.GetFieldData<realT>( "mu" );
  Array1dT<realT>& lambda = objectManager.GetFieldData<realT>( "lambda" );

  // five supergrid boundaries (==1), and one free surface (==2)
  m_bctype[0] = 1;
  m_bctype[1] = 1;
  m_bctype[2] = 1;
  m_bctype[3] = 1;
  m_bctype[4] = 2;
  m_bctype[5] = 1;

  //   std::cout << " sizeof R1Tensor is " << sizeof(R1Tensor) << std::endl;
  // Try to back out the dimensions, assuming Cartesian grid.
  Array1dT<R1Tensor>& coords = objectManager.GetFieldData<FieldInfo::referencePosition>();
  std::vector<realT> xpos;
  std::vector<realT> ypos;
  std::vector<realT> zpos;
  for( unsigned int i = 0 ; i < coords.size() ; i++ )
  {
    std::vector<realT>::iterator it = std::find( xpos.begin(), xpos.end(), coords[i][0] );
    if( it == xpos.end() )
      xpos.push_back( coords[i][0] );
    it = std::find( ypos.begin(), ypos.end(), coords[i][1] );
    if( it == ypos.end() )
      ypos.push_back( coords[i][1] );
    it = std::find( zpos.begin(), zpos.end(), coords[i][2] );
    if( it == zpos.end() )
      zpos.push_back( coords[i][2] );
  }
  std::cout << "Local dimensions found " << xpos.size() << " x " << ypos.size() << " x " << zpos.size() << std::endl;
  std::sort( xpos.begin(), xpos.end() );
  std::sort( ypos.begin(), ypos.end() );
  std::sort( zpos.begin(), zpos.end() );
  m_dx = xpos.size() > 0 ? xpos[1] - xpos[0] : 1;

  R1Tensor mindim, maxdim;
  partition.getSizes( mindim, maxdim );

  m_xmax_global = maxdim[0];
  m_ymax_global = maxdim[1];
  m_zmax_global = maxdim[2];
  m_xmin_global = mindim[0];
  m_ymin_global = mindim[1];
  m_zmin_global = mindim[2];

  // Adjust for ghost points outside the domain. User needs to specify
  // the domain size including ghost points.
  m_xmin_global += m_nghost * m_dx;
  m_ymin_global += m_nghost * m_dx;
  m_zmin_global += m_nghost * m_dx;
  m_xmax_global -= m_nghost * m_dx;
  m_ymax_global -= m_nghost * m_dx;
  m_zmax_global -= m_nghost * m_dx;

  m_ni_global = static_cast<int>( round( ( m_xmax_global - m_xmin_global ) / m_dx + 1 ) );
  m_nj_global = static_cast<int>( round( ( m_ymax_global - m_ymin_global ) / m_dx + 1 ) );
  m_nk_global = static_cast<int>( round( ( m_zmax_global - m_zmin_global ) / m_dx + 1 ) );

  std::cout << "Global dimensions " << m_xmin_global << " <= x <= " << m_xmax_global <<
            ", "
            << m_ymin_global << " <= y <= " << m_ymax_global <<
            ", "
            << m_zmin_global << " <= z <= " << m_zmax_global << std::endl;

  std::cout << "Grid size 1..N  Nx= " << m_ni_global << " Ny= " << m_nj_global << " Nz= " << m_nk_global << endl;
  // Define grid coordinate mapping s x[i] = (i-1)*dx + xmin  , i=-g+1,-g+2,..0,1,...,N,N+1,..,N+g
  //  where -g+1,..,0 and N+1,..,N+g are ghost points, g=number of ghost points
  //      i=1 is lower boundary and i=N is upper boundary.
  //  We will consider m_xmin_global and m_xmax_global as the locations of i=1 and i=N respectively.
  //  and  xmin = m_xmin_global. g=m_nghost, N= m_ni_global

  m_ifirst = static_cast<int>( round( ( xpos[0] - m_xmin_global ) / m_dx + 1 ) );
  m_ilast = static_cast<int>( round( ( xpos[xpos.size() - 1] - m_xmin_global ) / m_dx + 1 ) );
  std::cout << " local x-index range= " << m_ifirst << " " << m_ilast << std::endl;
  m_jfirst = static_cast<int>( round( ( ypos[0] - m_ymin_global ) / m_dx + 1 ) );
  m_jlast = static_cast<int>( round( ( ypos[ypos.size() - 1] - m_ymin_global ) / m_dx + 1 ) );
  std::cout << " local y-index range= " << m_jfirst << " " << m_jlast << std::endl;
  m_kfirst = static_cast<int>( round( ( zpos[0] - m_zmin_global ) / m_dx + 1 ) );
  m_klast = static_cast<int>( round( ( zpos[zpos.size() - 1] - m_zmin_global ) / m_dx + 1 ) );
  std::cout << " local z-index range= " << m_kfirst << " " << m_klast << std::endl;

  // Check that boundary stencil fits in processor
  if( !( ( m_kfirst <= 0 && m_klast >= 8 ) || ( m_kfirst >= 7 ) ) )
    std::cout << "ERROR: boundary stencil does not fit into processor" << std::endl;

  m_onesided = m_kfirst <= 0;

  // Set up interior dimensions
  if( m_ifirst == 1 - m_nghost )
    m_ifirst_int = 1;
  else
    m_ifirst_int = m_ifirst + m_npadding;

  if( m_ilast == m_ni_global + m_nghost )
    m_ilast_int = m_ni_global;
  else
    m_ilast_int = m_ilast - m_npadding;

  if( m_jfirst == 1 - m_nghost )
    m_jfirst_int = 1;
  else
    m_jfirst_int = m_jfirst + m_npadding;

  if( m_jlast == m_nj_global + m_nghost )
    m_jlast_int = m_nj_global;
  else
    m_jlast_int = m_jlast - m_npadding;

  if( m_kfirst == 1 - m_nghost )
    m_kfirst_int = 1;
  else
    m_kfirst_int = m_kfirst + m_npadding;

  if( m_klast == m_nk_global + m_nghost )
    m_klast_int = m_nk_global;
  else
    m_klast_int = m_klast - m_npadding;

  // No boundary conditions at processor boundaries
  if( m_ifirst > 1 - m_nghost )
    m_bctype[0] = 0;
  if( m_ilast < m_ni_global + m_nghost )
    m_bctype[1] = 0;
  if( m_jfirst > 1 - m_nghost )
    m_bctype[2] = 0;
  if( m_jlast < m_nj_global + m_nghost )
    m_bctype[3] = 0;
  if( m_kfirst > 1 - m_nghost )
    m_bctype[4] = 0;
  if( m_klast < m_nk_global + m_nghost )
    m_bctype[5] = 0;

  // Process SW4 input file
  parseInputFile( m_input_file );

  // Assign material
  for( unsigned int bl = 0 ; bl < m_material.size() ; bl++ )
    m_material[bl]->set_material_properties( rho, mu, lambda, coords, 0.0 );

  //   int imid = (m_ilast+m_ifirst)/2, jmid = (m_jlast+m_jfirst)/2;
  //   int nk=(m_klast-m_kfirst+1), nj=m_jlast-m_jfirst+1;
  //   size_t ind = -m_kfirst + nk*(jmid-m_jfirst) + nk*nj*(imid-m_ifirst);
  //   for( int k=m_kfirst ; k <= m_klast ; k++ )
  //      cout << "z= " << m_zmin_global + (k-1)*m_dx  << " k=" << k << " rho, cp, cs = " << rho[ind+k] << " " << mu[ind+k]<< " " << lambda[ind+k] << endl;

  convert_to_mulambda( rho, mu, lambda );

  // Define metric and grid
  setupCartesianMetric( metric, jac );

  const sArray1d dummy;
  SetMaxStableTimeStep( 0.0, domain, dummy, partition );

  // Supergrid damping layers at non-reflecting boundaries
  m_strx = new realT[m_ilast - m_ifirst + 1];
  m_stry = new realT[m_jlast - m_jfirst + 1];
  m_strz = new realT[m_klast - m_kfirst + 1];
  m_dcx = new realT[m_ilast - m_ifirst + 1];
  m_dcy = new realT[m_jlast - m_jfirst + 1];
  m_dcz = new realT[m_klast - m_kfirst + 1];
  m_cox = new realT[m_ilast - m_ifirst + 1];
  m_coy = new realT[m_jlast - m_jfirst + 1];
  m_coz = new realT[m_klast - m_kfirst + 1];
  double sg_width = m_sg_ptsinlayer * m_dx;

  m_supergrid_taper_x.define_taper( true, m_xmin_global, true, m_xmax_global, sg_width );
  m_supergrid_taper_y.define_taper( true, m_ymin_global, true, m_ymax_global, sg_width );
  m_supergrid_taper_z.define_taper( false, m_zmin_global, true, m_zmax_global, sg_width );

  setupSupergrid( m_strx, m_stry, m_strz, m_dcx, m_dcy, m_dcz, m_cox, m_coy, m_coz );

  realT dt = m_stabledt.m_maxdt;

  // CHANGE: Find out how to get nsteps
  //   int nsteps;
  //   realT tstart=0;
  //   preProcessSources( dt, nsteps, tstart );

  int nsteps = 3000;

  int interior[6] =
  { m_ifirst_int, m_ilast_int, m_jfirst_int, m_jlast_int, m_kfirst_int, m_klast_int };
  for( unsigned int s = 0 ; s < m_sources.size() ; s++ )
    m_sources[s]->set_grid_point_sources4( m_dx, m_xmin_global, m_ymin_global, m_zmin_global,
                                           m_ni_global,
                                           m_nj_global, m_nk_global, m_point_sources,
                                           interior );

  for( unsigned int s = 0 ; s < m_stations.size() ; s++ )
    m_stations[s]->allocateRecordingArrays( nsteps, 0.0, dt );

  xpos.clear();
  ypos.clear();
  zpos.clear();

  m_npts = ( static_cast<size_t>( m_ilast - m_ifirst + 1 ) ) * ( m_jlast - m_jfirst + 1 ) * ( m_klast - m_kfirst + 1 );
  // Give initial data
  Array1dT<R1Tensor>& U = objectManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>& Um = objectManager.GetFieldData<R1Tensor>( "Um" );

  U = 0;
  Um = 0;

  // DEBUG
  //   init_smooth( &U[0][0], &Um[0][0], &coords[0][0], dt );

}

//-----------------------------------------------------------------------
void SW4Solver::init_smooth( realT* u, realT* um, realT* coords, realT dt )
{
//  int ni = m_ilast - m_ifirst + 1,
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  double t = 0;
  double om = 1, cv = 1.2, ph = 0.1;
  for( int k = m_kfirst ; k <= m_klast ; k++ )
    for( int j = m_jfirst ; j <= m_jlast ; j++ )
      for( int i = m_ifirst ; i <= m_ilast ; i++ )
      {
        size_t ind = k - m_kfirst + nk * ( j - m_jfirst ) + nk * nj * ( i - m_ifirst );
        realT x = coords[3 * ind];
        realT y = coords[3 * ind + 1];
        realT z = coords[3 * ind + 2];
        u[3 * ind] = sin( om * ( x - cv * t ) ) * sin( om * y + ph ) * sin( om * z + ph );
        u[3 * ind + 1] = sin( om * x + ph ) * sin( om * ( y - cv * t ) ) * sin( om * z + ph );
        u[3 * ind + 2] = sin( om * x + ph ) * sin( om * y + ph ) * sin( om * ( z - cv * t ) );
        um[3 * ind] = sin( om * ( x - cv * ( t - dt ) ) ) * sin( om * y + ph ) * sin( om * z + ph );
        um[3 * ind + 1] = sin( om * x + ph ) * sin( om * ( y - cv * ( t - dt ) ) ) * sin( om * z + ph );
        um[3 * ind + 2] = sin( om * x + ph ) * sin( om * y + ph ) * sin( om * ( z - cv * ( t - dt ) ) );
      }
}

//-----------------------------------------------------------------------
void SW4Solver::debug_output( realT* Lu )
{
//  int ni = m_ilast - m_ifirst + 1,
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  int idbg = 12, jdbg = 14;
  for( int k = m_kfirst ; k <= m_klast ; k++ )
  {
    size_t ind = k - m_kfirst + nk * ( jdbg - m_jfirst ) + nk * nj * ( idbg - m_ifirst );
    cout << "Lu at k= " << k << " " << Lu[3 * ind] << " " << Lu[3 * ind + 1] << " " << Lu[3 * ind + 2] << endl;
  }
}

//-----------------------------------------------------------------------
void SW4Solver::InitializeCommunications( PartitionBase& partition )
{
  std::cout << "SW4 initialized communications  " << std::endl;
}

//-----------------------------------------------------------------------
void SW4Solver::RegisterFields( PhysicalDomainT& domain )
{
  std::cout << "SW4 Registered fields  " << std::endl;
  domain.m_feNodeManager.AddKeylessDataField<realT>( "rho", true, false );
  domain.m_feNodeManager.AddKeylessDataField<realT>( "mu", true, false );
  domain.m_feNodeManager.AddKeylessDataField<realT>( "lambda", true, false );
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>( "Src", false, false );
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>( "Lu", false, false );
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>( "Up", false, false );
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>( "Um", false, false );
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>( "Uacc", false, false );
  domain.m_feNodeManager.AddKeylessDataField<R2SymTensor>( "metric", false, false );
  domain.m_feNodeManager.AddKeylessDataField<realT>( "jacobian", false, false );
}

//-----------------------------------------------------------------------
void SW4Solver::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );
  m_input_file = hdn->GetAttributeString( "inputfile" );
  m_output_path = hdn->GetAttributeString( "outputpath" );
  if( m_output_path[m_output_path.size() - 1] != '/' )
    m_output_path.append( "/" );
  create_output_directory();
}

//-----------------------------------------------------------------------
void SW4Solver::create_output_directory()
{
  int myrank;
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  if( myrank == 0 )
  {

    cout << "----------------------------------------------------" << endl
         << " Making Output Directory: "
         << m_output_path << endl
         << "\t\t"
         << endl;
    // Create directory where all files will be written.
    int err = mkdirs( m_output_path );

    if( err == 0 )
      cout << "... Done!" << endl
           << "----------------------------------------------------"
           << endl;
    else
    {
// fatal error
      cerr << endl << "******** Failed to create the output directory *******" << endl << endl;
      //	MPI_Abort(MPI_COMM_WORLD,1);
    }

// check that we have write permission on the directory
    if( access( m_output_path.c_str(), W_OK ) != 0 )
    {
// fatal error
      cerr << endl << "Error: No write permission on output directory: " << m_output_path << endl;
      //	MPI_Abort(MPI_COMM_WORLD,1);
    }

  }
  // Let processor 0 finish first!

  cout.flush();
  cerr.flush();
  MPI_Barrier( MPI_COMM_WORLD );

// Check that the mPath directory exists from all processes
  struct stat statBuf;
  int statErr = stat( m_output_path.c_str(), &statBuf );
  if( !( statErr == 0 && S_ISDIR( statBuf.st_mode ) ) )
    cout << "Error: " << m_output_path << " is not a directory" << endl;

// check that all processes have write permission on the directory
  if( !access( m_output_path.c_str(), W_OK ) == 0 )
    cout << "Error: No write permission on output directory: " << m_output_path << endl;
}

//-----------------------------------------------------------------------
int SW4Solver::mkdirs( const string& path )
{
  int mVerbose = 1;
  string pathTemp( path.begin(), path.end() );
  //-----------------------------------------------------------------
  // Recursively call stat and then mkdir on each sub-directory in 'path'
  //-----------------------------------------------------------------
  string sep = "/";
  char* token = strtok( const_cast<char*>( pathTemp.c_str() ), sep.c_str() );

  stringstream pathsofar;

// for checking the status:
  struct stat statBuf;
  int statErr;

  // If there's a leading slash, put it back on...
  if( strncmp( pathTemp.c_str(), sep.c_str(), 1 ) == 0 )
    pathsofar << sep;

  while( token != NULL )
  {
    pathsofar << token << sep;

// test: check the status of the path so far...
//      cout << "Calling stat() on path: " << pathsofar.str() << endl;
    statErr = stat( pathsofar.str().c_str(), &statBuf );
    if( statErr == 0 )
    {
//	cout << "stat() returned successfully." << endl;
      if( S_ISDIR( statBuf.st_mode ) )
      {
//	  cout << "stat() says: '" << pathsofar.str() << "' is a directory." << endl;
// it already exists, this is okay, let's get the next directory in the string and skip to the while statement
        token = strtok( NULL, sep.c_str() );
        continue;
      }
      else
      {
        cerr << "stat() says: '" << pathsofar.str() << "' is not a directory." << endl;
        // real error, let's bail...
        return -1;
      }

    }
    else
    {
//	cerr << "stat() returned an error code." << endl;
      if( errno == EACCES )
      {
        cerr << "Error: **Search permission is denied for one of the directories in the path prefix of "
             << pathsofar.str() << endl;
        return -1;
      }
      else if( errno == ENOTDIR )
      {
        cerr << "Error: **A component of the path '" << pathsofar.str() << "' is not a directory. " << endl;
        return -1;
      }
      else if( errno == ENOENT )
      {
// this means that we need to call mkdir to create the directory
        if( mVerbose >= 2 )
          cout << "Info: **stat returned ENOENT (the path does not exist, or the path " << endl
               << "      is an empty string) "
               << pathsofar.str() << endl;
      }
      else
      {
        if( mVerbose >= 2 )
          cout << "Info: **stat returned other error code for path: " << pathsofar.str() << endl;
      }
    }

// if we got this far, then 'pathsofar' does not exists

// tmp
    if( mVerbose >= 2 )
      cout << "Calling mkdir() on path: " << pathsofar.str() << endl;
// old code for recursively making the output directory
    if( mkdir( pathsofar.str().c_str(),
    S_IWUSR | S_IXUSR | S_IRUSR | S_IRGRP | S_IXGRP ) // why do we need group permissions?
    == -1 )
    {
      if( mVerbose >= 2 )
        cout << "mkdir() returned an error code." << endl;
      // check error conditions
      if( errno == EEXIST )
      {
// can this ever happen since we called stat(), which said that the directory did not exist ???
        if( mVerbose >= 2 )
          cout << "Info: ** The directory already exists:" << pathsofar.str() << endl;

        // it already exists, this is okay!
        token = strtok( NULL, sep.c_str() );
        continue;
      }
      else if( errno == EACCES )
        cerr
            << "Error: **Write permission is denied for the parent directory in which the new directory is to be added."
            << pathsofar.str() << endl;
      else if( errno == EMLINK )
        cerr << "Error: **The parent directory has too many links (entries)." <<
             pathsofar.str()
             << endl;
      else if( errno == ENOSPC )
        cerr << "Error: **The file system doesn't have enough room to create the new directory." <<
             pathsofar.str()
             << endl;
      else if( errno == EROFS )
        cerr
            << "Error: **  The parent directory of the directory being created is on a read-only file system and cannot be modified."
            << pathsofar.str() << endl;
      else if( errno == ENOSPC )
        cerr << "Error: ** The new directory cannot be created because the user's disk quota is exhausted."
             << pathsofar.str() << endl;
      // real error, let's bail...
      return -1;
    }
    else
    {
      if( mVerbose >= 2 )
        cout << "mkdir() returned successfully." << endl;

// are there more directories to be made?
      token = strtok( NULL, sep.c_str() );
    }
  }
  return 0;
}

//-----------------------------------------------------------------------
void SW4Solver::SetMaxStableTimeStep( const realT& time,
                                      PhysicalDomainT& domain,
                                      const sArray1d& namesOfSolverRegions,
                                      SpatialPartition& partition )
{
  std::cout << "Calling SetMaxStableTimeStep" << std::endl;
  ObjectDataStructureBaseT& objectManager = domain.m_feNodeManager;

  Array1dT<realT>& rho = objectManager.GetFieldData<realT>( "rho" );
  Array1dT<realT>& mu = objectManager.GetFieldData<realT>( "mu" );
  Array1dT<realT>& lambda = objectManager.GetFieldData<realT>( "lambda" );
  Array1dT<R2SymTensor>& metric = objectManager.GetFieldData<R2SymTensor>( "metric" );
  Array1dT<realT>& jac = objectManager.GetFieldData<realT>( "jacobian" );

#define SQR(x) ((x)*(x))
  realT dtloc = 1e38;
  for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
  {

    realT m12 = SQR( metric[ind].Data()[0] );
    realT a2 = SQR(metric[ind].Data()[1]) + SQR( metric[ind].Data()[2] ) + SQR( metric[ind].Data()[3] );
    realT d = sqrt( SQR(m12+a2) - 4 * m12 * SQR( metric[ind].Data()[3] ) );
    realT eig;
    if( a2 + d > m12 )
      eig = ( 2.5 * m12 + 1.5 * a2 + 0.5 * d ) * mu[ind] + ( 0.5 * m12 + 0.5 * a2 + 0.5 * d ) * lambda[ind];
    else
      eig = ( 3 * m12 + a2 ) * mu[ind] + m12 * lambda[ind];

    realT dtgp = m_courant * sqrt( rho[ind] * jac[ind] / eig );
    if( dtgp < dtloc )
      dtloc = dtgp;
  }
  realT dt = dtloc;
  //   MPI_Allreduce( &dtloc, &dt, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
  m_stabledt.m_maxdt = dt;
  std::cout << "SetMaxStableTimeStep done, dt = " << dt << std::endl;
#undef SQR
}

//-----------------------------------------------------------------------
void SW4Solver::Forcing( realT time, Array1dT<R1Tensor>& F )
{
  F = 0;
  //int ni = m_ilast - m_ifirst + 1;
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  for( unsigned int s = 0 ; s < m_point_sources.size() ; s++ )
  {
    double fxyz[3];
    m_point_sources[s]->getFxyz( time, fxyz );
    size_t ind = nj * nk * ( m_point_sources[s]->m_i0 - m_ifirst ) + nk * ( m_point_sources[s]->m_j0 - m_jfirst )
        + m_point_sources[s]->m_k0 - m_kfirst;
    F[ind][0] += fxyz[0];
    F[ind][1] += fxyz[1];
    F[ind][2] += fxyz[2];
  }
}

//-----------------------------------------------------------------------
void SW4Solver::Forcing_tt( realT time, Array1dT<R1Tensor>& F )
{
  F = 0;
  //int ni = m_ilast - m_ifirst + 1;
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  for( unsigned int s = 0 ; s < m_point_sources.size() ; s++ )
  {
    double fxyz[3];
    m_point_sources[s]->getFxyztt( time, fxyz );
    size_t ind = nj * nk * ( m_point_sources[s]->m_i0 - m_ifirst ) + nk * ( m_point_sources[s]->m_j0 - m_jfirst )
        + m_point_sources[s]->m_k0 - m_kfirst;
    F[ind][0] += fxyz[0];
    F[ind][1] += fxyz[1];
    F[ind][2] += fxyz[2];
  }
}

//-----------------------------------------------------------------------
void SW4Solver::Forcing( realT time, realT* F )
{
  //int ni = m_ilast - m_ifirst + 1,
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  for( localIndex ind = 0 ; ind < 3 * m_npts ; ind++ )
    F[ind] = 0;
  for( unsigned int s = 0 ; s < m_point_sources.size() ; s++ )
  {
    double fxyz[3];
    m_point_sources[s]->getFxyz( time, fxyz );
    size_t ind = nj * nk * ( m_point_sources[s]->m_i0 - m_ifirst ) + nk * ( m_point_sources[s]->m_j0 - m_jfirst )
        + m_point_sources[s]->m_k0 - m_kfirst;
    F[3 * ind] += fxyz[0];
    F[3 * ind + 1] += fxyz[1];
    F[3 * ind + 2] += fxyz[2];
    //      if( fxyz[0] != 0 || fxyz[1] != 0 || fxyz[2] != 0 )
    //       cout << "source " << s << " pos= " << m_point_sources[s]->m_i0 << " " <<
    //        m_point_sources[s]->m_j0 << " " << m_point_sources[s]->m_k0 << " F= " << fxyz[0] << " " << fxyz[1] << " " << fxyz[2] << endl;
  }
}

//-----------------------------------------------------------------------
void SW4Solver::Forcing_tt( realT time, realT* F )
{
//  int ni = m_ilast - m_ifirst + 1,
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  for( localIndex ind = 0 ; ind < 3 * m_npts ; ind++ )
    F[ind] = 0;
  for( unsigned int s = 0 ; s < m_point_sources.size() ; s++ )
  {
    double fxyz[3];
    m_point_sources[s]->getFxyztt( time, fxyz );
    size_t ind = nj * nk * ( m_point_sources[s]->m_i0 - m_ifirst ) + nk * ( m_point_sources[s]->m_j0 - m_jfirst )
        + m_point_sources[s]->m_k0 - m_kfirst;
    F[3 * ind] += fxyz[0];
    F[3 * ind + 1] += fxyz[1];
    F[3 * ind + 2] += fxyz[2];
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalPredictor( Array1dT<R1Tensor>& Up, Array1dT<R1Tensor>& U, Array1dT<R1Tensor>& Um,
                               Array1dT<realT>& rho,
                               Array1dT<R1Tensor>& Lu, Array1dT<R1Tensor>& F, realT dt )
{
  for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
  {
    realT coef = dt * dt / rho[ind];
    Up[ind][0] = 2 * U[ind][0] - Um[ind][0] + coef * ( Lu[ind][0] + F[ind][0] );
    Up[ind][1] = 2 * U[ind][1] - Um[ind][1] + coef * ( Lu[ind][1] + F[ind][1] );
    Up[ind][2] = 2 * U[ind][2] - Um[ind][2] + coef * ( Lu[ind][2] + F[ind][2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalPredictor( realT* Up, realT* U, realT* Um, realT* rho,
                               realT* Lu,
                               realT* F, realT dt )
{
  for( localIndex ind = 0 ; ind < m_npts ; ind++ )
  {
    realT coef = dt * dt / rho[ind];
    Up[3 * ind] = 2 * U[3 * ind] - Um[3 * ind] + coef * ( Lu[3 * ind] + F[3 * ind] );
    Up[3 * ind + 1] = 2 * U[3 * ind + 1] - Um[3 * ind + 1] + coef * ( Lu[3 * ind + 1] + F[3 * ind + 1] );
    Up[3 * ind + 2] = 2 * U[3 * ind + 2] - Um[3 * ind + 2] + coef * ( Lu[3 * ind + 2] + F[3 * ind + 2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalDpDmInTime( Array1dT<R1Tensor>& Up, Array1dT<R1Tensor>& U, Array1dT<R1Tensor>& Um,
                                Array1dT<R1Tensor>& Uacc,
                                realT dt )
{
  realT idt2 = 1 / ( dt * dt );
  for( localIndex ind = 0 ; ind < Up.size() ; ind++ )
  {
    Uacc[ind][0] = idt2 * ( Up[ind][0] - 2 * U[ind][0] + Um[ind][0] );
    Uacc[ind][1] = idt2 * ( Up[ind][1] - 2 * U[ind][1] + Um[ind][1] );
    Uacc[ind][2] = idt2 * ( Up[ind][2] - 2 * U[ind][2] + Um[ind][2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalDpDmInTime( realT* Up, realT* U, realT* Um,
                                realT* Uacc,
                                realT dt )
{
  realT idt2 = 1 / ( dt * dt );
  for( localIndex ind = 0 ; ind < m_npts ; ind++ )
  {
    Uacc[3 * ind] = idt2 * ( Up[3 * ind] - 2 * U[3 * ind] + Um[3 * ind] );
    Uacc[3 * ind + 1] = idt2 * ( Up[3 * ind + 1] - 2 * U[3 * ind + 1] + Um[3 * ind + 1] );
    Uacc[3 * ind + 2] = idt2 * ( Up[3 * ind + 2] - 2 * U[3 * ind + 2] + Um[3 * ind + 2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalCorrector( Array1dT<R1Tensor>& Up, Array1dT<realT>& rho, Array1dT<R1Tensor>& Lu,
                               Array1dT<R1Tensor>& F,
                               realT dt )
{
  realT dt4 = dt * dt * dt * dt / 12;
  for( localIndex ind = 0 ; ind < Up.size() ; ind++ )
  {
    realT coef = dt4 / rho[ind];
    Up[ind][0] = Up[ind][0] + coef * ( Lu[ind][0] + F[ind][0] );
    Up[ind][1] = Up[ind][1] + coef * ( Lu[ind][1] + F[ind][1] );
    Up[ind][2] = Up[ind][2] + coef * ( Lu[ind][2] + F[ind][2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::evalCorrector( realT* Up, realT* rho, realT* Lu,
                               realT* F,
                               realT dt )
{
  realT dt4 = dt * dt * dt * dt / 12;
  for( localIndex ind = 0 ; ind < m_npts ; ind++ )
  {
    realT coef = dt4 / rho[ind];
    Up[3 * ind] = Up[3 * ind] + coef * ( Lu[3 * ind] + F[3 * ind] );
    Up[3 * ind + 1] = Up[3 * ind + 1] + coef * ( Lu[3 * ind + 1] + F[3 * ind + 1] );
    Up[3 * ind + 2] = Up[3 * ind + 2] + coef * ( Lu[3 * ind + 2] + F[3 * ind + 2] );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::enforceBC( Array1dT<R1Tensor>& Up, Array1dT<R2SymTensor>& metric,
                           Array1dT<realT>& lame_mu,
                           Array1dT<realT>& lame_lambda )
{
//  int ni = m_ilast - m_ifirst + 1,
  int nj = m_jlast - m_jfirst + 1, nk = m_klast - m_kfirst + 1;
  for( int side = 0 ; side < 6 ; side++ )
  {
    if( m_bctype[side] == 1 )
    {
      // supergrid boundary, enforce homogeneous Dirichlet conditions
      int kmin = m_kfirst, kmax = m_klast, jmin = m_jfirst, jmax = m_jlast, imin = m_ifirst, imax = m_ilast;
      int nghm = m_nghost - 1;
      if( side == 0 )
        imax = imin + nghm;
      else if( side == 1 )
        imin = imax - nghm;
      else if( side == 2 )
        jmax = jmin + nghm;
      else if( side == 3 )
        jmin = jmax - nghm;
      else if( side == 4 )
        kmax = kmin + nghm;
      else if( side == 5 )
        kmin = kmax - nghm;
      for( int i = imin ; i <= imax ; i++ )
        for( int j = jmin ; j <= jmax ; j++ )
          for( int k = kmin ; k <= kmax ; k++ )
          {
            localIndex ind = nj * nk * ( i - m_ifirst ) + nk * ( j - m_jfirst ) + ( k - m_kfirst );
            Up[ind][0] = Up[ind][1] = Up[ind][2] = 0;
          }
    }
    else if( m_bctype[side] == 2 )
    {
      // Free surface
      if( side != 4 )
        std::cout << "SW4Solver::enforceBC, error: Only the top side can be a free surface" << std::endl;
      else
      {
#define SQR(x) ((x)*(x))
#define u(c,i,j,k) (Up[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
#define mu(i,j,k)  (lame_mu[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define la(i,j,k) (lame_lambda[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define met(c,i,j,k) (metric[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)].Data()[c-1])
#define strx(i) (m_strx[i-m_ifirst])
#define stry(j) (m_stry[j-m_jfirst])
        const realT c1 = 2.0 / 3;
        const realT c2 = -1.0 / 12;
        const realT s0i = 1 / m_sbop[0];
        int k = 1;
        int kl = 1;
        for( int i = m_ifirst + 2 ; i <= m_ilast - 2 ; i++ )
        {
          realT istrx = 1 / strx( i );
          for( int j = m_jfirst + 2 ; j <= m_jlast - 2 ; j++ )
          {
            realT istry = 1 / stry( j );
            // Tangential derivatives
            realT rhs1 =
                // pr
                ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 2, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(1,i+2,j,k) - u( 1, i - 2, j, k ) ) +
                        c1 * ( u(1,i+1,j,k) - u( 1, i - 1, j, k ) ) ) * strx( i ) * istry
                    + mu(i,j,k) * met( 3, i, j, k ) * met( 1, i, j, k ) * (
                        c2 * ( u(2,i+2,j,k) - u( 2, i - 2, j, k ) ) +
                            c1 * ( u(2,i+1,j,k) - u( 2, i - 1, j, k ) ) )
                    + mu(i,j,k) * met( 4, i, j, k ) * met( 1, i, j, k ) * (
                        c2 * ( u(3,i+2,j,k) - u( 3, i - 2, j, k ) ) +
                            c1 * ( u(3,i+1,j,k) - u( 3, i - 1, j, k ) ) ) * istry
                    // qr
                    + mu(i,j,k) * met( 3, i, j, k ) * met( 1, i, j, k ) * (
                        c2 * ( u(1,i,j+2,k) - u( 1, i, j - 2, k ) ) +
                            c1 * ( u(1,i,j+1,k) - u( 1, i, j - 1, k ) ) ) * istrx * stry( j )
                    + la(i,j,k) * met( 2, i, j, k ) * met( 1, i, j, k ) * (
                        c2 * ( u(2,i,j+2,k) - u( 2, i, j - 2, k ) ) +
                            c1 * ( u(2,i,j+1,k) - u( 2, i, j - 1, k ) ) );
            //              - forcing(1,i,j)

            // (v-eq)
            realT rhs2 =
            // pr
            la(i,j,k) * met( 3, i, j, k ) * met( 1, i, j, k ) * (
                c2 * ( u(1,i+2,j,k) - u( 1, i - 2, j, k ) ) +
                    c1 * ( u(1,i+1,j,k) - u( 1, i - 1, j, k ) ) )
                + mu(i,j,k) * met( 2, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(2,i+2,j,k) - u( 2, i - 2, j, k ) ) +
                        c1 * ( u(2,i+1,j,k) - u( 2, i - 1, j, k ) ) ) * strx( i ) * istry
                // qr
                + mu(i,j,k) * met( 2, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(1,i,j+2,k) - u( 1, i, j - 2, k ) ) +
                        c1 * ( u(1,i,j+1,k) - u( 1, i, j - 1, k ) ) )
                + ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 3, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(2,i,j+2,k) - u( 2, i, j - 2, k ) ) +
                        c1 * ( u(2,i,j+1,k) - u( 2, i, j - 1, k ) ) ) * stry( j ) * istrx
                + mu(i,j,k) * met( 4, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(3,i,j+2,k) - u( 3, i, j - 2, k ) ) +
                        c1 * ( u(3,i,j+1,k) - u( 3, i, j - 1, k ) ) ) * istrx;
            //                       - forcing(2,i,j)

            // (w-eq)
            realT rhs3 =
            // pr
            la(i,j,k) * met( 4, i, j, k ) * met( 1, i, j, k ) * (
                c2 * ( u(1,i+2,j,k) - u( 1, i - 2, j, k ) ) +
                    c1 * ( u(1,i+1,j,k) - u( 1, i - 1, j, k ) ) ) * istry
                + mu(i,j,k) * met( 2, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(3,i+2,j,k) - u( 3, i - 2, j, k ) ) +
                        c1 * ( u(3,i+1,j,k) - u( 3, i - 1, j, k ) ) ) * strx( i ) * istry
                // qr
                + mu(i,j,k) * met( 3, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(3,i,j+2,k) - u( 3, i, j - 2, k ) ) +
                        c1 * ( u(3,i,j+1,k) - u( 3, i, j - 1, k ) ) ) * stry( j ) * istrx
                + la(i,j,k) * met( 4, i, j, k ) * met( 1, i, j, k ) * (
                    c2 * ( u(2,i,j+2,k) - u( 2, i, j - 2, k ) ) +
                        c1 * ( u(2,i,j+1,k) - u( 2, i, j - 1, k ) ) ) * istrx;
            //  - forcing(3,i,j)

            // Normal derivatives
            realT ac = strx(i) * istry * SQR( met(2,i,j,k) ) +
            stry(j) * istrx * SQR( met(3,i,j,k) ) +
                istry * istrx * SQR( met(4,i,j,k) );
            realT bc = 1 / ( mu(i,j,k) * ac );
            realT cc = ( mu(i,j,k) + la( i, j, k ) ) / ( 2 * mu( i, j, k ) + la( i, j, k ) ) * bc / ac;

            realT xoysqrt = sqrt( strx(i) * istry );
            realT yoxsqrt = 1 / xoysqrt;
            realT isqrtxy = istrx * xoysqrt;
            realT dc = cc * ( xoysqrt * met( 2, i, j, k ) * rhs1 +
                yoxsqrt * met( 3, i, j, k ) * rhs2 +
                isqrtxy * met( 4, i, j, k ) * rhs3 );

            u(1,i,j,k-kl) = -s0i * ( m_sbop[1] * u( 1, i, j, k ) + m_sbop[2] * u( 1, i, j, k + kl ) +
                m_sbop[3] * u( 1, i, j, k + 2 * kl ) + m_sbop[4] * u( 1, i, j, k + 3 * kl ) + bc * rhs1 -
                dc * met( 2, i, j, k ) * xoysqrt );
            u(2,i,j,k-kl) = -s0i * ( m_sbop[1] * u( 2, i, j, k ) + m_sbop[2] * u( 2, i, j, k + kl ) +
                m_sbop[3] * u( 2, i, j, k + 2 * kl ) + m_sbop[4] * u( 2, i, j, k + 3 * kl ) + bc * rhs2 -
                dc * met( 3, i, j, k ) * yoxsqrt );
            u(3,i,j,k-kl) = -s0i * ( m_sbop[1] * u( 3, i, j, k ) + m_sbop[2] * u( 3, i, j, k + kl ) +
                m_sbop[3] * u( 3, i, j, k + 2 * kl ) + m_sbop[4] * u( 3, i, j, k + 3 * kl ) + bc * rhs3 -
                dc * met( 4, i, j, k ) * isqrtxy );
#undef SQR
#undef u
#undef mu
#undef la
#undef met
#undef strx
#undef stry
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------
//void SW4Solver::addSuperGridDamping( Array1dT<R1Tensor>& Up, Array1dT<R1Tensor>& U,
//				     Array1dT<R1Tensor>& Um, Array1dT<realT>& _rho,
//				     realT* _dcx, realT* _dcy, realT* _strx, realT* _stry,
//				     Array1dT<realT>& _jac, realT* _cox, realT* _coy, realT coeff )
void SW4Solver::addSuperGridDamping( realT* Up, realT* U, realT* Um, realT* _rho,
                                     realT* _dcx,
                                     realT* _dcy, realT* _dcz,
                                     realT* _strx,
                                     realT* _stry, realT* _strz,
                                     realT* _jac,
                                     realT* _cox, realT* _coy, realT* _coz,
                                     realT coeff )
{
//  int ni = m_ilast - m_ifirst + 1;
  int nj = m_jlast - m_jfirst + 1;
  int nk = m_klast - m_kfirst + 1;

//#define up(c,i,j,k) (Up[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
//#define u(c,i,j,k) (U[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
//#define um(c,i,j,k) (Um[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
//#define rho(i,j,k) (_rho[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
//#define jac(i,j,k) (_jac[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define up(c,i,j,k) (Up[c-1+3*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define u(c,i,j,k)   (U[c-1+3*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define um(c,i,j,k) (Um[c-1+3*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define rho(i,j,k) (_rho[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define jac(i,j,k) (_jac[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])

#define strx(i) (_strx[i-m_ifirst])
#define stry(j) (_stry[j-m_jfirst])
#define strz(k) (_strz[k-m_kfirst])
#define dcx(i) (_dcx[i-m_ifirst])
#define dcy(j) (_dcy[j-m_jfirst])
#define dcz(k) (_dcz[k-m_kfirst])
#define cox(i) (_cox[i-m_ifirst])
#define coy(j) (_coy[j-m_jfirst])
#define coz(k) (_coz[k-m_kfirst])

  for( int i = m_ifirst + 2 ; i <= m_ilast - 2 ; i++ )
    for( int j = m_jfirst + 2 ; j <= m_jlast - 2 ; j++ )
      for( int k = m_kfirst + 2 ; k <= m_klast - 2 ; k++ )
      {
        realT irhoj = 1 / ( rho(i,j,k) * jac( i, j, k ) );
        for( int c = 1 ; c <= 3 ; c++ )
        {
          up(c,i,j,k) = up(c,i,j,k) - irhoj * coeff * ( strx(i) * coy( j ) * coz( k ) * (
          // x-differences
          rho(i+1,j,k) * dcx( i + 1 ) * jac( i + 1, j, k ) * (
          u(c,i+2,j,k) - 2 * u( c, i + 1, j, k ) + u( c, i, j, k )
              - ( um(c,i+2,j,k) - 2 * um( c, i + 1, j, k ) + um( c, i, j, k ) ) )
              - 2 * rho( i, j, k ) * dcx( i ) * jac( i, j, k ) * (
              u(c,i+1,j,k) - 2 * u( c, i, j, k ) + u( c, i - 1, j, k )
                  - ( um(c,i+1,j,k) - 2 * um( c, i, j, k ) + um( c, i - 1, j, k ) ) )
              + rho(i-1,j,k) * dcx( i - 1 ) * jac( i - 1, j, k ) * (
              u(c,i, j,k) - 2 * u( c, i - 1, j, k ) + u( c, i - 2, j, k )
                  - ( um(c,i, j,k) - 2 * um( c, i - 1, j, k ) + um( c, i - 2, j, k ) ) ) ) +
          stry(j) * cox( i ) * coz( k ) * (
          // y-differences
          rho(i,j+1,k) * dcy( j + 1 ) * jac( i, j + 1, k ) * (
          u(c,i,j+2,k) - 2 * u( c, i, j + 1, k ) + u( c, i, j, k )
              - ( um(c,i,j+2,k) - 2 * um( c, i, j + 1, k ) + um( c, i, j, k ) ) )
              - 2 * rho( i, j, k ) * dcy( j ) * jac( i, j, k ) * (
              u(c,i,j+1,k) - 2 * u( c, i, j, k ) + u( c, i, j - 1, k )
                  - ( um(c,i,j+1,k) - 2 * um( c, i, j, k ) + um( c, i, j - 1, k ) ) )
              + rho(i,j-1,k) * dcy( j - 1 ) * jac( i, j - 1, k ) * (
              u(c,i,j, k) - 2 * u( c, i, j - 1, k ) + u( c, i, j - 2, k )
                  - ( um(c,i,j, k) - 2 * um( c, i, j - 1, k ) + um( c, i, j - 2, k ) ) ) ) +
          cox(i) * coy( j ) * strz( k ) * (
          //    z-differences
          rho(i,j,k+1) * dcz( k + 1 ) * jac( i, j, k + 1 ) * (
          u(c,i,j,k+2) - 2 * u( c, i, j, k + 1 ) + u( c, i, j, k )
              - ( um(c,i,j,k+2) - 2 * um( c, i, j, k + 1 ) + um( c, i, j, k ) ) )
              - 2 * rho( i, j, k ) * dcz( k ) * jac( i, j, k ) * (
              u(c,i,j,k+1) - 2 * u( c, i, j, k ) + u( c, i, j, k - 1 )
                  - ( um(c,i,j,k+1) - 2 * um( c, i, j, k ) + um( c, i, j, k - 1 ) ) )
              + rho(i,j,k-1) * dcz( k - 1 ) * jac( i, j, k - 1 ) * (
              u(c,i,j, k) - 2 * u( c, i, j, k - 1 ) + u( c, i, j, k - 2 )
                  - ( um(c,i,j, k) - 2 * um( c, i, j, k - 1 ) + um( c, i, j, k - 2 ) ) ) ) );
        }
      }
#undef up
#undef u
#undef um
#undef rho
#undef jac
#undef strx
#undef stry
#undef strz
#undef dcx
#undef dcy
#undef dcz
#undef cox
#undef coy
#undef coz
}

//-----------------------------------------------------------------------
void SW4Solver::setupSupergrid( realT* _strx, realT* _stry, realT* _strz,
                                realT* _dcx,
                                realT* _dcy, realT* _dcz,
                                realT* _cox,
                                realT* _coy, realT* _coz )
{
#define strx(i) (_strx[i-m_ifirst])
#define stry(j) (_stry[j-m_jfirst])
#define strz(k) (_strz[k-m_kfirst])
#define dcx(i) (_dcx[i-m_ifirst])
#define dcy(j) (_dcy[j-m_jfirst])
#define dcz(k) (_dcz[k-m_kfirst])
#define cox(i) (_cox[i-m_ifirst])
#define coy(j) (_coy[j-m_jfirst])
#define coz(k) (_coz[k-m_kfirst])
  for( int i = m_ifirst ; i <= m_ilast ; i++ )
  {
    realT x = ( i - 1 ) * m_dx + m_xmin_global;
    dcx(i) = m_supergrid_taper_x.dampingCoeff( x );
    strx(i) = m_supergrid_taper_x.stretching( x );
    cox(i) = m_supergrid_taper_x.cornerTaper( x );
    //      cout << "sg x: i= " << i << " "  << strx(i) << " " << dcx(i) << " " << cox(i) << endl;
  }
  for( int j = m_jfirst ; j <= m_jlast ; j++ )
  {
    realT y = ( j - 1 ) * m_dx + m_ymin_global;
    dcy(j) = m_supergrid_taper_y.dampingCoeff( y );
    stry(j) = m_supergrid_taper_y.stretching( y );
    coy(j) = m_supergrid_taper_y.cornerTaper( y );
    //      cout << "sg y: j= " << j << " "  << stry(j) << " " << dcy(j) << " " << coy(j) << endl;
  }
  for( int k = m_kfirst ; k <= m_klast ; k++ )
  {
    realT z = ( k - 1 ) * m_dx + m_zmin_global;
    dcz(k) = m_supergrid_taper_z.dampingCoeff( z );
    strz(k) = m_supergrid_taper_z.stretching( z );
    coz(k) = m_supergrid_taper_z.cornerTaper( z );
    //      cout << "sg z: k= " << k << " "  << strz(k) << " " << dcz(k) << " " << coz(k) << endl;
  }
#undef strx
#undef stry
#undef strz
#undef dcx
#undef dcy
#undef dcz
#undef cox
#undef coy
#undef coz   
}

//-----------------------------------------------------------------------
void SW4Solver::cycleSolutionArrays( Array1dT<R1Tensor>& Up, Array1dT<R1Tensor>& U, Array1dT<R1Tensor>& Um )
{
  // There must be a better way to do this, copying pointers instead of moving the data.

  // Slow:
  //   for( localIndex ind = 0 ; ind < Up.size() ; ind++ )
  //      for( int c=0 ; c < 3 ; c++ )
  //      {
  //	 Um[ind][c] = U[ind][c];
  //	 U[ind][c]  = Up[ind][c];
  //      }

  // Slower:
  //   Um = U;
  //   U  = Up;

  // Fastest:
  realT* up_ptr = &Up[0][0];
  realT* u_ptr = &U[0][0];
  realT* um_ptr = &Um[0][0];
  for( localIndex ind = 0 ; ind < 3 * m_npts ; ind++ )
  {
    um_ptr[ind] = u_ptr[ind];
    u_ptr[ind] = up_ptr[ind];
  }
}

//-----------------------------------------------------------------------
void SW4Solver::preProcessSources( realT dt, int numberOfTimeSteps, realT tstart )
{

// Need to set the frequency to 1/dt for Dirac source
  for( unsigned int s = 0 ; s < m_sources.size() ; s++ )
    if( m_sources[s]->getTfunc() == iDirac )
      m_sources[s]->setFrequency( 1.0 / dt );

// Filter the sources
  if( m_source_filter != 0 )
  {
    // tell the filter about the time step and compute the second order sections
    m_source_filter->computeSOS( dt );

// 1. Make sure the smallest time offset is at least t0_min + (timeFcn dependent offset for centered fcn's)
    realT dt0 = 0;
    realT dt0loc, dt0max, t0_min;
    t0_min = m_source_filter->estimatePrecursor();
    for( unsigned int s = 0 ; s < m_sources.size() ; s++ )
    {
      dt0loc = m_sources[s]->compute_t0_increase( t0_min );
      if( dt0loc > dt0 )
        dt0 = dt0loc;
    }
    dt0max = dt0;
// If dt0max is positive, the t0 field in all source commands should be incremented 
// by at least this amount. Otherwise, there might be significant artifacts from 
// a sudden start of some source.
    if( dt0max > 0. )
    {
// Warn the user of potential transients due to unsmooth start
      std::cout << "\n*** WARNING: the 2 pass prefilter has an estimated precursor of length " << t0_min << "s\n" <<
                "*** To avoid artifacts due to sudden startup, increase t0 in all source commands by at least "
                << dt0max << endl;
    }
// Do the filtering
    for( unsigned int s = 0 ; s < m_sources.size() ; s++ )
      m_sources[s]->filter_timefunc( m_source_filter, tstart, dt, numberOfTimeSteps );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::setupCartesianMetric( Array1dT<R2SymTensor>& metric, Array1dT<realT>& jacobian )
{
  realT dx3 = m_dx * m_dx * m_dx;
  realT sqrtdx = sqrt( m_dx );
  for( localIndex ind = 0 ; ind < jacobian.size() ; ind++ )
  {
    jacobian[ind] = dx3;
    metric[ind].Data()[0] = sqrtdx;
    metric[ind].Data()[1] = 0;
    metric[ind].Data()[2] = 0;
    metric[ind].Data()[3] = sqrtdx;
  }
}

//-----------------------------------------------------------------------
bool SW4Solver::startswith( const char begin[], char *line )
{
  int lenb = strlen( begin );

  // We ignore any preceeding whitespace
  while( strncmp( line, " ", 1 ) == 0 || strncmp( line, "\t", 1 ) == 0 )
    line++;

  if( strncmp( begin, line, lenb ) == 0 )
    return true;
  else
    return false;
}

//-----------------------------------------------------------------------
void SW4Solver::parseInputFile( string fname )
{
  std::ifstream inputFile( fname.c_str() );
  if( inputFile.is_open() )
  {
    char buffer[1024];
    while( !inputFile.eof() )
    {
      inputFile.getline( buffer, 1024 );
      if( strlen( buffer ) > 0 ) // empty lines produce this
      {
        if( startswith( "source", buffer ) )
          processSource( buffer );
        else if( startswith( "block", buffer ) )
          processMaterialBlock( buffer );
        else if( startswith( "pfile", buffer ) )
          processMaterialPfile( buffer );
        else if( startswith( "supergrid", buffer ) )
          processSupergrid( buffer );
        else if( startswith( "prefilter", buffer ) )
          processPrefilter( buffer );
        else if( startswith( "rec", buffer ) )
          processTimeSeries( buffer );
        else if( !( startswith( "#", buffer ) ||
            startswith( "\n", buffer ) || startswith( "\r", buffer ) ) )
        {
          std::cout << "ERROR, command " << buffer << " not recognized " << std::endl;
        }
      }
    }
    inputFile.close();
  }
  else
    std::cout << "ERROR, could not open file " << fname << std::endl;
}

//-----------------------------------------------------------------------
void SW4Solver::processSource( char* buffer )
{

  // Start time and frequency
  double t0 = 0.0, freq = 1.0;

  // Position
  double x = 0.0, y = 0.0, z = 0.0;

  // Moment tensor
  double m0 = 1.0, mxx = 0.0, mxy = 0.0, mxz = 0.0, myy = 0.0, myz = 0.0, mzz = 0.0;

  // Point force
  double f0 = 1.0, fx = 0.0, fy = 0.0, fz = 0.0;

  // Type of source
  //   bool strikeDipRake = false;
  bool isMomentType = false;
  bool isForceType = false;

  // Time function
  timeDep tDep = iRickerInt;
  char formstring[100];
  char dfile[100];
  strcpy( formstring, "Ricker" );
  int ncyc = 0;
  bool ncyc_set = false;
  bool dfileset = false;

  char* token = strtok( buffer, " \t" );
  token = strtok( NULL, " \t" );

  while( token != NULL )
  {
    // while there are tokens in the string still
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      // Ignore commented lines and lines with just a space.
      break;
    if( startswith( "m0=", token ) )
    {
      token += 3; // skip m0=
      m0 = atof( token );
    }
    else if( startswith( "x=", token ) )
    {
      token += 2; // skip x=
      x = atof( token );
    }
    else if( startswith( "y=", token ) )
    {
      token += 2; // skip y=
      y = atof( token );
    }
    else if( startswith( "z=", token ) )
    {
      token += 2; // skip z=
      z = atof( token );
    }
    else if( startswith( "Mxx=", token ) || startswith( "mxx=", token ) )
    {
      token += 4; // skip Mxx=
      mxx = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Mxy=", token ) || startswith( "mxy=", token ) )
    {
      token += 4; // skip Mxy=
      mxy = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Mxz=", token ) || startswith( "mxz=", token ) )
    {
      token += 4; // skip Mxz=
      mxz = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Myy=", token ) || startswith( "myy=", token ) )
    {
      token += 4; // skip Myy=
      myy = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Myz=", token ) || startswith( "myz=", token ) )
    {
      token += 4; // skip Myz=
      myz = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Mzz=", token ) || startswith( "mzz=", token ) )
    {
      token += 4; // skip Mzz=
      mzz = atof( token );
      isMomentType = true;
    }
    else if( startswith( "Fz=", token ) || startswith( "fz=", token ) )
    {
      token += 3; // skip Fz=
      fz = atof( token );
      isForceType = true;
    }
    else if( startswith( "Fx=", token ) || startswith( "fx=", token ) )
    {
      token += 3; // skip Fx=
      fx = atof( token );
      isForceType = true;
    }
    else if( startswith( "Fy=", token ) || startswith( "fy=", token ) )
    {
      token += 3; // skip Fy=
      fy = atof( token );
      isForceType = true;
    }
    //      else if (startswith("Rake=", token) || startswith("rake=", token))
    //      {
    //         token += 5; // skip Rake=
    //         rake = atof(token);
    //	 strikeDipRake = true;
    //         isMomentType = true;
    //      }
    //      else if (startswith("Strike=", token) || startswith("strike=", token))
    //      {
    //         token += 7; // skip Strike=
    //         strike = atof(token);
    //	 strikeDipRake = true;
    //         isMomentType = true;
    //      }
    //      else if (startswith("Dip=", token) || startswith("dip=", token))
    //      {
    //         token += 4; // skip Dip=
    //         dip = atof(token);
    //	 strikeDipRake = true;
    //         isMomentType = true;
    //      }
    else if( startswith( "t0=", token ) )
    {
      token += 3; // skip t0=
      t0 = atof( token );
    }
    else if( startswith( "freq=", token ) )
    {
      token += 5; // skip freq=
      freq = atof( token );
    }
    else if( startswith( "f0=", token ) )
    {
      token += strlen( "f0=" );
      f0 = atof( token );
    }
    else if( startswith( "type=", token ) )
    {
      token += 5;
      strncpy( formstring, token, 100 );
      if( !strcmp( "Ricker", formstring ) )
        tDep = iRicker;
      else if( !strcmp( "Gaussian", formstring ) )
        tDep = iGaussian;
      else if( !strcmp( "Ramp", formstring ) )
        tDep = iRamp;
      else if( !strcmp( "Triangle", formstring ) )
        tDep = iTriangle;
      else if( !strcmp( "Sawtooth", formstring ) )
        tDep = iSawtooth;
      else if( !strcmp( "SmoothWave", formstring ) )
        tDep = iSmoothWave;
      else if( !strcmp( "Erf", formstring ) || !strcmp( "GaussianInt", formstring ) )
        tDep = iErf;
      else if( !strcmp( "VerySmoothBump", formstring ) )
        tDep = iVerySmoothBump;
      else if( !strcmp( "RickerInt", formstring ) )
        tDep = iRickerInt;
      else if( !strcmp( "Brune", formstring ) )
        tDep = iBrune;
      else if( !strcmp( "BruneSmoothed", formstring ) )
        tDep = iBruneSmoothed;
      else if( !strcmp( "DBrune", formstring ) )
        tDep = iDBrune;
      else if( !strcmp( "GaussianWindow", formstring ) )
        tDep = iGaussianWindow;
      else if( !strcmp( "Liu", formstring ) )
        tDep = iLiu;
      else if( !strcmp( "Dirac", formstring ) )
        tDep = iDirac;
      else if( !strcmp( "C6SmoothBump", formstring ) )
        tDep = iC6SmoothBump;
      else
        std::cout << "unknown time function: " << formstring << endl << " using default RickerInt function."
                  << std::endl;
    }
    else if( startswith( "ncyc=", token ) )
    {
      token += 5; // skip ncyc=
      ncyc = atoi( token );
      ncyc_set = true;
    }
    else if( startswith( "dfile=", token ) )
    {
      token += 6;
      strncpy( dfile, token, 100 );
      dfileset = true;
    }
    else
    {
      std::cout << "SW4Solver: unknown source command option " << token << std::endl;
    }
    token = strtok( NULL, " \t" );
  }

  if( tDep == iGaussianWindow && !ncyc_set )
    std::cout << "SW4Solver: source command: ncyc must be set for Gaussian Window function" << std::endl;

  // Discrete source time function
  double* par = NULL;
  int* ipar = NULL;
  int npar = 0, nipar = 0;
  if( dfileset )
  {
    tDep = iDiscrete;
    //  g(t) defined by spline points on a uniform grid, read from file.
    //  Format: t0, dt, npts
    //          g_1
    //          g_2
    //         ....
    FILE* fd = fopen( dfile, "r" );
    if( fd == NULL )
      std::cout << "SW4Solver: ERROR could not open file " << dfile << std::endl;
    else
    {
      double t0file, dt;
      int npts;
      fscanf( fd, " %lg %lg %i", &t0file, &dt, &npts );
      par = new double[npts + 1];
      par[0] = t0file;
      freq = 1 / dt;
      ipar = new int[1];
      ipar[0] = npts;
      for( int i = 0 ; i < npts ; i++ )
        fscanf( fd, "%lg", &par[i + 1] );
      npar = npts + 1;
      nipar = 1;
      fclose( fd );
    }
  }

  // Useful test, do not know how to find global max now, uncomment later
  //  if ( x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || 
  //       z < zmin || z > m_global_zmax  ) 
  //  {
  //     stringstream sourceposerr;
  //     sourceposerr << endl
  //		  << "***************************************************" << endl
  //		  << " FATAL ERROR:  Source positioned outside grid!  " << endl
  //		  << endl
  //		  << " Source Type: " << formstring << endl
  //		  << "              @ x=" << x 
  //		  << " y=" << y << " z=" << z << endl 
  //		  << endl;

  //     if ( x < xmin )
  //	sourceposerr << " x is " << xmin - x << 
  //	  " meters away from min x (" << xmin << ")" << endl;
  //     else if ( x > m_global_xmax)
  //	sourceposerr << " x is " << x - m_global_xmax << 
  //	  " meters away from max x (" << m_global_xmax << ")" << endl;
  //     if ( y < ymin )
  //	sourceposerr << " y is " << ymin - y << 
  //	  " meters away from min y (" << ymin << ")" << endl;
  //     else if ( y > m_global_ymax)
  //	sourceposerr << " y is " << y - m_global_ymax << 
  //	  " meters away from max y (" << m_global_ymax << ")" << endl;
  //     if ( z < zmin )
  //	sourceposerr << " z is " << zmin - z << 
  //	  " meters away from min z (" << zmin << ")" << endl;
  //     else if ( z > m_global_zmax)
  //	sourceposerr << " z is " << z - m_global_zmax << 
  //	  " meters away from max z (" << m_global_zmax << ")" << endl;
  //     sourceposerr << "***************************************************" << endl;
  //     if (m_myRank == 0)
  //	cout << sourceposerr.str();
  //  }

  // if strike, dip and rake have been given we need to convert into M_{ij} form
  //   if ( strikeDipRake )
  //   {
  //      double radconv = M_PI / 180.;
  //      double S, D, R;
  //      strike -= mGeoAz; // subtract off the grid azimuth
  //      S = strike*radconv; D = dip*radconv; R = rake*radconv;

  //      mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) );
  //      myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) );
  //      mzz = -1.0 * ( mxx + myy );
  //      mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );
  //      mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );
  //      myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );
  //   }

  if( isMomentType && isForceType )
  {
    std::cout << "SW4Solver: can not set both moment and point force type source. Using moment source. " << std::endl;
  }

  if( isMomentType )
  {
    // Remove amplitude variable
    mxx *= m0;
    mxy *= m0;
    mxz *= m0;
    myy *= m0;
    myz *= m0;
    mzz *= m0;
    bool topodepth = false;
    // these have global location since they will be used by all processors
    Source *sourcePtr = new Source( freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
                                    tDep,
                                    formstring, topodepth, ncyc, par, npar, ipar, nipar );

    if( sourcePtr->ignore() )
      delete sourcePtr;
    else
      m_sources.push_back( sourcePtr );
  }
  else // point force
  {
    // Remove amplitude variable
    fx *= f0;
    fy *= f0;
    fz *= f0;
    // global version (gets real coordinates)
    bool topodepth = false;
    Source *sourcePtr = new Source( freq, t0, x, y, z, fx, fy, fz, tDep, formstring, topodepth, ncyc,
                                    par,
                                    npar, ipar, nipar );
    if( sourcePtr->ignore() )
      delete sourcePtr;
    else
      m_sources.push_back( sourcePtr );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::processMaterialBlock( char* buffer )
{
  double vpgrad = 0.0, vsgrad = 0.0, rhograd = 0.0;

  double x1 = -m_xmax_global, x2 = 2 * m_xmax_global;
  double y1 = -m_ymax_global, y2 = 2 * m_ymax_global;
  double z1 = -m_zmax_global, z2 = 2 * m_zmax_global;

  double vp = -1, vs = -1, rho = -1, qp = -1, qs = -1, freq = 1;

  char* token = strtok( buffer, " \t" );
  token = strtok( NULL, " \t" );

  while( token != NULL )
  {
    // while there are tokens in the string still
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      // Ignore commented lines and lines with just a space.
      break;
// the xygrad keywords must occur before the corresponding xy keywords
    if( startswith( "rhograd=", token ) )
    {
      token += 8; // skip rhograd=
      rhograd = atof( token );
    }
    else if( startswith( "vpgrad=", token ) )
    {
      token += 7; // skip vpgrad=
      vpgrad = atof( token );
    }
    else if( startswith( "vsgrad=", token ) )
    {
      token += 7; // skip vsgrad=
      vsgrad = atof( token );
    }
    else if( startswith( "vp=", token ) )
    {
      token += 3; // skip vp=
      vp = atof( token );
    }
    else if( startswith( "vs=", token ) )
    {
      token += 3; // skip vs=
      vs = atof( token );
    }
    else if( startswith( "rho=", token ) )
    {
      token += 4; // skip rho=
      rho = atof( token );
    }
    //      else if (startswith("Qs=", token) || startswith("qs=",token) )
    //      {
    //         token += 3; // skip qs=
    //         qs = atof(token);
    //      }
    //      else if (startswith("Qp=", token) || startswith("qp=",token) )
    //      {
    //         token += 3; // skip qp=
    //         qp = atof(token);
    //      }
    else if( startswith( "x1=", token ) )
    {
      token += 3; // skip x1=
      x1 = atof( token );
    }
    else if( startswith( "x2=", token ) )
    {
      token += 3; // skip x2=
      x2 = atof( token );
    }
    else if( startswith( "y1=", token ) )
    {
      token += 3; // skip y1=
      y1 = atof( token );
    }
    else if( startswith( "y2=", token ) )
    {
      token += 3; // skip y2=
      y2 = atof( token );
    }
    else if( startswith( "z1=", token ) )
    {
      token += 3; // skip z1=
      z1 = atof( token );
    }
    else if( startswith( "z2=", token ) )
    {
      token += 3; // skip z2=
      z2 = atof( token );
    }
    else
    {
      std::cout << "SW4Solver: unknown block command option " << token << std::endl;
    }
    token = strtok( NULL, " \t" );
  }
  if( x1 > m_xmax_global || x2 < m_xmin_global || y1 > m_ymax_global || y2 < m_ymin_global
      || z1 > m_zmax_global || z2 < m_zmin_global )
    std::cout << "SW4Solver: specified material block is outside the domain. Dimensions " <<
              x1
              << " < x < " << x2 << " , " << y1 << " < y < " << y2 << " , " << z1 << " < z < " << z2 << std::endl;
  if( x1 > x2 || y1 > y2 || z1 > z2 )
    std::cout << "SW4Solver: specified material block is empty. Dimensions " <<
              x1
              << " < x < " << x2 << " , " << y1 << " < y < " << y2 << " , " << z1 << " < z < " << z2 << std::endl;
  if( vs <= 0 || vp <= 0 || rho <= 0 )
    std::cout << "SW4Solver: error in material block vp vs rho are   "
              << vp
              << " " << vs << " " << rho << std::endl;

  MaterialBlock* bl = new MaterialBlock( rho, vs, vp, x1, x2, y1, y2, z1, z2, qs, qp, freq );
  bl->set_gradients( rhograd, vsgrad, vpgrad );
  m_material.push_back( bl );
}

//-----------------------------------------------------------------------
void SW4Solver::processMaterialPfile( char * buffer )
{
  string name = "pfile";

  // Used for pfiles
  string filename = "NONE";
  string directory = "NONE";
//  double a_ppm = 0.,
  double vpmin_ppm = 0., vsmin_ppm = 0, rhomin_ppm = 0.;
  //   string cflatten = "NONE";
  //   bool flatten = false;
  //   bool coords_geographic = true;
  int nstenc = 5;

  char* token = strtok( buffer, " \t" );
  //   CHECK_INPUT(strcmp("pfile", token) == 0,
  //	      "ERROR: material data can only be set by an pfile line, not: " << token);

  string err = token;
  err += " Error: ";
  token = strtok( NULL, " \t" );
  while( token != NULL )
  {
    // while there are tokens in the string still
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      // Ignore commented lines and lines with just a space.
      break;
    //      else if (startswith("a=", token))
    //      {
    //         token += 2; // skip a=
    //         a_ppm = atof(token);
    //      }
    else if( startswith( "smoothingsize=", token ) )
    {
      token += 14;
      nstenc = atoi( token );
      if( nstenc < 1 )
        cout << "processMaterialPfile Error: nstenc is " << nstenc << "but should be >= 1\n";
    }
    else if( startswith( "vpmin=", token ) )
    {
      token += 6; // skip vpmin=
      vpmin_ppm = atof( token );
    }
    else if( startswith( "vsmin=", token ) )
    {
      token += 6; // skip vsmin=
      vsmin_ppm = atof( token );
    }
    else if( startswith( "rhomin=", token ) )
    {
      token += 7; // skip rhomin=
      rhomin_ppm = atof( token );
    }
    //      else if (startswith("flatten=", token))
    //      {
    //	token += 8; // skip flatten=
    //	cflatten = token;
    //	VERIFY2( (int)cflatten.find('T')>=0 || (int)cflatten.find('t')>=0 ||
    //		 (int)cflatten.find('F')>=0 || (int)cflatten.find('f')>=0,
    //		 "processMaterialPfile Error: value of flatten unclear\n" );
    //	if ((int)cflatten.find('T')>=0||(int)cflatten.find('t')>=0)
    //	  flatten=true;
    //	else if ((int)cflatten.find('F')>=0||(int)cflatten.find('f')>=0)
    //	  flatten=false;
    //	else
    //	  flatten=false;
    //      }
    else if( startswith( "filename=", token ) )
    {
      token += 9; // skip filename=
      filename = token;
    }
    else if( startswith( "directory=", token ) )
    {
      token += 10; // skip directory=
      directory = token;
    }
    //      else if (startswith("style=", token))
    //      {
    //	token += 6; // skip style=
    //	if( strcmp(token,"geographic") == 0 || strcmp(token,"Geographic")==0 )
    //	  coords_geographic = true;
    //	else if( strcmp(token,"cartesian") == 0 || strcmp(token,"Cartesian")==0 )
    //	  coords_geographic = false;
    //      else
    //	  CHECK_INPUT( false, "processMaterialPfile Error: style= " << token << " not recognized\n" );
    //      }
    else
    {
      cout << token << " is not a pfile option " << endl;
    }
    token = strtok( NULL, " \t" );
  }
  // End parsing...

  //----------------------------------------------------------------
  // Check parameters
  //----------------------------------------------------------------
  if( strcmp( directory.c_str(), "NONE" ) == 0 )
    directory = string( "./" );

  //   if (m_myRank == 0)
  //      cout << "*** Reading data from Pfile " << filename << " in directory " << directory << endl;
  MaterialPfile* pf = new MaterialPfile( filename, directory, nstenc, vpmin_ppm, vsmin_ppm,
                                         rhomin_ppm );
  m_material.push_back( pf );
}

//-----------------------------------------------------------------------
void SW4Solver::processSupergrid( char * buffer )
{
  char* token = strtok( buffer, " \t" );
  token = strtok( NULL, " \t" );
  while( token != NULL )
  {
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      // Ignore commented lines and lines with just a space.
      break;
    if( startswith( "gp=", token ) ) // in number of grid points
    {
      token += 3;
      m_sg_ptsinlayer = atoi( token );
    }
    else if( startswith( "dc=", token ) )
    {
      token += 3;
      m_sg_damping_coefficient = atof( token );
    }
    else
    {
      std::cout << "SW4Solver: unknown supergrid command option " << token << std::endl;
    }
    token = strtok( NULL, " \t" );
  }
}

//-----------------------------------------------------------------------
void SW4Solver::processPrefilter( char* buffer )
{
  char* token = strtok( buffer, " \t" );
  token = strtok( NULL, " \t" );
  double fc1 = 0.1, fc2 = 1.0; // corner frequencies, only fc2 is used for low-pass
  FilterType passband = bandPass; //
  int passes = 2; // forwards and backwards gives a zero-phase filter
  int order = 2;

  while( token != NULL )
  {
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      break;
    if( startswith( "fc1=", token ) )
    {
      token += 4;
      fc1 = atof( token );
    }
    else if( startswith( "fc2=", token ) )
    {
      token += 4;
      fc2 = atof( token );
    }
    else if( startswith( "type=", token ) )
    {
      token += 5;
      if( strcmp( token, "lowpass" ) == 0 )
        passband = lowPass;
      else if( strcmp( token, "bandpass" ) == 0 )
        passband = bandPass;
      else
        std::cout << "SW4Solver, processPrefilter: Error: type= " << token <<
                  " Only lowpass or bandpass are recognized"
                  << std::endl;
    }
    else if( startswith( "passes=", token ) )
    {
      token += 7;
      passes = atoi( token );
      if( !( passes == 1 || passes == 2 ) )
        std::cout << "SW4Solver, processPrefilter: Error: passes must be 1 or 2, not = "
                  << token
                  << std::endl;
    }
    else if( startswith( "order=", token ) )
    {
      token += 6;
      order = atoi( token );
      if( !( order > 0 && order <= 10 ) )
        std::cout << "SW4Solver, processPrefilter: Error: order = "
                  << token
                  << " out of bounds" << std::endl;
    }
    else
    {
      std::cout << "SW4Solver: unknown prefilter command option " << token << std::endl;
    }
    token = strtok( NULL, " \t" );
  }
  m_source_filter = new Filter( passband, order, passes, fc1, fc2 );
}

//-----------------------------------------------------------------------
void SW4Solver::convert_to_mulambda( Array1dT<realT>& rho, Array1dT<realT>& mu,
                                     Array1dT<realT>& lambda )
{
  // On input, we have stored cs in MU, cp in Lambda
  // use mu = rho*cs*cs and lambda = rho*cp*cp  - 2*mu
  for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
  {
    mu[ind] = rho[ind] * mu[ind] * mu[ind];
    lambda[ind] = rho[ind] * lambda[ind] * lambda[ind] - 2 * mu[ind];
  }
}

//-----------------------------------------------------------------------
void SW4Solver::processTimeSeries( char* buffer )
{
  double xs = 0.0, ys = 0.0, zs = 0.0;
  string fileName = "station";
  string staName = "station";
  bool staNameGiven = false;
  int writeEvery = 1000;
  TimeSeries::receiverMode mode = TimeSeries::Displacement;

  char* token = strtok( buffer, " \t" );
  token = strtok( NULL, " \t" );
  while( token != NULL )
  {
    if( startswith( "#", token ) || startswith( " ", buffer ) )
      // Ignore commented lines and lines with just a space.
      break;
    if( startswith( "x=", token ) )
    {
      token += 2; // skip x=
      xs = atof( token );
    }
    else if( startswith( "y=", token ) )
    {
      token += 2; // skip y=
      ys = atof( token );
    }
    else if( startswith( "z=", token ) )
    {
      token += 2; // skip z=
      zs = atof( token );
    }
    else if( startswith( "file=", token ) )
    {
      token += 5; // skip file=
      fileName = token;
    }
    else if( startswith( "sta=", token ) )
    {
      token += strlen( "sta=" );
      staName = token;
      staNameGiven = true;
    }
    else if( startswith( "writeEvery=", token ) )
    {
      token += strlen( "writeEvery=" );
      writeEvery = atoi( token );
    }
    else if( startswith( "variables=", token ) )
    {
      token += strlen( "variables=" );
      if( strcmp( "displacement", token ) == 0 )
      {
        mode = TimeSeries::Displacement;
      }
      else if( strcmp( "velocity", token ) == 0 )
      {
        mode = TimeSeries::Velocity;
      }
      else if( strcmp( "div", token ) == 0 )
      {
        mode = TimeSeries::Div;
      }
      else if( strcmp( "curl", token ) == 0 )
      {
        mode = TimeSeries::Curl;
      }
      else if( strcmp( "strains", token ) == 0 )
      {
        mode = TimeSeries::Strains;
      }
      else
      {
        std::cout << "receiver command: variables=" << token << " not understood" << std::endl
                  << "using default mode (displacement)"
                  << std::endl;
        mode = TimeSeries::Displacement;
      }
    }
    else
    {
      std::cout << "SW4Solver: unknown rec command option " << token << std::endl;
    }
    token = strtok( NULL, " \t" );
  }
  if( !staNameGiven )
    staName = fileName;

  if( xs < m_xmin_global || xs > m_xmax_global || ys < m_ymin_global || ys > m_ymax_global || zs < m_zmin_global || zs > m_zmax_global )
  {
    std::cout << "SW4Solver: warning receiver " << fileName << " is outside the computational domain" << std::endl;
    std::cout << "(x,y,z) = " << xs << " " << ys << " " << zs << std::endl;
  }
  else
  {
    // Nearest grid point
    int i = static_cast<int>( round( ( xs - m_xmin_global ) / m_dx + 1 ) );
    int j = static_cast<int>( round( ( ys - m_ymin_global ) / m_dx + 1 ) );
    int k = static_cast<int>( round( ( zs - m_zmin_global ) / m_dx + 1 ) );
    // Point in processor?
    if( m_ifirst_int <= i && i <= m_ilast_int &&
        m_jfirst_int <= j && j <= m_jlast_int &&
        m_kfirst_int <= k && k <= m_klast_int )
    {
      // coordinates of nearest point
      double xg = ( i - 1 ) * m_dx + m_xmin_global;
      double yg = ( j - 1 ) * m_dx + m_ymin_global;
      double zg = ( k - 1 ) * m_dx + m_zmin_global;
      TimeSeries* time_series = new TimeSeries( fileName, staName, mode, xs, ys, zs,
                                                i,
                                                j, k, xg, yg, zg, writeEvery, m_output_path );
      m_stations.push_back( time_series );
    }
  }
}

//-----------------------------------------------------------------------
void SW4Solver::sw4discretization( realT* displacement, realT* lame_mu, realT* lame_lambda,
                                   realT* metric,
                                   realT* jacobian, realT* lhs,
                                   char op,
                                   realT* _strx, realT* _stry, realT* _strz,
                                   realT* _acof,
                                   realT* _ghcof, realT* _bope )
{
//  int ni = m_ilast - m_ifirst + 1;
  int nj = m_jlast - m_jfirst + 1;
  int nk = m_klast - m_kfirst + 1;

  const double c1 = 2.0 / 3;
  const double c2 = -1.0 / 12;
  const double tf = 3.0 / 4;
  const double i6 = 1.0 / 6;
  const double i144 = 1.0 / 144;

  //*** Routine with supergrid stretchings strx and stry. No stretching
  //*** in z, since top is always topography, and bottom always interface
  //*** to a deeper Cartesian grid.

  //   int ifirst, ilast, jfirst, jlast, kfirst, klast;

  // met(1) is sqrt(J)*px = sqrt(J)*qy
  // met(2) is sqrt(J)*rx
  // met(3) is sqrt(J)*ry
  // met(4) is sqrt(J)*rz
  int kcarteff = m_kcart;
  // Array1dT  arrays
  //#define u(c,i,j,k) (displacement[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
  //#define mu(i,j,k) (lame_mu[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
  //#define la(i,j,k) (lame_lambda[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
  //#define met(c,i,j,k) (metric[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)].Data()[c-1])
  //#define jac(i,j,k) (jacobian[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
  //#define lu(c,i,j,k) (lhs[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)][c-1])
  // double* or vector<double>
#define u(c,i,j,k) (displacement[c-1+3*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define mu(i,j,k) (lame_mu[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define la(i,j,k) (lame_lambda[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
  //#define met(c,i,j,k) (metric[c-1+4*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define met(c,i,j,k) (metric[c-1+6*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])
#define jac(i,j,k) (jacobian[k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst)])
#define lu(c,i,j,k) (lhs[c-1+3*(k-m_kfirst+nk*(j-m_jfirst)+nk*nj*(i-m_ifirst))])

#define acof(q,k,m) (_acof[q-1+6*(k-1)+48*(m-1)])
#define ghcof(k) (_ghcof[k-1])
#define bope(q,k) (_bope[q-1+6*(k-1)])
#define strx(i) (_strx[i-m_ifirst])
#define stry(j) (_stry[j-m_jfirst])
#define strz(k) (_strz[k-m_kfirst])
#define SQR(x) ((x)*(x))

  int a1 = -1;
  int sgn = 0;
  if( op == '=' )
  {
    a1 = 0;
    sgn = 1;
  }
  else if( op == '+' )
  {
    a1 = 1;
    sgn = 1;
  }
  else if( op == '-' )
  {
    a1 = 1;
    sgn = -1;
  }
  if( kcarteff >= m_klast - 3 )
    kcarteff = m_klast - 3;
  if( kcarteff <= 4 )
    kcarteff = 5;

// SBP Boundary closure terms
  double cof1, cof2, cof3, cof4, cof5, mux1, mux2, mux3, mux4;
  int kb = m_kfirst + 2;
  if( m_onesided )
  {
    kb = 7;
    for( int i = m_ifirst + 2 ; i <= m_ilast - 2 ; i++ )
      for( int j = m_jfirst + 2 ; j <= m_jlast - 2 ; j++ )
        for( int k = 1 ; k <= 6 ; k++ )
        {
          double ijac = strx(i) * stry( j ) / jac( i, j, k );
          double istry = 1 / ( stry( j ) );
          double istrx = 1 / ( strx( i ) );
          double istrxy = istry * istrx;
          double r1 = 0;
          double r2 = 0;
          double r3 = 0;

// pp derivative (u) (u-eq)
          cof1 = ( 2 * mu( i - 2, j, k ) + la( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx(
              i - 2 );
          cof2 = ( 2 * mu( i - 1, j, k ) + la( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx(
              i - 1 );
          cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
          cof4 = ( 2 * mu( i + 1, j, k ) + la( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx(
              i + 1 );
          cof5 = ( 2 * mu( i + 2, j, k ) + la( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx(
              i + 2 );

          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r1 = r1 + i6 * (
              mux1 * ( u(1,i-2,j,k) - u( 1, i, j, k ) ) +
                  mux2 * ( u(1,i-1,j,k) - u( 1, i, j, k ) ) +
                  mux3 * ( u(1,i+1,j,k) - u( 1, i, j, k ) ) +
                  mux4 * ( u(1,i+2,j,k) - u( 1, i, j, k ) ) ) * istry;

// qq derivative (u) (u-eq)
          cof1 = ( mu( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry( j - 2 );
          cof2 = ( mu( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry( j - 1 );
          cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
          cof4 = ( mu( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry( j + 1 );
          cof5 = ( mu( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry( j + 2 );

          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r1 = r1 + i6 * (
              mux1 * ( u(1,i,j-2,k) - u( 1, i, j, k ) ) +
                  mux2 * ( u(1,i,j-1,k) - u( 1, i, j, k ) ) +
                  mux3 * ( u(1,i,j+1,k) - u( 1, i, j, k ) ) +
                  mux4 * ( u(1,i,j+2,k) - u( 1, i, j, k ) ) ) * istrx;

// pp derivative (v) (v-eq)
          cof1 = ( mu( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx( i - 2 );
          cof2 = ( mu( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx( i - 1 );
          cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
          cof4 = ( mu( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx( i + 1 );
          cof5 = ( mu( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx( i + 2 );

          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r2 = r2 + i6 * (
              mux1 * ( u(2,i-2,j,k) - u( 2, i, j, k ) ) +
                  mux2 * ( u(2,i-1,j,k) - u( 2, i, j, k ) ) +
                  mux3 * ( u(2,i+1,j,k) - u( 2, i, j, k ) ) +
                  mux4 * ( u(2,i+2,j,k) - u( 2, i, j, k ) ) ) * istry;

// qq derivative (v) (v-eq)
          cof1 = ( 2 * mu( i, j - 2, k ) + la( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry(
              j - 2 );
          cof2 = ( 2 * mu( i, j - 1, k ) + la( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry(
              j - 1 );
          cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
          cof4 = ( 2 * mu( i, j + 1, k ) + la( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry(
              j + 1 );
          cof5 = ( 2 * mu( i, j + 2, k ) + la( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry(
              j + 2 );
          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r2 = r2 + i6 * (
              mux1 * ( u(2,i,j-2,k) - u( 2, i, j, k ) ) +
                  mux2 * ( u(2,i,j-1,k) - u( 2, i, j, k ) ) +
                  mux3 * ( u(2,i,j+1,k) - u( 2, i, j, k ) ) +
                  mux4 * ( u(2,i,j+2,k) - u( 2, i, j, k ) ) ) * istrx;

// pp derivative (w) (w-eq)
          cof1 = ( mu( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx( i - 2 );
          cof2 = ( mu( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx( i - 1 );
          cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
          cof4 = ( mu( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx( i + 1 );
          cof5 = ( mu( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx( i + 2 );

          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r3 = r3 + i6 * (
              mux1 * ( u(3,i-2,j,k) - u( 3, i, j, k ) ) +
                  mux2 * ( u(3,i-1,j,k) - u( 3, i, j, k ) ) +
                  mux3 * ( u(3,i+1,j,k) - u( 3, i, j, k ) ) +
                  mux4 * ( u(3,i+2,j,k) - u( 3, i, j, k ) ) ) * istry;

// qq derivative (w) (w-eq)
          cof1 = ( mu( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry( j - 2 );
          cof2 = ( mu( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry( j - 1 );
          cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
          cof4 = ( mu( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry( j + 1 );
          cof5 = ( mu( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry( j + 2 );
          mux1 = cof2 - tf * ( cof3 + cof1 );
          mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
          mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
          mux4 = cof4 - tf * ( cof3 + cof5 );

          r3 = r3 + i6 * (
              mux1 * ( u(3,i,j-2,k) - u( 3, i, j, k ) ) +
                  mux2 * ( u(3,i,j-1,k) - u( 3, i, j, k ) ) +
                  mux3 * ( u(3,i,j+1,k) - u( 3, i, j, k ) ) +
                  mux4 * ( u(3,i,j+2,k) - u( 3, i, j, k ) ) ) * istrx;

// All rr-derivatives at once
// averaging the coefficient
          double mucofu2, mucofuv, mucofuw, mucofvw, mucofv2, mucofw2;
          for( int q = 1 ; q <= 8 ; q++ )
          {
            mucofu2 = 0;
            mucofuv = 0;
            mucofuw = 0;
            mucofvw = 0;
            mucofv2 = 0;
            mucofw2 = 0;
            for( int m = 1 ; m <= 8 ; m++ )
            {
              mucofu2 = mucofu2 +
              acof(k,q,m) * ( ( 2 * mu( i, j, m ) + la( i, j, m ) ) *
                  SQR( met(2,i,j,m)*strx(i) )
                  + mu(i,j,m) *
                      ( SQR(met(3,i,j,m)*stry(j)) + SQR( met(4,i,j,m) ) ) );
              mucofv2 = mucofv2 +
              acof(k,q,m) * ( ( 2 * mu( i, j, m ) + la( i, j, m ) ) *
                  SQR( met(3,i,j,m)*stry(j) )
                  + mu(i,j,m) *
                      ( SQR(met(2,i,j,m)*strx(i)) + SQR( met(4,i,j,m) ) ) );
              mucofw2 = mucofw2 +
              acof(k,q,m) * ( ( 2 * mu( i, j, m ) + la( i, j, m ) ) * SQR( met(4,i,j,m) )
                  + mu(i,j,m) *
                      ( SQR(met(2,i,j,m)*strx(i)) + SQR( met(3,i,j,m)*stry(j) ) ) );
              mucofuv = mucofuv + acof(k,q,m) * ( mu(i,j,m) + la( i, j, m ) ) *
                  met( 2, i, j, m ) * met( 3, i, j, m );
              mucofuw = mucofuw + acof(k,q,m) * ( mu(i,j,m) + la( i, j, m ) ) *
                  met( 2, i, j, m ) * met( 4, i, j, m );
              mucofvw = mucofvw + acof(k,q,m) * ( mu(i,j,m) + la( i, j, m ) ) *
                  met( 3, i, j, m ) * met( 4, i, j, m );
            }
// Computing the second derivative,
            r1 = r1 + istrxy * mucofu2 * u( 1, i, j, q ) + mucofuv * u( 2, i, j, q ) +
                istry * mucofuw * u( 3, i, j, q );
            r2 = r2 + mucofuv * u( 1, i, j, q ) + istrxy * mucofv2 * u( 2, i, j, q ) +
                istrx * mucofvw * u( 3, i, j, q );
            r3 = r3 + istry * mucofuw * u( 1, i, j, q ) +
                istrx * mucofvw * u( 2, i, j, q ) + istrxy * mucofw2 * u( 3, i, j, q );
          }

// Ghost point values, only nonzero for k=1.
          mucofu2 = ghcof(k) * ( ( 2 * mu( i, j, 1 ) + la( i, j, 1 ) ) *
              SQR( met(2,i,j,1)*strx(i) )
              + mu(i,j,1) * ( SQR(met(3,i,j,1)*stry(j)) + SQR( met(4,i,j,1) ) ) );
          mucofv2 = ghcof(k) * ( ( 2 * mu( i, j, 1 ) + la( i, j, 1 ) ) *
              SQR( met(3,i,j,1)*stry(j) )
              + mu(i,j,1) * ( SQR(met(2,i,j,1)*strx(i)) + SQR( met(4,i,j,1) ) ) );
          mucofw2 = ghcof(k) * ( ( 2 * mu( i, j, 1 ) + la( i, j, 1 ) ) * SQR( met(4,i,j,1) )
              + mu(i,j,1) *
                  ( SQR(met(2,i,j,1)*strx(i)) + SQR( met(3,i,j,1)*stry(j) ) ) );
          mucofuv = ghcof(k) * ( mu(i,j,1) + la( i, j, 1 ) ) *
              met( 2, i, j, 1 ) * met( 3, i, j, 1 );
          mucofuw = ghcof(k) * ( mu(i,j,1) + la( i, j, 1 ) ) *
              met( 2, i, j, 1 ) * met( 4, i, j, 1 );
          mucofvw = ghcof(k) * ( mu(i,j,1) + la( i, j, 1 ) ) *
              met( 3, i, j, 1 ) * met( 4, i, j, 1 );
          r1 = r1 + istrxy * mucofu2 * u( 1, i, j, 0 ) + mucofuv * u( 2, i, j, 0 ) +
              istry * mucofuw * u( 3, i, j, 0 );
          r2 = r2 + mucofuv * u( 1, i, j, 0 ) + istrxy * mucofv2 * u( 2, i, j, 0 ) +
              istrx * mucofvw * u( 3, i, j, 0 );
          r3 = r3 + istry * mucofuw * u( 1, i, j, 0 ) +
              istrx * mucofvw * u( 2, i, j, 0 ) + istrxy * mucofw2 * u( 3, i, j, 0 );

// pq-derivatives (u-eq)
          r1 = r1 +
              c2 * ( mu(i,j+2,k) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                  c2 * ( u(2,i+2,j+2,k) - u( 2, i - 2, j + 2, k ) ) +
                      c1 * ( u(2,i+1,j+2,k) - u( 2, i - 1, j + 2, k ) ) )
                  - mu(i,j-2,k) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                      c2 * ( u(2,i+2,j-2,k) - u( 2, i - 2, j - 2, k ) ) +
                          c1 * ( u(2,i+1,j-2,k) - u( 2, i - 1, j - 2, k ) ) )
                  ) +
              c1 * ( mu(i,j+1,k) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                  c2 * ( u(2,i+2,j+1,k) - u( 2, i - 2, j + 1, k ) ) +
                      c1 * ( u(2,i+1,j+1,k) - u( 2, i - 1, j + 1, k ) ) )
                  - mu(i,j-1,k) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                      c2 * ( u(2,i+2,j-1,k) - u( 2, i - 2, j - 1, k ) ) +
                          c1 * ( u(2,i+1,j-1,k) - u( 2, i - 1, j - 1, k ) ) ) );

// qp-derivatives (u-eq)
          r1 = r1 +
              c2 * ( la(i+2,j,k) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                  c2 * ( u(2,i+2,j+2,k) - u( 2, i + 2, j - 2, k ) ) +
                      c1 * ( u(2,i+2,j+1,k) - u( 2, i + 2, j - 1, k ) ) )
                  - la(i-2,j,k) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                      c2 * ( u(2,i-2,j+2,k) - u( 2, i - 2, j - 2, k ) ) +
                          c1 * ( u(2,i-2,j+1,k) - u( 2, i - 2, j - 1, k ) ) )
                  ) +
              c1 * ( la(i+1,j,k) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                  c2 * ( u(2,i+1,j+2,k) - u( 2, i + 1, j - 2, k ) ) +
                      c1 * ( u(2,i+1,j+1,k) - u( 2, i + 1, j - 1, k ) ) )
                  - la(i-1,j,k) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                      c2 * ( u(2,i-1,j+2,k) - u( 2, i - 1, j - 2, k ) ) +
                          c1 * ( u(2,i-1,j+1,k) - u( 2, i - 1, j - 1, k ) ) ) );

// pq-derivatives (v-eq)
          r2 = r2 +
              c2 * ( la(i,j+2,k) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                  c2 * ( u(1,i+2,j+2,k) - u( 1, i - 2, j + 2, k ) ) +
                      c1 * ( u(1,i+1,j+2,k) - u( 1, i - 1, j + 2, k ) ) )
                  - la(i,j-2,k) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                      c2 * ( u(1,i+2,j-2,k) - u( 1, i - 2, j - 2, k ) ) +
                          c1 * ( u(1,i+1,j-2,k) - u( 1, i - 1, j - 2, k ) ) )
                  ) +
              c1 * ( la(i,j+1,k) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                  c2 * ( u(1,i+2,j+1,k) - u( 1, i - 2, j + 1, k ) ) +
                      c1 * ( u(1,i+1,j+1,k) - u( 1, i - 1, j + 1, k ) ) )
                  - la(i,j-1,k) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                      c2 * ( u(1,i+2,j-1,k) - u( 1, i - 2, j - 1, k ) ) +
                          c1 * ( u(1,i+1,j-1,k) - u( 1, i - 1, j - 1, k ) ) ) );

// qp-derivatives (v-eq)
          r2 = r2 +
              c2 * ( mu(i+2,j,k) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                  c2 * ( u(1,i+2,j+2,k) - u( 1, i + 2, j - 2, k ) ) +
                      c1 * ( u(1,i+2,j+1,k) - u( 1, i + 2, j - 1, k ) ) )
                  - mu(i-2,j,k) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                      c2 * ( u(1,i-2,j+2,k) - u( 1, i - 2, j - 2, k ) ) +
                          c1 * ( u(1,i-2,j+1,k) - u( 1, i - 2, j - 1, k ) ) )
                  ) +
              c1 * ( mu(i+1,j,k) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                  c2 * ( u(1,i+1,j+2,k) - u( 1, i + 1, j - 2, k ) ) +
                      c1 * ( u(1,i+1,j+1,k) - u( 1, i + 1, j - 1, k ) ) )
                  - mu(i-1,j,k) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                      c2 * ( u(1,i-1,j+2,k) - u( 1, i - 1, j - 2, k ) ) +
                          c1 * ( u(1,i-1,j+1,k) - u( 1, i - 1, j - 1, k ) ) ) );

// rp - derivatives
          double dudrm2 = 0, dudrm1 = 0, dudrp1 = 0, dudrp2 = 0;
          double dvdrm2 = 0, dvdrm1 = 0, dvdrp1 = 0, dvdrp2 = 0;
          double dwdrm2 = 0, dwdrm1 = 0, dwdrp1 = 0, dwdrp2 = 0;
          for( int q = 1 ; q <= 8 ; q++ )
          {
            dudrm2 = dudrm2 + bope(k,q) * u( 1, i - 2, j, q );
            dvdrm2 = dvdrm2 + bope(k,q) * u( 2, i - 2, j, q );
            dwdrm2 = dwdrm2 + bope(k,q) * u( 3, i - 2, j, q );
            dudrm1 = dudrm1 + bope(k,q) * u( 1, i - 1, j, q );
            dvdrm1 = dvdrm1 + bope(k,q) * u( 2, i - 1, j, q );
            dwdrm1 = dwdrm1 + bope(k,q) * u( 3, i - 1, j, q );
            dudrp2 = dudrp2 + bope(k,q) * u( 1, i + 2, j, q );
            dvdrp2 = dvdrp2 + bope(k,q) * u( 2, i + 2, j, q );
            dwdrp2 = dwdrp2 + bope(k,q) * u( 3, i + 2, j, q );
            dudrp1 = dudrp1 + bope(k,q) * u( 1, i + 1, j, q );
            dvdrp1 = dvdrp1 + bope(k,q) * u( 2, i + 1, j, q );
            dwdrp1 = dwdrp1 + bope(k,q) * u( 3, i + 1, j, q );
          }

// rp derivatives (u-eq)
          r1 = r1 + ( c2 * (
              ( 2 * mu( i + 2, j, k ) + la( i + 2, j, k ) ) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) *
                  strx( i + 2 ) * dudrp2
                  + la(i+2,j,k) * met( 3, i + 2, j, k ) * met( 1, i + 2, j, k ) * dvdrp2 * stry( j )
                  + la(i+2,j,k) * met( 4, i + 2, j, k ) * met( 1, i + 2, j, k ) * dwdrp2
                  - ( ( 2 * mu( i - 2, j, k ) + la( i - 2, j, k ) ) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) *
                      strx( i - 2 ) * dudrm2
                      + la(i-2,j,k) * met( 3, i - 2, j, k ) * met( 1, i - 2, j, k ) * dvdrm2 * stry( j )
                      + la(i-2,j,k) * met( 4, i - 2, j, k ) * met( 1, i - 2, j, k ) * dwdrm2 )
              ) + c1 * (
              ( 2 * mu( i + 1, j, k ) + la( i + 1, j, k ) ) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) *
                  strx( i + 1 ) * dudrp1
                  + la(i+1,j,k) * met( 3, i + 1, j, k ) * met( 1, i + 1, j, k ) * dvdrp1 * stry( j )
                  + la(i+1,j,k) * met( 4, i + 1, j, k ) * met( 1, i + 1, j, k ) * dwdrp1
                  - ( ( 2 * mu( i - 1, j, k ) + la( i - 1, j, k ) ) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) *
                      strx( i - 1 ) * dudrm1
                      + la(i-1,j,k) * met( 3, i - 1, j, k ) * met( 1, i - 1, j, k ) * dvdrm1 * stry( j )
                      + la(i-1,j,k) * met( 4, i - 1, j, k ) * met( 1, i - 1, j, k ) * dwdrm1 ) ) ) * istry;

// rp derivatives (v-eq)
          r2 = r2 + c2 * (
          mu(i+2,j,k) * met( 3, i + 2, j, k ) * met( 1, i + 2, j, k ) * dudrp2
              + mu(i+2,j,k) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) * dvdrp2 *
                  strx( i + 2 ) * istry
              - ( mu(i-2,j,k) * met( 3, i - 2, j, k ) * met( 1, i - 2, j, k ) * dudrm2
                  + mu(i-2,j,k) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) * dvdrm2 *
                      strx( i - 2 ) * istry )
              ) + c1 * (
          mu(i+1,j,k) * met( 3, i + 1, j, k ) * met( 1, i + 1, j, k ) * dudrp1
              + mu(i+1,j,k) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) * dvdrp1 *
                  strx( i + 1 ) * istry
              - ( mu(i-1,j,k) * met( 3, i - 1, j, k ) * met( 1, i - 1, j, k ) * dudrm1
                  + mu(i-1,j,k) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) * dvdrm1 *
                      strx( i - 1 ) * istry )
              );

// rp derivatives (w-eq)
          r3 = r3 + istry * ( c2 * (
          mu(i+2,j,k) * met( 4, i + 2, j, k ) * met( 1, i + 2, j, k ) * dudrp2
              + mu(i+2,j,k) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) * dwdrp2 * strx( i + 2 )
              - ( mu(i-2,j,k) * met( 4, i - 2, j, k ) * met( 1, i - 2, j, k ) * dudrm2
                  + mu(i-2,j,k) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) * dwdrm2 * strx( i - 2 ) )
              ) + c1 * (
          mu(i+1,j,k) * met( 4, i + 1, j, k ) * met( 1, i + 1, j, k ) * dudrp1
              + mu(i+1,j,k) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) * dwdrp1 * strx( i + 1 )
              - ( mu(i-1,j,k) * met( 4, i - 1, j, k ) * met( 1, i - 1, j, k ) * dudrm1
                  + mu(i-1,j,k) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) * dwdrm1 * strx( i - 1 ) )
              ) );

// rq - derivatives
          dudrm2 = dudrm1 = dudrp1 = dudrp2 = 0;
          dvdrm2 = dvdrm1 = dvdrp1 = dvdrp2 = 0;
          dwdrm2 = dwdrm1 = dwdrp1 = dwdrp2 = 0;
          for( int q = 1 ; q <= 8 ; q++ )
          {
            dudrm2 = dudrm2 + bope(k,q) * u( 1, i, j - 2, q );
            dvdrm2 = dvdrm2 + bope(k,q) * u( 2, i, j - 2, q );
            dwdrm2 = dwdrm2 + bope(k,q) * u( 3, i, j - 2, q );
            dudrm1 = dudrm1 + bope(k,q) * u( 1, i, j - 1, q );
            dvdrm1 = dvdrm1 + bope(k,q) * u( 2, i, j - 1, q );
            dwdrm1 = dwdrm1 + bope(k,q) * u( 3, i, j - 1, q );
            dudrp2 = dudrp2 + bope(k,q) * u( 1, i, j + 2, q );
            dvdrp2 = dvdrp2 + bope(k,q) * u( 2, i, j + 2, q );
            dwdrp2 = dwdrp2 + bope(k,q) * u( 3, i, j + 2, q );
            dudrp1 = dudrp1 + bope(k,q) * u( 1, i, j + 1, q );
            dvdrp1 = dvdrp1 + bope(k,q) * u( 2, i, j + 1, q );
            dwdrp1 = dwdrp1 + bope(k,q) * u( 3, i, j + 1, q );
          }

// rq derivatives (u-eq)
          r1 = r1 + c2 * (
          mu(i,j+2,k) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * dudrp2 *
              stry( j + 2 ) * istrx
              + mu(i,j+2,k) * met( 2, i, j + 2, k ) * met( 1, i, j + 2, k ) * dvdrp2
              - ( mu(i,j-2,k) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * dudrm2 *
                  stry( j - 2 ) * istrx
                  + mu(i,j-2,k) * met( 2, i, j - 2, k ) * met( 1, i, j - 2, k ) * dvdrm2 )
              ) + c1 * (
          mu(i,j+1,k) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * dudrp1 *
              stry( j + 1 ) * istrx
              + mu(i,j+1,k) * met( 2, i, j + 1, k ) * met( 1, i, j + 1, k ) * dvdrp1
              - ( mu(i,j-1,k) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * dudrm1 *
                  stry( j - 1 ) * istrx
                  + mu(i,j-1,k) * met( 2, i, j - 1, k ) * met( 1, i, j - 1, k ) * dvdrm1 )
              );

// rq derivatives (v-eq)
          r2 = r2 + c2 * (
          la(i,j+2,k) * met( 2, i, j + 2, k ) * met( 1, i, j + 2, k ) * dudrp2
              + ( 2 * mu( i, j + 2, k ) + la( i, j + 2, k ) ) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * dvdrp2
                  * stry( j + 2 ) * istrx
              + la(i,j+2,k) * met( 4, i, j + 2, k ) * met( 1, i, j + 2, k ) * dwdrp2 * istrx
              - ( la(i,j-2,k) * met( 2, i, j - 2, k ) * met( 1, i, j - 2, k ) * dudrm2
                  + ( 2 * mu( i, j - 2, k ) + la( i, j - 2, k ) ) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * dvdrm2
                      * stry( j - 2 ) * istrx
                  + la(i,j-2,k) * met( 4, i, j - 2, k ) * met( 1, i, j - 2, k ) * dwdrm2 * istrx )
              ) + c1 * (
          la(i,j+1,k) * met( 2, i, j + 1, k ) * met( 1, i, j + 1, k ) * dudrp1
              + ( 2 * mu( i, j + 1, k ) + la( i, j + 1, k ) ) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * dvdrp1
                  * stry( j + 1 ) * istrx
              + la(i,j+1,k) * met( 4, i, j + 1, k ) * met( 1, i, j + 1, k ) * dwdrp1 * istrx
              - ( la(i,j-1,k) * met( 2, i, j - 1, k ) * met( 1, i, j - 1, k ) * dudrm1
                  + ( 2 * mu( i, j - 1, k ) + la( i, j - 1, k ) ) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * dvdrm1
                      * stry( j - 1 ) * istrx
                  + la(i,j-1,k) * met( 4, i, j - 1, k ) * met( 1, i, j - 1, k ) * dwdrm1 * istrx )
              );

// rq derivatives (w-eq)
          r3 = r3 + ( c2 * (
          mu(i,j+2,k) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * dwdrp2 * stry( j + 2 )
              + mu(i,j+2,k) * met( 4, i, j + 2, k ) * met( 1, i, j + 2, k ) * dvdrp2
              - ( mu(i,j-2,k) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * dwdrm2 * stry( j - 2 )
                  + mu(i,j-2,k) * met( 4, i, j - 2, k ) * met( 1, i, j - 2, k ) * dvdrm2 )
              ) + c1 * (
          mu(i,j+1,k) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * dwdrp1 * stry( j + 1 )
              + mu(i,j+1,k) * met( 4, i, j + 1, k ) * met( 1, i, j + 1, k ) * dvdrp1
              - ( mu(i,j-1,k) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * dwdrm1 * stry( j - 1 )
                  + mu(i,j-1,k) * met( 4, i, j - 1, k ) * met( 1, i, j - 1, k ) * dvdrm1 )
              ) ) * istrx;

// pr and qr derivatives at once

          for( int q = 1 ; q <= 8 ; q++ )
          {
// (u-eq)
            r1 = r1 + bope(k,q) * (
                // pr
                ( 2 * mu( i, j, q ) + la( i, j, q ) ) * met( 2, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(1,i+2,j,q) - u( 1, i - 2, j, q ) ) +
                        c1 * ( u(1,i+1,j,q) - u( 1, i - 1, j, q ) ) ) * strx( i ) * istry
                    + mu(i,j,q) * met( 3, i, j, q ) * met( 1, i, j, q ) * (
                        c2 * ( u(2,i+2,j,q) - u( 2, i - 2, j, q ) ) +
                            c1 * ( u(2,i+1,j,q) - u( 2, i - 1, j, q ) ) )
                    + mu(i,j,q) * met( 4, i, j, q ) * met( 1, i, j, q ) * (
                        c2 * ( u(3,i+2,j,q) - u( 3, i - 2, j, q ) ) +
                            c1 * ( u(3,i+1,j,q) - u( 3, i - 1, j, q ) ) ) * istry
                    // qr
                    + mu(i,j,q) * met( 3, i, j, q ) * met( 1, i, j, q ) * (
                        c2 * ( u(1,i,j+2,q) - u( 1, i, j - 2, q ) ) +
                            c1 * ( u(1,i,j+1,q) - u( 1, i, j - 1, q ) ) ) * stry( j ) * istrx
                    + la(i,j,q) * met( 2, i, j, q ) * met( 1, i, j, q ) * (
                        c2 * ( u(2,i,j+2,q) - u( 2, i, j - 2, q ) ) +
                            c1 * ( u(2,i,j+1,q) - u( 2, i, j - 1, q ) ) ) );

// (v-eq)
            r2 = r2 + bope(k,q) * (
            // pr
            la(i,j,q) * met( 3, i, j, q ) * met( 1, i, j, q ) * (
                c2 * ( u(1,i+2,j,q) - u( 1, i - 2, j, q ) ) +
                    c1 * ( u(1,i+1,j,q) - u( 1, i - 1, j, q ) ) )
                + mu(i,j,q) * met( 2, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(2,i+2,j,q) - u( 2, i - 2, j, q ) ) +
                        c1 * ( u(2,i+1,j,q) - u( 2, i - 1, j, q ) ) ) * strx( i ) * istry
                // qr
                + mu(i,j,q) * met( 2, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(1,i,j+2,q) - u( 1, i, j - 2, q ) ) +
                        c1 * ( u(1,i,j+1,q) - u( 1, i, j - 1, q ) ) )
                + ( 2 * mu( i, j, q ) + la( i, j, q ) ) * met( 3, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(2,i,j+2,q) - u( 2, i, j - 2, q ) ) +
                        c1 * ( u(2,i,j+1,q) - u( 2, i, j - 1, q ) ) ) * stry( j ) * istrx
                + mu(i,j,q) * met( 4, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(3,i,j+2,q) - u( 3, i, j - 2, q ) ) +
                        c1 * ( u(3,i,j+1,q) - u( 3, i, j - 1, q ) ) ) * istrx );

// (w-eq)
            r3 = r3 + bope(k,q) * (
            // pr
            la(i,j,q) * met( 4, i, j, q ) * met( 1, i, j, q ) * (
                c2 * ( u(1,i+2,j,q) - u( 1, i - 2, j, q ) ) +
                    c1 * ( u(1,i+1,j,q) - u( 1, i - 1, j, q ) ) ) * istry
                + mu(i,j,q) * met( 2, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(3,i+2,j,q) - u( 3, i - 2, j, q ) ) +
                        c1 * ( u(3,i+1,j,q) - u( 3, i - 1, j, q ) ) ) * strx( i ) * istry
                // qr
                + mu(i,j,q) * met( 3, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(3,i,j+2,q) - u( 3, i, j - 2, q ) ) +
                        c1 * ( u(3,i,j+1,q) - u( 3, i, j - 1, q ) ) ) * stry( j ) * istrx
                + la(i,j,q) * met( 4, i, j, q ) * met( 1, i, j, q ) * (
                    c2 * ( u(2,i,j+2,q) - u( 2, i, j - 2, q ) ) +
                        c1 * ( u(2,i,j+1,q) - u( 2, i, j - 1, q ) ) ) * istrx );

          }
          //         lu(1,i,j,k) = r1*ijac
          lu(1,i,j,k) = a1 * lu( 1, i, j, k ) + sgn * r1 * ijac;
          //         lu(2,i,j,k) = r2*ijac
          lu(2,i,j,k) = a1 * lu( 2, i, j, k ) + sgn * r2 * ijac;
          //          lu(3,i,j,k) = r3*ijac
          lu(3,i,j,k) = a1 * lu( 3, i, j, k ) + sgn * r3 * ijac;
        }
  }

  // Go to kcart+1, since Cartesian computation should
  // start at kcart+2, to ensure all five points in the
  // stencil are in the Cartesian region.
  for( int i = m_ifirst + 2 ; i <= m_ilast - 2 ; i++ )
    for( int j = m_jfirst + 2 ; j <= m_jlast - 2 ; j++ )
      for( int k = kb ; k <= kcarteff + 1 ; k++ )
      {
        double ijac = strx(i) * stry( j ) / jac( i, j, k );
        double istry = 1 / ( stry( j ) );
        double istrx = 1 / ( strx( i ) );
        double istrxy = istry * istrx;

        double r1 = 0;

        // pp derivative (u)
        cof1 = ( 2 * mu( i - 2, j, k ) + la( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx(
            i - 2 );
        cof2 = ( 2 * mu( i - 1, j, k ) + la( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx(
            i - 1 );
        cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
        cof4 = ( 2 * mu( i + 1, j, k ) + la( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx(
            i + 1 );
        cof5 = ( 2 * mu( i + 2, j, k ) + la( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx(
            i + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(1,i-2,j,k) - u( 1, i, j, k ) ) +
                mux2 * ( u(1,i-1,j,k) - u( 1, i, j, k ) ) +
                mux3 * ( u(1,i+1,j,k) - u( 1, i, j, k ) ) +
                mux4 * ( u(1,i+2,j,k) - u( 1, i, j, k ) ) ) * istry;
// qq derivative (u)
        cof1 = ( mu( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry( j - 2 );
        cof2 = ( mu( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry( j - 1 );
        cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
        cof4 = ( mu( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry( j + 1 );
        cof5 = ( mu( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry( j + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(1,i,j-2,k) - u( 1, i, j, k ) ) +
                mux2 * ( u(1,i,j-1,k) - u( 1, i, j, k ) ) +
                mux3 * ( u(1,i,j+1,k) - u( 1, i, j, k ) ) +
                mux4 * ( u(1,i,j+2,k) - u( 1, i, j, k ) ) ) * istrx;
// rr derivative (u)
        cof1 = ( 2 * mu( i, j, k - 2 ) + la( i, j, k - 2 ) ) * SQR( met(2,i,j,k-2)*strx(i) )
            + mu(i,j,k-2) * ( SQR(met(3,i,j,k-2)*stry(j)) + SQR( met(4,i,j,k-2) ) );
        cof2 = ( 2 * mu( i, j, k - 1 ) + la( i, j, k - 1 ) ) * SQR( met(2,i,j,k-1)*strx(i) )
            + mu(i,j,k-1) * ( SQR(met(3,i,j,k-1)*stry(j)) + SQR( met(4,i,j,k-1) ) );
        cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * SQR( met(2,i,j,k)*strx(i) ) +
        mu(i,j,k) * ( SQR(met(3,i,j,k)*stry(j)) + SQR( met(4,i,j,k) ) );
        cof4 = ( 2 * mu( i, j, k + 1 ) + la( i, j, k + 1 ) ) * SQR( met(2,i,j,k+1)*strx(i) )
            + mu(i,j,k+1) * ( SQR(met(3,i,j,k+1)*stry(j)) + SQR( met(4,i,j,k+1) ) );
        cof5 = ( 2 * mu( i, j, k + 2 ) + la( i, j, k + 2 ) ) * SQR( met(2,i,j,k+2)*strx(i) )
            + mu(i,j,k+2) * ( SQR(met(3,i,j,k+2)*stry(j)) + SQR( met(4,i,j,k+2) ) );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(1,i,j,k-2) - u( 1, i, j, k ) ) +
                mux2 * ( u(1,i,j,k-1) - u( 1, i, j, k ) ) +
                mux3 * ( u(1,i,j,k+1) - u( 1, i, j, k ) ) +
                mux4 * ( u(1,i,j,k+2) - u( 1, i, j, k ) ) ) * istrxy;

        // rr derivative (v)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 2, i, j, k - 2 ) * met( 3, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 2, i, j, k - 1 ) * met( 3, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 2, i, j, k ) * met( 3, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 2, i, j, k + 1 ) * met( 3, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 2, i, j, k + 2 ) * met( 3, i, j, k + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(2,i,j,k-2) - u( 2, i, j, k ) ) +
                mux2 * ( u(2,i,j,k-1) - u( 2, i, j, k ) ) +
                mux3 * ( u(2,i,j,k+1) - u( 2, i, j, k ) ) +
                mux4 * ( u(2,i,j,k+2) - u( 2, i, j, k ) ) );

// rr derivative (w)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 2, i, j, k - 2 ) * met( 4, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 2, i, j, k - 1 ) * met( 4, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 2, i, j, k ) * met( 4, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 2, i, j, k + 1 ) * met( 4, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 2, i, j, k + 2 ) * met( 4, i, j, k + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(3,i,j,k-2) - u( 3, i, j, k ) ) +
                mux2 * ( u(3,i,j,k-1) - u( 3, i, j, k ) ) +
                mux3 * ( u(3,i,j,k+1) - u( 3, i, j, k ) ) +
                mux4 * ( u(3,i,j,k+2) - u( 3, i, j, k ) ) ) * istry;

        // pq-derivatives
        r1 = r1 +
            c2 * ( mu(i,j+2,k) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(2,i+2,j+2,k) - u( 2, i - 2, j + 2, k ) ) +
                    c1 * ( u(2,i+1,j+2,k) - u( 2, i - 1, j + 2, k ) ) )
                - mu(i,j-2,k) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(2,i+2,j-2,k) - u( 2, i - 2, j - 2, k ) ) +
                        c1 * ( u(2,i+1,j-2,k) - u( 2, i - 1, j - 2, k ) ) )
                ) +
            c1 * ( mu(i,j+1,k) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(2,i+2,j+1,k) - u( 2, i - 2, j + 1, k ) ) +
                    c1 * ( u(2,i+1,j+1,k) - u( 2, i - 1, j + 1, k ) ) )
                - mu(i,j-1,k) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(2,i+2,j-1,k) - u( 2, i - 2, j - 1, k ) ) +
                        c1 * ( u(2,i+1,j-1,k) - u( 2, i - 1, j - 1, k ) ) ) );

        // qp-derivatives
        r1 = r1 +
            c2 * ( la(i+2,j,k) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                c2 * ( u(2,i+2,j+2,k) - u( 2, i + 2, j - 2, k ) ) +
                    c1 * ( u(2,i+2,j+1,k) - u( 2, i + 2, j - 1, k ) ) )
                - la(i-2,j,k) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                    c2 * ( u(2,i-2,j+2,k) - u( 2, i - 2, j - 2, k ) ) +
                        c1 * ( u(2,i-2,j+1,k) - u( 2, i - 2, j - 1, k ) ) )
                ) +
            c1 * ( la(i+1,j,k) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                c2 * ( u(2,i+1,j+2,k) - u( 2, i + 1, j - 2, k ) ) +
                    c1 * ( u(2,i+1,j+1,k) - u( 2, i + 1, j - 1, k ) ) )
                - la(i-1,j,k) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                    c2 * ( u(2,i-1,j+2,k) - u( 2, i - 1, j - 2, k ) ) +
                        c1 * ( u(2,i-1,j+1,k) - u( 2, i - 1, j - 1, k ) ) ) );

        // pr-derivatives
        r1 = r1 + c2 * (
            ( 2 * mu( i, j, k + 2 ) + la( i, j, k + 2 ) ) * met( 2, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i - 2, j, k + 2 ) ) +
                    c1 * ( u(1,i+1,j,k+2) - u( 1, i - 1, j, k + 2 ) ) ) * strx( i ) * istry
                + mu(i,j,k+2) * met( 3, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                    c2 * ( u(2,i+2,j,k+2) - u( 2, i - 2, j, k + 2 ) ) +
                        c1 * ( u(2,i+1,j,k+2) - u( 2, i - 1, j, k + 2 ) ) )
                + mu(i,j,k+2) * met( 4, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                    c2 * ( u(3,i+2,j,k+2) - u( 3, i - 2, j, k + 2 ) ) +
                        c1 * ( u(3,i+1,j,k+2) - u( 3, i - 1, j, k + 2 ) ) ) * istry
                - ( ( 2 * mu( i, j, k - 2 ) + la( i, j, k - 2 ) ) * met( 2, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(1,i+2,j,k-2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i+1,j,k-2) - u( 1, i - 1, j, k - 2 ) ) ) * strx( i ) * istry
                    + mu(i,j,k-2) * met( 3, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                        c2 * ( u(2,i+2,j,k-2) - u( 2, i - 2, j, k - 2 ) ) +
                            c1 * ( u(2,i+1,j,k-2) - u( 2, i - 1, j, k - 2 ) ) )
                    + mu(i,j,k-2) * met( 4, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                        c2 * ( u(3,i+2,j,k-2) - u( 3, i - 2, j, k - 2 ) ) +
                            c1 * ( u(3,i+1,j,k-2) - u( 3, i - 1, j, k - 2 ) ) ) * istry )
            ) + c1 * (
            ( 2 * mu( i, j, k + 1 ) + la( i, j, k + 1 ) ) * met( 2, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(1,i+2,j,k+1) - u( 1, i - 2, j, k + 1 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i - 1, j, k + 1 ) ) ) * strx( i ) * istry
                + mu(i,j,k+1) * met( 3, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                    c2 * ( u(2,i+2,j,k+1) - u( 2, i - 2, j, k + 1 ) ) +
                        c1 * ( u(2,i+1,j,k+1) - u( 2, i - 1, j, k + 1 ) ) )
                + mu(i,j,k+1) * met( 4, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                    c2 * ( u(3,i+2,j,k+1) - u( 3, i - 2, j, k + 1 ) ) +
                        c1 * ( u(3,i+1,j,k+1) - u( 3, i - 1, j, k + 1 ) ) ) * istry
                - ( ( 2 * mu( i, j, k - 1 ) + la( i, j, k - 1 ) ) * met( 2, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(1,i+2,j,k-1) - u( 1, i - 2, j, k - 1 ) ) +
                        c1 * ( u(1,i+1,j,k-1) - u( 1, i - 1, j, k - 1 ) ) ) * strx( i ) * istry
                    + mu(i,j,k-1) * met( 3, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                        c2 * ( u(2,i+2,j,k-1) - u( 2, i - 2, j, k - 1 ) ) +
                            c1 * ( u(2,i+1,j,k-1) - u( 2, i - 1, j, k - 1 ) ) )
                    + mu(i,j,k-1) * met( 4, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                        c2 * ( u(3,i+2,j,k-1) - u( 3, i - 2, j, k - 1 ) ) +
                            c1 * ( u(3,i+1,j,k-1) - u( 3, i - 1, j, k - 1 ) ) ) * istry ) );

// rp derivatives
        r1 = r1 + ( c2 * (
            ( 2 * mu( i + 2, j, k ) + la( i + 2, j, k ) ) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i + 2, j, k - 2 ) ) +
                    c1 * ( u(1,i+2,j,k+1) - u( 1, i + 2, j, k - 1 ) ) ) * strx( i + 2 )
                + la(i+2,j,k) * met( 3, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                    c2 * ( u(2,i+2,j,k+2) - u( 2, i + 2, j, k - 2 ) ) +
                        c1 * ( u(2,i+2,j,k+1) - u( 2, i + 2, j, k - 1 ) ) ) * stry( j )
                + la(i+2,j,k) * met( 4, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                    c2 * ( u(3,i+2,j,k+2) - u( 3, i + 2, j, k - 2 ) ) +
                        c1 * ( u(3,i+2,j,k+1) - u( 3, i + 2, j, k - 1 ) ) )
                - ( ( 2 * mu( i - 2, j, k ) + la( i - 2, j, k ) ) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                    c2 * ( u(1,i-2,j,k+2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i-2,j,k+1) - u( 1, i - 2, j, k - 1 ) ) ) * strx( i - 2 )
                    + la(i-2,j,k) * met( 3, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                        c2 * ( u(2,i-2,j,k+2) - u( 2, i - 2, j, k - 2 ) ) +
                            c1 * ( u(2,i-2,j,k+1) - u( 2, i - 2, j, k - 1 ) ) ) * stry( j )
                    + la(i-2,j,k) * met( 4, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                        c2 * ( u(3,i-2,j,k+2) - u( 3, i - 2, j, k - 2 ) ) +
                            c1 * ( u(3,i-2,j,k+1) - u( 3, i - 2, j, k - 1 ) ) ) )
            ) + c1 * (
            ( 2 * mu( i + 1, j, k ) + la( i + 1, j, k ) ) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                c2 * ( u(1,i+1,j,k+2) - u( 1, i + 1, j, k - 2 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i + 1, j, k - 1 ) ) ) * strx( i + 1 )
                + la(i+1,j,k) * met( 3, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                    c2 * ( u(2,i+1,j,k+2) - u( 2, i + 1, j, k - 2 ) ) +
                        c1 * ( u(2,i+1,j,k+1) - u( 2, i + 1, j, k - 1 ) ) ) * stry( j )
                + la(i+1,j,k) * met( 4, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                    c2 * ( u(3,i+1,j,k+2) - u( 3, i + 1, j, k - 2 ) ) +
                        c1 * ( u(3,i+1,j,k+1) - u( 3, i + 1, j, k - 1 ) ) )
                - ( ( 2 * mu( i - 1, j, k ) + la( i - 1, j, k ) ) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                    c2 * ( u(1,i-1,j,k+2) - u( 1, i - 1, j, k - 2 ) ) +
                        c1 * ( u(1,i-1,j,k+1) - u( 1, i - 1, j, k - 1 ) ) ) * strx( i - 1 )
                    + la(i-1,j,k) * met( 3, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                        c2 * ( u(2,i-1,j,k+2) - u( 2, i - 1, j, k - 2 ) ) +
                            c1 * ( u(2,i-1,j,k+1) - u( 2, i - 1, j, k - 1 ) ) ) * stry( j )
                    + la(i-1,j,k) * met( 4, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                        c2 * ( u(3,i-1,j,k+2) - u( 3, i - 1, j, k - 2 ) ) +
                            c1 * ( u(3,i-1,j,k+1) - u( 3, i - 1, j, k - 1 ) ) ) ) ) ) * istry;

        // qr derivatives
        r1 = r1 + c2 * (
        mu(i,j,k+2) * met( 3, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
            c2 * ( u(1,i,j+2,k+2) - u( 1, i, j - 2, k + 2 ) ) +
                c1 * ( u(1,i,j+1,k+2) - u( 1, i, j - 1, k + 2 ) ) ) * stry( j ) * istrx
            + la(i,j,k+2) * met( 2, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j - 2, k + 2 ) ) +
                    c1 * ( u(2,i,j+1,k+2) - u( 2, i, j - 1, k + 2 ) ) )
            - ( mu(i,j,k-2) * met( 3, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                c2 * ( u(1,i,j+2,k-2) - u( 1, i, j - 2, k - 2 ) ) +
                    c1 * ( u(1,i,j+1,k-2) - u( 1, i, j - 1, k - 2 ) ) ) * stry( j ) * istrx
                + la(i,j,k-2) * met( 2, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(2,i,j+2,k-2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j+1,k-2) - u( 2, i, j - 1, k - 2 ) ) ) )
            ) + c1 * (
        mu(i,j,k+1) * met( 3, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
            c2 * ( u(1,i,j+2,k+1) - u( 1, i, j - 2, k + 1 ) ) +
                c1 * ( u(1,i,j+1,k+1) - u( 1, i, j - 1, k + 1 ) ) ) * stry( j ) * istrx
            + la(i,j,k+1) * met( 2, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(2,i,j+2,k+1) - u( 2, i, j - 2, k + 1 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j - 1, k + 1 ) ) )
            - ( mu(i,j,k-1) * met( 3, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                c2 * ( u(1,i,j+2,k-1) - u( 1, i, j - 2, k - 1 ) ) +
                    c1 * ( u(1,i,j+1,k-1) - u( 1, i, j - 1, k - 1 ) ) ) * stry( j ) * istrx
                + la(i,j,k-1) * met( 2, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(2,i,j+2,k-1) - u( 2, i, j - 2, k - 1 ) ) +
                        c1 * ( u(2,i,j+1,k-1) - u( 2, i, j - 1, k - 1 ) ) ) ) );

// rq derivatives
        r1 = r1 + c2 * (
        mu(i,j+2,k) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
            c2 * ( u(1,i,j+2,k+2) - u( 1, i, j + 2, k - 2 ) ) +
                c1 * ( u(1,i,j+2,k+1) - u( 1, i, j + 2, k - 1 ) ) ) * stry( j + 2 ) * istrx
            + mu(i,j+2,k) * met( 2, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j + 2, k - 2 ) ) +
                    c1 * ( u(2,i,j+2,k+1) - u( 2, i, j + 2, k - 1 ) ) )
            - ( mu(i,j-2,k) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                c2 * ( u(1,i,j-2,k+2) - u( 1, i, j - 2, k - 2 ) ) +
                    c1 * ( u(1,i,j-2,k+1) - u( 1, i, j - 2, k - 1 ) ) ) * stry( j - 2 ) * istrx
                + mu(i,j-2,k) * met( 2, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(2,i,j-2,k+2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j-2,k+1) - u( 2, i, j - 2, k - 1 ) ) ) )
            ) + c1 * (
        mu(i,j+1,k) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
            c2 * ( u(1,i,j+1,k+2) - u( 1, i, j + 1, k - 2 ) ) +
                c1 * ( u(1,i,j+1,k+1) - u( 1, i, j + 1, k - 1 ) ) ) * stry( j + 1 ) * istrx
            + mu(i,j+1,k) * met( 2, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(2,i,j+1,k+2) - u( 2, i, j + 1, k - 2 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j + 1, k - 1 ) ) )
            - ( mu(i,j-1,k) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                c2 * ( u(1,i,j-1,k+2) - u( 1, i, j - 1, k - 2 ) ) +
                    c1 * ( u(1,i,j-1,k+1) - u( 1, i, j - 1, k - 1 ) ) ) * stry( j - 1 ) * istrx
                + mu(i,j-1,k) * met( 2, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(2,i,j-1,k+2) - u( 2, i, j - 1, k - 2 ) ) +
                        c1 * ( u(2,i,j-1,k+1) - u( 2, i, j - 1, k - 1 ) ) ) ) );

        //          lu(1,i,j,k) = r1/jac(i,j,k)
        //          lu(1,i,j,k) = r1*ijac
        lu(1,i,j,k) = a1 * lu( 1, i, j, k ) + sgn * r1 * ijac;
// v-equation

        r1 = 0;
        // pp derivative (v)
        cof1 = ( mu( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx( i - 2 );
        cof2 = ( mu( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx( i - 1 );
        cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
        cof4 = ( mu( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx( i + 1 );
        cof5 = ( mu( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx( i + 2 );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(2,i-2,j,k) - u( 2, i, j, k ) ) +
                mux2 * ( u(2,i-1,j,k) - u( 2, i, j, k ) ) +
                mux3 * ( u(2,i+1,j,k) - u( 2, i, j, k ) ) +
                mux4 * ( u(2,i+2,j,k) - u( 2, i, j, k ) ) ) * istry;
        // qq derivative (v)
        cof1 = ( 2 * mu( i, j - 2, k ) + la( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry(
            j - 2 );
        cof2 = ( 2 * mu( i, j - 1, k ) + la( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry(
            j - 1 );
        cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
        cof4 = ( 2 * mu( i, j + 1, k ) + la( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry(
            j + 1 );
        cof5 = ( 2 * mu( i, j + 2, k ) + la( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry(
            j + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(2,i,j-2,k) - u( 2, i, j, k ) ) +
                mux2 * ( u(2,i,j-1,k) - u( 2, i, j, k ) ) +
                mux3 * ( u(2,i,j+1,k) - u( 2, i, j, k ) ) +
                mux4 * ( u(2,i,j+2,k) - u( 2, i, j, k ) ) ) * istrx;

        // rr derivative (u)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 2, i, j, k - 2 ) * met( 3, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 2, i, j, k - 1 ) * met( 3, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 2, i, j, k ) * met( 3, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 2, i, j, k + 1 ) * met( 3, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 2, i, j, k + 2 ) * met( 3, i, j, k + 2 );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(1,i,j,k-2) - u( 1, i, j, k ) ) +
                mux2 * ( u(1,i,j,k-1) - u( 1, i, j, k ) ) +
                mux3 * ( u(1,i,j,k+1) - u( 1, i, j, k ) ) +
                mux4 * ( u(1,i,j,k+2) - u( 1, i, j, k ) ) );

        // rr derivative (v)
        cof1 = ( 2 * mu( i, j, k - 2 ) + la( i, j, k - 2 ) ) * SQR( met(3,i,j,k-2)*stry(j) )
            + mu(i,j,k-2) * ( SQR(met(2,i,j,k-2)*strx(i)) + SQR( met(4,i,j,k-2) ) );
        cof2 = ( 2 * mu( i, j, k - 1 ) + la( i, j, k - 1 ) ) * SQR( met(3,i,j,k-1)*stry(j) )
            + mu(i,j,k-1) * ( SQR(met(2,i,j,k-1)*strx(i)) + SQR( met(4,i,j,k-1) ) );
        cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * SQR( met(3,i,j,k)*stry(j) ) +
        mu(i,j,k) * ( SQR(met(2,i,j,k)*strx(i)) + SQR( met(4,i,j,k) ) );
        cof4 = ( 2 * mu( i, j, k + 1 ) + la( i, j, k + 1 ) ) * SQR( met(3,i,j,k+1)*stry(j) )
            + mu(i,j,k+1) * ( SQR(met(2,i,j,k+1)*strx(i)) + SQR( met(4,i,j,k+1) ) );
        cof5 = ( 2 * mu( i, j, k + 2 ) + la( i, j, k + 2 ) ) * SQR( met(3,i,j,k+2)*stry(j) )
            + mu(i,j,k+2) * ( SQR(met(2,i,j,k+2)*strx(i)) + SQR( met(4,i,j,k+2) ) );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(2,i,j,k-2) - u( 2, i, j, k ) ) +
                mux2 * ( u(2,i,j,k-1) - u( 2, i, j, k ) ) +
                mux3 * ( u(2,i,j,k+1) - u( 2, i, j, k ) ) +
                mux4 * ( u(2,i,j,k+2) - u( 2, i, j, k ) ) ) * istrxy;

// rr derivative (w)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 3, i, j, k - 2 ) * met( 4, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 3, i, j, k - 1 ) * met( 4, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 3, i, j, k ) * met( 4, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 3, i, j, k + 1 ) * met( 4, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 3, i, j, k + 2 ) * met( 4, i, j, k + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(3,i,j,k-2) - u( 3, i, j, k ) ) +
                mux2 * ( u(3,i,j,k-1) - u( 3, i, j, k ) ) +
                mux3 * ( u(3,i,j,k+1) - u( 3, i, j, k ) ) +
                mux4 * ( u(3,i,j,k+2) - u( 3, i, j, k ) ) ) * istrx;

// pq-derivatives
        r1 = r1 +
            c2 * ( la(i,j+2,k) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(1,i+2,j+2,k) - u( 1, i - 2, j + 2, k ) ) +
                    c1 * ( u(1,i+1,j+2,k) - u( 1, i - 1, j + 2, k ) ) )
                - la(i,j-2,k) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(1,i+2,j-2,k) - u( 1, i - 2, j - 2, k ) ) +
                        c1 * ( u(1,i+1,j-2,k) - u( 1, i - 1, j - 2, k ) ) )
                ) +
            c1 * ( la(i,j+1,k) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(1,i+2,j+1,k) - u( 1, i - 2, j + 1, k ) ) +
                    c1 * ( u(1,i+1,j+1,k) - u( 1, i - 1, j + 1, k ) ) )
                - la(i,j-1,k) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(1,i+2,j-1,k) - u( 1, i - 2, j - 1, k ) ) +
                        c1 * ( u(1,i+1,j-1,k) - u( 1, i - 1, j - 1, k ) ) ) );

        // qp-derivatives
        r1 = r1 +
            c2 * ( mu(i+2,j,k) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                c2 * ( u(1,i+2,j+2,k) - u( 1, i + 2, j - 2, k ) ) +
                    c1 * ( u(1,i+2,j+1,k) - u( 1, i + 2, j - 1, k ) ) )
                - mu(i-2,j,k) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                    c2 * ( u(1,i-2,j+2,k) - u( 1, i - 2, j - 2, k ) ) +
                        c1 * ( u(1,i-2,j+1,k) - u( 1, i - 2, j - 1, k ) ) )
                ) +
            c1 * ( mu(i+1,j,k) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                c2 * ( u(1,i+1,j+2,k) - u( 1, i + 1, j - 2, k ) ) +
                    c1 * ( u(1,i+1,j+1,k) - u( 1, i + 1, j - 1, k ) ) )
                - mu(i-1,j,k) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                    c2 * ( u(1,i-1,j+2,k) - u( 1, i - 1, j - 2, k ) ) +
                        c1 * ( u(1,i-1,j+1,k) - u( 1, i - 1, j - 1, k ) ) ) );

        // pr-derivatives
        r1 = r1 + c2 * (
            ( la( i, j, k + 2 ) ) * met( 3, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i - 2, j, k + 2 ) ) +
                    c1 * ( u(1,i+1,j,k+2) - u( 1, i - 1, j, k + 2 ) ) )
                + mu(i,j,k+2) * met( 2, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                    c2 * ( u(2,i+2,j,k+2) - u( 2, i - 2, j, k + 2 ) ) +
                        c1 * ( u(2,i+1,j,k+2) - u( 2, i - 1, j, k + 2 ) ) ) * strx( i ) * istry
                - ( ( la( i, j, k - 2 ) ) * met( 3, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(1,i+2,j,k-2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i+1,j,k-2) - u( 1, i - 1, j, k - 2 ) ) )
                    + mu(i,j,k-2) * met( 2, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                        c2 * ( u(2,i+2,j,k-2) - u( 2, i - 2, j, k - 2 ) ) +
                            c1 * ( u(2,i+1,j,k-2) - u( 2, i - 1, j, k - 2 ) ) ) * strx( i ) * istry )
            ) + c1 * (
            ( la( i, j, k + 1 ) ) * met( 3, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(1,i+2,j,k+1) - u( 1, i - 2, j, k + 1 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i - 1, j, k + 1 ) ) )
                + mu(i,j,k+1) * met( 2, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                    c2 * ( u(2,i+2,j,k+1) - u( 2, i - 2, j, k + 1 ) ) +
                        c1 * ( u(2,i+1,j,k+1) - u( 2, i - 1, j, k + 1 ) ) ) * strx( i ) * istry
                - ( la(i,j,k-1) * met( 3, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(1,i+2,j,k-1) - u( 1, i - 2, j, k - 1 ) ) +
                        c1 * ( u(1,i+1,j,k-1) - u( 1, i - 1, j, k - 1 ) ) )
                    + mu(i,j,k-1) * met( 2, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                        c2 * ( u(2,i+2,j,k-1) - u( 2, i - 2, j, k - 1 ) ) +
                            c1 * ( u(2,i+1,j,k-1) - u( 2, i - 1, j, k - 1 ) ) ) * strx( i ) * istry ) );

        // rp derivatives
        r1 = r1 + c2 * (
            ( mu( i + 2, j, k ) ) * met( 3, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i + 2, j, k - 2 ) ) +
                    c1 * ( u(1,i+2,j,k+1) - u( 1, i + 2, j, k - 1 ) ) )
                + mu(i+2,j,k) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                    c2 * ( u(2,i+2,j,k+2) - u( 2, i + 2, j, k - 2 ) ) +
                        c1 * ( u(2,i+2,j,k+1) - u( 2, i + 2, j, k - 1 ) ) ) * strx( i + 2 ) * istry
                - ( mu(i-2,j,k) * met( 3, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                    c2 * ( u(1,i-2,j,k+2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i-2,j,k+1) - u( 1, i - 2, j, k - 1 ) ) )
                    + mu(i-2,j,k) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                        c2 * ( u(2,i-2,j,k+2) - u( 2, i - 2, j, k - 2 ) ) +
                            c1 * ( u(2,i-2,j,k+1) - u( 2, i - 2, j, k - 1 ) ) ) * strx( i - 2 ) * istry )
            ) + c1 * (
            ( mu( i + 1, j, k ) ) * met( 3, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                c2 * ( u(1,i+1,j,k+2) - u( 1, i + 1, j, k - 2 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i + 1, j, k - 1 ) ) )
                + mu(i+1,j,k) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                    c2 * ( u(2,i+1,j,k+2) - u( 2, i + 1, j, k - 2 ) ) +
                        c1 * ( u(2,i+1,j,k+1) - u( 2, i + 1, j, k - 1 ) ) ) * strx( i + 1 ) * istry
                - ( mu(i-1,j,k) * met( 3, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                    c2 * ( u(1,i-1,j,k+2) - u( 1, i - 1, j, k - 2 ) ) +
                        c1 * ( u(1,i-1,j,k+1) - u( 1, i - 1, j, k - 1 ) ) )
                    + mu(i-1,j,k) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                        c2 * ( u(2,i-1,j,k+2) - u( 2, i - 1, j, k - 2 ) ) +
                            c1 * ( u(2,i-1,j,k+1) - u( 2, i - 1, j, k - 1 ) ) ) * strx( i - 1 ) * istry ) );

        // qr derivatives
        r1 = r1 + c2 * (
        mu(i,j,k+2) * met( 2, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
            c2 * ( u(1,i,j+2,k+2) - u( 1, i, j - 2, k + 2 ) ) +
                c1 * ( u(1,i,j+1,k+2) - u( 1, i, j - 1, k + 2 ) ) )
            + ( 2 * mu( i, j, k + 2 ) + la( i, j, k + 2 ) ) * met( 3, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j - 2, k + 2 ) ) +
                    c1 * ( u(2,i,j+1,k+2) - u( 2, i, j - 1, k + 2 ) ) ) * stry( j ) * istrx
            + mu(i,j,k+2) * met( 4, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(3,i,j+2,k+2) - u( 3, i, j - 2, k + 2 ) ) +
                    c1 * ( u(3,i,j+1,k+2) - u( 3, i, j - 1, k + 2 ) ) ) * istrx
            - ( mu(i,j,k-2) * met( 2, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                c2 * ( u(1,i,j+2,k-2) - u( 1, i, j - 2, k - 2 ) ) +
                    c1 * ( u(1,i,j+1,k-2) - u( 1, i, j - 1, k - 2 ) ) )
                + ( 2 * mu( i, j, k - 2 ) + la( i, j, k - 2 ) ) * met( 3, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(2,i,j+2,k-2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j+1,k-2) - u( 2, i, j - 1, k - 2 ) ) ) * stry( j ) * istrx +
                mu(i,j,k-2) * met( 4, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(3,i,j+2,k-2) - u( 3, i, j - 2, k - 2 ) ) +
                        c1 * ( u(3,i,j+1,k-2) - u( 3, i, j - 1, k - 2 ) ) ) * istrx )
            ) + c1 * (
        mu(i,j,k+1) * met( 2, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
            c2 * ( u(1,i,j+2,k+1) - u( 1, i, j - 2, k + 1 ) ) +
                c1 * ( u(1,i,j+1,k+1) - u( 1, i, j - 1, k + 1 ) ) )
            + ( 2 * mu( i, j, k + 1 ) + la( i, j, k + 1 ) ) * met( 3, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(2,i,j+2,k+1) - u( 2, i, j - 2, k + 1 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j - 1, k + 1 ) ) ) * stry( j ) * istrx
            + mu(i,j,k+1) * met( 4, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(3,i,j+2,k+1) - u( 3, i, j - 2, k + 1 ) ) +
                    c1 * ( u(3,i,j+1,k+1) - u( 3, i, j - 1, k + 1 ) ) ) * istrx
            - ( mu(i,j,k-1) * met( 2, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                c2 * ( u(1,i,j+2,k-1) - u( 1, i, j - 2, k - 1 ) ) +
                    c1 * ( u(1,i,j+1,k-1) - u( 1, i, j - 1, k - 1 ) ) )
                + ( 2 * mu( i, j, k - 1 ) + la( i, j, k - 1 ) ) * met( 3, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(2,i,j+2,k-1) - u( 2, i, j - 2, k - 1 ) ) +
                        c1 * ( u(2,i,j+1,k-1) - u( 2, i, j - 1, k - 1 ) ) ) * stry( j ) * istrx
                + mu(i,j,k-1) * met( 4, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(3,i,j+2,k-1) - u( 3, i, j - 2, k - 1 ) ) +
                        c1 * ( u(3,i,j+1,k-1) - u( 3, i, j - 1, k - 1 ) ) ) * istrx ) );

        // rq derivatives
        r1 = r1 + c2 * (
        la(i,j+2,k) * met( 2, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
            c2 * ( u(1,i,j+2,k+2) - u( 1, i, j + 2, k - 2 ) ) +
                c1 * ( u(1,i,j+2,k+1) - u( 1, i, j + 2, k - 1 ) ) )
            + ( 2 * mu( i, j + 2, k ) + la( i, j + 2, k ) ) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j + 2, k - 2 ) ) +
                    c1 * ( u(2,i,j+2,k+1) - u( 2, i, j + 2, k - 1 ) ) ) * stry( j + 2 ) * istrx
            + la(i,j+2,k) * met( 4, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(3,i,j+2,k+2) - u( 3, i, j + 2, k - 2 ) ) +
                    c1 * ( u(3,i,j+2,k+1) - u( 3, i, j + 2, k - 1 ) ) ) * istrx
            - ( la(i,j-2,k) * met( 2, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                c2 * ( u(1,i,j-2,k+2) - u( 1, i, j - 2, k - 2 ) ) +
                    c1 * ( u(1,i,j-2,k+1) - u( 1, i, j - 2, k - 1 ) ) )
                + ( 2 * mu( i, j - 2, k ) + la( i, j - 2, k ) ) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(2,i,j-2,k+2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j-2,k+1) - u( 2, i, j - 2, k - 1 ) ) ) * stry( j - 2 ) * istrx
                + la(i,j-2,k) * met( 4, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(3,i,j-2,k+2) - u( 3, i, j - 2, k - 2 ) ) +
                        c1 * ( u(3,i,j-2,k+1) - u( 3, i, j - 2, k - 1 ) ) ) * istrx )
            ) + c1 * (
        la(i,j+1,k) * met( 2, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
            c2 * ( u(1,i,j+1,k+2) - u( 1, i, j + 1, k - 2 ) ) +
                c1 * ( u(1,i,j+1,k+1) - u( 1, i, j + 1, k - 1 ) ) )
            + ( 2 * mu( i, j + 1, k ) + la( i, j + 1, k ) ) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(2,i,j+1,k+2) - u( 2, i, j + 1, k - 2 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j + 1, k - 1 ) ) ) * stry( j + 1 ) * istrx
            + la(i,j+1,k) * met( 4, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(3,i,j+1,k+2) - u( 3, i, j + 1, k - 2 ) ) +
                    c1 * ( u(3,i,j+1,k+1) - u( 3, i, j + 1, k - 1 ) ) ) * istrx
            - ( la(i,j-1,k) * met( 2, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                c2 * ( u(1,i,j-1,k+2) - u( 1, i, j - 1, k - 2 ) ) +
                    c1 * ( u(1,i,j-1,k+1) - u( 1, i, j - 1, k - 1 ) ) )
                + ( 2 * mu( i, j - 1, k ) + la( i, j - 1, k ) ) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(2,i,j-1,k+2) - u( 2, i, j - 1, k - 2 ) ) +
                        c1 * ( u(2,i,j-1,k+1) - u( 2, i, j - 1, k - 1 ) ) ) * stry( j - 1 ) * istrx
                + la(i,j-1,k) * met( 4, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(3,i,j-1,k+2) - u( 3, i, j - 1, k - 2 ) ) +
                        c1 * ( u(3,i,j-1,k+1) - u( 3, i, j - 1, k - 1 ) ) ) * istrx ) );

        //         lu(2,i,j,k) = r1/jac(i,j,k)
        //         lu(2,i,j,k) = r1*ijac
        lu(2,i,j,k) = a1 * lu( 2, i, j, k ) + sgn * r1 * ijac;
        // w-equation

        r1 = 0;
        // pp derivative (w)
        cof1 = ( mu( i - 2, j, k ) ) * met( 1, i - 2, j, k ) * met( 1, i - 2, j, k ) * strx( i - 2 );
        cof2 = ( mu( i - 1, j, k ) ) * met( 1, i - 1, j, k ) * met( 1, i - 1, j, k ) * strx( i - 1 );
        cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * strx( i );
        cof4 = ( mu( i + 1, j, k ) ) * met( 1, i + 1, j, k ) * met( 1, i + 1, j, k ) * strx( i + 1 );
        cof5 = ( mu( i + 2, j, k ) ) * met( 1, i + 2, j, k ) * met( 1, i + 2, j, k ) * strx( i + 2 );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(3,i-2,j,k) - u( 3, i, j, k ) ) +
                mux2 * ( u(3,i-1,j,k) - u( 3, i, j, k ) ) +
                mux3 * ( u(3,i+1,j,k) - u( 3, i, j, k ) ) +
                mux4 * ( u(3,i+2,j,k) - u( 3, i, j, k ) ) ) * istry;

// qq derivative (w)
        cof1 = ( mu( i, j - 2, k ) ) * met( 1, i, j - 2, k ) * met( 1, i, j - 2, k ) * stry( j - 2 );
        cof2 = ( mu( i, j - 1, k ) ) * met( 1, i, j - 1, k ) * met( 1, i, j - 1, k ) * stry( j - 1 );
        cof3 = ( mu( i, j, k ) ) * met( 1, i, j, k ) * met( 1, i, j, k ) * stry( j );
        cof4 = ( mu( i, j + 1, k ) ) * met( 1, i, j + 1, k ) * met( 1, i, j + 1, k ) * stry( j + 1 );
        cof5 = ( mu( i, j + 2, k ) ) * met( 1, i, j + 2, k ) * met( 1, i, j + 2, k ) * stry( j + 2 );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(3,i,j-2,k) - u( 3, i, j, k ) ) +
                mux2 * ( u(3,i,j-1,k) - u( 3, i, j, k ) ) +
                mux3 * ( u(3,i,j+1,k) - u( 3, i, j, k ) ) +
                mux4 * ( u(3,i,j+2,k) - u( 3, i, j, k ) ) ) * istrx;
        // rr derivative (u)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 2, i, j, k - 2 ) * met( 4, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 2, i, j, k - 1 ) * met( 4, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 2, i, j, k ) * met( 4, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 2, i, j, k + 1 ) * met( 4, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 2, i, j, k + 2 ) * met( 4, i, j, k + 2 );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(1,i,j,k-2) - u( 1, i, j, k ) ) +
                mux2 * ( u(1,i,j,k-1) - u( 1, i, j, k ) ) +
                mux3 * ( u(1,i,j,k+1) - u( 1, i, j, k ) ) +
                mux4 * ( u(1,i,j,k+2) - u( 1, i, j, k ) ) ) * istry;

        // rr derivative (v)
        cof1 = ( mu(i,j,k-2) + la( i, j, k - 2 ) ) * met( 3, i, j, k - 2 ) * met( 4, i, j, k - 2 );
        cof2 = ( mu(i,j,k-1) + la( i, j, k - 1 ) ) * met( 3, i, j, k - 1 ) * met( 4, i, j, k - 1 );
        cof3 = ( mu(i,j,k) + la( i, j, k ) ) * met( 3, i, j, k ) * met( 4, i, j, k );
        cof4 = ( mu(i,j,k+1) + la( i, j, k + 1 ) ) * met( 3, i, j, k + 1 ) * met( 4, i, j, k + 1 );
        cof5 = ( mu(i,j,k+2) + la( i, j, k + 2 ) ) * met( 3, i, j, k + 2 ) * met( 4, i, j, k + 2 );

        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(2,i,j,k-2) - u( 2, i, j, k ) ) +
                mux2 * ( u(2,i,j,k-1) - u( 2, i, j, k ) ) +
                mux3 * ( u(2,i,j,k+1) - u( 2, i, j, k ) ) +
                mux4 * ( u(2,i,j,k+2) - u( 2, i, j, k ) ) ) * istrx;

// rr derivative (w)
        cof1 = ( 2 * mu( i, j, k - 2 ) + la( i, j, k - 2 ) ) * SQR( met(4,i,j,k-2) ) +
        mu(i,j,k-2) * ( SQR(met(2,i,j,k-2)*strx(i)) +
            SQR( met(3,i,j,k-2)*stry(j) ) );
        cof2 = ( 2 * mu( i, j, k - 1 ) + la( i, j, k - 1 ) ) * SQR( met(4,i,j,k-1) ) +
        mu(i,j,k-1) * ( SQR(met(2,i,j,k-1)*strx(i)) +
            SQR( met(3,i,j,k-1)*stry(j) ) );
        cof3 = ( 2 * mu( i, j, k ) + la( i, j, k ) ) * SQR( met(4,i,j,k) ) +
        mu(i,j,k) * ( SQR(met(2,i,j,k)*strx(i)) +
            SQR( met(3,i,j,k)*stry(j) ) );
        cof4 = ( 2 * mu( i, j, k + 1 ) + la( i, j, k + 1 ) ) * SQR( met(4,i,j,k+1) ) +
        mu(i,j,k+1) * ( SQR(met(2,i,j,k+1)*strx(i)) +
            SQR( met(3,i,j,k+1)*stry(j) ) );
        cof5 = ( 2 * mu( i, j, k + 2 ) + la( i, j, k + 2 ) ) * SQR( met(4,i,j,k+2) ) +
        mu(i,j,k+2) * ( SQR(met(2,i,j,k+2)*strx(i)) +
            SQR( met(3,i,j,k+2)*stry(j) ) );
        mux1 = cof2 - tf * ( cof3 + cof1 );
        mux2 = cof1 + cof4 + 3 * ( cof3 + cof2 );
        mux3 = cof2 + cof5 + 3 * ( cof4 + cof3 );
        mux4 = cof4 - tf * ( cof3 + cof5 );

        r1 = r1 + i6 * (
            mux1 * ( u(3,i,j,k-2) - u( 3, i, j, k ) ) +
                mux2 * ( u(3,i,j,k-1) - u( 3, i, j, k ) ) +
                mux3 * ( u(3,i,j,k+1) - u( 3, i, j, k ) ) +
                mux4 * ( u(3,i,j,k+2) - u( 3, i, j, k ) ) ) * istrxy

        // pr-derivatives
        //	    r1 = r1

        + c2 * (
            ( la( i, j, k + 2 ) ) * met( 4, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i - 2, j, k + 2 ) ) +
                    c1 * ( u(1,i+1,j,k+2) - u( 1, i - 1, j, k + 2 ) ) ) * istry
                + mu(i,j,k+2) * met( 2, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                    c2 * ( u(3,i+2,j,k+2) - u( 3, i - 2, j, k + 2 ) ) +
                        c1 * ( u(3,i+1,j,k+2) - u( 3, i - 1, j, k + 2 ) ) ) * strx( i ) * istry
                - ( ( la( i, j, k - 2 ) ) * met( 4, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(1,i+2,j,k-2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i+1,j,k-2) - u( 1, i - 1, j, k - 2 ) ) ) * istry
                    + mu(i,j,k-2) * met( 2, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                        c2 * ( u(3,i+2,j,k-2) - u( 3, i - 2, j, k - 2 ) ) +
                            c1 * ( u(3,i+1,j,k-2) - u( 3, i - 1, j, k - 2 ) ) ) * strx( i ) * istry )
            ) + c1 * (
            ( la( i, j, k + 1 ) ) * met( 4, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(1,i+2,j,k+1) - u( 1, i - 2, j, k + 1 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i - 1, j, k + 1 ) ) ) * istry
                + mu(i,j,k+1) * met( 2, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                    c2 * ( u(3,i+2,j,k+1) - u( 3, i - 2, j, k + 1 ) ) +
                        c1 * ( u(3,i+1,j,k+1) - u( 3, i - 1, j, k + 1 ) ) ) * strx( i ) * istry
                - ( la(i,j,k-1) * met( 4, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(1,i+2,j,k-1) - u( 1, i - 2, j, k - 1 ) ) +
                        c1 * ( u(1,i+1,j,k-1) - u( 1, i - 1, j, k - 1 ) ) ) * istry
                    + mu(i,j,k-1) * met( 2, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                        c2 * ( u(3,i+2,j,k-1) - u( 3, i - 2, j, k - 1 ) ) +
                            c1 * ( u(3,i+1,j,k-1) - u( 3, i - 1, j, k - 1 ) ) ) * strx( i ) * istry ) )

        // rp derivatives
        //      r1 = r1

        + istry * ( c2 * (
            ( mu( i + 2, j, k ) ) * met( 4, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                c2 * ( u(1,i+2,j,k+2) - u( 1, i + 2, j, k - 2 ) ) +
                    c1 * ( u(1,i+2,j,k+1) - u( 1, i + 2, j, k - 1 ) ) )
                + mu(i+2,j,k) * met( 2, i + 2, j, k ) * met( 1, i + 2, j, k ) * (
                    c2 * ( u(3,i+2,j,k+2) - u( 3, i + 2, j, k - 2 ) ) +
                        c1 * ( u(3,i+2,j,k+1) - u( 3, i + 2, j, k - 1 ) ) ) * strx( i + 2 )
                - ( mu(i-2,j,k) * met( 4, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                    c2 * ( u(1,i-2,j,k+2) - u( 1, i - 2, j, k - 2 ) ) +
                        c1 * ( u(1,i-2,j,k+1) - u( 1, i - 2, j, k - 1 ) ) )
                    + mu(i-2,j,k) * met( 2, i - 2, j, k ) * met( 1, i - 2, j, k ) * (
                        c2 * ( u(3,i-2,j,k+2) - u( 3, i - 2, j, k - 2 ) ) +
                            c1 * ( u(3,i-2,j,k+1) - u( 3, i - 2, j, k - 1 ) ) ) * strx( i - 2 ) )
            ) + c1 * (
            ( mu( i + 1, j, k ) ) * met( 4, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                c2 * ( u(1,i+1,j,k+2) - u( 1, i + 1, j, k - 2 ) ) +
                    c1 * ( u(1,i+1,j,k+1) - u( 1, i + 1, j, k - 1 ) ) )
                + mu(i+1,j,k) * met( 2, i + 1, j, k ) * met( 1, i + 1, j, k ) * (
                    c2 * ( u(3,i+1,j,k+2) - u( 3, i + 1, j, k - 2 ) ) +
                        c1 * ( u(3,i+1,j,k+1) - u( 3, i + 1, j, k - 1 ) ) ) * strx( i + 1 )
                - ( mu(i-1,j,k) * met( 4, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                    c2 * ( u(1,i-1,j,k+2) - u( 1, i - 1, j, k - 2 ) ) +
                        c1 * ( u(1,i-1,j,k+1) - u( 1, i - 1, j, k - 1 ) ) )
                    + mu(i-1,j,k) * met( 2, i - 1, j, k ) * met( 1, i - 1, j, k ) * (
                        c2 * ( u(3,i-1,j,k+2) - u( 3, i - 1, j, k - 2 ) ) +
                            c1 * ( u(3,i-1,j,k+1) - u( 3, i - 1, j, k - 1 ) ) ) * strx( i - 1 ) ) ) )

        // qr derivatives

        //      r1 = r1
        + c2 * (
        mu(i,j,k+2) * met( 3, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
            c2 * ( u(3,i,j+2,k+2) - u( 3, i, j - 2, k + 2 ) ) +
                c1 * ( u(3,i,j+1,k+2) - u( 3, i, j - 1, k + 2 ) ) ) * stry( j ) * istrx
            + la(i,j,k+2) * met( 4, i, j, k + 2 ) * met( 1, i, j, k + 2 ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j - 2, k + 2 ) ) +
                    c1 * ( u(2,i,j+1,k+2) - u( 2, i, j - 1, k + 2 ) ) ) * istrx
            - ( mu(i,j,k-2) * met( 3, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                c2 * ( u(3,i,j+2,k-2) - u( 3, i, j - 2, k - 2 ) ) +
                    c1 * ( u(3,i,j+1,k-2) - u( 3, i, j - 1, k - 2 ) ) ) * stry( j ) * istrx
                + la(i,j,k-2) * met( 4, i, j, k - 2 ) * met( 1, i, j, k - 2 ) * (
                    c2 * ( u(2,i,j+2,k-2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j+1,k-2) - u( 2, i, j - 1, k - 2 ) ) ) * istrx )
            ) + c1 * (
        mu(i,j,k+1) * met( 3, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
            c2 * ( u(3,i,j+2,k+1) - u( 3, i, j - 2, k + 1 ) ) +
                c1 * ( u(3,i,j+1,k+1) - u( 3, i, j - 1, k + 1 ) ) ) * stry( j ) * istrx
            + la(i,j,k+1) * met( 4, i, j, k + 1 ) * met( 1, i, j, k + 1 ) * (
                c2 * ( u(2,i,j+2,k+1) - u( 2, i, j - 2, k + 1 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j - 1, k + 1 ) ) ) * istrx
            - ( mu(i,j,k-1) * met( 3, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                c2 * ( u(3,i,j+2,k-1) - u( 3, i, j - 2, k - 1 ) ) +
                    c1 * ( u(3,i,j+1,k-1) - u( 3, i, j - 1, k - 1 ) ) ) * stry( j ) * istrx
                + la(i,j,k-1) * met( 4, i, j, k - 1 ) * met( 1, i, j, k - 1 ) * (
                    c2 * ( u(2,i,j+2,k-1) - u( 2, i, j - 2, k - 1 ) ) +
                        c1 * ( u(2,i,j+1,k-1) - u( 2, i, j - 1, k - 1 ) ) ) * istrx ) )

        // rq derivatives

        //      r1 = r1
        + istrx * ( c2 * (
        mu(i,j+2,k) * met( 3, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
            c2 * ( u(3,i,j+2,k+2) - u( 3, i, j + 2, k - 2 ) ) +
                c1 * ( u(3,i,j+2,k+1) - u( 3, i, j + 2, k - 1 ) ) ) * stry( j + 2 )
            + mu(i,j+2,k) * met( 4, i, j + 2, k ) * met( 1, i, j + 2, k ) * (
                c2 * ( u(2,i,j+2,k+2) - u( 2, i, j + 2, k - 2 ) ) +
                    c1 * ( u(2,i,j+2,k+1) - u( 2, i, j + 2, k - 1 ) ) )
            - ( mu(i,j-2,k) * met( 3, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                c2 * ( u(3,i,j-2,k+2) - u( 3, i, j - 2, k - 2 ) ) +
                    c1 * ( u(3,i,j-2,k+1) - u( 3, i, j - 2, k - 1 ) ) ) * stry( j - 2 )
                + mu(i,j-2,k) * met( 4, i, j - 2, k ) * met( 1, i, j - 2, k ) * (
                    c2 * ( u(2,i,j-2,k+2) - u( 2, i, j - 2, k - 2 ) ) +
                        c1 * ( u(2,i,j-2,k+1) - u( 2, i, j - 2, k - 1 ) ) ) )
            ) + c1 * (
        mu(i,j+1,k) * met( 3, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
            c2 * ( u(3,i,j+1,k+2) - u( 3, i, j + 1, k - 2 ) ) +
                c1 * ( u(3,i,j+1,k+1) - u( 3, i, j + 1, k - 1 ) ) ) * stry( j + 1 )
            + mu(i,j+1,k) * met( 4, i, j + 1, k ) * met( 1, i, j + 1, k ) * (
                c2 * ( u(2,i,j+1,k+2) - u( 2, i, j + 1, k - 2 ) ) +
                    c1 * ( u(2,i,j+1,k+1) - u( 2, i, j + 1, k - 1 ) ) )
            - ( mu(i,j-1,k) * met( 3, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                c2 * ( u(3,i,j-1,k+2) - u( 3, i, j - 1, k - 2 ) ) +
                    c1 * ( u(3,i,j-1,k+1) - u( 3, i, j - 1, k - 1 ) ) ) * stry( j - 1 )
                + mu(i,j-1,k) * met( 4, i, j - 1, k ) * met( 1, i, j - 1, k ) * (
                    c2 * ( u(2,i,j-1,k+2) - u( 2, i, j - 1, k - 2 ) ) +
                        c1 * ( u(2,i,j-1,k+1) - u( 2, i, j - 1, k - 1 ) ) ) ) ) );
        //         lu(3,i,j,k) = r1/jac(i,j,k)
        //         lu(3,i,j,k) = r1*ijac
        lu(3,i,j,k) = a1 * lu( 3, i, j, k ) + sgn * r1 * ijac;
      }

// Cartesian grid discretization, use for k >= kcart + 2
  double h;
  double cof = 0.0;
  if( kcarteff + 2 <= m_klast - 2 )
  {
    h = met(1,m_ifirst+2,m_jfirst+2,kcarteff+2) * met( 1, m_ifirst + 2, m_jfirst + 2, kcarteff + 2 );
    cof = 1.0 / ( h * h );
    if( op == '-' )
      cof = -cof;
  }
  double muy1, muy2, muy3, muy4, muz1, muz2, muz3, muz4;
  // Use Cartesian grid assumption
  for( int i = m_ifirst + 2 ; i <= m_ilast - 2 ; i++ )
    for( int j = m_jfirst + 2 ; j <= m_jlast - 2 ; j++ )
      for( int k = kcarteff + 2 ; k <= m_klast - 2 ; k++ )
      {
        double r1, r2, r3;
        mux1 = mu(i-1,j,k) * strx( i - 1 ) -
            tf * ( mu(i,j,k) * strx( i ) + mu(i-2,j,k) * strx( i - 2 ) );
        mux2 = mu(i-2,j,k) * strx( i - 2 ) + mu(i+1,j,k) * strx( i + 1 ) +
            3 * ( mu(i,j,k) * strx( i ) + mu(i-1,j,k) * strx( i - 1 ) );
        mux3 = mu(i-1,j,k) * strx( i - 1 ) + mu(i+2,j,k) * strx( i + 2 ) +
            3 * ( mu(i+1,j,k) * strx( i + 1 ) + mu(i,j,k) * strx( i ) );
        mux4 = mu(i+1,j,k) * strx( i + 1 ) -
            tf * ( mu(i,j,k) * strx( i ) + mu(i+2,j,k) * strx( i + 2 ) );

        muy1 = mu(i,j-1,k) * stry( j - 1 ) -
            tf * ( mu(i,j,k) * stry( j ) + mu(i,j-2,k) * stry( j - 2 ) );
        muy2 = mu(i,j-2,k) * stry( j - 2 ) + mu(i,j+1,k) * stry( j + 1 ) +
            3 * ( mu(i,j,k) * stry( j ) + mu(i,j-1,k) * stry( j - 1 ) );
        muy3 = mu(i,j-1,k) * stry( j - 1 ) + mu(i,j+2,k) * stry( j + 2 ) +
            3 * ( mu(i,j+1,k) * stry( j + 1 ) + mu(i,j,k) * stry( j ) );
        muy4 = mu(i,j+1,k) * stry( j + 1 ) -
            tf * ( mu(i,j,k) * stry( j ) + mu(i,j+2,k) * stry( j + 2 ) );

        muz1 = mu(i,j,k-1) * strz( k - 1 ) -
            tf * ( mu(i,j,k) * strz( k ) + mu(i,j,k-2) * strz( k - 2 ) );
        muz2 = mu(i,j,k-2) * strz( k - 2 ) + mu(i,j,k+1) * strz( k + 1 ) +
            3 * ( mu(i,j,k) * strz( k ) + mu(i,j,k-1) * strz( k - 1 ) );
        muz3 = mu(i,j,k-1) * strz( k - 1 ) + mu(i,j,k+2) * strz( k + 2 ) +
            3 * ( mu(i,j,k+1) * strz( k + 1 ) + mu(i,j,k) * strz( k ) );
        muz4 = mu(i,j,k+1) * strz( k + 1 ) -
            tf * ( mu(i,j,k) * strz( k ) + mu(i,j,k+2) * strz( k + 2 ) );

        // xx, yy, and zz derivatives:
        r1 = i6 * ( strx(i) * ( ( 2 * mux1 + la(i-1,j,k) * strx( i - 1 ) -
            tf * ( la(i,j,k) * strx( i ) + la(i-2,j,k) * strx( i - 2 ) ) ) *
            ( u(1,i-2,j,k) - u( 1, i, j, k ) ) +
            ( 2 * mux2 + la(i-2,j,k) * strx( i - 2 ) + la(i+1,j,k) * strx( i + 1 ) +
                3 * ( la(i,j,k) * strx( i ) + la(i-1,j,k) * strx( i - 1 ) ) ) *
                ( u(1,i-1,j,k) - u( 1, i, j, k ) ) +
            ( 2 * mux3 + la(i-1,j,k) * strx( i - 1 ) + la(i+2,j,k) * strx( i + 2 ) +
                3 * ( la(i+1,j,k) * strx( i + 1 ) + la(i,j,k) * strx( i ) ) ) *
                ( u(1,i+1,j,k) - u( 1, i, j, k ) ) +
            ( 2 * mux4 + la(i+1,j,k) * strx( i + 1 ) -
                tf * ( la(i,j,k) * strx( i ) + la(i+2,j,k) * strx( i + 2 ) ) ) *
                ( u(1,i+2,j,k) - u( 1, i, j, k ) ) ) + stry(j) * (
            muy1 * ( u(1,i,j-2,k) - u( 1, i, j, k ) ) +
                muy2 * ( u(1,i,j-1,k) - u( 1, i, j, k ) ) +
                muy3 * ( u(1,i,j+1,k) - u( 1, i, j, k ) ) +
                muy4 * ( u(1,i,j+2,k) - u( 1, i, j, k ) ) ) + strz(k) * (
            muz1 * ( u(1,i,j,k-2) - u( 1, i, j, k ) ) +
                muz2 * ( u(1,i,j,k-1) - u( 1, i, j, k ) ) +
                muz3 * ( u(1,i,j,k+1) - u( 1, i, j, k ) ) +
                muz4 * ( u(1,i,j,k+2) - u( 1, i, j, k ) ) ) );

        r2 = i6 * ( strx(i) * ( mux1 * ( u(2,i-2,j,k) - u( 2, i, j, k ) ) +
            mux2 * ( u(2,i-1,j,k) - u( 2, i, j, k ) ) +
            mux3 * ( u(2,i+1,j,k) - u( 2, i, j, k ) ) +
            mux4 * ( u(2,i+2,j,k) - u( 2, i, j, k ) ) ) + stry(j) * (
            ( 2 * muy1 + la(i,j-1,k) * stry( j - 1 ) -
                tf * ( la(i,j,k) * stry( j ) + la(i,j-2,k) * stry( j - 2 ) ) ) *
                ( u(2,i,j-2,k) - u( 2, i, j, k ) ) +
                ( 2 * muy2 + la(i,j-2,k) * stry( j - 2 ) + la(i,j+1,k) * stry( j + 1 ) +
                    3 * ( la(i,j,k) * stry( j ) + la(i,j-1,k) * stry( j - 1 ) ) ) *
                    ( u(2,i,j-1,k) - u( 2, i, j, k ) ) +
                ( 2 * muy3 + la(i,j-1,k) * stry( j - 1 ) + la(i,j+2,k) * stry( j + 2 ) +
                    3 * ( la(i,j+1,k) * stry( j + 1 ) + la(i,j,k) * stry( j ) ) ) *
                    ( u(2,i,j+1,k) - u( 2, i, j, k ) ) +
                ( 2 * muy4 + la(i,j+1,k) * stry( j + 1 ) -
                    tf * ( la(i,j,k) * stry( j ) + la(i,j+2,k) * stry( j + 2 ) ) ) *
                    ( u(2,i,j+2,k) - u( 2, i, j, k ) ) ) + strz(k) * (
            muz1 * ( u(2,i,j,k-2) - u( 2, i, j, k ) ) +
                muz2 * ( u(2,i,j,k-1) - u( 2, i, j, k ) ) +
                muz3 * ( u(2,i,j,k+1) - u( 2, i, j, k ) ) +
                muz4 * ( u(2,i,j,k+2) - u( 2, i, j, k ) ) ) );

        r3 = i6 * ( strx(i) * ( mux1 * ( u(3,i-2,j,k) - u( 3, i, j, k ) ) +
            mux2 * ( u(3,i-1,j,k) - u( 3, i, j, k ) ) +
            mux3 * ( u(3,i+1,j,k) - u( 3, i, j, k ) ) +
            mux4 * ( u(3,i+2,j,k) - u( 3, i, j, k ) ) ) + stry(j) * (
            muy1 * ( u(3,i,j-2,k) - u( 3, i, j, k ) ) +
                muy2 * ( u(3,i,j-1,k) - u( 3, i, j, k ) ) +
                muy3 * ( u(3,i,j+1,k) - u( 3, i, j, k ) ) +
                muy4 * ( u(3,i,j+2,k) - u( 3, i, j, k ) ) ) + strz(k) * (
            ( 2 * muz1 + la(i,j,k-1) * strz( k - 1 ) -
                tf * ( la(i,j,k) * strz( k ) + la(i,j,k-2) * strz( k - 2 ) ) ) *
                ( u(3,i,j,k-2) - u( 3, i, j, k ) ) +
                ( 2 * muz2 + la(i,j,k-2) * strz( k - 2 ) + la(i,j,k+1) * strz( k + 1 ) +
                    3 * ( la(i,j,k) * strz( k ) + la(i,j,k-1) * strz( k - 1 ) ) ) *
                    ( u(3,i,j,k-1) - u( 3, i, j, k ) ) +
                ( 2 * muz3 + la(i,j,k-1) * strz( k - 1 ) + la(i,j,k+2) * strz( k + 2 ) +
                    3 * ( la(i,j,k+1) * strz( k + 1 ) + la(i,j,k) * strz( k ) ) ) *
                    ( u(3,i,j,k+1) - u( 3, i, j, k ) ) +
                ( 2 * muz4 + la(i,j,k+1) * strz( k + 1 ) -
                    tf * ( la(i,j,k) * strz( k ) + la(i,j,k+2) * strz( k + 2 ) ) ) *
                    ( u(3,i,j,k+2) - u( 3, i, j, k ) ) ) );

// Mixed derivatives:
//   (la*v_y)_x
        r1 = r1 + strx(i) * stry( j ) *
            i144 * ( la(i-2,j,k) * ( u(2,i-2,j-2,k) - u( 2, i - 2, j + 2, k ) +
                8 * ( -u( 2, i - 2, j - 1, k ) + u( 2, i - 2, j + 1, k ) ) ) - 8 * (
            la(i-1,j,k) * ( u(2,i-1,j-2,k) - u( 2, i - 1, j + 2, k ) +
                8 * ( -u( 2, i - 1, j - 1, k ) + u( 2, i - 1, j + 1, k ) ) ) ) + 8 * (
            la(i+1,j,k) * ( u(2,i+1,j-2,k) - u( 2, i + 1, j + 2, k ) +
                8 * ( -u( 2, i + 1, j - 1, k ) + u( 2, i + 1, j + 1, k ) ) ) ) - (
            la(i+2,j,k) * ( u(2,i+2,j-2,k) - u( 2, i + 2, j + 2, k ) +
                8 * ( -u( 2, i + 2, j - 1, k ) + u( 2, i + 2, j + 1, k ) ) ) ) )
            //   (la*w_z)_x
            + strx(i) * strz( k ) *
                i144 * ( la(i-2,j,k) * ( u(3,i-2,j,k-2) - u( 3, i - 2, j, k + 2 ) +
                    8 * ( -u( 3, i - 2, j, k - 1 ) + u( 3, i - 2, j, k + 1 ) ) ) - 8 * (
                la(i-1,j,k) * ( u(3,i-1,j,k-2) - u( 3, i - 1, j, k + 2 ) +
                    8 * ( -u( 3, i - 1, j, k - 1 ) + u( 3, i - 1, j, k + 1 ) ) ) ) + 8 * (
                la(i+1,j,k) * ( u(3,i+1,j,k-2) - u( 3, i + 1, j, k + 2 ) +
                    8 * ( -u( 3, i + 1, j, k - 1 ) + u( 3, i + 1, j, k + 1 ) ) ) ) - (
                la(i+2,j,k) * ( u(3,i+2,j,k-2) - u( 3, i + 2, j, k + 2 ) +
                    8 * ( -u( 3, i + 2, j, k - 1 ) + u( 3, i + 2, j, k + 1 ) ) ) ) )
            //  (mu*v_x)_y
            + strx(i) * stry( j ) *
                i144 * ( mu(i,j-2,k) * ( u(2,i-2,j-2,k) - u( 2, i + 2, j - 2, k ) +
                    8 * ( -u( 2, i - 1, j - 2, k ) + u( 2, i + 1, j - 2, k ) ) ) - 8 * (
                mu(i,j-1,k) * ( u(2,i-2,j-1,k) - u( 2, i + 2, j - 1, k ) +
                    8 * ( -u( 2, i - 1, j - 1, k ) + u( 2, i + 1, j - 1, k ) ) ) ) + 8 * (
                mu(i,j+1,k) * ( u(2,i-2,j+1,k) - u( 2, i + 2, j + 1, k ) +
                    8 * ( -u( 2, i - 1, j + 1, k ) + u( 2, i + 1, j + 1, k ) ) ) ) - (
                mu(i,j+2,k) * ( u(2,i-2,j+2,k) - u( 2, i + 2, j + 2, k ) +
                    8 * ( -u( 2, i - 1, j + 2, k ) + u( 2, i + 1, j + 2, k ) ) ) ) )
            //   (mu*w_x)_z
            + strx(i) * strz( k ) *
                i144 * ( mu(i,j,k-2) * ( u(3,i-2,j,k-2) - u( 3, i + 2, j, k - 2 ) +
                    8 * ( -u( 3, i - 1, j, k - 2 ) + u( 3, i + 1, j, k - 2 ) ) ) - 8 * (
                mu(i,j,k-1) * ( u(3,i-2,j,k-1) - u( 3, i + 2, j, k - 1 ) +
                    8 * ( -u( 3, i - 1, j, k - 1 ) + u( 3, i + 1, j, k - 1 ) ) ) ) + 8 * (
                mu(i,j,k+1) * ( u(3,i-2,j,k+1) - u( 3, i + 2, j, k + 1 ) +
                    8 * ( -u( 3, i - 1, j, k + 1 ) + u( 3, i + 1, j, k + 1 ) ) ) ) - (
                mu(i,j,k+2) * ( u(3,i-2,j,k+2) - u( 3, i + 2, j, k + 2 ) +
                    8 * ( -u( 3, i - 1, j, k + 2 ) + u( 3, i + 1, j, k + 2 ) ) ) ) );

//  (mu*u_y)_x
        r2 = r2 + strx(i) * stry( j ) *
            i144 * ( mu(i-2,j,k) * ( u(1,i-2,j-2,k) - u( 1, i - 2, j + 2, k ) +
                8 * ( -u( 1, i - 2, j - 1, k ) + u( 1, i - 2, j + 1, k ) ) ) - 8 * (
            mu(i-1,j,k) * ( u(1,i-1,j-2,k) - u( 1, i - 1, j + 2, k ) +
                8 * ( -u( 1, i - 1, j - 1, k ) + u( 1, i - 1, j + 1, k ) ) ) ) + 8 * (
            mu(i+1,j,k) * ( u(1,i+1,j-2,k) - u( 1, i + 1, j + 2, k ) +
                8 * ( -u( 1, i + 1, j - 1, k ) + u( 1, i + 1, j + 1, k ) ) ) ) - (
            mu(i+2,j,k) * ( u(1,i+2,j-2,k) - u( 1, i + 2, j + 2, k ) +
                8 * ( -u( 1, i + 2, j - 1, k ) + u( 1, i + 2, j + 1, k ) ) ) ) )
            // (la*u_x)_y
            + strx(i) * stry( j ) *
                i144 * ( la(i,j-2,k) * ( u(1,i-2,j-2,k) - u( 1, i + 2, j - 2, k ) +
                    8 * ( -u( 1, i - 1, j - 2, k ) + u( 1, i + 1, j - 2, k ) ) ) - 8 * (
                la(i,j-1,k) * ( u(1,i-2,j-1,k) - u( 1, i + 2, j - 1, k ) +
                    8 * ( -u( 1, i - 1, j - 1, k ) + u( 1, i + 1, j - 1, k ) ) ) ) + 8 * (
                la(i,j+1,k) * ( u(1,i-2,j+1,k) - u( 1, i + 2, j + 1, k ) +
                    8 * ( -u( 1, i - 1, j + 1, k ) + u( 1, i + 1, j + 1, k ) ) ) ) - (
                la(i,j+2,k) * ( u(1,i-2,j+2,k) - u( 1, i + 2, j + 2, k ) +
                    8 * ( -u( 1, i - 1, j + 2, k ) + u( 1, i + 1, j + 2, k ) ) ) ) )
            // (la*w_z)_y
            + stry(j) * strz( k ) *
                i144 * ( la(i,j-2,k) * ( u(3,i,j-2,k-2) - u( 3, i, j - 2, k + 2 ) +
                    8 * ( -u( 3, i, j - 2, k - 1 ) + u( 3, i, j - 2, k + 1 ) ) ) - 8 * (
                la(i,j-1,k) * ( u(3,i,j-1,k-2) - u( 3, i, j - 1, k + 2 ) +
                    8 * ( -u( 3, i, j - 1, k - 1 ) + u( 3, i, j - 1, k + 1 ) ) ) ) + 8 * (
                la(i,j+1,k) * ( u(3,i,j+1,k-2) - u( 3, i, j + 1, k + 2 ) +
                    8 * ( -u( 3, i, j + 1, k - 1 ) + u( 3, i, j + 1, k + 1 ) ) ) ) - (
                la(i,j+2,k) * ( u(3,i,j+2,k-2) - u( 3, i, j + 2, k + 2 ) +
                    8 * ( -u( 3, i, j + 2, k - 1 ) + u( 3, i, j + 2, k + 1 ) ) ) ) )
            // (mu*w_y)_z
            + stry(j) * strz( k ) *
                i144 * ( mu(i,j,k-2) * ( u(3,i,j-2,k-2) - u( 3, i, j + 2, k - 2 ) +
                    8 * ( -u( 3, i, j - 1, k - 2 ) + u( 3, i, j + 1, k - 2 ) ) ) - 8 * (
                mu(i,j,k-1) * ( u(3,i,j-2,k-1) - u( 3, i, j + 2, k - 1 ) +
                    8 * ( -u( 3, i, j - 1, k - 1 ) + u( 3, i, j + 1, k - 1 ) ) ) ) + 8 * (
                mu(i,j,k+1) * ( u(3,i,j-2,k+1) - u( 3, i, j + 2, k + 1 ) +
                    8 * ( -u( 3, i, j - 1, k + 1 ) + u( 3, i, j + 1, k + 1 ) ) ) ) - (
                mu(i,j,k+2) * ( u(3,i,j-2,k+2) - u( 3, i, j + 2, k + 2 ) +
                    8 * ( -u( 3, i, j - 1, k + 2 ) + u( 3, i, j + 1, k + 2 ) ) ) ) );
        //  (mu*u_z)_x
        r3 = r3 + strx(i) * strz( k ) *
            i144 * ( mu(i-2,j,k) * ( u(1,i-2,j,k-2) - u( 1, i - 2, j, k + 2 ) +
                8 * ( -u( 1, i - 2, j, k - 1 ) + u( 1, i - 2, j, k + 1 ) ) ) - 8 * (
            mu(i-1,j,k) * ( u(1,i-1,j,k-2) - u( 1, i - 1, j, k + 2 ) +
                8 * ( -u( 1, i - 1, j, k - 1 ) + u( 1, i - 1, j, k + 1 ) ) ) ) + 8 * (
            mu(i+1,j,k) * ( u(1,i+1,j,k-2) - u( 1, i + 1, j, k + 2 ) +
                8 * ( -u( 1, i + 1, j, k - 1 ) + u( 1, i + 1, j, k + 1 ) ) ) ) - (
            mu(i+2,j,k) * ( u(1,i+2,j,k-2) - u( 1, i + 2, j, k + 2 ) +
                8 * ( -u( 1, i + 2, j, k - 1 ) + u( 1, i + 2, j, k + 1 ) ) ) ) )
            //  (mu*v_z)_y
            + stry(j) * strz( k ) *
                i144 * ( mu(i,j-2,k) * ( u(2,i,j-2,k-2) - u( 2, i, j - 2, k + 2 ) +
                    8 * ( -u( 2, i, j - 2, k - 1 ) + u( 2, i, j - 2, k + 1 ) ) ) - 8 * (
                mu(i,j-1,k) * ( u(2,i,j-1,k-2) - u( 2, i, j - 1, k + 2 ) +
                    8 * ( -u( 2, i, j - 1, k - 1 ) + u( 2, i, j - 1, k + 1 ) ) ) ) + 8 * (
                mu(i,j+1,k) * ( u(2,i,j+1,k-2) - u( 2, i, j + 1, k + 2 ) +
                    8 * ( -u( 2, i, j + 1, k - 1 ) + u( 2, i, j + 1, k + 1 ) ) ) ) - (
                mu(i,j+2,k) * ( u(2,i,j+2,k-2) - u( 2, i, j + 2, k + 2 ) +
                    8 * ( -u( 2, i, j + 2, k - 1 ) + u( 2, i, j + 2, k + 1 ) ) ) ) )
            //   (la*u_x)_z
            + strx(i) * strz( k ) *
                i144 * ( la(i,j,k-2) * ( u(1,i-2,j,k-2) - u( 1, i + 2, j, k - 2 ) +
                    8 * ( -u( 1, i - 1, j, k - 2 ) + u( 1, i + 1, j, k - 2 ) ) ) - 8 * (
                la(i,j,k-1) * ( u(1,i-2,j,k-1) - u( 1, i + 2, j, k - 1 ) +
                    8 * ( -u( 1, i - 1, j, k - 1 ) + u( 1, i + 1, j, k - 1 ) ) ) ) + 8 * (
                la(i,j,k+1) * ( u(1,i-2,j,k+1) - u( 1, i + 2, j, k + 1 ) +
                    8 * ( -u( 1, i - 1, j, k + 1 ) + u( 1, i + 1, j, k + 1 ) ) ) ) - (
                la(i,j,k+2) * ( u(1,i-2,j,k+2) - u( 1, i + 2, j, k + 2 ) +
                    8 * ( -u( 1, i - 1, j, k + 2 ) + u( 1, i + 1, j, k + 2 ) ) ) ) )
            //  (la*v_y)_z
            + stry(j) * strz( k ) *
                i144 * ( la(i,j,k-2) * ( u(2,i,j-2,k-2) - u( 2, i, j + 2, k - 2 ) +
                    8 * ( -u( 2, i, j - 1, k - 2 ) + u( 2, i, j + 1, k - 2 ) ) ) - 8 * (
                la(i,j,k-1) * ( u(2,i,j-2,k-1) - u( 2, i, j + 2, k - 1 ) +
                    8 * ( -u( 2, i, j - 1, k - 1 ) + u( 2, i, j + 1, k - 1 ) ) ) ) + 8 * (
                la(i,j,k+1) * ( u(2,i,j-2,k+1) - u( 2, i, j + 2, k + 1 ) +
                    8 * ( -u( 2, i, j - 1, k + 1 ) + u( 2, i, j + 1, k + 1 ) ) ) ) - (
                la(i,j,k+2) * ( u(2,i,j-2,k+2) - u( 2, i, j + 2, k + 2 ) +
                    8 * ( -u( 2, i, j - 1, k + 2 ) + u( 2, i, j + 1, k + 2 ) ) ) ) );
        lu(1,i,j,k) = a1 * lu( 1, i, j, k ) + cof * r1;
        lu(2,i,j,k) = a1 * lu( 2, i, j, k ) + cof * r2;
        lu(3,i,j,k) = a1 * lu( 3, i, j, k ) + cof * r3;
      }
#undef u
#undef mu
#undef la
#undef met
#undef jac
#undef lu
#undef acof
#undef ghcof
#undef bope
#undef SQR
#undef strx
#undef stry
#undef strz
}

//-----------------------------------------------------------------------
void SW4Solver::GetStencilCoefficients( realT* _acof, realT* _ghcof,
                                        realT* _bope,
                                        realT* _sbop )
{
#define acof(q,k,m) (_acof[q-1+6*(k-1)+48*(m-1)])
#define ghcof(k) (_ghcof[k-1])
#define bope(q,k) (_bope[q-1+6*(k-1)])

  ghcof(1) = 12.0 / 17;
  ghcof(2) = 0;
  ghcof(3) = 0;
  ghcof(4) = 0;
  ghcof(5) = 0;
  ghcof(6) = 0;

  acof(1,1,1) = 104.0 / 289.0;
  acof(1,1,2) = -2476335.0 / 2435692.0;
  acof(1,1,3) = -16189.0 / 84966.0;
  acof(1,1,4) = -9.0 / 3332.0;
  acof(1,1,5) = 0;
  acof(1,1,6) = 0;
  acof(1,1,7) = 0;
  acof(1,1,8) = 0;
  acof(1,2,1) = -516.0 / 289.0;
  acof(1,2,2) = 544521.0 / 1217846.0;
  acof(1,2,3) = 2509879.0 / 3653538.0;
  acof(1,2,4) = 0;
  acof(1,2,5) = 0;
  acof(1,2,6) = 0;
  acof(1,2,7) = 0;
  acof(1,2,8) = 0;
  acof(1,3,1) = 312.0 / 289.0;
  acof(1,3,2) = 1024279.0 / 2435692.0;
  acof(1,3,3) = -687797.0 / 1217846.0;
  acof(1,3,4) = 177.0 / 3332.0;
  acof(1,3,5) = 0;
  acof(1,3,6) = 0;
  acof(1,3,7) = 0;
  acof(1,3,8) = 0;
  acof(1,4,1) = -104.0 / 289.0;
  acof(1,4,2) = 181507.0 / 1217846.0;
  acof(1,4,3) = 241309.0 / 3653538.0;
  acof(1,4,4) = 0;
  acof(1,4,5) = 0;
  acof(1,4,6) = 0;
  acof(1,4,7) = 0;
  acof(1,4,8) = 0;
  acof(1,5,1) = 0;
  acof(1,5,2) = 0;
  acof(1,5,3) = 5.0 / 2193.0;
  acof(1,5,4) = -48.0 / 833.0;
  acof(1,5,5) = 0;
  acof(1,5,6) = 0;
  acof(1,5,7) = 0;
  acof(1,5,8) = 0;
  acof(1,6,1) = 0;
  acof(1,6,2) = 0;
  acof(1,6,3) = 0;
  acof(1,6,4) = 6.0 / 833.0;
  acof(1,6,5) = 0;
  acof(1,6,6) = 0;
  acof(1,6,7) = 0;
  acof(1,6,8) = 0;
  acof(1,7,1) = 0;
  acof(1,7,2) = 0;
  acof(1,7,3) = 0;
  acof(1,7,4) = 0;
  acof(1,7,5) = 0;
  acof(1,7,6) = 0;
  acof(1,7,7) = 0;
  acof(1,7,8) = 0;
  acof(1,8,1) = 0;
  acof(1,8,2) = 0;
  acof(1,8,3) = 0;
  acof(1,8,4) = 0;
  acof(1,8,5) = 0;
  acof(1,8,6) = 0;
  acof(1,8,7) = 0;
  acof(1,8,8) = 0;
  acof(2,1,1) = 12.0 / 17.0;
  acof(2,1,2) = 544521.0 / 4226642.0;
  acof(2,1,3) = 2509879.0 / 12679926.0;
  acof(2,1,4) = 0;
  acof(2,1,5) = 0;
  acof(2,1,6) = 0;
  acof(2,1,7) = 0;
  acof(2,1,8) = 0;
  acof(2,2,1) = -59.0 / 68.0;
  acof(2,2,2) = -1633563.0 / 4226642.0;
  acof(2,2,3) = -21510077.0 / 25359852.0;
  acof(2,2,4) = -12655.0 / 372939.0;
  acof(2,2,5) = 0;
  acof(2,2,6) = 0;
  acof(2,2,7) = 0;
  acof(2,2,8) = 0;
  acof(2,3,1) = 2.0 / 17.0;
  acof(2,3,2) = 1633563.0 / 4226642.0;
  acof(2,3,3) = 2565299.0 / 4226642.0;
  acof(2,3,4) = 40072.0 / 372939.0;
  acof(2,3,5) = 0;
  acof(2,3,6) = 0;
  acof(2,3,7) = 0;
  acof(2,3,8) = 0;
  acof(2,4,1) = 3.0 / 68.0;
  acof(2,4,2) = -544521.0 / 4226642.0;
  acof(2,4,3) = 987685.0 / 25359852.0;
  acof(2,4,4) = -14762.0 / 124313.0;
  acof(2,4,5) = 0;
  acof(2,4,6) = 0;
  acof(2,4,7) = 0;
  acof(2,4,8) = 0;
  acof(2,5,1) = 0;
  acof(2,5,2) = 0;
  acof(2,5,3) = 1630.0 / 372939.0;
  acof(2,5,4) = 18976.0 / 372939.0;
  acof(2,5,5) = 0;
  acof(2,5,6) = 0;
  acof(2,5,7) = 0;
  acof(2,5,8) = 0;
  acof(2,6,1) = 0;
  acof(2,6,2) = 0;
  acof(2,6,3) = 0;
  acof(2,6,4) = -1.0 / 177.0;
  acof(2,6,5) = 0;
  acof(2,6,6) = 0;
  acof(2,6,7) = 0;
  acof(2,6,8) = 0;
  acof(2,7,1) = 0;
  acof(2,7,2) = 0;
  acof(2,7,3) = 0;
  acof(2,7,4) = 0;
  acof(2,7,5) = 0;
  acof(2,7,6) = 0;
  acof(2,7,7) = 0;
  acof(2,7,8) = 0;
  acof(2,8,1) = 0;
  acof(2,8,2) = 0;
  acof(2,8,3) = 0;
  acof(2,8,4) = 0;
  acof(2,8,5) = 0;
  acof(2,8,6) = 0;
  acof(2,8,7) = 0;
  acof(2,8,8) = 0;
  acof(3,1,1) = -96.0 / 731.0;
  acof(3,1,2) = 1024279.0 / 6160868.0;
  acof(3,1,3) = -687797.0 / 3080434.0;
  acof(3,1,4) = 177.0 / 8428.0;
  acof(3,1,5) = 0;
  acof(3,1,6) = 0;
  acof(3,1,7) = 0;
  acof(3,1,8) = 0;
  acof(3,2,1) = 118.0 / 731.0;
  acof(3,2,2) = 1633563.0 / 3080434.0;
  acof(3,2,3) = 2565299.0 / 3080434.0;
  acof(3,2,4) = 40072.0 / 271803.0;
  acof(3,2,5) = 0;
  acof(3,2,6) = 0;
  acof(3,2,7) = 0;
  acof(3,2,8) = 0;
  acof(3,3,1) = -16.0 / 731.0;
  acof(3,3,2) = -5380447.0 / 6160868.0;
  acof(3,3,3) = -3569115.0 / 3080434.0;
  acof(3,3,4) = -331815.0 / 362404.0;
  acof(3,3,5) = -283.0 / 6321.0;
  acof(3,3,6) = 0;
  acof(3,3,7) = 0;
  acof(3,3,8) = 0;
  acof(3,4,1) = -6.0 / 731.0;
  acof(3,4,2) = 544521.0 / 3080434.0;
  acof(3,4,3) = 2193521.0 / 3080434.0;
  acof(3,4,4) = 8065.0 / 12943.0;
  acof(3,4,5) = 381.0 / 2107.0;
  acof(3,4,6) = 0;
  acof(3,4,7) = 0;
  acof(3,4,8) = 0;
  acof(3,5,1) = 0;
  acof(3,5,2) = 0;
  acof(3,5,3) = -14762.0 / 90601.0;
  acof(3,5,4) = 32555.0 / 271803.0;
  acof(3,5,5) = -283.0 / 2107.0;
  acof(3,5,6) = 0;
  acof(3,5,7) = 0;
  acof(3,5,8) = 0;
  acof(3,6,1) = 0;
  acof(3,6,2) = 0;
  acof(3,6,3) = 0;
  acof(3,6,4) = 9.0 / 2107.0;
  acof(3,6,5) = -11.0 / 6321.0;
  acof(3,6,6) = 0;
  acof(3,6,7) = 0;
  acof(3,6,8) = 0;
  acof(3,7,1) = 0;
  acof(3,7,2) = 0;
  acof(3,7,3) = 0;
  acof(3,7,4) = 0;
  acof(3,7,5) = 0;
  acof(3,7,6) = 0;
  acof(3,7,7) = 0;
  acof(3,7,8) = 0;
  acof(3,8,1) = 0;
  acof(3,8,2) = 0;
  acof(3,8,3) = 0;
  acof(3,8,4) = 0;
  acof(3,8,5) = 0;
  acof(3,8,6) = 0;
  acof(3,8,7) = 0;
  acof(3,8,8) = 0;
  acof(4,1,1) = -36.0 / 833.0;
  acof(4,1,2) = 181507.0 / 3510262.0;
  acof(4,1,3) = 241309.0 / 10530786.0;
  acof(4,1,4) = 0;
  acof(4,1,5) = 0;
  acof(4,1,6) = 0;
  acof(4,1,7) = 0;
  acof(4,1,8) = 0;
  acof(4,2,1) = 177.0 / 3332.0;
  acof(4,2,2) = -544521.0 / 3510262.0;
  acof(4,2,3) = 987685.0 / 21061572.0;
  acof(4,2,4) = -14762.0 / 103243.0;
  acof(4,2,5) = 0;
  acof(4,2,6) = 0;
  acof(4,2,7) = 0;
  acof(4,2,8) = 0;
  acof(4,3,1) = -6.0 / 833.0;
  acof(4,3,2) = 544521.0 / 3510262.0;
  acof(4,3,3) = 2193521.0 / 3510262.0;
  acof(4,3,4) = 8065.0 / 14749.0;
  acof(4,3,5) = 381.0 / 2401.0;
  acof(4,3,6) = 0;
  acof(4,3,7) = 0;
  acof(4,3,8) = 0;
  acof(4,4,1) = -9.0 / 3332.0;
  acof(4,4,2) = -181507.0 / 3510262.0;
  acof(4,4,3) = -2647979.0 / 3008796.0;
  acof(4,4,4) = -80793.0 / 103243.0;
  acof(4,4,5) = -1927.0 / 2401.0;
  acof(4,4,6) = -2.0 / 49.0;
  acof(4,4,7) = 0;
  acof(4,4,8) = 0;
  acof(4,5,1) = 0;
  acof(4,5,2) = 0;
  acof(4,5,3) = 57418.0 / 309729.0;
  acof(4,5,4) = 51269.0 / 103243.0;
  acof(4,5,5) = 1143.0 / 2401.0;
  acof(4,5,6) = 8.0 / 49.0;
  acof(4,5,7) = 0;
  acof(4,5,8) = 0;
  acof(4,6,1) = 0;
  acof(4,6,2) = 0;
  acof(4,6,3) = 0;
  acof(4,6,4) = -283.0 / 2401.0;
  acof(4,6,5) = 403.0 / 2401.0;
  acof(4,6,6) = -6.0 / 49.0;
  acof(4,6,7) = 0;
  acof(4,6,8) = 0;
  acof(4,7,1) = 0;
  acof(4,7,2) = 0;
  acof(4,7,3) = 0;
  acof(4,7,4) = 0;
  acof(4,7,5) = 0;
  acof(4,7,6) = 0;
  acof(4,7,7) = 0;
  acof(4,7,8) = 0;
  acof(4,8,1) = 0;
  acof(4,8,2) = 0;
  acof(4,8,3) = 0;
  acof(4,8,4) = 0;
  acof(4,8,5) = 0;
  acof(4,8,6) = 0;
  acof(4,8,7) = 0;
  acof(4,8,8) = 0;
  acof(5,1,1) = 0;
  acof(5,1,2) = 0;
  acof(5,1,3) = 5.0 / 6192.0;
  acof(5,1,4) = -1.0 / 49.0;
  acof(5,1,5) = 0;
  acof(5,1,6) = 0;
  acof(5,1,7) = 0;
  acof(5,1,8) = 0;
  acof(5,2,1) = 0;
  acof(5,2,2) = 0;
  acof(5,2,3) = 815.0 / 151704.0;
  acof(5,2,4) = 1186.0 / 18963.0;
  acof(5,2,5) = 0;
  acof(5,2,6) = 0;
  acof(5,2,7) = 0;
  acof(5,2,8) = 0;
  acof(5,3,1) = 0;
  acof(5,3,2) = 0;
  acof(5,3,3) = -7381.0 / 50568.0;
  acof(5,3,4) = 32555.0 / 303408.0;
  acof(5,3,5) = -283.0 / 2352.0;
  acof(5,3,6) = 0;
  acof(5,3,7) = 0;
  acof(5,3,8) = 0;
  acof(5,4,1) = 0;
  acof(5,4,2) = 0;
  acof(5,4,3) = 28709.0 / 151704.0;
  acof(5,4,4) = 51269.0 / 101136.0;
  acof(5,4,5) = 381.0 / 784.0;
  acof(5,4,6) = 1.0 / 6.0;
  acof(5,4,7) = 0;
  acof(5,4,8) = 0;
  acof(5,5,1) = 0;
  acof(5,5,2) = 0;
  acof(5,5,3) = -349.0 / 7056.0;
  acof(5,5,4) = -247951.0 / 303408.0;
  acof(5,5,5) = -577.0 / 784.0;
  acof(5,5,6) = -5.0 / 6.0;
  acof(5,5,7) = -1.0 / 24.0;
  acof(5,5,8) = 0;
  acof(5,6,1) = 0;
  acof(5,6,2) = 0;
  acof(5,6,3) = 0;
  acof(5,6,4) = 1135.0 / 7056.0;
  acof(5,6,5) = 1165.0 / 2352.0;
  acof(5,6,6) = 1.0 / 2.0;
  acof(5,6,7) = 1.0 / 6.0;
  acof(5,6,8) = 0;
  acof(5,7,1) = 0;
  acof(5,7,2) = 0;
  acof(5,7,3) = 0;
  acof(5,7,4) = 0;
  acof(5,7,5) = -1.0 / 8.0;
  acof(5,7,6) = 1.0 / 6.0;
  acof(5,7,7) = -1.0 / 8.0;
  acof(5,7,8) = 0;
  acof(5,8,1) = 0;
  acof(5,8,2) = 0;
  acof(5,8,3) = 0;
  acof(5,8,4) = 0;
  acof(5,8,5) = 0;
  acof(5,8,6) = 0;
  acof(5,8,7) = 0;
  acof(5,8,8) = 0;
  acof(6,1,1) = 0;
  acof(6,1,2) = 0;
  acof(6,1,3) = 0;
  acof(6,1,4) = 1.0 / 392.0;
  acof(6,1,5) = 0;
  acof(6,1,6) = 0;
  acof(6,1,7) = 0;
  acof(6,1,8) = 0;
  acof(6,2,1) = 0;
  acof(6,2,2) = 0;
  acof(6,2,3) = 0;
  acof(6,2,4) = -1.0 / 144.0;
  acof(6,2,5) = 0;
  acof(6,2,6) = 0;
  acof(6,2,7) = 0;
  acof(6,2,8) = 0;
  acof(6,3,1) = 0;
  acof(6,3,2) = 0;
  acof(6,3,3) = 0;
  acof(6,3,4) = 3.0 / 784.0;
  acof(6,3,5) = -11.0 / 7056.0;
  acof(6,3,6) = 0;
  acof(6,3,7) = 0;
  acof(6,3,8) = 0;
  acof(6,4,1) = 0;
  acof(6,4,2) = 0;
  acof(6,4,3) = 0;
  acof(6,4,4) = -283.0 / 2352.0;
  acof(6,4,5) = 403.0 / 2352.0;
  acof(6,4,6) = -1.0 / 8.0;
  acof(6,4,7) = 0;
  acof(6,4,8) = 0;
  acof(6,5,1) = 0;
  acof(6,5,2) = 0;
  acof(6,5,3) = 0;
  acof(6,5,4) = 1135.0 / 7056.0;
  acof(6,5,5) = 1165.0 / 2352.0;
  acof(6,5,6) = 1.0 / 2.0;
  acof(6,5,7) = 1.0 / 6.0;
  acof(6,5,8) = 0;
  acof(6,6,1) = 0;
  acof(6,6,2) = 0;
  acof(6,6,3) = 0;
  acof(6,6,4) = -47.0 / 1176.0;
  acof(6,6,5) = -5869.0 / 7056.0;
  acof(6,6,6) = -3.0 / 4.0;
  acof(6,6,7) = -5.0 / 6.0;
  acof(6,6,8) = -1.0 / 24.0;
  acof(6,7,1) = 0;
  acof(6,7,2) = 0;
  acof(6,7,3) = 0;
  acof(6,7,4) = 0;
  acof(6,7,5) = 1.0 / 6.0;
  acof(6,7,6) = 1.0 / 2.0;
  acof(6,7,7) = 1.0 / 2.0;
  acof(6,7,8) = 1.0 / 6.0;
  acof(6,8,1) = 0;
  acof(6,8,2) = 0;
  acof(6,8,3) = 0;
  acof(6,8,4) = 0;
  acof(6,8,5) = 0;
  acof(6,8,6) = -1.0 / 8.0;
  acof(6,8,7) = 1.0 / 6.0;
  acof(6,8,8) = -1.0 / 8.0;

  bope(1,1) = -24.0 / 17.0;
  bope(1,2) = 59.0 / 34.0;
  bope(1,3) = -4.0 / 17.0;
  bope(1,4) = -3.0 / 34.0;
  bope(1,5) = 0;
  bope(1,6) = 0;
  bope(1,7) = 0;
  bope(1,8) = 0;
  bope(2,1) = -1.0 / 2.0;
  bope(2,2) = 0;
  bope(2,3) = 1.0 / 2.0;
  bope(2,4) = 0;
  bope(2,5) = 0;
  bope(2,6) = 0;
  bope(2,7) = 0;
  bope(2,8) = 0;
  bope(3,1) = 4.0 / 43.0;
  bope(3,2) = -59.0 / 86.0;
  bope(3,3) = 0;
  bope(3,4) = 59.0 / 86.0;
  bope(3,5) = -4.0 / 43.0;
  bope(3,6) = 0;
  bope(3,7) = 0;
  bope(3,8) = 0;
  bope(4,1) = 3.0 / 98.0;
  bope(4,2) = 0;
  bope(4,3) = -59.0 / 98.0;
  bope(4,4) = 0;
  bope(4,5) = 32.0 / 49.0;
  bope(4,6) = -4.0 / 49.0;
  bope(4,7) = 0;
  bope(4,8) = 0;
  bope(5,1) = 0;
  bope(5,2) = 0;
  realT d4a = 2.0 / 3;
  realT d4b = -1.0 / 12;
  bope(5,3) = -d4b;
  bope(5,4) = -d4a;
  bope(5,5) = 0;
  bope(5,6) = d4a;
  bope(5,7) = d4b;
  bope(5,8) = 0;
  bope(6,1) = 0;
  bope(6,2) = 0;
  bope(6,3) = 0;
  bope(6,4) = -d4b;
  bope(6,5) = -d4a;
  bope(6,6) = 0;
  bope(6,7) = d4a;
  bope(6,8) = d4b;
#undef acof
#undef ghcof
#undef bope
  _sbop[0] = -1.0 / 4;
  _sbop[1] = -5.0 / 6;
  _sbop[2] = 3.0 / 2;
  _sbop[3] = -1.0 / 2;
  _sbop[4] = 1.0 / 12;
}

/// Register solver in the solver factory
REGISTER_SOLVER( SW4Solver )

