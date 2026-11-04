// Tests use the production mapper/neighbor/cleanup/event bodies via CPU adapters.
static void close(double actual,double expected,double tolerance=3.e-11)
{
  if(std::abs(actual-expected)>tolerance*std::max({1.,std::abs(actual),std::abs(expected)}))
  {std::cerr<<"Mismatch "<<actual<<" != "<<expected<<"\n";std::abort();}
}
static void grid(int (&map)[9][9][9])
{for(int i=0;i<9;++i)for(int j=0;j<9;++j)for(int k=0;k<9;++k)map[i][j][k]=(i*9+j)*9+k;}
static void coordinates(int node,double const (&origin)[3],double const (&h)[3],double (&x)[3])
{x[0]=origin[0]+(node/81)*h[0];x[1]=origin[1]+((node/9)%9)*h[1];x[2]=origin[2]+(node%9)*h[2];}
static void surfaceTests()
{
  constexpr double dr=.02,R=1.;
  int count=0;
  // Uniform fill, including queries within one sample width of the axis.
  for(double radius:{.01,.15,.6,1.2,3.2})
  {
    std::vector<SurfaceSample> samples;
    for(int i=0;i<250;++i)for(int j=-60;j<=60;++j)
    {
      double const r=(i+.5)*dr,y=j*dr;
      if(std::hypot(r-radius,y)<R) samples.push_back({r,y,axisymmetric::volume(dr*dr,r)});
    }
    assert(detect(samples,radius,0,R)==0);
    // The old ring-weighted centroid has a spurious outward displacement.
    if(std::abs(radius-1.2)<1.e-12)
    {
      double w=0,offset=0;
      for(auto const &p:samples) {double const q=std::hypot(p.r-radius,p.y)/R;
        double const a=p.volume*(1-q)*(1-q)*(1+2*q);w+=a;offset+=a*(p.r-radius);}
      assert(offset/w>.08*R);close(offset/w,5./(42*radius),3.e-5);
    }
    // The real axial free surface remains detectable even at the axis.
    std::vector<SurfaceSample> top;
    for(auto const &p:samples)if(p.y<=.05)top.push_back(p);
    assert(detect(top,radius,0,R)==1);
    ++count;
  }
  // Cylindrical exterior, a cavity wall and mixed sample areas.
  for(bool outer:{false,true})
  {
    std::vector<SurfaceSample> samples;
    for(int i=0;i<250;++i)for(int j=-60;j<=60;++j)
    {
      double const r=(i+.5)*dr,y=j*dr;
      if(std::hypot(r-2.,y)<R && (outer?r<2.05:r>1.95))
      {double const v=axisymmetric::volume(dr*dr,r);samples.push_back({r,y,v*.3});samples.push_back({r,y,v*.7});}
    }
    assert(detect(samples,2,0,R)==1);++count;
  }
  assert(detect({},.1,0,R)==0);assert(detect({},.1,0,R,3)==3);
  std::cout<<count<<" production surface neighborhoods: uniform-fill bias removed, axis hidden, real free surfaces retained.\n";
}
static void cleanupTests()
{
  int map[9][9][9];grid(map);double const origin[3]={-1,-1,0},h[3]={1,1,1};int const cells[3]={8,8,0};
  Kernel k;k.m_axisymmetric=true;
  int count=0;
  for(auto type:{ParticleType::CPDI,ParticleType::CPTI})for(int mode=0;mode<3;++mode)
  {
    k.m_cpdiDomainScaling=mode>0;k.m_cpdiDomainScalingType=mode==1?mpm::CPDIDomainScalingTypeOption::Homel:mpm::CPDIDomainScalingTypeOption::VPHencky;
    double center[3]={.1,3,0},g[3][3]={{3.6,.1,0},{.2,1,0},{0,0,1}},original[3][3];
    LvArray::tensorOps::copy<3,3>(original,g);
    int nodes[16];double w[16],grad[16][3];
    assert(k.computeParticleMapping<2>(type,center,g,map,origin,h,cells,nodes,w,grad));
    assert(!k.cleanup<2>(type,center,g));
    for(int a=0;a<3;++a)for(int d=0;d<3;++d)close(g[a][d],original[a][d]);
    ++count;
  }
  // True escapes are still removed; Cartesian cleanup still uses its own projection.
  k.m_cpdiDomainScaling=0;
  double g[3][3]={{3.6,0,0},{0,1,0},{0,0,1}},center[3]={.1,3,0};
  k.m_axisymmetric=false;assert(k.cleanup<2>(ParticleType::CPDI,center,g));
  k.m_axisymmetric=true;center[0]=6.5;assert(k.cleanup<2>(ParticleType::CPDI,center,g));
  center[0]=3;center[1]=6.9;assert(k.cleanup<2>(ParticleType::CPDI,center,g));
  center[1]=6.5;assert(!k.cleanup<2>(ParticleType::CPDI,center,g)); // corner exactly at buffer edge
  std::cout<<count<<" production cleanup cases: mapping fit agrees, stored geometry unchanged, true escapes retained.\n";
}
static void basisTests()
{
  int map[9][9][9];grid(map);double const origin[3]={-1,-1,0},h[3]={1,1,1};int const cells[3]={8,8,0};
  Kernel k;k.m_axisymmetric=true;
  double const geometry[3][3]={{.1,0,0},{0,.1,0},{0,0,1}};
  int checked=0;
  for(double r:{1.e-20,1.e-14,1.e-8,.001,.02,.2,.8,.999999,1.,1.000001,2.3,5.99})
  for(double y:{0.,.01,.9,2.5,5.99,6.})
  {
    double const p[3]={r,y,0};int nodes[16];double w[16],grad[16][3];
    assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,nodes,w,grad));
    double sum=0,ds[2]={},m[2]={},moment[2][2]={},hoop=0,work=0;
    for(int a=0;a<16;++a)
    {
      double x[3];coordinates(nodes[a],origin,h,x);assert(w[a]>=-1.e-15);sum+=w[a];
      if(x[0]<0) {close(w[a],0,0);for(double v:grad[a])close(v,0,0);}
      for(int d=0;d<2;++d){m[d]+=x[d]*w[a];ds[d]+=grad[a][d];
        for(int e=0;e<2;++e)moment[d][e]+=x[d]*grad[a][e];}
      hoop+=x[0]*grad[a][2];work+=2*x[0]*(-7*grad[a][0]-11*grad[a][2])+3*x[1]*(-13*grad[a][1]);
    }
    close(sum,1);close(m[0],r);close(m[1],y);close(ds[0],0);close(ds[1],0);
    for(int d=0;d<2;++d)for(int e=0;e<2;++e)close(moment[d][e],d==e?1.:0.);
    close(hoop,1);close(work,-75);++checked;
  }
  // The BC acts on the interpolated radial field: only axis controls have
  // nonzero weight at r=0. Test the basis directly; particles require r>0.
  mpmMapping::CubicBoundary closed;closed.lower[0]=1;
  double p[3]={0,2.3,0};int nodes[16];double w[16],grad[16][3];
  assert((mpmMapping::point<2,true>(p,map,origin,h,cells,nodes,w,grad,closed)));
  for(int a=0;a<16;++a){double x[3];coordinates(nodes[a],origin,h,x);if(x[0]>0)close(w[a],0,0);}
  // Coordinate finite differences, including both sides of the transition at h.
  for(double r:{.05,.4,.99,1.,1.01,2.3})
  {
    p[0]=r;p[1]=2.3;k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,nodes,w,grad);
    std::map<int,std::array<double,3>> deriv;
    for(int a=0;a<16;++a)for(int d=0;d<3;++d)deriv[nodes[a]][d]=grad[a][d];
    for(int d=0;d<2;++d)
    {
      double const eps=1.e-5;std::map<int,double> diff;
      for(int sign:{-1,1}){p[d]+=sign*eps;int ni[16];double wi[16],gi[16][3];
        assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,ni,wi,gi));
        for(int a=0;a<16;++a)diff[ni[a]]+=sign*wi[a]/(2*eps);
        p[d]-=sign*eps;}
      for(auto const &entry:diff)close(entry.second,deriv[entry.first][d],1.e-9);
    }
  }
  // A genuinely exterior center is rejected rather than silently clamped.
  p[0]=-1.e-4;assert(!k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,nodes,w,grad));
  for(int a=0;a<16;++a){assert(nodes[a]==-1);close(w[a],0);}
  std::cout<<checked<<" axis B-spline rows: positivity, moments, hoop work, tiny-radius regularity and boundary/derivative checks passed.\n";
}
static void partitionTests()
{
  int map[9][9][9];grid(map);double const geometry[3][3]={{.1,0,0},{0,.1,0},{0,0,1}},h[3]={1,1,1};
  int const cells[3]={8,8,0};Kernel k;k.m_axisymmetric=true;k.m_xGlobalMax[0]=10;
  // Two overlapping rank stencils represent the same global controls. Neither
  // rank may accidentally close the internal patch edge.
  for(double r:{4.1,4.9,5.,5.1,5.9})
  {
    std::map<std::pair<int,int>,std::array<double,4>> rows[2];double const p[3]={r,2.3,0};
    for(int rank=0;rank<2;++rank)
    {
      double const origin[3]={rank==0?-1.:3.,-1,0};int n[16];double w[16],g[16][3];
      assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,n,w,g));
      for(int a=0;a<16;++a){double x[3];coordinates(n[a],origin,h,x);rows[rank][{int(x[0]),int(x[1])}]={w[a],g[a][0],g[a][1],g[a][2]};}
    }
    for(auto const &entry:rows[0])for(int d=0;d<4;++d)close(entry.second[d],rows[1][entry.first][d]);
  }
  // Cartesian and annular maps still use the exact original cardinal path.
  for(bool annular:{false,true})
  {
    k.m_axisymmetric=annular;k.m_xGlobalMin[0]=annular?2.:0.;
    double const origin[3]={annular?1.:-1.,-1,0},p[3]={annular?2.01:.01,2.3,0};
    int n[16],nr[16];double w[16],wr[16],g[16][3],gr[16][3];
    assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,geometry,map,origin,h,cells,n,w,g));
    assert((mpmMapping::point<2,true>(p,map,origin,h,cells,nr,wr,gr)));
    for(int a=0;a<16;++a){assert(n[a]==nr[a]);close(w[a],wr[a],0);close(g[a][0],gr[a][0],0);close(g[a][1],gr[a][1],0);}
  }
  std::cout<<"Partition-overlap mappings agree; Cartesian/F-table and annular cardinal weights are unchanged.\n";
}
static void pressureTests()
{
  int map[9][9][9];grid(map);double const origin[3]={-1,-1,0},h[3]={1,1,1},g0[3][3]={{.1,0,0},{0,.1,0},{0,0,1}};
  int const cells[3]={8,8,0};Kernel k;k.m_axisymmetric=true;
  double force[9][9]={};
  // Two-point Gauss quadrature integrates r*N'+N exactly in each cell.
  for(int i=0;i<6;++i)for(int j=0;j<6;++j)for(int sr:{-1,1})for(int sy:{-1,1})
  {
    double const p[3]={i+.5+sr/(2*std::sqrt(3.)),j+.5+sy/(2*std::sqrt(3.)),0};
    int n[16];double w[16],g[16][3];assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,p,g0,map,origin,h,cells,n,w,g));
    for(int a=0;a<16;++a)force[n[a]/81][(n[a]/9)%9]-=axisymmetric::volume(.25,p[0])*(g[a][0]+g[a][2]);
  }
  // Radial interior supports, including the axis and its first neighbor.
  for(int i=0;i<=5;++i)for(int j=2;j<=5;++j)close(force[i][j],0,2.e-11);
  double reaction=0;for(auto const &row:force)for(double value:row)reaction+=value;
  close(reaction,-axisymmetric::twoPi*6*6);
  std::cout<<"Hydrostatic cylinder: no spurious axis/interior radial force; correct full-ring resultant.\n";
}
// A scalar, full-ring transfer check using the actual generated B-spline rows.
// It checks the common transfer algebra, not the complete GEOS FMPM timestep.
static void transferTests()
{
  int map[9][9][9];grid(map);
  double const origin[3]={-1,-1,0},h[3]={1,1,1},geometry[3][3]={{.1,0,0},{0,.1,0},{0,0,1}};
  int const cells[3]={8,8,0};Kernel k;k.m_axisymmetric=true;
  constexpr int np=20*20,nn=729;
  std::vector<std::array<int,16>> nodes(np);
  std::vector<std::array<double,16>> weights(np);
  std::vector<double> mass(np),v(np),nodalMass(nn),constant(np,3.5);
  for(int i=0;i<20;++i)for(int j=0;j<20;++j)
  {
    int const p=i*20+j;double const position[3]={(i+.5)/4,(j+.5)/4,0};double gradients[16][3];
    assert(k.computeParticleMapping<2>(ParticleType::SinglePointBSpline,position,geometry,map,origin,h,cells,
      nodes[p].data(),weights[p].data(),gradients));
    mass[p]=axisymmetric::volume(1./16,position[0]);v[p]=std::sin(position[0])*.7+std::cos(position[1]*1.3);
    for(int a=0;a<16;++a)nodalMass[nodes[p][a]]+=mass[p]*weights[p][a];
  }
  auto project=[&](std::vector<double> const &input) {
    std::vector<double> gridValue(nn),out(np);
    for(int p=0;p<np;++p)for(int a=0;a<16;++a)gridValue[nodes[p][a]]+=mass[p]*input[p]*weights[p][a];
    for(int n=0;n<nn;++n)if(nodalMass[n]>0)gridValue[n]/=nodalMass[n];
    for(int p=0;p<np;++p)for(int a=0;a<16;++a)out[p]+=weights[p][a]*gridValue[nodes[p][a]];
    return out;
  };
  auto pic=project(v),twice=project(pic),translated=project(constant);
  double momentum=0,picMomentum=0,fmpmMomentum=0,energy=0,picEnergy=0,fmpmEnergy=0;
  for(int p=0;p<np;++p)
  {
    close(translated[p],3.5);
    double const fmpm2=2*pic[p]-twice[p];
    momentum+=mass[p]*v[p];picMomentum+=mass[p]*pic[p];fmpmMomentum+=mass[p]*fmpm2;
    energy+=mass[p]*v[p]*v[p];picEnergy+=mass[p]*pic[p]*pic[p];fmpmEnergy+=mass[p]*fmpm2*fmpm2;
  }
  close(momentum,picMomentum);close(momentum,fmpmMomentum);
  assert(picEnergy<=fmpmEnergy && fmpmEnergy<=energy);
  std::cout<<"400 full-ring particles: axial translation/momentum preserved; PIC and FMPM(2) transfer algebra is non-amplifying.\n";
}
static void eventTests()
{
  FrictionHarness k;FrictionCoefficientSwapMPMEvent e;k.trigger(e);
  close(k.m_frictionCoefficient,.25);assert(k.m_frictionCoefficientTable==e.table);assert(e.complete&&k.initialized==1);
  e.coefficient=.9;k.trigger(e);close(k.m_frictionCoefficient,.25);assert(k.initialized==1);
  e.complete=false;e.table.clear();k.trigger(e);close(k.m_frictionCoefficient,.9);assert(k.m_frictionCoefficientTable.empty());
  bool rejected=false;e.catalog="TransformParticles";try{k.trigger(e);}catch(std::runtime_error const &){rejected=true;}assert(rejected);
  std::cout<<"Actual event guard/swap handler: scalar/table updates accepted, completion retained, unsupported event rejected.\n";
}
int main(){surfaceTests();cleanupTests();basisTests();partitionTests();pressureTests();transferTests();eventTests();}
