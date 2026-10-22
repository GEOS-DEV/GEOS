
static bool near(double a,double b) {return std::abs(a-b)<2e-10*std::max({1.,std::abs(a),std::abs(b)});}
int main()
{
  double xyz[729][3]; int map[9][9][9];
  double origin[3]={-1,-1,0},h[3]={1,1,1};int cells[3]={8,8,0};
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k)
  {int n=(i*9+j)*9+k;map[i][j][k]=n;xyz[n][0]=i-1;xyz[n][1]=j-1;xyz[n][2]=0;}
  std::mt19937 rng(1486);std::uniform_real_distribution<double> jitter(-.15,.15);
  int checked=0;
  for(auto type : {ParticleType::SinglePoint,ParticleType::CPDI,ParticleType::CPTI})
  for(int mode=0;mode<3;++mode) for(int trial=0;trial<200;++trial)
  {
    Kernel kernel;kernel.m_axisymmetric=true;kernel.m_cpdiDomainScaling=mode>0;
    kernel.m_cpdiDomainScalingType=mode==1?mpm::CPDIDomainScalingTypeOption::Homel:mpm::CPDIDomainScalingTypeOption::VPHencky;
    double point[3]={trial%2 ? 2.3 : .02,2.3,0};
    double g[3][3]={{.8+jitter(rng),jitter(rng),0},{jitter(rng),.7+jitter(rng),0},{0,0,1}};
    double copy[3][3];LvArray::tensorOps::copy<3,3>(copy,g);
    int nodes[16];double weights[16],grad[16][3];
    assert(kernel.computeParticleMapping<2>(type,point,g,map,origin,h,cells,nodes,weights,grad));
    int count=type==ParticleType::SinglePoint?4:(type==ParticleType::CPTI?12:16);
    double sum=0,ds[3]={},moment[2][3]={},virtualWork=0,hoopWork=0;
    for(int a=0;a<count;++a)
    {
      sum+=weights[a];assert(weights[a]>=-1e-14);
      double x=xyz[nodes[a]][0],y=xyz[nodes[a]][1];
      // No physical support at negative r, even when an unscaled domain crosses the axis.
      if(x<0) assert(std::abs(weights[a])+std::abs(grad[a][0])+std::abs(grad[a][1])+std::abs(grad[a][2])<1e-12);
      for(int d=0;d<3;++d){ds[d]+=grad[a][d];moment[0][d]+=x*grad[a][d];moment[1][d]+=y*grad[a][d];}
      // v_r=2r, v_z=3z: stress power = 2*s_rr+3*s_zz+2*s_tt.
      double fr=-(7*grad[a][0]+11*grad[a][2]);double fz=-13*grad[a][1];
      virtualWork+=2*x*fr+3*y*fz;
      hoopWork+=2*x*grad[a][2];
    }
    assert(near(sum,1));assert(near(ds[0],0));assert(near(ds[1],0));assert(near(ds[2],1/point[0]));
    assert(near(moment[0][0],1)&&near(moment[1][1],1));
    assert(near(moment[0][1],0)&&near(moment[1][0],0));
    assert(near(moment[0][2],1));assert(near(hoopWork,2));assert(near(virtualWork,-75));
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) assert(near(copy[i][j],g[i][j]));
    ++checked;
  }
  // Uniform hydrostatic stress in a tiled cylinder: all interior radial
  // nodal forces cancel, including axis nodes. Only the exterior traction remains.
  double radialForce[9][9]={};
  Kernel equilibrium;equilibrium.m_axisymmetric=true;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j)
  {
    double center[3]={i+.5,j+.5,0},edges[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    int nodes[16];double w[16],grad[16][3];
    assert(equilibrium.computeParticleMapping<2>(ParticleType::CPDI,center,edges,map,origin,h,cells,nodes,w,grad));
    for(int a=0;a<16;++a)
    {
      int const ni=nodes[a]/81,nj=(nodes[a]/9)%9;
      radialForce[ni][nj]-=axisymmetric::volume(1,center[0])*(grad[a][0]+grad[a][2]);
    }
  }
  for(int i=1;i<5;++i) for(int j=1;j<=5;++j) assert(near(radialForce[i][j],0));
  double outerForce=0;for(int j=0;j<9;++j) outerForce+=radialForce[5][j];
  assert(near(outerForce,-axisymmetric::twoPi*4*4));
  // Compare radial box moments to an independently integrated linear field.

  // N=x*y on a single background cell is reproduced exactly by Q1.
  Kernel kernel;kernel.m_axisymmetric=true;
  double point[3]={2.5,2.5,0},g[3][3]={{.4,.1,0},{.15,.3,0},{0,0,1}};
  int nodes[16];double w[16],grad[16][3];
  assert(kernel.computeParticleMapping<2>(ParticleType::CPDI,point,g,map,origin,h,cells,nodes,w,grad));
  double actual=0;for(int a=0;a<16;++a) actual+=w[a]*xyz[nodes[a]][0];
  double expected=point[0]+(g[0][0]*g[0][0]+g[1][0]*g[1][0])/(12*point[0]);assert(near(actual,expected));
  // Real box/triangle subdivision plans conserve full-ring mass and affine momentum.
  for(auto shape : {particleGeometry::Shape::Parallelogram,particleGeometry::Shape::Triangle})
  {
    auto parent=particleGeometry::affine::Domain::fromGenerators(shape,g);
    auto plan=shape==particleGeometry::Shape::Parallelogram
      ? particleGeometry::subdivision::boxPlan(3,2)
      : particleGeometry::subdivision::sizePlan(parent,2,.2);
    assert(plan.count>1);
    double sum=0,mean[3]={};
    for(int c=0;c<plan.count;++c)
    { double dx[3];plan.children[c].offset(parent,dx);double f=axisymmetric::splitFraction(plan.children[c].fraction,point[0],dx[0]);sum+=f;for(int d=0;d<3;++d) mean[d]+=f*dx[d]; }
    assert(near(sum,1));double momentum=0;
    for(int c=0;c<plan.count;++c)
    { double dx[3];plan.children[c].offset(parent,dx);double f=axisymmetric::splitFraction(plan.children[c].fraction,point[0],dx[0]);momentum+=f*(4+2*(dx[0]-mean[0])+3*(dx[1]-mean[1])); }
    assert(near(momentum,4));
  }
  assert(near(axisymmetric::boxVolume(0,2,3),12*3.14159265358979323846));
  assert(near(axisymmetric::volume(2,3),12*3.14159265358979323846));
  point[0]=0;assert(!kernel.computeParticleMapping<2>(ParticleType::CPDI,point,g,map,origin,h,cells,nodes,w,grad));
  for(int a=0;a<16;++a){assert(nodes[a]==-1);assert(near(w[a],0));}
  std::cout<<checked<<" production axisymmetric mappings: moments, hoop strain, stress work, axis restriction, geometry immutability and conservative splits passed.\n";
}
