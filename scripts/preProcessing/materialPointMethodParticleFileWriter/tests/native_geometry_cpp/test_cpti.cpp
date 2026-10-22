// Production CPTI mapping checked against independent Gaussian elimination and
// a separately accumulated trilinear stencil, with analytic moment identities.
template<class T> struct View { T *data; T &operator[](int i) const {return data[i];} };
struct GradientView { double (*data)[3]; double *operator[](int i) const {return data[i];} };
static bool close(double a,double b,double tol=3e-10)
{return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});}

static void inverse(double const (&g)[3][3],double (&inv)[3][3])
{
  double a[3][6]={};
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) {a[i][j]=g[i][j];a[i][j+3]=i==j?1.:0.;}
  for(int j=0;j<3;++j)
  {
    int pivot=j;
    for(int i=j+1;i<3;++i) if(std::abs(a[i][j])>std::abs(a[pivot][j])) pivot=i;
    for(int k=0;k<6;++k) std::swap(a[j][k],a[pivot][k]);
    double d=a[j][j];assert(std::abs(d)>0);
    for(int k=0;k<6;++k) a[j][k]/=d;
    for(int i=0;i<3;++i) if(i!=j) {d=a[i][j];for(int k=0;k<6;++k) a[i][k]-=d*a[j][k];}
  }
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) inv[i][j]=a[i][j+3];
}

static int checked=0;
static void check(Kernel const &kernel,double const (&point)[3],double const (&g)[3][3])
{
  int map[9][9][9]; double grid[729][3],h[3]={.8,1.1,1.3},xmin[3]={-.2,.1,-.4}; int nel[3]={8,8,8};
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k)
  {int n=(i*9+j)*9+k;map[i][j][k]=n;grid[n][0]=xmin[0]+i*h[0];grid[n][1]=xmin[1]+j*h[1];grid[n][2]=xmin[2]+k*h[2];}
  int nodes[34],nodes2[34];double weights[34],weights2[34],grad[34][3],grad2[34][3];
  std::fill_n(nodes,34,-123); std::fill_n(nodes2,34,-123);
  assert(kernel.computeCPTIShapeFunctions(point,g,map,xmin,h,nel,nodes+1,weights+1,grad+1));
  assert(kernel.computeCPTIShapeFunctions(point,g,map,xmin,h,nel,View<int>{nodes2+1},View<double>{weights2+1},GradientView{grad2+1}));
  assert(nodes[0]==-123 && nodes[33]==-123 && nodes2[0]==-123 && nodes2[33]==-123);
  double projected[3][3],inv[3][3],alpha[4][3]={};
  assert(!(kernel.projectCPTIIntegrationGeometry(g,projected)&4));inverse(projected,inv);
  for(int d=0;d<3;++d) for(int a=0;a<3;++a) {alpha[a+1][d]=inv[d][a];alpha[0][d]-=inv[d][a];}
  double reference[729]={},referenceGrad[729][3]={},actual[729]={},actualGrad[729][3]={};
  for(int a=0;a<4;++a)
  {
    double vertex[3];int cell[3];double xi[3];
    for(int d=0;d<3;++d)
    {
      vertex[d]=point[d];for(int e=0;e<3;++e) vertex[d]+=((a==e+1?1.:0.)-.25)*projected[e][d];
      double u=(vertex[d]-xmin[d])/h[d];cell[d]=static_cast<int>(std::floor(u));xi[d]=u-cell[d];
    }
    for(int i=0;i<2;++i) for(int j=0;j<2;++j) for(int k=0;k<2;++k)
    {
      int n=map[cell[0]+i][cell[1]+j][cell[2]+k];
      double w=(i?xi[0]:1-xi[0])*(j?xi[1]:1-xi[1])*(k?xi[2]:1-xi[2]);
      reference[n]+=.25*w;for(int d=0;d<3;++d) referenceGrad[n][d]+=w*alpha[a][d];
    }
  }
  double sum=0,position[3]={},gradient[3]={},moment[3][3]={};
  for(int n=1;n<=32;++n)
  {
    assert(nodes[n]==nodes2[n]);assert(close(weights[n],weights2[n]));assert(weights[n]>=-1e-15);
    int id=nodes[n];actual[id]+=weights[n];sum+=weights[n];
    for(int d=0;d<3;++d)
    {
      assert(close(grad[n][d],grad2[n][d]));actualGrad[id][d]+=grad[n][d];
      position[d]+=weights[n]*grid[id][d];gradient[d]+=grad[n][d];
      for(int j=0;j<3;++j) moment[d][j]+=grid[id][d]*grad[n][j];
    }
  }
  assert(close(sum,1));
  for(int d=0;d<3;++d)
  {assert(close(position[d],point[d]));assert(close(gradient[d],0));for(int j=0;j<3;++j) assert(close(moment[d][j],d==j?1:0));}
  for(int n=0;n<729;++n)
  {assert(close(actual[n],reference[n]));for(int d=0;d<3;++d) assert(close(actualGrad[n][d],referenceGrad[n][d]));}
  ++checked;
}

int main()
{
  std::mt19937 rng(2131);std::uniform_real_distribution<double> jitter(-.2,.2);
  for(int mode=0;mode<3;++mode) for(int trial=0;trial<1000;++trial)
  {
    Kernel kernel;kernel.m_cpdiDomainScaling=mode!=0;
    kernel.m_cpdiDomainScalingType=mode==1?mpm::CPDIDomainScalingTypeOption::Homel:mpm::CPDIDomainScalingTypeOption::VPHencky;
    double g[3][3],original[3][3],point[3]={3.13+jitter(rng),4.21+jitter(rng),4.39+jitter(rng)};
    for(int a=0;a<3;++a) for(int d=0;d<3;++d) g[a][d]=(a==d?.4+(trial%9)*.12:0)+jitter(rng);
    // Proper rotations cover non-axis-aligned and negative edge components.
    double theta=trial*.1,c=std::cos(theta),s=std::sin(theta);
    for(int a=0;a<3;++a) {double x=g[a][0],y=g[a][1];g[a][0]=c*x-s*y;g[a][1]=s*x+c*y;}
    LvArray::tensorOps::copy<3,3>(original,g);check(kernel,point,g);
    double projected[3][3];int status=kernel.projectCPTIIntegrationGeometry(g,projected);
    assert(!(status&4));
    auto before=particleGeometry::affine::Domain::fromGenerators(particleGeometry::Shape::Tetrahedron,g);
    auto after=particleGeometry::affine::Domain::fromGenerators(particleGeometry::Shape::Tetrahedron,projected);
    double limit=mode==1?.99998:.99999*kernel.m_neighborRadius;
    double scale=mode==0?1.:std::min(1.,limit/before.diameter());
    assert(close(after.measure(),before.measure()*scale*scale*scale));
    for(int a=0;a<3;++a) for(int d=0;d<3;++d) {assert(close(projected[a][d],scale*g[a][d]));assert(close(g[a][d],original[a][d],0));}
    if(trial%20==0)
    {
      auto plan=particleGeometry::subdivision::sizePlan(before,3,.65);
      double weighted[3]={},fraction=0;
      for(int child=0;child<plan.count;++child)
      {
        auto domain=plan.children[child].geometry(before);double offset[3],center[3];
        plan.children[child].offset(before,offset);
        for(int d=0;d<3;++d) {center[d]=point[d]+offset[d];weighted[d]+=plan.children[child].fraction*center[d];}
        fraction+=plan.children[child].fraction;check(kernel,center,domain.generators);
      }
      assert(close(fraction,1));for(int d=0;d<3;++d) assert(close(weighted[d],point[d]));
      double newF[3][3]={{1.1,.1,0},{0,.9,.1},{0,0,1}},current[3][3];
      LvArray::tensorOps::copy<3,3>(current,g);particleGeometry::affine::reset(original,current,newF,false);
      check(kernel,point,current);
      for(int a=0;a<3;++a) for(int d=0;d<3;++d) assert(close(current[a][d],g[a][d],0));
      particleGeometry::affine::reset(original,current,newF,true);check(kernel,point,current);
    }
  }
  // Reject bad geometry/positions before any map access; every output is cleared.
  Kernel kernel;int map[9][9][9]={},nodes[32],nel[3]={8,8,8};double h[3]={1,1,1},xmin[3]={};
  double weights[32],grad[32][3];
  for(int bad=0;bad<7;++bad)
  {
    double g[3][3]={{1,0,0},{0,1,0},{0,0,1}},point[3]={4,4,4};
    if(bad==0) g[0][0]=0;
    if(bad==1) g[0][0]=-1;
    if(bad==2) g[1][1]=std::numeric_limits<double>::quiet_NaN();
    if(bad==3) point[0]=std::numeric_limits<double>::infinity();
    if(bad==4) point[1]=1e100;
    if(bad==5) point[2]=-.01;
    if(bad==6) point[0]=8.;
    std::fill_n(nodes,32,123);
    assert(!kernel.computeCPTIShapeFunctions(point,g,map,xmin,h,nel,nodes,weights,grad));
    for(int n=0;n<32;++n) {assert(nodes[n]==-1);assert(close(weights[n],0,0));for(int d=0;d<3;++d) assert(close(grad[n][d],0,0));}
  }
  // Exactly on either global boundary, including one-ulp outward roundoff.
  for(bool upper : {false,true})
  {
    double g[3][3]={{1,0,0},{0,1,0},{0,0,1}},point[3]={.25,.25,.25};
    if(upper) point[0]=7.25;
    point[0]=std::nextafter(point[0],upper?9.:-1.);
    assert(kernel.computeCPTIShapeFunctions(point,g,map,xmin,h,nel,nodes,weights,grad));
  }
  // Geometrically similar particles at very different length scales.
  for(double scale : {1e-90,1e90})
  {
    double g[3][3]={{scale,0,0},{0,scale,0},{0,0,scale}},gradient[4][3];
    assert(particleGeometry::affine::tetrahedronGradients(g,gradient));
    assert(close(gradient[1][0]*scale,1));
  }
  std::cout<<"CPTI: "<<checked<<" mapping cases passed; reference stencils, moments, views/raw arrays, projection, splitting, resets, invalid inputs and boundaries.\n";
}
