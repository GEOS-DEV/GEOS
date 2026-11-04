// Repository/array adapters surrounding the extracted production force lambda.
struct Particle { double x[3]={},geometry[3][3]={},volume=0; int group=0; ParticleType type=ParticleType::CPDI; };
int m_ijkMap[65][65][1]={},m_nEl[3]={};
double m_xLocalMin[3]={};
double gridPosition[65*65][3]={},gridMaterialVolume[65*65][3]={},gridSurfaceTensionForce[65*65][2][3]={};
double boundary[65*65][2]={};
std::vector<Particle> particles;
SurfaceKernel(int n)
{
  m_axisymmetric=true;m_hEl[0]=m_hEl[1]=1./n;m_nEl[0]=m_nEl[1]=n+6;
  m_xLocalMin[0]=m_xLocalMin[1]=-3./n;m_xGlobalMax[0]=m_xGlobalMax[1]=1;
  for(int d=0;d<4;++d)m_boundaryConditionTypes[d]=1;
  for(int i=0;i<=n+6;++i)for(int j=0;j<=n+6;++j)
  {int g=i*65+j;m_ijkMap[i][j][0]=g;gridPosition[g][0]=(i-3.)/n;gridPosition[g][1]=(j-3.)/n;}
}
void addRectangle(int n,double radius,bool split,ParticleType type)
{
  double dx=1./(4*n);
  for(int i=0;i<4*n;++i)for(int j=0;j<4*n;++j)
  {
    Particle p;p.x[0]=(i+.5)*dx;p.x[1]=(j+.5)*dx;
    if(p.x[0]>radius)continue;
    p.geometry[0][0]=p.geometry[1][1]=dx;
    p.volume=axisymmetricSurface::twoPi*p.x[0]*dx*dx;p.type=type;
    p.group=split&&p.x[0]>.5?1:0;
    if(type==ParticleType::CPTI)
    {
      for(int triangle=0;triangle<2;++triangle)
      {
        Particle t=p;
        if(triangle==0){t.x[0]-=dx/6;t.x[1]-=dx/6;}
        else{t.x[0]+=dx/6;t.x[1]+=dx/6;t.geometry[0][0]=-dx;t.geometry[1][1]=-dx;}
        t.volume=axisymmetricSurface::twoPi*t.x[0]*dx*dx/2;particles.push_back(t);
      }
    }
    else particles.push_back(p);
  }
}
void project()
{
  for(auto const &p:particles)
  {
    bool valid=axisymmetricSurface::projectGeometry(p.type==ParticleType::CPTI,p.x,p.geometry,m_xLocalMin,m_hEl,m_nEl,
      [&](int i,int j,double weight,double br,double by)
      {int g=m_ijkMap[i][j][0];gridMaterialVolume[g][p.group]+=p.volume*weight;boundary[g][0]+=p.volume*br;boundary[g][1]+=p.volume*by;});
    assert(valid);
  }
}
double work(bool axial)
{
  double w=0;
  for(int i=0;i<=m_nEl[0];++i)for(int j=0;j<=m_nEl[1];++j)
  {int g=m_ijkMap[i][j][0];for(int f=0;f<2;++f)w+=gridSurfaceTensionForce[g][f][0]*gridPosition[g][0]+(axial?gridSurfaceTensionForce[g][f][1]*gridPosition[g][1]:0);}
  return w;
}
