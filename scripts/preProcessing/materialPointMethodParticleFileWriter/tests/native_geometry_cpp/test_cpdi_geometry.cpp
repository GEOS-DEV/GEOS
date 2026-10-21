static bool close(double a,double b) {return std::abs(a-b)<2e-10*std::max({1.,std::abs(a),std::abs(b)});}
template<int DIM> int check()
{
  double grid[729][3]; int map[9][9][9];
  for(int i=0;i<9;++i) for(int j=0;j<9;++j) for(int k=0;k<9;++k)
  {int n=(i*9+j)*9+k;map[i][j][k]=n;grid[n][0]=i;grid[n][1]=j;grid[n][2]=k;}
  double xmin[3]={},h[3]={1,1,1};int nel[3]={8,8,8};
  std::mt19937 random(859);std::uniform_real_distribution<double> jitter(-.15,.15);
  int checked=0;
  constexpr int dim=DIM;
  if constexpr(DIM==2) nel[2]=0;
  for(int mode : {0,1,2}) for(int trial=0;trial<500;++trial)
  {
    Kernel kernel;
    kernel.m_cpdiDomainScaling=mode!=0;
    kernel.m_cpdiDomainScalingType=mode==1 ? mpm::CPDIDomainScalingTypeOption::Homel : mpm::CPDIDomainScalingTypeOption::VPHencky;
    double g[3][3]={},original[3][3],projected[3][3];
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) g[i][j]=(i==j ? .2+(trial%11)*.19 : 0)+jitter(random);
    if(dim==2) {g[0][2]=g[1][2]=g[2][0]=g[2][1]=0;g[2][2]=.2;}
    // A material domain must have positive orientation.
    if(particleGeometry::affine::determinant(g)<0) for(int d=0;d<3;++d) std::swap(g[0][d],g[1][d]);
    LvArray::tensorOps::copy<3,3>(original,g);
    int status=kernel.projectIntegrationGeometry<DIM>(g,projected);assert(!(status&4));
    double inputMeasure=particleDomainMeasure(g,dim),outputMeasure=particleDomainMeasure(projected,dim);
    assert(outputMeasure>0);
    if(mode==2)
    {
      double cap=particleDomainMaximumFeasibleMeasure(dim,.99999*kernel.m_neighborRadius);
      assert(close(outputMeasure,std::min(inputMeasure,cap)));
      assert(particleDomainMaximumDiagonal(projected,dim)<=.99999*kernel.m_neighborRadius*(1+1e-9));
    }
    if(mode==0) for(int i=0;i<DIM;++i) for(int j=0;j<DIM;++j) assert(close(g[i][j],projected[i][j]));
    double point[3]={4.13+jitter(random),4.26+jitter(random),4.39+jitter(random)};
    if constexpr(DIM==2) point[2]=0;
    constexpr int count=1<<(2*DIM);
    int nodes[count],nodes2[count];double weights[count],weights2[count],grad[count][3],grad2[count][3];
    if(!kernel.computeParticleMapping<DIM>(ParticleType::CPDI,point,g,map,xmin,h,nel,nodes,weights,grad)) {
      std::cerr<<"mapping failed DIM="<<DIM<<" mode="<<mode<<" trial="<<trial<<" measure="<<inputMeasure<<"\n";
      for(auto & row:g) {for(auto v:row) std::cerr<<v<<" "; std::cerr<<"\n";}
      std::abort();
    }
    assert(kernel.computeParticleMapping<DIM>(ParticleType::CPDI,point,g,map,xmin,h,nel,nodes2,weights2,grad2));
    double sum=0,position[3]={},gradient[3]={},moment[3][3]={};
    for(int n=0;n<count;++n)
    {
      assert(nodes[n]==nodes2[n]);assert(close(weights[n],weights2[n]));sum+=weights[n];
      for(int i=0;i<3;++i)
      {
        assert(close(grad[n][i],grad2[n][i]));position[i]+=weights[n]*grid[nodes[n]][i];gradient[i]+=grad[n][i];
        for(int j=0;j<3;++j) moment[i][j]+=grid[nodes[n]][i]*grad[n][j];
      }
    }
    assert(close(sum,1));
    for(int i=0;i<3;++i) {assert(close(position[i],point[i]));assert(close(gradient[i],0));
      for(int j=0;j<3;++j) {assert(close(moment[i][j],i==j && i<DIM?1:0));assert(close(g[i][j],original[i][j]));}}
    // The projection and mapping inputs remain valid after repeated calls.
    ++checked;
  }
  // Invalid projection is reported and preserves the original geometry.
  Kernel kernel;kernel.m_cpdiDomainScaling=1;double bad[3][3]={},out[3][3];
  assert(kernel.projectIntegrationGeometry<DIM>(bad,out)==4);
  // Infeasible reset projection preserves requested measure, unlike mapping.
  double big[3][3]={{3,0,0},{0,2,0},{0,0,1}};
  auto result=projectVPHenckyParticleDomain(big,3,1,6,true,out);
  assert(result.valid && !result.volumeFeasible && !result.volumeRelaxed);
  assert(close(particleDomainMeasure(out,3),6));
  return checked;
}
int main()
{
  int checked=check<2>()+check<3>();
  std::cout<<"Production CPDI: "<<checked<<" native 2D/3D unscaled/Homel/VP-Hencky cases; projection measure, interpolation moments and input immutability passed.\n";
}
