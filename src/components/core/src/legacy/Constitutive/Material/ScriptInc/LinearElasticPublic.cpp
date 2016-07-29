//FUNCTION_BEGIN_PARSE
void LinearElastic::InitializeStates(const localIndex index)
{

  for (localIndex a = 0; a < m_stateData.Dimension(1); ++a)
  {
    const localIndex paramIndex = m_parameterData.size() > 1 ? index : 0;

    m_stateData[index][a].density = m_parameterData[paramIndex].init_density;
    m_stateData[index][a].BulkModulus = m_parameterData[paramIndex].init_bulkModulus;
    m_stateData[index][a].ShearModulus = m_parameterData[paramIndex].init_shearModulus;
    m_stateData[index][a].ElasticBulkModulus = m_parameterData[paramIndex].init_bulkModulus;
    m_stateData[index][a].ElasticShearModulus = m_parameterData[paramIndex].init_shearModulus;
  }

}

//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::StrainDrivenUpdateMember(const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT<3>& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot, const realT dt)
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0;
  const LinearElasticParameterData& matParams = m_parameterData[paramIndex];
  LinearElasticStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;
  const realT& K = matParams.init_bulkModulus;
  const realT& G = matParams.init_shearModulus;

  const realT trDdt = Ddt.Trace();

  R2SymTensorT<3> temp;

  pressure += trDdt * K;

  temp = Ddt;
  temp.PlusIdentity(-trDdt / 3.0);
  temp *= 2.0 * G;
  devStress += temp;

  matState.RotateState(Rot);
  return;
}

//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::StrainDrivenUpdateMember(const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT<3>& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot, const realT& volume_n,
                                        const realT& volume_np1, const realT dt)
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0;
  const LinearElasticParameterData& matParams = m_parameterData[paramIndex];
  LinearElasticStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;
  const realT& K = matParams.init_bulkModulus;
  const realT& G = matParams.init_shearModulus;

  const realT trDdt = Ddt.Trace();

  realT StressPowerIncrement = 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_n;
  //  realT strainEnergyIncrement = - 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_n;

  pressure += trDdt * K;

  {
    R2SymTensorT<3> temp;
    temp = Ddt;
    temp.PlusIdentity(-trDdt / 3.0);
    temp *= 2.0 * G;
    devStress += temp;
  }

  StressPowerIncrement += 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_np1;
  //  strainEnergyIncrement += 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_np1;

  matState.ElasticStrainEnergy += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1);
  matState.StressPower += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1 );

  matState.RotateState(Rot);
  return;
}

//FUNCTION_BEGIN_PARSE
void LinearElastic::MeanPressureDevStress(const localIndex index, realT& pressure,
                                          R2SymTensor& devStress) const
{
  MeanPressureDevStressFromDerived<LinearElastic>(index, pressure, devStress);
}
