//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::PreSetValues(const sArray1d& names)
{
  //reset mechanical constants
  for (localIndex a = 0; a < m_parameterData.size(); ++a)
  {
    m_parameterData[a].E = 0;
    m_parameterData[a].Lame = 0;
    m_parameterData[a].Nu = 0;
    m_parameterData[a].init_bulkModulus = 0;
    m_parameterData[a].init_shearModulus = 0;
  }
}

//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::PostSetValues(const sArray1d& names)
{
  //recalculate mechanical parameters that were not already explicitly set
  for(localIndex a = 0; a < m_parameterData.size(); a++)
  {
    realT M = 0;
    FillLinearElasticModuli(m_parameterData[a].init_bulkModulus,
                            m_parameterData[a].init_shearModulus,
                            m_parameterData[a].E,
                            m_parameterData[a].Nu,
                            m_parameterData[a].Lame,
                            M);
  }
}
