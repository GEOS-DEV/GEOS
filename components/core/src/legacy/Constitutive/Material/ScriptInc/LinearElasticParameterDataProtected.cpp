//FUNCTION_BEGIN_PARSE
virtual_void
LinearElasticParameterData::PostReadXML( const TICPP::HierarchicalDataNode& node )
{
  realT M = 0;
  FillLinearElasticModuli(init_bulkModulus, init_shearModulus, E, Nu, Lame, M);
}
