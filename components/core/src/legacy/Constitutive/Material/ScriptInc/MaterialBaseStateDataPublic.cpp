//FUNCTION_BEGIN_PARSE
void
MaterialBaseStateData::TotalStress(R2SymTensor& totalStress) const
{
  totalStress = devStress;
  totalStress.PlusIdentity(pressure);
}

//FUNCTION_BEGIN_PARSE
void
MaterialBaseStateData::RotateState( const R2Tensor& Rot )
{
  R2SymTensor temp;

  devStress.PlusIdentity( pressure );
  temp.QijAjkQlk(devStress,Rot);
  pressure = temp.Trace() / 3.0;
  devStress = temp;
  devStress.PlusIdentity(-pressure);
}
