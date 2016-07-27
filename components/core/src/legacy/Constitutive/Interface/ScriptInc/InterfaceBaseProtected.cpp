//FUNCTION_BEGIN_PARSE
virtual_realT InterfaceBase::ShearStrength(const InterfaceBaseParameterData& matParams,
                                               InterfaceBaseStateData& matState) const
{
  return matState.stress * matParams.mu0;
}

//FUNCTION_BEGIN_PARSE
virtual_void InterfaceBase::UpdateFriction(const InterfaceBaseParameterData&,
                                               InterfaceBaseStateData&) const
{
}

//FUNCTION_BEGIN_PARSE
virtual_realT
InterfaceBase::NormalStiffness(const InterfaceBaseParameterData&,
                                   InterfaceBaseStateData&,
                                   const realT ,
                                   const bool ) const
{
  return 0.0;
}

//FUNCTION_BEGIN_PARSE
virtual_void
InterfaceBase::ThresholdToFailureSurface(const InterfaceBaseParameterData& matParams,
                                         InterfaceBaseStateData& matState,
                                         const realT dxs)
{
  //determine whether strength exceeded and return to failure surface if necessary
  const realT strength = ShearStrength( matParams, matState);
  if(strength < fabs(matState.stressShear))
  {
    matState.DissipatedEnergy += fabs(dxs) * strength;
    matState.stressShear = matState.stressShear < 0 ? -strength : strength;
  }
}
