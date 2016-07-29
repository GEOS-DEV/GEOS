//FUNCTION_BEGIN_PARSE
virtual_void
InterfaceBase::Initialize( const localIndex index,
                           const realT stressNormal,
                           const realT stressShear)
{
  throw GPException("Cannot call InterfaceBase::Initialize\n");
}

//FUNCTION_BEGIN_PARSE
virtual_void
InterfaceBase::UpdateProperties(const localIndex , std::map<std::string, realT>&, std::map<std::string, realT>& )
{
}

//FUNCTION_BEGIN_PARSE
virtual_realT
InterfaceBase::SetPermeabilityTerm(const localIndex index)
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  if(matState.normalApproach <= 0)
  {
    matState.kappa = matParams.kappa0 - (matState.normalApproach * matState.normalApproach * matState.normalApproach);
  }
  else
  {
    matState.kappa = matParams.kappa0 +
        matParams.dkappadnFct * (isEqual(matParams.dkappadnExp, 1.0) ? matState.normalApproach : pow(matState.normalApproach, matParams.dkappadnExp)) +
        matParams.dkappadsFct * (isEqual(matParams.dkappadsExp, 1.0) ? matState.xs : pow(matState.xs, matParams.dkappadsExp));
  }
  return matState.kappa;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
InterfaceBase::StiffnessProjected( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  const realT normalApproachEstimate = matState.normalApproach +
      (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);

  return NormalStiffness(matParams, matState, normalApproachEstimate, false);
}

//FUNCTION_BEGIN_PARSE
virtual_void
InterfaceBase::StrainDrivenUpdate( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  //evolve the normal force ("stress" ... bad terminology, but I like the inheritance)
  const realT stiffness = NormalStiffness(matParams, matState, matState.normalApproach, true);

  //evolve the friction
  UpdateFriction(matParams, matState);

  //evolve the shear stress
  {
    matState.stressShear += -stiffness * matParams.ston * matState.dxsdt * matState.dt;
    matState.xs += matState.dxsdt * matState.dt;
  }

  //determine whether strength exceeded and return to failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  //update vector representation of shear stress
  //NOTE: JointBaseStateData::UpdateOrientation sets the direction as -(shear slip direction)
  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}
