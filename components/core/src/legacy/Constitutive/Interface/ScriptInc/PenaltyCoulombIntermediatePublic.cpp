//FUNCTION_BEGIN_PARSE
virtual_void
PenaltyCoulombIntermediate::StrainDrivenUpdate( const localIndex index )
{
  //get temporary references
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);

  //(1) evolve normal stress
  matState.stress = matParams.Stress(matState.normalApproach);

  //(2) evolve friction
  UpdateFriction(matParams, matState);

  //(3) evolve shear stress
  const realT dssdt = matParams.kshear * matState.dxsdt;
  matState.stressShear += dssdt * matState.dt;

  //(4) return shear stress to the failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
PenaltyCoulombIntermediate::StiffnessProjected(const localIndex index)
{
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  const PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);
  const realT normalApproachEstimate = matState.normalApproach +
      (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);
  return matParams.Stiffness(normalApproachEstimate);
}
