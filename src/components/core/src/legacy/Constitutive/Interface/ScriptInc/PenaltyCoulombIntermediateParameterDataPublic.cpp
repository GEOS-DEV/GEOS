//FUNCTION_BEGIN_PARSE
realT
PenaltyCoulombIntermediateParameterData::Stiffness(const realT normalApproach) const
{
  realT arealStiffness = 0.0;

  if(normalApproach <= 0)
    arealStiffness = 0.0;
  else if(normalApproach < normalApproachYield)
    arealStiffness = ktildeAperture / (aperture - normalApproach);
  else if(normalApproach < normalApproachSoften)
    arealStiffness = kyield;
  else
    arealStiffness = ksoften;

  return arealStiffness;
}

//FUNCTION_BEGIN_PARSE
realT
PenaltyCoulombIntermediateParameterData::Stress(const realT normalApproach) const
{
  realT normalStress = 0.0;

  if(normalApproach <= 0)
    normalStress = 0;
  else if(normalApproach < normalApproachYield)
    normalStress = ktildeAperture * log(aperture / (aperture - normalApproach));
  else if(normalApproach < normalApproachSoften)
    normalStress = stressYield + kyield * (normalApproach - normalApproachYield);
  else
    normalStress = stressSoften + ksoften * (normalApproach - normalApproachSoften);

  return normalStress;
}
