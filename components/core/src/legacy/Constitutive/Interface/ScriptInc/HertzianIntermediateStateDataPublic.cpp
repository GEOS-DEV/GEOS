//FUNCTION_BEGIN_PARSE
virtual_void
HertzianIntermediateStateData::Update(const realT curvature1,
                                  const realT curvature2)
{
  radius = EffectiveRadius(curvature1, curvature2);
  //set other derived values:
  //get the Hertzian coefficient, k, where f = k*a^1.5, k = Hertzian coefficient
  hertzCf = youngs * sqrt(radius) * 4.0 / 3.0;
}

//FUNCTION_BEGIN_PARSE
virtual_void
HertzianIntermediateStateData::Initialize(const realT curvature1, const realT curvature2,
                                          const realT poissons1, const realT poissons2,
                                          const realT youngs1, const realT youngs2,
                                          const realT mass1, const realT mass2,
                                          const realT , const realT ,
                                          const realT , const realT ,
                                          const realT , const realT ,
                                          const realT , const realT ,
                                          const realT , const realT )
{
  youngs = EffectiveYoungsModulus(poissons1, poissons2, youngs1, youngs2);
  mass = EffectiveMass(mass1, mass2);

  //set mechanical properties of the contact
  Update(curvature1, curvature2);
}

//FUNCTION_BEGIN_PARSE
realT
HertzianIntermediateStateData::EffectiveRadius(const realT curvature1,
                                               const realT curvature2) const
{
  realT r = curvature1 + curvature2;
  r = 1. / r;
  return r;
}

//FUNCTION_BEGIN_PARSE
realT
HertzianIntermediateStateData::EffectiveYoungsModulus(const realT poissons1,
                                                      const realT poissons2,
                                                      const realT youngs1,
                                                      const realT youngs2) const
{
  realT Eeff = (1 - poissons1 * poissons1) / youngs1;
  Eeff += (1 - poissons2 * poissons2) / youngs2;
  Eeff = 1. / Eeff;
  return Eeff;
}

//FUNCTION_BEGIN_PARSE
realT
HertzianIntermediateStateData::EffectiveMass(const realT mass1,
                                             const realT mass2) const
{
  realT em1 = 1/mass1 + 1/mass2;
  em1 = 1. / em1;
  return em1;
// EBH comment: if either mass1 or mass2 is large due to boundary conditions, etc. (e.g. 1.e100)
// EBH comment: then c++ will produce garbage for the line below.  The fix above accounts for this.
//  return mass1 * mass2 / (mass1 + mass2);
}
