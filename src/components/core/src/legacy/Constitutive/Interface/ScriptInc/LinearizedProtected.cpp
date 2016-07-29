//FUNCTION_BEGIN_PARSE
virtual_realT
Linearized::NormalStiffness(const InterfaceBaseParameterData& ,
                            InterfaceBaseStateData& matStateBase,
                            const realT normalApproach,
                            const bool setForces) const
{
  LinearizedStateData& matState = static_cast < LinearizedStateData& > ( matStateBase );

  const realT youngs2 = matState.youngs * matState.youngs;

  realT fnmag_coh = 0.0;
  realT fnmag_mat = 0.0;

  //------COHESION------
  //if < a0, then we are (potentially) cohesive
  if (normalApproach < matState.a0)
  {
    //this is now out of contact with the other body ... unless it is cohesive
    realT coh_stiffness = std::numeric_limits<realT>::min();
    if (matState.isCohesive > 0)
    {
      //cohesive stiffness
      coh_stiffness = matState.at <= matState.ac ?
          matState.linearStiffnessElastic :
          PlasticUnloadingSlope(matState.youngs, youngs2,
                                matState.mass, matState.vark2,
                                matState.linearStiffnessElastic,
                                matState.at, matState.ac, matState.halfRestVel);
      //get the cohesion distance
      if (setForces)
      {
        realT dcoh0 = matState.pullOffForce / (-1.0 * coh_stiffness);
        if(normalApproach >= (matState.a0 - dcoh0))
        {
          fnmag_coh = coh_stiffness * (matState.a0 - dcoh0 - normalApproach);
          matState.ElasticStrainEnergy += 0.5 * fnmag_coh * fnmag_coh / coh_stiffness;
          matState.stress += fnmag_coh;
        }
      }
    }
    return coh_stiffness;
  }

  //------GRAIN LOADING/UNLOADING------
  //otherwise, contact is in a different regime ... note dcoh0 = 0 for !is_cohesive
  realT at; //at >= ac ... if > ac, then it is yielding
  if(setForces)
  {
    matState.at = matState.at < matState.ac ? matState.ac : matState.at;
    at = matState.at;
  }
  else
  {
    at = matState.at < matState.ac ? matState.ac : matState.at;
  }

  //previous update stored in the "store" state
  if (normalApproach > at)
  {
    //------POST-YIELD LOADING------
    at = normalApproach;
    const realT k2 = PlasticUnloadingSlope(matState.youngs, youngs2,
                                           matState.mass, matState.vark2,
                                           matState.linearStiffnessElastic,
                                           at, matState.ac, matState.halfRestVel);
    if(setForces)
    {
      matState.at = at;
      if ( !isEqual( k2, matState.linearStiffnessElastic ) )
        matState.a0 = (k2 - matState.linearStiffnessElastic) * at / k2;

      fnmag_mat = k2 * (normalApproach - matState.a0);
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / k2;
      //why was I adding in a cohesive force during post-yield loading????
      //fnmag_coh = matState.pullOffForce; //just a constant value here ... could alter this in the future
      matState.stress += fnmag_mat;//+fnmag_coh
    }
    return k2;
  }
  else if (at > matState.ac)
  {
    //------POST-YIELD UNLOADING------
    const realT k2 = PlasticUnloadingSlope(matState.youngs, youngs2,
                                           matState.mass, matState.vark2,
                                           matState.linearStiffnessElastic,
                                           at, matState.ac,
                                           matState.halfRestVel);
    if (setForces)
    {
      fnmag_mat = k2 * (normalApproach - matState.a0);
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / k2;
      matState.stress += fnmag_mat;
    }

    return k2;
  }
  else
  {
    //------ELASTIC LOADING/UNLOADING------
    if (setForces)
    {
      //still elastic
      fnmag_mat = matState.linearStiffnessElastic * normalApproach;
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / matState.linearStiffnessElastic;
      matState.stress += fnmag_mat;
    }
    return matState.linearStiffnessElastic;
  }
}

//FUNCTION_BEGIN_PARSE
realT
Linearized::PlasticUnloadingSlope(const realT Eeff, const realT Eeff2, const realT tmass,
                                  const int vark2, const realT linearStiffnessElastic,
                                  const realT at, const realT ac, const realT velHalf) const
{
  if (vark2!=0)
  {
    //Get variable latching term
    const realT slope = 3 * sqrt(linearStiffnessElastic / tmass) / velHalf;

    //Calculate k2 with saturation at e = 0.1
    realT k2 = at > ac ? linearStiffnessElastic + slope * linearStiffnessElastic * (at - ac) : linearStiffnessElastic;

    if (k2 > 100.0 * linearStiffnessElastic) //saturate at 100*K1, e = sqrt(K1/K2) = 0.1
      k2 = 100.0 * linearStiffnessElastic;
    return k2;
  }
  else
  {
    return linearStiffnessElastic / Eeff2;
  }
}
