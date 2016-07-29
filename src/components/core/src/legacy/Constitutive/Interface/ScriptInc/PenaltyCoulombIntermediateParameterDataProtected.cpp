//FUNCTION_BEGIN_PARSE
void
PenaltyCoulombIntermediateParameterData::Initialize()
{
  //check user settings
  const realT tol = 1e-6;
  if(normalApproachYield > aperture)
    normalApproachYield = (1.0 - tol) * aperture;
  if(stressYield > stressSoften)
    stressYield = stressSoften;

  //set derived values
  const realT factor = 1.0/(aperture - normalApproachYield);
  ktildeAperture = stressYield / log(aperture * factor);
  kyield = ktildeAperture * factor;
  if(isEqual(stressSoften, std::numeric_limits<realT>::max()))
    normalApproachSoften = std::numeric_limits<realT>::max();
  else
    normalApproachSoften = normalApproachYield + (stressSoften - stressYield) / kyield;
}

//FUNCTION_BEGIN_PARSE
virtual_void
PenaltyCoulombIntermediateParameterData::PostReadXML( const TICPP::HierarchicalDataNode& )
{
  if(aperture <= 0)
    aperture = 1e-2;

  //normal approach at the onset of yielding
  if(normalApproachYield <= 0)
    throw GPException("You must set a normal approach for the onset of yielding > 0");
  if(normalApproachYield >= aperture)
    throw GPException("You must set a normal approach for the onset of yielding < aperture");
  if(stressYield <= 0)
    throw GPException("You must set a stress for the onset of yielding > 0");
  if(stressSoften <= 0)
    stressSoften = std::numeric_limits<realT>::max();
  if(stressSoften < stressYield)
    stressSoften = stressYield;

  //set all derived values
  Initialize();

  //set the rest ... arealStiffnessSoften is not calculated or used in SimpleInitialize
  if(ksoften <= 0)
    ksoften = (1e-5) * kyield;
  if(kshear <= 0)
    kshear = 0.7 * kyield;
}

//FUNCTION_BEGIN_PARSE
virtual_void
PenaltyCoulombIntermediateParameterData::PostSetValue()
{
  const TICPP::HierarchicalDataNode node;
  PostReadXML(node);
}
