//FUNCTION_BEGIN_PARSE
virtual_void
HertzianIntermediate::UpdateProperties(const localIndex index, std::map<std::string, realT>& p0, std::map<std::string, realT>& p1)
{
  if(index >= NumStateIndex0())
    throw GPException("index of out range: HertzianIntermediate::UpdateProperties\n");

  //radius(0), poissons(1), youngs(2), mass(3), rest(4),
  //yield(5), velHalf(6), surfaceEnergy(7), cement(8)

  HertzianIntermediateStateData& hstate = *StateData(index, 0);
  if (hstate.youngs > 0)
    hstate.Update     ( p0["radius"], p1["radius"]);
  else
    hstate.Initialize ( p0["radius"], p1["radius"],
                        p0["Nu"],     p1["Nu"],
                        p0["E"],      p1["E"],
                        p0["mass"],   p1["mass"],
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}
