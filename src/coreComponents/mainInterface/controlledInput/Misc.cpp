#include "Misc.hpp"

namespace geos::input
{

real64 convertTime( string const & time )
{
  std::vector< string > tokens = stringutilities::tokenize< std::vector >( time, " " );
  std::size_t const numTokens = tokens.size();
  GEOS_ASSERT( numTokens == 1 || numTokens == 2 );

  real64 const t = std::stod( tokens.front() ); // TODO check cast.

  real64 scaling = 1.;
  if( numTokens == 2 )
  {
    GEOS_ASSERT_EQ( tokens.size(), 2 );
    int constexpr sec = 1;
    int constexpr min = 60;
    int constexpr hour = 60 * min;
    int constexpr day = 24 * hour;
    int constexpr week = 7 * day;
    real64 constexpr year = 365.25 * day;
    real64 constexpr month = year / 12.;

    std::map< string, real64 > const conv{  // TODO To real64 because of the 365.25
      { "s",         sec },
      { "sec",       sec },
      { "second",    sec },
      { "second(s)", sec },
      { "seconds",   sec },
      { "min",       min },
      { "minute",    min },
      { "minute(s)", min },
      { "minutes",   min },
      { "h",         hour },
      { "hour",      hour },
      { "hour(s)",   hour },
      { "hours",     hour },
      { "d",         day },
      { "day",       day },
      { "day(s)",    day },
      { "days",      day },
      { "w",         week },
      { "week",      week },
      { "week(s)",   week },
      { "weeks",     week },
      { "m",         month },
      { "month",     month },
      { "month(s)",  month },
      { "months",    month },
      { "y",         year },
      { "year",      year },
      { "year(s)",   year },
      { "years",     year },
    };
    scaling = conv.at( tokens.back() );// TODO check if value is found
  }
  return t * scaling;
}

string convertYamlElementTypeToGeosElementType( string const yamlElementType )
{
  std::map< string, string > m{
    { "tetrahedra",         "C3D4" },
    { "pyramids",           "C3D5" },
    { "wedges",             "C3D6" },
    { "hexahedra",          "C3D8" },
    { "pentagonal_prism",   "PentagonalPrism" },
    { "hexagonal_prism",    "HexagonalPrism" },
    { "heptagonal_prism",   "HeptagonalPrism" },
    { "octagonal_prism",    "OctagonalPrism" },
    { "nonagonal_prism",    "NonagonalPrism" },
    { "decagonal_prism",    "DecagonalPrism" },
    { "hendecagonal_prism", "HendecagonalPrism" },
    { "polyhedron",         "Polyhedron" },
  };

  auto const geosIt = m.find( yamlElementType );
  if( geosIt == m.cend() )
  {
    GEOS_ERROR( "Could not find element type " << yamlElementType << " in supported element list." );
  }
  return geosIt->second;
}

} // end of namespace