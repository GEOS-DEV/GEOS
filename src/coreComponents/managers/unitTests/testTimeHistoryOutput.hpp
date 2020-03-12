#include "managers/Outputs/TimeHistoryOutput.hpp"
#include <gtest/gtest.h>

using namespace geosx;

TEST( testTimeHistoryOutput, ScalarTimeHistory )
{
  TimeHistory hist( );

  TimeHistoryUpdate update_event;
  TimeHistoryOutput output_event;
}