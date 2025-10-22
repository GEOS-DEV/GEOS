#include <gtest/gtest.h>
#include "common/DataTypes.hpp"
#include "prescribedStressPathGeomechanics/OedometricStressPath.hpp"

using namespace geos;

TEST(HelloWorldTest, BasicAssertions) {
    EXPECT_EQ(1, 1);
    std::cout << "Hello, World!" << std::endl;
}

TEST( OedometricTest, StressAssertion )
{
  //OedometricStressPath name OedometricStressPath parent SinglePhaseFlow

  real64 const pressure = 100041.0;
  array1d< real64 > const normal(3);
  normal[0] = 0.0;
  normal[1] = 1.0;
  normal[2] = 0.0;
  EXPECT_DOUBLE_EQ( pressure, 100041.0 );
  std::cout << "pressure " << pressure << std::endl;

  OedometricStressPath oedometric("OedometricStressPath", nullptr );
  real64 stress = oedometric.computeFractureStress(pressure, normal);
  EXPECT_DOUBLE_EQ( stress, 8.489998e7 );
  std::cout << "stress " << stress << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}