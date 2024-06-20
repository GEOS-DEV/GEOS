/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>

#include "common/Format.hpp"

// Test the GEOS_FMT macro for simple formatting
TEST(CommonFormatTests, FormatString) {
    std::string name = "World";
    int year = 2024;
    std::string formatted = GEOS_FMT("Hello, {}! See you in {}", name, year);
    EXPECT_EQ(formatted, "Hello, World! See you in 2024");
}

// Test the GEOS_FMT_TO macro for output iterator formatting
TEST(CommonFormatTests, FormatTo) {
    char buffer[50];
    std::string greeting = "Hello";
    double number = 123.456;
    GEOS_FMT_TO(buffer, sizeof(buffer), "{} number is {:.2f}", greeting, number);
    std::string result(buffer);
    EXPECT_EQ(result, "Hello number is 123.46");
}

// Test the dynamic format argument store for named arguments
TEST(CommonFormatTests, FormatWithArgMap) {
    std::map<std::string, double> namedArgs = {{"pi", 3.14159}, {"e", 2.71828}};
    auto argStore = GEOS_FMT_ARG_MAP(namedArgs);
    std::string formatted = GEOS_VFMT("pi = {pi}, e = {e}", argStore);
    EXPECT_EQ(formatted, "pi = 3.14159, e = 2.71828");
}
