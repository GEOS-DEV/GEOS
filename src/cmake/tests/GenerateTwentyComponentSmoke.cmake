# SPDX-License-Identifier: LGPL-2.1-only

# Reuse the nine-component test topology, splitting C10 into 17 identical
# pseudocomponents to exercise the supported upper limit without duplicating it.
file( READ "${INPUT_FILE}" smoke )
foreach( attribute IN ITEMS componentNames componentCriticalPressure componentCriticalTemperature
                            componentAcentricFactor componentMolarWeight componentVolumeShift )
  string( REGEX MATCH "${attribute}=\"([^\"]*)\"" matched "${smoke}" )
  set( values "${CMAKE_MATCH_1}" )
  string( REGEX REPLACE "[{} ]" "" values "${values}" )
  string( REPLACE "," ";" values "${values}" )
  list( SUBLIST values 0 4 extended )
  list( GET values 1 c10Value )
  foreach( component RANGE 4 19 )
    if( attribute STREQUAL "componentNames" )
      list( APPEND extended "C10_${component}" )
    else()
      list( APPEND extended "${c10Value}" )
    endif()
  endforeach()
  list( JOIN extended ", " extended )
  string( REPLACE "${matched}" "${attribute}=\"{ ${extended} }\"" smoke "${smoke}" )
endforeach()

# Inject C10 through its first pseudocomponent. Keeping the original four
# nonzero fractions avoids roundoff in the well's strict sum-to-one check.
set( stream 0.1 0.1 0.1 0.7 )
foreach( component RANGE 4 19 )
  list( APPEND stream 0.0 )
endforeach()
list( JOIN stream ", " stream )
string( REGEX REPLACE "injectionStream=\"[^\"]*\"" "injectionStream=\"{ ${stream} }\"" smoke "${smoke}" )
string( REPLACE "scale=\"0.05\"" "scale=\"0.017647058823529411\"" smoke "${smoke}" )
string( REPLACE "six identical" "seventeen identical" smoke "${smoke}" )
string( REPLACE "nine-component" "twenty-component" smoke "${smoke}" )
set( additionalFields "" )
foreach( component RANGE 9 19 )
  set( fraction 0.017647058823529411 )
  if( component EQUAL 19 )
    set( fraction 0.01764705882352946 )
  endif()
  string( APPEND additionalFields "    <FieldSpecification name=\"initialComposition_C10_${component}\" initialCondition=\"1\"
      setNames=\"{ all }\" objectPath=\"ElementRegions/Region1/cb1\"
      fieldName=\"globalCompFraction\" component=\"${component}\" scale=\"${fraction}\"/>
" )
endforeach()
string( REPLACE "  </FieldSpecifications>" "${additionalFields}  </FieldSpecifications>" smoke "${smoke}" )
file( WRITE "${OUTPUT_FILE}" "${smoke}" )

# Apply the same composition to the boundary-condition regression.
file( READ "${DIRICHLET_INPUT_FILE}" dirichlet )
get_filename_component( smokeName "${OUTPUT_FILE}" NAME )
string( REPLACE "compositional_multiphase_wells_9comp_smoke.xml" "${smokeName}" dirichlet "${dirichlet}" )
string( REPLACE "scale=\"0.05\"" "scale=\"0.017647058823529411\"" dirichlet "${dirichlet}" )
set( additionalFields "" )
foreach( component RANGE 9 19 )
  set( fraction 0.017647058823529411 )
  if( component EQUAL 19 )
    set( fraction 0.01764705882352946 )
  endif()
  string( APPEND additionalFields "    <FieldSpecification name=\"boundaryComponent${component}\" setNames=\"{ rightEnd }\"
      objectPath=\"faceManager\" fieldName=\"globalCompFraction\" component=\"${component}\" scale=\"${fraction}\"/>
" )
endforeach()
string( REPLACE "  </FieldSpecifications>" "${additionalFields}  </FieldSpecifications>" dirichlet "${dirichlet}" )
file( WRITE "${DIRICHLET_OUTPUT_FILE}" "${dirichlet}" )
