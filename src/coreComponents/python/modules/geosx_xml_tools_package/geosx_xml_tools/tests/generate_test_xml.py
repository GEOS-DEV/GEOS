"""Tool for generating test xml files for processing."""

from lxml import etree as ElementTree
import os
from geosx_xml_tools import xml_formatter


def generate_test_xml_files(root_dir):
    """Build example input/output xml files, which can be used to test the parser.
     These are derived from a GEOSX integrated test xml.

     @param root_dir The folder to write the example xml files.
    """

    # Build segments of an xml file that can be compiled to form a test
    # File header/footer
    xml_header = """<Problem>"""

    xml_footer = """</Problem>"""

    # Parameters
    xml_parameters = """
<Parameters>
    <Parameter name="permeability" value="2.0e-16" />
    <Parameter name="porosity" value="0.05" />
    <Parameter name="pressure" value="5 [MPa]" />
  </Parameters>"""

    # Includes
    xml_includes = """
<Included>
    <File name="./included/included_a.xml"/>
    <File name="./included/included_b.xml"/>
    <File name="./included/included_c.xml"/>
</Included>"""

    # Base segments
    xml_base_a = """
<Solvers
    gravityVector="{ 0.0, 0.0, -9.81 }">

    <SinglePhaseFlow name="SinglePhaseFlow"
                     logLevel="0"
                     gravityFlag="1"
                     discretization="singlePhaseTPFA"
                     fluidName="water"
                     solidName="rock"
                     targetRegions="{Region1}">
      <LinearSolverParameters krylovTol="1.0e-10"
                              newtonTol="1.0e-6"
                              maxIterNewton="8"/>
    </SinglePhaseFlow>
  </Solvers>

  <Mesh>
    <InternalMesh name="mesh1"
                  elementTypes="{C3D8}"
                  xCoords="{0, 10}"
                  yCoords="{0, 1}"
                  zCoords="{0, 1}"
                  nx="{10}"
                  ny="{1}"
                  nz="{1}"
                  cellBlockNames="{block1}"/>
  </Mesh>

  <Geometry>
    <Box name="source" xMin="{-0.01, -0.01, -0.01}" xMax="{ 1.01, 1.01, 1.01}"/>
    <Box name="sink"   xMin="{ 8.99, -0.01, -0.01}" xMax="{10.01, 1.01, 1.01}"/>
  </Geometry>
"""

    xml_base_b = """
<Events maxTime="2e4">

    <PeriodicEvent name="outputs"
                   timeFrequency="1000.0"
                   targetExactTimestep="1"
                   target="/Outputs/siloOutput" />

    <PeriodicEvent name="solverApplications"
                   forceDt="1e3"
                   target="/Solvers/SinglePhaseFlow" />

    <PeriodicEvent name="restarts"
                   timeFrequency="1e4"
                   targetExactTimestep="0"
                   target="/Outputs/restartOutput"/>

  </Events>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA"
                                 fieldName="pressure"
                                 coefficientName="permeability"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region1"
                   cellBlocks="{block1}"
                   materialList="{water, rock}"/>
  </ElementRegions>

  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  defaultDensity="1000"
                                  defaultViscosity="0.001"
                                  referencePressure="0.0"
                                  referenceDensity="1000"
                                  compressibility="5e-10"
                                  referenceViscosity="0.001"
                                  viscosibility="0.0"/>

    <PoreVolumeCompressibleSolid name="rock"
                                 referencePressure="0.0"
                                 compressibility="1e-9"/>
  </Constitutive>

  <Outputs>
    <Silo name="siloOutput"/>
    <Restart name="restartOutput"/>
  </Outputs>
"""

    # Field specifications with parameters, symbolic math, and their compiled equivalents
    field_string_with_parameters = """
<FieldSpecifications>
    <FieldSpecification name="permx"
                        component="0"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="$permeability$"/>

    <FieldSpecification name="permy"
                        component="1"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="$permeability$"/>

    <FieldSpecification name="permz"
                        component="2"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="$permeability$"/>

    <FieldSpecification name="referencePorosity"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="referencePorosity"
                        scale="$porosity$"/>

    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"/>

    <SourceFlux name="sourceTerm"
                objectPath="ElementRegions/Region1/block1"
                scale="-0.00001"
                setNames="{source}"/>

    <FieldSpecification name="sinkTerm"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"
                        setNames="{sink}"/>
  </FieldSpecifications>
"""

    field_string_with_symbolic = """
<FieldSpecifications>
    <FieldSpecification name="permx"
                        component="0"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="`2*$permeability$`"/>

    <FieldSpecification name="permy"
                        component="1"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="`0.5*$permeability$`"/>

    <FieldSpecification name="permz"
                        component="2"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="`1*$permeability$`"/>

    <FieldSpecification name="referencePorosity"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="referencePorosity"
                        scale="$porosity$"/>

    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"/>

    <SourceFlux name="sourceTerm"
                objectPath="ElementRegions/Region1/block1"
                scale="-0.00001"
                setNames="{source}"/>

    <FieldSpecification name="sinkTerm"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"
                        setNames="{sink}"/>
  </FieldSpecifications>
"""

    field_string_base = """
<FieldSpecifications>
    <FieldSpecification name="permx"
                        component="0"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="2.0e-16"/>

    <FieldSpecification name="permy"
                        component="1"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="2.0e-16"/>

    <FieldSpecification name="permz"
                        component="2"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="2.0e-16"/>

    <FieldSpecification name="referencePorosity"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="referencePorosity"
                        scale="0.05"/>

    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"/>

    <SourceFlux name="sourceTerm"
                objectPath="ElementRegions/Region1/block1"
                scale="-0.00001"
                setNames="{source}"/>

    <FieldSpecification name="sinkTerm"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"
                        setNames="{sink}"/>
  </FieldSpecifications>
"""

    field_string_alt = """
<FieldSpecifications>
    <FieldSpecification name="permx"
                        component="0"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="4e-16"/>

    <FieldSpecification name="permy"
                        component="1"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="1e-16"/>

    <FieldSpecification name="permz"
                        component="2"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="permeability"
                        scale="2e-16"/>

    <FieldSpecification name="referencePorosity"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="referencePorosity"
                        scale="0.05"/>

    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"/>

    <SourceFlux name="sourceTerm"
                objectPath="ElementRegions/Region1/block1"
                scale="-0.00001"
                setNames="{source}"/>

    <FieldSpecification name="sinkTerm"
                        objectPath="ElementRegions/Region1/block1"
                        fieldName="pressure"
                        scale="5e6"
                        setNames="{sink}"/>
  </FieldSpecifications>
"""

    # Write the files, and apply pretty_print to targets for easy matches
    # No advanced features case
    with open("%s/no_advanced_features_input.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_base_b + field_string_base + xml_footer)
    with open("%s/no_advanced_features_target.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_base_b + field_string_base + xml_footer)
    xml_formatter.format_file("%s/no_advanced_features_target.xml" % (root_dir))

    # Parameters case
    with open("%s/parameters_input.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_parameters + xml_base_a + xml_base_b + field_string_with_parameters + xml_footer)
    with open("%s/parameters_target.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_base_b + field_string_base + xml_footer)
    xml_formatter.format_file("%s/parameters_target.xml" % (root_dir))

    # Symbolic + parameters case
    with open("%s/symbolic_parameters_input.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_parameters + xml_base_a + xml_base_b + field_string_with_symbolic + xml_footer)
    with open("%s/symbolic_parameters_target.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_base_b + field_string_alt + xml_footer)
    xml_formatter.format_file("%s/symbolic_parameters_target.xml" % (root_dir))

    # Included case
    os.makedirs("%s/included" % (root_dir), exist_ok=True)
    with open("%s/included_input.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_includes + xml_footer)
    with open("%s/included/included_a.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_footer)
    with open("%s/included/included_b.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_b + xml_footer)
    with open("%s/included/included_c.xml" % (root_dir), "w") as f:
        f.write(xml_header + field_string_base + xml_footer)
    with open("%s/included_target.xml" % (root_dir), "w") as f:
        f.write(xml_header + xml_base_a + xml_base_b + field_string_base + xml_footer)
    xml_formatter.format_file("%s/included_target.xml" % (root_dir))
