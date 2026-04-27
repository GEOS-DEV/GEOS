import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.normpath(os.path.join(_this_dir, '../../../../../../../inputFiles/singlePhaseFlow/scripts')))
import wellboreTemperatureFigure

wellboreTemperatureFigure.main(xmlFilePrefix='thermalCompressible_temperatureDependentSinglePhaseThermalConductivity', outputDir=_this_dir)
