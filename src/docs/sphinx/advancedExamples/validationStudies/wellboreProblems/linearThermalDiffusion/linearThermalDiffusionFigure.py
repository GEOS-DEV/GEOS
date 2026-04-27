import os
import sys

try:
	_this_dir = os.path.dirname(os.path.abspath(__file__))
	_scripts_dir = os.path.normpath(os.path.join(_this_dir, '../../../../../../../inputFiles/singlePhaseFlow/scripts'))
except NameError:
	_src_dir = os.getcwd()
	_this_dir = os.path.join(_src_dir, 'docs', 'sphinx', 'advancedExamples', 'validationStudies', 'wellboreProblems', 'linearThermalDiffusion')
	_scripts_dir = os.path.join(_src_dir, 'inputFiles', 'singlePhaseFlow', 'scripts')
sys.path.append(_scripts_dir)
import wellboreTemperatureFigure

wellboreTemperatureFigure.main(xmlFilePrefix='thermalCompressible_2d', outputDir=_this_dir)
