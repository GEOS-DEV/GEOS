import os
try:
    _this_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _this_dir = os.getcwd()

import wellboreTemperatureFigure
wellboreTemperatureFigure.main(xmlFilePrefix='thermalCompressible_temperatureDependentSinglePhaseThermalConductivity', outputDir=_this_dir)
