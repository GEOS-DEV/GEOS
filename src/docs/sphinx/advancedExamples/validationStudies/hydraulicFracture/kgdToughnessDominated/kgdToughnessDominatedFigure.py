import sys
import os
import argparse


sys.path.append('../../../../../../../inputFiles/hydraulicFracturing/scripts')
import hydrofractureFigure

hydrofractureFigure.main( geosPath='../../../../../../..',xmlFilePrefix='kgdToughnessDominated')
