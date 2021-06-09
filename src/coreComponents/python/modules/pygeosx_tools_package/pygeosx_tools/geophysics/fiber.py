import os
import numpy as np
from pygeosx_tools import wrapper, parallel_io, plot_tools
import matplotlib.pyplot as plt
from matplotlib import cm
from pyevtk.hl import gridToVTK


class Fiber():
    def __init__(self):
        self.time = []
        self.channel_position = []
        self.gage_length = -1.0

        self.x = []
        self.dudx = []
        self.dudt = []
        self.dudxdt = []


class FiberAnalysis():
    def __init__(self):
        """
        InSAR Analysis class
        """
        self.set_names = []

