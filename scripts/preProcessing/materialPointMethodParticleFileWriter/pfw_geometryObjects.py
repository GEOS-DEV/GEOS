# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
Updated on Wed Apr 17 10:30:00 2024
@author: homel1, crook5
Geometry object functions for the particle file writer.
"""

import numpy as np   
# from typing import Annotated, TypeVar, Literal, Union, List
# import numpy.typing as npt                
from sklearn.neighbors import KDTree         
import random
import math
from scipy.optimize import minimize_scalar
from scipy.spatial import Voronoi, ConvexHull
from abc import ABCMeta, abstractmethod
import scipy
import itertools
from datetime import datetime
np.random.seed(1)

from mpi4py import MPI
comm = MPI.COMM_WORLD # gets communication pool 
g_rank = comm.Get_rank()  # gets rank of current process 
num_ranks = comm.Get_size() # total number of processes

# ===========================================
# TYPES
# ===========================================
floatType = (int, float, np.byte, np.ubyte, np.short, np.ushort, np.intc, np.uintc, np.long, np.longlong, np.ulonglong, np.half, np.single, np.double, np.longdouble)

def type_check_scalar(s, varName, scalarType, canBeFunc=False):
  if not isinstance(s, scalarType):
    if canBeFunc:
      if not callable(s):
        raise TypeError( varName + " must be callable or have type: " + ", ".join([t.__name__ for t in scalarType]))
    else:
      raise TypeError( varName + " must have type: " + ", ".join([t.__name__ for t in scalarType]))

def type_check_array(arr, varName, size, scalarType, canBeFunc=False):
  if isinstance(arr, (list, np.ndarray)):
    found = False
    if type(size) is int:
          found = True
    elif type(size) is list or type(size) is tuple:
      for s in size:
        if len(arr) == s:
            found = True
            break
    assert found, varName + " must have dimensions from: " + ", ".join([str(s) for s in size])
    for a in arr:
      assert isinstance(a, scalarType), varName + " elements must be of type: " + ", ".join([t.__name__ for t in scalarType])
  elif not canBeFunc or not callable(arr):
      raise TypeError(varName + " must be callable if not a list")


# DType = TypeVar("DType", bound=np.generic)
# Array2 = Annotated[npt.NDArray[DType], Literal[2]]
# fArray3 = Annotated[list[float], 3]
# Vector = Union[Array2, Array3]

# ===========================================
# END TYPES
# ===========================================


# ===========================================
# TYPES
# ===========================================
floatType = (int, float, np.byte, np.ubyte, np.short, np.ushort, np.intc, np.uintc, np.long, np.longlong, np.ulonglong, np.half, np.single, np.double, np.longdouble)

def type_check_scalar(s, varName, scalarType, canBeFunc=False):
  if not isinstance(s, scalarType):
    if canBeFunc:
      if not callable(s):
        raise TypeError( varName + " must be callable or have type: " + ", ".join([t.__name__ for t in scalarType]))
    else:
      raise TypeError( varName + " must have type: " + ", ".join([t.__name__ for t in scalarType]))

def type_check_array(arr, varName, size, scalarType, canBeFunc=False):
  if isinstance(arr, (list, np.ndarray)):
    found = False
    if type(size) is int:
          found = True
    elif type(size) is list or type(size) is tuple:
      for s in size:
        if len(arr) == s:
            found = True
            break
    assert found, varName + " must have dimensions from: " + ", ".join([str(s) for s in size])
    for a in arr:
      assert isinstance(a, scalarType), varName + " elements must be of type: " + ", ".join([t.__name__ for t in scalarType])
  elif not canBeFunc or not callable(arr):
      raise TypeError(varName + " must be callable if not a list")


# DType = TypeVar("DType", bound=np.generic)
# Array2 = Annotated[npt.NDArray[DType], Literal[2]]
# fArray3 = Annotated[list[float], 3]
# Vector = Union[Array2, Array3]

# ===========================================
# END TYPES
# ===========================================


# ===========================================
# DEFAULTS
# ===========================================

_defaultMat = 0
_defaultGroup = 0
_defaultVelocity = np.array([0.0, 0.0, 0.0])
_defaultDamage = 0.0
_defaultPorosity = 0.0
_defaultStrengthScale = 1.0
_defaultTemperature = 300.0
_defaultSurfaceFlag = 2 # Basic surface flag (0 = internal, 1=fully damaged, 2=normal surface flag, 3=cohesive)
_defaultSurfaceNormal = np.array([0.0, 0.0, 0.0]) # 0 vector turns off surface normal use in solver (defaults to implicit surface normals)
_defaultSurfacePosition = np.array([0.0, 0.0, 0.0])
_defaultMatDir = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) 
_defaultSurfaceTraction = np.array([0.0,0.0,0.0])
_defaultParticleType = 2 # 0 (Single point), 1 (Single Point Bsplines), 2 (CPDI), 3 (CPTI), 4 (CPDI2)
_defaultCZTag = 0

# ===========================================
# END DEFAULTS
# ===========================================


# ===========================================
_gyroid_rhoVsC = np.array([[0.0000000000000000, 0.0000000000000000],
                          [0.0097600000000000, 0.0158666295635848],
                          [0.0204160000000000, 0.0317332591271696],
                          [0.0304960000000000, 0.0475998886907544],
                          [0.0399040000000000, 0.0634665182543393],
                          [0.0499840000000000, 0.0793331478179241],
                          [0.0604000000000000, 0.0951997773815089],
                          [0.0710240000000000, 0.1110664069450940],
                          [0.0810080000000000, 0.1269330365086790],
                          [0.0908480000000000, 0.1427996660722630],
                          [0.1011200000000000, 0.1586662956358480],
                          [0.1116800000000000, 0.1745329251994330],
                          [0.1209920000000000, 0.1903995547630180],
                          [0.1322400000000000, 0.2062661843266030],
                          [0.1411680000000000, 0.2221328138901870],
                          [0.1519200000000000, 0.2379994434537720],
                          [0.1621920000000000, 0.2538660730173570],
                          [0.1721760000000000, 0.2697327025809420],
                          [0.1821120000000000, 0.2855993321445270],
                          [0.1924480000000000, 0.3014659617081110],
                          [0.2032960000000000, 0.3173325912716960],
                          [0.2130880000000000, 0.3331992208352810],
                          [0.2232160000000000, 0.3490658503988660],
                          [0.2331040000000000, 0.3649324799624510],
                          [0.2442080000000000, 0.3807991095260360],
                          [0.2541920000000000, 0.3966657390896200],
                          [0.2642720000000000, 0.4125323686532050],
                          [0.2740640000000000, 0.4283989982167900],
                          [0.2852000000000000, 0.4442656277803750],
                          [0.2943680000000000, 0.4601322573439600],
                          [0.3052800000000000, 0.4759988869075440],
                          [0.3150240000000000, 0.4918655164711290],
                          [0.3260160000000000, 0.5077321460347140],
                          [0.3362400000000000, 0.5235987755982990],
                          [0.3458880000000000, 0.5394654051618840],
                          [0.3569280000000000, 0.5553320347254680],
                          [0.3675520000000000, 0.5711986642890530],
                          [0.3776800000000000, 0.5870652938526380],
                          [0.3878560000000000, 0.6029319234162230],
                          [0.3992320000000000, 0.6187985529798080],
                          [0.4091680000000000, 0.6346651825433930],
                          [0.4196800000000000, 0.6505318121069770],
                          [0.4294880000000000, 0.6663984416705620],
                          [0.4405280000000000, 0.6822650712341470],
                          [0.4508480000000000, 0.6981317007977320],
                          [0.4613600000000000, 0.7139983303613170],
                          [0.4727360000000000, 0.7298649599249010],
                          [0.4826400000000000, 0.7457315894884860],
                          [0.4941600000000000, 0.7615982190520710],
                          [0.5037600000000000, 0.7774648486156560],
                          [0.5157120000000000, 0.7933314781792410],
                          [0.5259360000000000, 0.8091981077428260],
                          [0.5365120000000000, 0.8250647373064100],
                          [0.5471200000000000, 0.8409313668699950],
                          [0.5584480000000000, 0.8567979964335800],
                          [0.5691040000000000, 0.8726646259971650],
                          [0.5805760000000000, 0.8885312555607500],
                          [0.5908640000000000, 0.9043978851243340],
                          [0.6026240000000000, 0.9202645146879190],
                          [0.6137600000000000, 0.9361311442515040],
                          [0.6247040000000000, 0.9519977738150890],
                          [0.6359360000000000, 0.9678644033786740],
                          [0.6479040000000000, 0.9837310329422580],
                          [0.6594720000000000, 0.9995976625058430],
                          [0.6711360000000000, 1.0154642920694300],
                          [0.6830400000000000, 1.0313309216330100],
                          [0.6943360000000000, 1.0471975511966000],
                          [0.7055680000000000, 1.0630641807601800],
                          [0.7171840000000000, 1.0789308103237700],
                          [0.7285120000000000, 1.0947974398873500],
                          [0.7403840000000000, 1.1106640694509400],
                          [0.7519040000000000, 1.1265306990145200],
                          [0.7631360000000000, 1.1423973285781100],
                          [0.7747520000000000, 1.1582639581416900],
                          [0.7862400000000000, 1.1741305877052800],
                          [0.7979520000000000, 1.1899972172688600],
                          [0.8096160000000000, 1.2058638468324500],
                          [0.8206560000000000, 1.2217304763960300],
                          [0.8325280000000000, 1.2375971059596200],
                          [0.8447680000000000, 1.2534637355232000],
                          [0.8568640000000000, 1.2693303650867900],
                          [0.8688320000000000, 1.2851969946503700],
                          [0.8807840000000000, 1.3010636242139500],
                          [0.8925920000000000, 1.3169302537775400],
                          [0.9049440000000000, 1.3327968833411200],
                          [0.9166560000000000, 1.3486635129047100],
                          [0.9291840000000000, 1.3645301424682900],
                          [0.9419200000000000, 1.3803967720318800],
                          [0.9544000000000000, 1.3962634015954600],
                          [0.9667520000000000, 1.4121300311590500],
                          [0.9769760000000000, 1.4279966607226300],
                          [0.9847200000000000, 1.4438632902862200],
                          [0.9913920000000000, 1.4597299198498000],
                          [0.9959680000000000, 1.4755965494133900],
                          [0.9992000000000000, 1.4914631789769700],
                          [1.0000000000000000, 1.5073298085405600]])

_schwarzDiamond_rhoVsC = np.array([[0.000000000000000, 0.0000000000000000],
                                  [0.019200000000000, 0.0158666295635848],
                                  [0.040416000000000, 0.0317332591271696],
                                  [0.057728000000000, 0.0475998886907544],
                                  [0.076832000000000, 0.0634665182543393],
                                  [0.095936000000000, 0.0793331478179241],
                                  [0.114272000000000, 0.0951997773815089],
                                  [0.133504000000000, 0.1110664069450940],
                                  [0.151264000000000, 0.1269330365086790],
                                  [0.169600000000000, 0.1427996660722630],
                                  [0.189088000000000, 0.1586662956358480],
                                  [0.208128000000000, 0.1745329251994330],
                                  [0.225792000000000, 0.1903995547630180],
                                  [0.244800000000000, 0.2062661843266030],
                                  [0.263328000000000, 0.2221328138901870],
                                  [0.281024000000000, 0.2379994434537720],
                                  [0.301184000000000, 0.2538660730173570],
                                  [0.317600000000000, 0.2697327025809420],
                                  [0.337376000000000, 0.2855993321445270],
                                  [0.355936000000000, 0.3014659617081110],
                                  [0.373216000000000, 0.3173325912716960],
                                  [0.391840000000000, 0.3331992208352810],
                                  [0.411136000000000, 0.3490658503988660],
                                  [0.429696000000000, 0.3649324799624510],
                                  [0.446880000000000, 0.3807991095260360],
                                  [0.465024000000000, 0.3966657390896200],
                                  [0.485088000000000, 0.4125323686532050],
                                  [0.501824000000000, 0.4283989982167900],
                                  [0.522368000000000, 0.4442656277803750],
                                  [0.540704000000000, 0.4601322573439600],
                                  [0.559232000000000, 0.4759988869075440],
                                  [0.577888000000000, 0.4918655164711290],
                                  [0.596512000000000, 0.5077321460347140],
                                  [0.616672000000000, 0.5235987755982990],
                                  [0.634912000000000, 0.5394654051618840],
                                  [0.651552000000000, 0.5553320347254680],
                                  [0.671808000000000, 0.5711986642890530],
                                  [0.690528000000000, 0.5870652938526380],
                                  [0.710144000000000, 0.6029319234162230],
                                  [0.728672000000000, 0.6187985529798080],
                                  [0.748064000000000, 0.6346651825433930],
                                  [0.768736000000000, 0.6505318121069770],
                                  [0.786112000000000, 0.6663984416705620],
                                  [0.805600000000000, 0.6822650712341470],
                                  [0.825088000000000, 0.6981317007977320],
                                  [0.843840000000000, 0.7139983303613170],
                                  [0.859456000000000, 0.7298649599249010],
                                  [0.874016000000000, 0.7457315894884860],
                                  [0.888448000000000, 0.7615982190520710],
                                  [0.900480000000000, 0.7774648486156560],
                                  [0.911104000000000, 0.7933314781792410],
                                  [0.922560000000000, 0.8091981077428260],
                                  [0.932800000000000, 0.8250647373064100],
                                  [0.941952000000000, 0.8409313668699950],
                                  [0.950720000000000, 0.8567979964335800],
                                  [0.959712000000000, 0.8726646259971650],
                                  [0.966944000000000, 0.8885312555607500],
                                  [0.973728000000000, 0.9043978851243340],
                                  [0.979616000000000, 0.9202645146879190],
                                  [0.985120000000000, 0.9361311442515040],
                                  [0.990368000000000, 0.9519977738150890],
                                  [0.994400000000000, 0.9678644033786740],
                                  [0.998016000000000, 0.9837310329422580],
                                  [0.999936000000000, 0.9995976625058430],
                                  [1.000000000000000, 1.0154642920694300]])

_schwarzPrimitive_rhoVsC = np.array([[0.000000000000000, 0.0000000000000000],
                                    [0.017808000000000, 0.0317332591271696],
                                    [0.036560000000000, 0.0634665182543393],
                                    [0.053792000000000, 0.0951997773815089],
                                    [0.072536000000000, 0.1269330365086790],
                                    [0.089680000000000, 0.1586662956358480],
                                    [0.108160000000000, 0.1903995547630180],
                                    [0.125728000000000, 0.2221328138901870],
                                    [0.143952000000000, 0.2538660730173570],
                                    [0.161592000000000, 0.2855993321445270],
                                    [0.180024000000000, 0.3173325912716960],
                                    [0.197624000000000, 0.3490658503988660],
                                    [0.215936000000000, 0.3807991095260360],
                                    [0.234152000000000, 0.4125323686532050],
                                    [0.251992000000000, 0.4442656277803750],
                                    [0.270160000000000, 0.4759988869075440],
                                    [0.287896000000000, 0.5077321460347140],
                                    [0.306312000000000, 0.5394654051618840],
                                    [0.324120000000000, 0.5711986642890530],
                                    [0.342648000000000, 0.6029319234162230],
                                    [0.360176000000000, 0.6346651825433930],
                                    [0.378824000000000, 0.6663984416705620],
                                    [0.396680000000000, 0.6981317007977320],
                                    [0.415192000000000, 0.7298649599249010],
                                    [0.433024000000000, 0.7615982190520710],
                                    [0.451384000000000, 0.7933314781792410],
                                    [0.469320000000000, 0.8250647373064100],
                                    [0.488088000000000, 0.8567979964335800],
                                    [0.505992000000000, 0.8885312555607500],
                                    [0.525176000000000, 0.9202645146879190],
                                    [0.542528000000000, 0.9519977738150890],
                                    [0.562360000000000, 0.9837310329422580],
                                    [0.579016000000000, 1.0154642920694300],
                                    [0.595840000000000, 1.0471975511966000],
                                    [0.610416000000000, 1.0789308103237700],
                                    [0.625440000000000, 1.1106640694509400],
                                    [0.638520000000000, 1.1423973285781100],
                                    [0.651656000000000, 1.1741305877052800],
                                    [0.664400000000000, 1.2058638468324500],
                                    [0.676448000000000, 1.2375971059596200],
                                    [0.687976000000000, 1.2693303650867900],
                                    [0.699256000000000, 1.3010636242139500],
                                    [0.710736000000000, 1.3327968833411200],
                                    [0.721248000000000, 1.3645301424682900],
                                    [0.731832000000000, 1.3962634015954600],
                                    [0.741632000000000, 1.4279966607226300],
                                    [0.751568000000000, 1.4597299198498000],
                                    [0.760976000000000, 1.4914631789769700],
                                    [0.770440000000000, 1.5231964381041400],
                                    [0.779320000000000, 1.5549296972313100],
                                    [0.788280000000000, 1.5866629563584800],
                                    [0.796848000000000, 1.6183962154856500],
                                    [0.805152000000000, 1.6501294746128200],
                                    [0.813176000000000, 1.6818627337399900],
                                    [0.821600000000000, 1.7135959928671600],
                                    [0.829120000000000, 1.7453292519943300],
                                    [0.836536000000000, 1.7770625111215000],
                                    [0.843928000000000, 1.8087957702486700],
                                    [0.851016000000000, 1.8405290293758400],
                                    [0.857976000000000, 1.8722622885030100],
                                    [0.864560000000000, 1.9039955476301800],
                                    [0.871400000000000, 1.9357288067573500],
                                    [0.877672000000000, 1.9674620658845200],
                                    [0.884200000000000, 1.9991953250116900],
                                    [0.890080000000000, 2.0309285841388600],
                                    [0.896184000000000, 2.0626618432660300],
                                    [0.901752000000000, 2.0943951023932000],
                                    [0.907424000000000, 2.1261283615203700],
                                    [0.912848000000000, 2.1578616206475300],
                                    [0.918208000000000, 2.1895948797747000],
                                    [0.923344000000000, 2.2213281389018700],
                                    [0.928224000000000, 2.2530613980290400],
                                    [0.932832000000000, 2.2847946571562100],
                                    [0.937784000000000, 2.3165279162833800],
                                    [0.942248000000000, 2.3482611754105500],
                                    [0.946624000000000, 2.3799944345377200],
                                    [0.950896000000000, 2.4117276936648900],
                                    [0.955056000000000, 2.4434609527920600],
                                    [0.958928000000000, 2.4751942119192300],
                                    [0.962624000000000, 2.5069274710464000],
                                    [0.966304000000000, 2.5386607301735700],
                                    [0.969856000000000, 2.5703939893007400],
                                    [0.973176000000000, 2.6021272484279100],
                                    [0.976304000000000, 2.6338605075550800],
                                    [0.979400000000000, 2.6655937666822500],
                                    [0.982144000000000, 2.6973270258094200],
                                    [0.984960000000000, 2.7290602849365900],
                                    [0.987488000000000, 2.7607935440637600],
                                    [0.989872000000000, 2.7925268031909300],
                                    [0.992208000000000, 2.8242600623181000],
                                    [0.994184000000000, 2.8559933214452700],
                                    [0.995872000000000, 2.8877265805724400],
                                    [0.997464000000000, 2.9194598396996100],
                                    [0.998800000000000, 2.9511930988267800],
                                    [0.999728000000000, 2.9829263579539500],
                                    [1.000000000000000, 3.0146596170811100]])

# ===========================================

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

log_file = "job_file_log"
def log2file(msg):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    msg = current_time + " Rank " + str(g_rank) + ": " + msg + "\n"

    # Prefer MPI-IO when available, but keep pfw importable in lightweight
    # preflight environments that provide only COMM_WORLD.
    try:
      mode = MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND
      file_handle = MPI.File.Open(MPI.COMM_SELF, log_file + "_" + str(g_rank), mode)
      file_handle.Set_atomicity(True)

      b = bytearray()
      b.extend(map(ord, msg))

      file_handle.Write(b)
      file_handle.Sync()
      file_handle.Close()
    except Exception:
      with open(log_file + "_" + str(g_rank), 'a') as f:
        f.write(msg)


def print2file(file_name, text):
    with open(file_name, 'a') as f:
        f.write(text + "\n")


def countFileLines(fname):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count

def fOffsetLineNum(fname, line_number):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    if line_number == 0:
        return 0

    count = 0
    with open(fname, "rb") as f:
        for buf in _make_gen(f.raw.read):
            count += buf.count(b"\n")
            if count >= line_number:
                chunk_pos = f.tell()
                break

        # Back track to line_number
        for ci in range(len(buf)-1, -1, -1):
            c = buf[ci]
            # print("Character ",c, )
            chunk_pos -= 1
            if c == 10: #b"\n":
                if count == line_number:
                    return chunk_pos + 1 # Add one char because first char of new line immediately follows new line char
                count -= 1


#Might be a faster way to do this
def fileOffsetFromLineNum(fname, line_number):
    char_pos = 0
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i == line_number:
                offset = len(line)
                return char_pos
            else:
                char_pos += len(line)


# used for periodicity
def mapToRange(x, xmin, xmax):
    # Make sure all inputs are np arrays
    xmin = np.asarray(xmin)
    xmax = np.asarray(xmax)
    dx = xmax-xmin
    x = np.asarray(x) - xmin

    x = x * (x < dx) + (x >= dx) * (np.fmod(x, dx))
    x = x * (x >= 0) + (x < 0) * (dx - np.fmod(np.absolute(x), dx))

    return x + xmin


# Smooth sigmoid interpolation between rf and rf0
def smoothStep3(r, rf, rf0):
  if r >= rf0:
    return 0

  if r <= rf:
    return 1

  r = (rf0-r)/(rf0-rf)

  return r*r*(3-2*r)


def inside_box(x, x0=np.array([0.0,0.0,0.0]), dx=np.array([1.0,1.0,1.0]), periodic=[False, False, False]):
  x = np.asarray(x, dtype=float)
  return np.all(np.logical_or(periodic, np.logical_and(x >= x0, x < x0+dx)))


# Checks if point x lies outside of box with minimum corner at x0 and dimensions dx
def outside(x, x0=np.array([0.0,0.0,0.0]), dx=np.array([1.0,1.0,1.0]), periodic=[False, False, False]):
    x = periodic * (np.fmod(x - x0, dx) + x0) + np.logical_not(periodic) * x
    return np.any(np.logical_or(x < x0, x > x0 + dx)), x
    

def generate_orthonormal_axes(v):
    """
    Given a 3D vector v, return two unit vectors u2, u3 that for an orthonormal basis
    """
    v = np.asarray(v, dtype=float)
    norm_v = np.linalg.norm(v)
    if norm_v == 0:
        raise ValueError("Input vector v must be nonzero")

    # First axis aligned with v
    u1 = v / norm_v

    # Pick a vector not parallel to u1 to start the construction
    if np.abs(u1[0]) < 0.9:
        # u1 is not too aligned with x axis, use x as helper
        helper = np.array([1.0, 0.0, 0.0])
    else:
        # otherwise use y as helper
        helper = np.array([0.0, 1.0, 0.0])

    # Make u2 orthogonal to u1 by subtracting the projection of helper onto u1
    proj = np.dot(helper, u1) * u1
    u2 = helper - proj
    u2 /= np.linalg.norm(u2)

    # u3 is orthogonal to both u1 and u2
    u3 = np.cross(u1, u2)

    return u1, u2, u3

def matDirs_from_matDir(v):
  dir1,dir2,dir3 = generate_orthonormal_axes(v)
  matDirs = []
  matDirs.append(np.vstack((dir1, dir2, dir3)))
  return matDirs[0]


class maxFlawRadius:
    def __init__(self,
                 maxFlawRadius_file,
                 simgrid,
                 grid,
                 x0,
                 x1,
                 pores,
                 dim: int=3,
                 periodic=[False, False, False],
                 readFromFile=False):
        self.dim = dim
        self.x0 = np.asarray(x0[:self.dim]) # [x0, y0] should define box bigger than subObject
        self.x1 = np.asarray(x1[:self.dim]) # [x1, y1] should define box bigger than subObject
        self.dx = self.x1-self.x0
        self.periodic = periodic[:self.dim]
        self.pores = pores
        self.simcells = np.asarray(simgrid[:self.dim])
        self.cells = np.asarray(grid[:self.dim])
        self.maxFlawRadius_file = maxFlawRadius_file
       
        print("Setting up max flaw radius field")
        if not readFromFile:
            self.ranks_per_side = int(np.floor(np.cbrt(g_num_ranks)))
            self.num_ranks_used = self.ranks_per_side**self.dim
            self.comm = g_comm_world.Create(g_comm_world.group.Incl([i for i in range(self.num_ranks_used)]))

            self.rank_coord = self.rankIndexToCoord(g_rank)

            self.ranks_in_x = self.ranks_per_side
            self.ranks_in_y = self.ranks_per_side
            self.ranks_in_z = self.ranks_per_side

            self.simcell_range_x = [math.ceil(self.simcells[0]*self.rank_coord[0]/self.ranks_in_x), math.ceil(self.simcells[0]*(self.rank_coord[0]+1)/self.ranks_in_x)]
            self.simcell_range_y = [math.ceil(self.simcells[1]*self.rank_coord[1]/self.ranks_in_y), math.ceil(self.simcells[1]*(self.rank_coord[1]+1)/self.ranks_in_y)]
            self.simcell_range_z = [math.ceil(self.simcells[2]*self.rank_coord[2]/self.ranks_in_z), math.ceil(self.simcells[2]*(self.rank_coord[2]+1)/self.ranks_in_z)]

            self.cell_range_x = [math.ceil(self.cells[0]*self.rank_coord[0]/self.ranks_in_x), math.ceil(self.cells[0]*(self.rank_coord[0]+1)/self.ranks_in_x)]
            self.cell_range_y = [math.ceil(self.cells[1]*self.rank_coord[1]/self.ranks_in_y), math.ceil(self.cells[1]*(self.rank_coord[1]+1)/self.ranks_in_y)]
            self.cell_range_z = [math.ceil(self.cells[2]*self.rank_coord[2]/self.ranks_in_z), math.ceil(self.cells[2]*(self.rank_coord[2]+1)/self.ranks_in_z)]

            start_time = time.perf_counter()
            if g_rank < self.num_ranks_used:
                print("Creating distance field")
                log2file("Creating distance field for max flaw radius object...")
                self.createDistanceField()
                log2file("Finding critical points for max flaw radius object...")
                print("Finding critical points")
                self.findCriticalPoints()
                log2file("Propagating max flaw radius from critical points for max flaw radius object...")
                print("Creating max flaw radius field")
                # self.createMaxFlawRadiusField()
                self.watershed_cf()

            g_comm_world.Barrier()
            log2file("Syncing max flaw radius field across ranks")
            self.syncField()
            g_comm_world.Barrier()

            end_time = time.perf_counter()
            print(f"Computation time was {(end_time-start_time)/60.0:.3f} min")
            log2file(f"Computation time was {(end_time-start_time)/60.0:.3f} min")
            
            if g_rank == 0:
                print("Writing max flaw radius to file")
                self.write2File()

            print("Max flaw radius field generation complete")
            log2file("Finished creating max flaw radius field")
        else:
            if g_rank == 0:
                print("Reading max flaw radius from file " + self.maxFlawRadius_file)
            log2file("Reading max flaw radius from file " + self.maxFlawRadius_file)
            
            self.readFromFile()
            
            if g_rank == 0:
                print("Max flaw radius read from " + maxFlawRadius_file)
            log2file("Max flaw radius read from " + maxFlawRadius_file)


    # For periodic boundaries list of pores should include their images across those boundaries as well
    def createDistanceField(self):
        df_rank = np.full(self.cells, np.inf)

        # Loop over every cell and check distance to nearest surface
        for i in range(self.cell_range_x[0], self.cell_range_x[1]):
          for j in range(self.cell_range_y[0], self.cell_range_y[1]):
            for k in range(self.cell_range_z[0], self.cell_range_z[1]):
              pos = self.dx*(np.array([i,j,k])+0.5*np.ones(self.dim))/self.cells+self.x0

              dist2pore = [None] * len(self.pores)
              for pi in range(0, len(self.pores)):
                  p = self.pores[pi]
                  dx = (p[1] - pos[0])
                  dy = (p[2] - pos[1])
                  dz = (p[3] - pos[2])
                  dist2pore[pi] = math.sqrt(dx*dx + dy*dy + dz*dz) - p[0]

              for dd in range(self.dim):
                  if not self.periodic[dd]:
                    dist2pore.append(abs(pos[dd] - self.x0[dd]))
                    dist2pore.append(abs(self.x1[dd] - pos[dd]))

              df_rank[i][j][k] = min(dist2pore)

        self.df = np.empty(self.cells, dtype=float)
        self.comm.Allreduce([df_rank, MPI.DOUBLE], [self.df, MPI.DOUBLE], op=MPI.MIN)


    def findCriticalPoints(self):
        #Find critical points of distance field
        crit_data = np.empty((0, self.dim+1), dtype=float)

        for i in range(self.cell_range_x[0], self.cell_range_x[1]) if self.periodic[0] else range(max(1, self.cell_range_x[0]), min(self.cells[0]-2, self.cell_range_x[1])):
          for j in range(self.cell_range_y[0], self.cell_range_y[1]) if self.periodic[1] else range(max(1, self.cell_range_y[0]), min(self.cells[1]-2, self.cell_range_y[1])):
            for k in range(self.cell_range_z[0], self.cell_range_z[1]) if self.periodic[2] else range(max(1, self.cell_range_z[0]), min(self.cells[2]-2, self.cell_range_z[1])):
              if (self.df[i][j][k] >= 0) or \
                 (self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][j][k] and self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][j][k]) or \
                 (self.df[i][j][k] > self.df[i][wrap(j-1, 0, self.cells[1]-1)][k] and self.df[i][j][k] > self.df[i][wrap(j+1, 0, self.cells[1]-1)][k]) or \
                 (self.df[i][j][k] > self.df[i][j][wrap(k-1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[i][j][wrap(k+1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][k] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][k]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][k] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][k]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][j][wrap(k-1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][j][wrap(k+1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][j][wrap(k+1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][j][wrap(k-1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[i][wrap(j+1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[i][wrap(j-1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[i][wrap(j+1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[i][wrap(j-1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)]) or \
                 (self.df[i][j][k] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)][wrap(k+1, 0, self.cells[2]-1)] and self.df[i][j][k] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)][wrap(k-1, 0, self.cells[2]-1)]):
                crit_data = np.append(crit_data, [[i, j, k, self.df[i][j][k]]], axis=0)

        # # Syncs cp data asynchronously (doesn't work in 3D but does in 2D, why?) # CC debug
        # # Async send
        # for r in range(self.num_ranks_used):
        #   if r  == g_rank:
        #     continue

        #   self.comm.Isend([np.array(crit_data.shape[0]), MPI.INT], dest=r, tag=0)
        #   self.comm.Isend([crit_data, MPI.DOUBLE], dest=r, tag=1)

        # # Recv
        # self.crit_pts = np.copy(crit_data)
        # for r in range(self.num_ranks_used):
        #   if r == g_rank:
        #     continue

        #   num_pts = np.empty([1,1], dtype=int)
        #   self.comm.Recv([num_pts, MPI.INT], source=r, tag=0)

        #   recv_pts = np.empty([num_pts[0][0],self.dim+1], dtype=float)
        #   self.comm.Recv([recv_pts, MPI.DOUBLE], source=r, tag=1)
          
        #   self.crit_pts = np.append(self.crit_pts, recv_pts, axis=0)

        # Performings syncs serially
        self.crit_pts = np.empty((0, self.dim+1)) # CC debug
        for r in range(self.num_ranks_used):
            cp_buffer = None
            if r == g_rank:
                cp_buffer = np.copy(crit_data)

            cp_buffer = self.comm.bcast(cp_buffer, root=r)
            self.crit_pts = np.append(self.crit_pts, cp_buffer, axis=0)

        self.crit_pts = self.crit_pts[np.argsort(-self.crit_pts[:, self.dim])]


    def createMaxFlawRadiusField(self):
        # Create max flaw radius field
        maxFlawRadius_rank = np.copy(self.df)
        
        # Compute grid search domain for each critical point once based on flaw size
        ds = self.crit_pts[:, self.dim].reshape(self.crit_pts.shape[0],1)
        deltas = np.ceil(np.divide(np.tile(ds,(1,self.dim)), np.tile(np.divide(self.dx, self.cells), (ds.shape[0], 1) ) ))
        
        for i in range(self.cell_range_x[0], self.cell_range_x[1]):
          print("Processing row ", i, " of ", self.cells[0], "...")
          for j in range(self.cell_range_y[0], self.cell_range_y[1]):
            for k in range(self.cell_range_z[0], self.cell_range_z[1]):
              if maxFlawRadius_rank[wrap(i, 0, self.cells[0]-1)][wrap(j,0,self.cells[1]-1)][wrap(k, 0, self.cells[2]-1)] < 0:
                  continue

              for m, c in enumerate(self.crit_pts):
                d = c[self.dim]
                delta = deltas[m]
                if np.any(delta == 0):
                  continue

                dgx = c[0]-i
                dgy = c[1]-j
                dgz = c[2]-k
                
                if self.periodic[0]:
                  dgx = min([abs(dgx), abs(dgx-self.cells[0]), abs(dgx+self.cells[0])])
                
                if self.periodic[1]:
                  dgy = min([abs(dgy), abs(dgy-self.cells[1]), abs(dgy+self.cells[1])])
                
                if self.periodic[2]:
                  dgz = min([abs(dgz), abs(dgz-self.cells[2]), abs(dgz+self.cells[2])])

                if dgx*dgx + dgy*dgy + dgz*dgz < delta[0]*delta[0]:
                    maxFlawRadius_rank[i][j][k] = d
                    break

        self.maxFlawRadius_complete = np.empty(self.cells, dtype=float)
        self.comm.Allreduce([maxFlawRadius_rank, MPI.DOUBLE], [self.maxFlawRadius_complete, MPI.DOUBLE], op=MPI.MAX)

    
    def rankIndexToCoord(self, index):
        i = wrap(index, 0, self.num_ranks_used-1)

        coord = np.array([(i % self.ranks_per_side**2) % self.ranks_per_side,
                          math.floor((i % self.ranks_per_side**2)/self.ranks_per_side),
                          math.floor(i/self.ranks_per_side**2)])

        return coord


    def rankCoordToIndex(self, coord):    
        for dd in range(self.dim):
            if self.periodic[dd]:
                coord[dd] = wrap(coord[dd], 0, self.ranks_per_side-1)

        return coord[0] + self.ranks_per_side*coord[1] + self.ranks_per_side**2*coord[2]


    def getRankNeighbors(self):
         # Get rank indices of neighbors
        neighbors = []
        for k in range(self.dim):
            for f in [-1, 1]:
                offset = np.array([0,0,0])
                offset[k] = f
                neighbor_coord = offset + self.rank_coord
                if self.periodic[k]:
                    neighbor_coord = wrap(neighbor_coord, np.zeros(self.dim), (self.ranks_per_side-1)*np.ones(self.dim)) 

                # Can't communicate with neighbor that is outside grid
                if np.any(neighbor_coord < 0):
                    continue

                neighbors.append(self.rankCoordToIndex(neighbor_coord))

        return neighbors


    def watershed_cf(self):
        neighbors = self.getRankNeighbors()

        print("Propagating critical points using water shed algorithm...")
        log2file("Propagating critical points using water shed algorithm...")
        cell_size = (self.dx[0]/self.cells[0])

        max_df = 0
        maxFlawRadius_rank = np.copy(self.df)
        for cp in self.crit_pts:
            maxFlawRadius_rank[int(cp[0])][int(cp[1])][int(cp[2])] = cp[3]
            if cp[3] > max_df:
                max_df = cp[3]
        max_n = int(max(math.ceil(max_df/cell_size), 1))

        cpIndex = np.empty(np.append(self.cells, self.dim))
        for i in range(0,self.cells[0]):
            for j in range(0,self.cells[1]):
                for k in range(0, self.cells[2]):
                    cpIndex[i][j][k] = np.array([i, j, k])

        i_range = np.arange(self.cell_range_x[0], self.cell_range_x[1], 1)
        j_range = np.arange(self.cell_range_y[0], self.cell_range_y[1], 1)
        k_range = np.arange(self.cell_range_z[0], self.cell_range_z[1], 1)

        for n in range(0,max_n+1):
            print("\tWatershed Iter: ", n, " of ", max_n)
            log2file("\tWatershed Iter: " + str(n) + " of " + str(max_n))
            np.random.shuffle(i_range)
            np.random.shuffle(j_range)
            np.random.shuffle(k_range)

            for i0 in i_range:
                for j0 in j_range:
                    for k0 in k_range:

                        for i in range(i0-1, i0+2):
                            i = wrap(i, 0, self.cells[0]-1)
                            for j in range(j0-1, j0+2):
                                j = wrap(j, 0, self.cells[1]-1)
                                for k in range(k0-1, k0+2):
                                    k = wrap(k, 0, self.cells[2]-1)
                                    if maxFlawRadius_rank[i][j][k] > maxFlawRadius_rank[i0][j0][k0]: 
                                        cpi = cpIndex[i][j][k]
                                        i00 = cpi[0]
                                        j00 = cpi[1]
                                        k00 = cpi[2]

                                        d = max(math.ceil(maxFlawRadius_rank[i][j][k]/cell_size), 1)
                                        
                                        di = i00-i0
                                        dj = j00-j0
                                        dk = k00-k0

                                        if self.periodic[0]:
                                          di = min([abs(di), abs(di-self.cells[0]), abs(di+self.cells[0])])
                                        
                                        if self.periodic[1]:
                                          dj = min([abs(dj), abs(dj-self.cells[1]), abs(dj+self.cells[1])])
                                        
                                        if self.periodic[2]:
                                          dk = min([abs(dk), abs(dk-self.cells[2]), abs(dk+self.cells[2])])

                                        if di*di + dj*dj + dk*dk < d*d and maxFlawRadius_rank[i0][j0][k0] > 0:
                                            maxFlawRadius_rank[i0][j0][k0] = maxFlawRadius_rank[i][j][k] 
                                            cpIndex[i0][j0][k0] = cpIndex[i][j][k]

            # Between iterations need to check cells on outer edge against neighbors
            for r in range(self.num_ranks_used):
                flaw_buffer = np.empty(self.cells, dtype=float) 
                cp_buffer = np.empty(np.append(self.cells, 2))
                if r == g_rank:
                    flaw_buffer = np.copy(maxFlawRadius_rank)
                    cp_buffer = np.copy(cpIndex)

                flaw_buffer = self.comm.bcast(flaw_buffer, root=r)
                cp_buffer = self.comm.bcast(cp_buffer, root=r)

                if r in neighbors:
                    neighbor_coord = self.rankIndexToCoord(r)
                    dr = neighbor_coord - self.rank_coord
                    for dd in range(self.dim):
                        # Bit of a hack for neighbors across periodic boundaries
                        # since neighbors must be within one cell distance
                        if dr[dd] > 1:
                            dr[dd] = -1

                        if dr[dd] < -1:
                            dr[dd] = 1

                    # Merge field and cpIndex (should be a way to condense the following)
                    # X lower boundary
                    xx = wrap(self.cell_range_x[0]-1,0,self.cells[0]-1)
                    if np.all(dr == np.array([-1, 0, 0])):
                        for v in j_range:
                            for w in k_range:
                                if flaw_buffer[xx][v][w] > maxFlawRadius_rank[xx][v][w]:
                                    maxFlawRadius_rank[xx][v][w] = flaw_buffer[xx][v][w]
                                    cpIndex[xx][v][w] = cp_buffer[xx][v][w]

                    # X upper boundary
                    xx = wrap(self.cell_range_x[1],0,self.cells[0]-1)
                    if np.all(dr == np.array([1, 0, 0])):
                        for v in j_range:
                            for w in k_range:
                                 if flaw_buffer[xx][v][w] > maxFlawRadius_rank[xx][v][w]:
                                    maxFlawRadius_rank[xx][v][w] = flaw_buffer[xx][v][w]
                                    cpIndex[xx][v][w] = cp_buffer[xx][v][w]

                    # Y lower boundary
                    yy = wrap(self.cell_range_y[0]-1,0,self.cells[1]-1)
                    if np.all(dr == np.array([0, -1, 0])):
                        for u in i_range:
                            for w in k_range:
                                 if flaw_buffer[u][yy][w] > maxFlawRadius_rank[u][yy][w]:
                                    maxFlawRadius_rank[u][yy][w] = flaw_buffer[u][yy][w]
                                    cpIndex[u][yy][w] = cp_buffer[u][yy][w]

                    # Y upper boundary
                    yy = wrap(self.cell_range_y[1],0,self.cells[1]-1)
                    if np.all(dr == np.array([0, 1, 0])):
                        for u in i_range:
                            for w in k_range:
                                 if flaw_buffer[u][yy][w] > maxFlawRadius_rank[u][yy][w]:
                                    maxFlawRadius_rank[u][yy][w] = flaw_buffer[u][yy][w]
                                    cpIndex[u][yy][w] = cp_buffer[u][yy][w]

                    # Z lower boundary
                    zz = wrap(self.cell_range_z[0]-1,0,self.cells[1]-1)
                    if np.all(dr == np.array([0, 0, -1])):
                        for u in i_range:
                            for v in j_range:
                                 if flaw_buffer[u][v][zz] > maxFlawRadius_rank[u][v][zz]:
                                    maxFlawRadius_rank[u][v][zz] = flaw_buffer[u][v][zz]
                                    cpIndex[u][v][zz] = cp_buffer[u][v][zz]

                    # Z upper boundary
                    zz = wrap(self.cell_range_z[1],0,self.cells[1]-1)
                    if np.all(dr == np.array([0, 0, 1])):
                        for u in i_range:
                            for v in j_range:
                                 if flaw_buffer[u][v][zz] > maxFlawRadius_rank[u][v][zz]:
                                    maxFlawRadius_rank[u][v][zz] = flaw_buffer[u][v][zz]
                                    cpIndex[u][v][zz] = cp_buffer[u][v][zz]

                self.comm.Barrier()

        self.maxFlawRadius_complete = np.empty(self.cells, dtype=float)
        self.comm.Allreduce([maxFlawRadius_rank, MPI.DOUBLE], [self.maxFlawRadius_complete, MPI.DOUBLE], op=MPI.MAX)


    def syncField(self):
        self.maxFlawRadius = np.empty(self.cells, dtype=float)
        if g_rank == 0:
            self.maxFlawRadius = self.maxFlawRadius_complete     

        self.maxFlawRadius = g_comm_world.bcast(self.maxFlawRadius,root=0) # does this need to return for field to be updated?
        

    def write2File(self):
        # Save max flaw radius field to file
        with open(self.maxFlawRadius_file, 'w') as f:
            f.write(str(self.cells[0]) + " " + str(self.cells[1]) + " " + str(self.cells[2]) + "\n")
            for i in range(self.cells[0]):
                for j in range(self.cells[1]):
                    for k in range(self.cells[2]): 
                        f.write(str(i) + " " + str(j) + " " + str(k) + " " + str(self.maxFlawRadius[i][j][k]) + "\n")


    def readFromFile(self):
        num_lines = countFileLines(self.maxFlawRadius_file)-1 # First line gives the number of cells
        line_range = np.ceil(num_lines*np.array([g_rank, g_rank+1])/g_num_ranks)
        file_pos = fileOffsetFromLineNum(self.maxFlawRadius_file, line_range[0]+1)

        count = 0
        with open(self.maxFlawRadius_file, 'r') as f:
            first_line = f.readline()
            self.cells = [int(t) for t in first_line.split()]            

            #Preallocating array sizes speeds things up
            self.pts = np.empty((np.prod(self.cells),self.dim), dtype=float)
            self.maxFlawRadius = np.empty(np.prod(self.cells), dtype=float)

            f.seek(file_pos)
            for i in range(int(line_range[0]), int(line_range[1])):
                line = f.readline()
                if i % 10000 == 0:
                    print('MaxFlaw: (' + str(i) + '/' + str(num_lines) + ')...')
                    log2file('MaxFlaw: (' + str(i) + '/' + str(num_lines) + ')...')

                terms = line.split()

                self.maxFlawRadius[i] = float(terms[self.dim]) 
                self.pts[i] = np.array([[int(t) for t in terms[:self.dim]]])

        # Sync max flaw radius field between processes
        for r in range(g_num_ranks):
            print("Syncing max flaw radius field on rank " + str(r) + "...")
            log2file("Syncing max flaw radius field on rank " + str(r) + "...")
            data_range = np.ceil(num_lines*np.array([r, r+1])/g_num_ranks)

            flaw_buffer = None
            pts_buffer = None
            if r == g_rank:
                flaw_buffer = np.copy(self.maxFlawRadius)
                pts_buffer = np.copy(self.pts)

            flaw_buffer = g_comm_world.bcast(flaw_buffer, root=r)
            pts_buffer = g_comm_world.bcast(pts_buffer, root=r)

            self.maxFlawRadius[int(data_range[0]):int(data_range[1])] = flaw_buffer[int(data_range[0]):int(data_range[1])]
            self.pts[int(data_range[0]):int(data_range[1]),:] = pts_buffer[int(data_range[0]):int(data_range[1]),:]

        keep = self.maxFlawRadius > 0 
        log2file(str(keep.shape) + " " + str(self.maxFlawRadius.shape) + " " + str(self.pts.shape))
        self.maxFlawRadius = self.maxFlawRadius[keep]
        self.pts = self.pts[keep,:]
        self.kdt = KDTree(self.pts, leaf_size=int(np.ceil(len(self.pts)/2)), metric='euclidean')


    def getFlawSize(self, pt):
        # Grid position of particle
        gx = math.floor(self.cells[0]*pt[0])
        gy = math.floor(self.cells[1]*pt[1])
        gz = math.floor(self.cells[2]*pt[2])
        _, index = self.kdt.query(np.array([gx,gy,gz]).reshape(1,-1),k=1)
        return self.maxFlawRadius[index[0][0]]


class maxFlawRadius2D:
    def __init__(self,
    maxFlawRadius_file,
    simgrid,
    grid,
    x0,
    x1,
    pores,
    periodic = [False, False],
    readFromFile=False):
        dim = 2
        self.dim = dim
        self.x0 = np.asarray(x0[:self.dim]) # [x0, y0] should define box bigger than subObject
        self.x1 = np.asarray(x1[:self.dim]) # [x1, y1] should define box bigger than subObject
        self.dx = self.x1-self.x0 
        self.periodic = periodic[:self.dim]
        self.pores = pores
        self.simcells = np.asarray(simgrid[:self.dim])
        self.cells = np.asarray(grid[:self.dim])
        self.maxFlawRadius_file = maxFlawRadius_file

        print("Setting up max flaw radius field")
        if not readFromFile:
            self.ranks_per_side = math.floor(g_num_ranks**(1/self.dim))
            self.num_ranks_used = self.ranks_per_side**self.dim
            self.comm = g_comm_world.Create(g_comm_world.group.Incl([i for i in range(self.num_ranks_used)]))

            self.rank_coord = self.rankIndexToCoord(g_rank)

            # CC: Logic for conversion to maxflawradius object that handles both 2d and 3d
            # self.ranks_in_dir = self.ranks_per_side*np.ones(self.dim)
            # self.simcell_range_low = mp.ceil(self.simcells*self.rank_coord/self.ranks_in_dir)
            # self.simcell_range_high = np.ceil(self.simcells*(self.rank_coord+np.ones(self.dim))/self.ranks_in_dir)
            # self.cell_range_low = np.ceil(self.cells*self.rank_coord/self.ranks_in_dir)
            # self.cell_range_high = np.ceil(self.cells*(self.rank_coord+np.ones(self.dim))/self.ranks_in_dir)

            self.ranks_in_x = self.ranks_per_side
            self.ranks_in_y = self.ranks_per_side

            self.simcell_range_x = [math.ceil(self.simcells[0]*self.rank_coord[0]/self.ranks_in_x), math.ceil(self.simcells[0]*(self.rank_coord[0]+1)/self.ranks_in_x)]
            self.simcell_range_y = [math.ceil(self.simcells[1]*self.rank_coord[1]/self.ranks_in_y), math.ceil(self.simcells[1]*(self.rank_coord[1]+1)/self.ranks_in_y)]

            self.cell_range_x = [math.ceil(self.cells[0]*self.rank_coord[0]/self.ranks_in_x), math.ceil(self.cells[0]*(self.rank_coord[0]+1)/self.ranks_in_x)]
            self.cell_range_y = [math.ceil(self.cells[1]*self.rank_coord[1]/self.ranks_in_y), math.ceil(self.cells[1]*(self.rank_coord[1]+1)/self.ranks_in_y)]
            
            start_time = time.perf_counter()
            if g_rank < self.num_ranks_used:
                print("Creating distance field")
                self.createDistanceField()
                print("Finding critical points")
                self.findCriticalPoints()
                print("Creating max flaw radius field")
                # self.createMaxFlawRadiusField()
                self.watershed_cf()

            g_comm_world.Barrier()
            log2file("Syncing max flaw radius field across ranks")
            self.syncField()
            g_comm_world.Barrier()
            
            end_time = time.perf_counter()
            print(f"Computation time was {(end_time-start_time)/60.0:.3f} min")
            log2file(f"Computation time was {(end_time-start_time)/60.0:.3f} min")
            
            if g_rank == 0:
                print("Writing max flaw radius to file")
                self.write2File()

            print("Max flaw radius field generation complete")
            log2file("Finished creating max flaw radius field")
        else:
            if g_rank == 0:
                print("Reading max flaw radius from file " + self.maxFlawRadius_file)
            log2file("Reading max flaw radius from file " + self.maxFlawRadius_file)
            
            self.readFromFile()
            
            if g_rank == 0:
                print("Max flaw radius read from " + maxFlawRadius_file)
            log2file("Max flaw radius read from " + maxFlawRadius_file)

    # For periodic boundaries list of pores should include their images across those boundaries as well
    def createDistanceField(self):
        df_rank = np.full(self.cells, np.inf)

        # Loop over every cell and check distance to nearest surface
        for i in range(self.cell_range_x[0], self.cell_range_x[1]):
          for j in range(self.cell_range_y[0], self.cell_range_y[1]):
              pos = self.dx*(np.array([i,j])+0.5*np.ones(self.dim))/self.cells+self.x0

              dist2pore = [None] * len(self.pores)
              for pi in range(0, len(self.pores)):
                  p = self.pores[pi]
                  dx = (p[1] - pos[0])
                  dy = (p[2] - pos[1])
                  dist2pore[pi] = math.sqrt(dx*dx + dy*dy) - p[0]

              for dd in range(self.dim):
                  if not self.periodic[dd]:
                    dist2pore.append(abs(pos[dd] - self.x0[dd]))
                    dist2pore.append(abs(self.x1[dd] - pos[dd]))

              df_rank[i][j] = min(dist2pore)

        self.df = np.empty(self.cells, dtype=float)
        self.comm.Allreduce([df_rank, MPI.DOUBLE], [self.df, MPI.DOUBLE], op=MPI.MIN)


    def findCriticalPoints(self):
        #Find critical points of distance field
        crit_pts_rank = np.full(self.cells, True)
        crit_data = np.empty((0, self.dim+1), dtype=float)

        for i in range(self.cell_range_x[0], self.cell_range_x[1]) if self.periodic[0] else range(max(1, self.cell_range_x[0]), min(self.cells[0]-2, self.cell_range_x[1])):
          for j in range(self.cell_range_y[0], self.cell_range_y[1]) if self.periodic[1] else range(max(1, self.cell_range_y[0]), min(self.cells[1]-2, self.cell_range_y[1])):
              if (self.df[i][j] >= 0) or \
                 (self.df[i][j] > self.df[wrap(i-1, 0, self.cells[0]-1)][j] and self.df[i][j] > self.df[wrap(i+1, 0, self.cells[0]-1)][j]) or \
                 (self.df[i][j] > self.df[i][wrap(j-1, 0, self.cells[1]-1)] and self.df[i][j] > self.df[i][wrap(j+1, 0, self.cells[1]-1)]) or \
                 (self.df[i][j] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)] and self.df[i][j] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)]) or \
                 (self.df[i][j] > self.df[wrap(i+1, 0, self.cells[0]-1)][wrap(j+1, 0, self.cells[1]-1)] and self.df[i][j] > self.df[wrap(i-1, 0, self.cells[0]-1)][wrap(j-1, 0, self.cells[1]-1)]):
                crit_pts_rank[i][j] = True
                crit_data = np.append(crit_data, [[i, j, self.df[i][j]]], axis=0)

        # # Syncs cp data asynchronously (doesn't work in 3D but does in 2D, why?) # CC debug
        # # Async send
        # for r in range(self.num_ranks_used):
        #   if r  == g_rank:
        #     continue

        #   self.comm.Isend([np.array(crit_data.shape[0]), MPI.INT], dest=r, tag=0)
        #   self.comm.Isend([crit_data, MPI.DOUBLE], dest=r, tag=1)

        # # Recv
        # self.crit_pts = np.copy(crit_data)
        # for r in range(self.num_ranks_used):
        #   if r == g_rank:
        #     continue

        #   num_pts = np.empty([1,1], dtype=int)
        #   self.comm.Recv([num_pts, MPI.INT], source=r, tag=0)

        #   recv_pts = np.empty([num_pts[0][0],self.dim+1], dtype=float)
        #   self.comm.Recv([recv_pts, MPI.DOUBLE], source=r, tag=1)
          
        #   self.crit_pts = np.append(self.crit_pts, recv_pts, axis=0)

        # Performings syncs serially
        self.crit_pts = np.empty((0, self.dim+1)) # CC debug
        for r in range(self.num_ranks_used):
            cp_buffer = None
            if r == g_rank:
                cp_buffer = np.copy(crit_data)

            cp_buffer = self.comm.bcast(cp_buffer, root=r)
            self.crit_pts = np.append(self.crit_pts, cp_buffer, axis=0)

        self.crit_pts = self.crit_pts[np.argsort(-self.crit_pts[:, self.dim])]


    def createMaxFlawRadiusField(self):
        # Create max flaw radius field
        maxFlawRadius_rank = np.copy(self.df)
        self.maxFlawRadius_complete = np.empty(self.cells, dtype=float) # copy.deepcopy(self.df)

        # Compute grid search domain for each critical point once based on flaw size
        ds = self.crit_pts[:, self.dim].reshape(self.crit_pts.shape[0],1)
        deltas = np.ceil(np.divide(np.tile(ds,(1,2)), np.tile(np.divide(self.dx, self.cells), (ds.shape[0], 1) ) ))
        
        for i in range(self.cell_range_x[0], self.cell_range_x[1]):
          print("Processing row ", i, " of ", self.cells[0], "...")
          for j in range(self.cell_range_y[0], self.cell_range_y[1]):
              if maxFlawRadius_rank[wrap(i, 0, self.cells[0]-1)][wrap(j,0,self.cells[1]-1)] < 0:
                  continue
              
              for m, c in enumerate(self.crit_pts):  
                d = c[2]
                delta = deltas[m]
                if delta[0] == 0 or delta[1] == 0:
                  continue

                dgx = c[0]-i
                dgy = c[1]-j
                
                if self.periodic[0]:
                  dgx = min([abs(dgx), abs(dgx-self.cells[0]), abs(dgx+self.cells[0])])
                
                if self.periodic[1]:
                  dgy = min([abs(dgy), abs(dgy-self.cells[1]), abs(dgy+self.cells[1])])

                if dgx*dgx + dgy*dgy < delta[0]*delta[0]:
                    maxFlawRadius_rank[i][j] = d
                    break

        self.comm.Allreduce([maxFlawRadius_rank, MPI.DOUBLE], [self.maxFlawRadius_complete, MPI.DOUBLE], op=MPI.MAX)


    def rankIndexToCoord(self, index):
        i = wrap(index, 0, self.num_ranks_used-1)

        coord = np.array([i % self.ranks_per_side,
                          math.floor(i / self.ranks_per_side)])

        return coord


    def rankCoordToIndex(self, coord):
        for dd in range(self.dim):
            if self.periodic[dd]:
                coord[dd] = wrap(coord[dd], 0, self.ranks_per_side-1)

        return coord[0] + self.ranks_per_side*coord[1]


    def getRankNeighbors(self):
         # Get rank indices of neighbors
        neighbors = []
        for k in range(self.dim):
            for f in [-1, 1]:
                offset = np.array([0,0])
                offset[k] = f
                neighbor_coord = offset + self.rank_coord
                if self.periodic[k]:
                    neighbor_coord = wrap(neighbor_coord, np.zeros(self.dim), (self.ranks_per_side-1)*np.ones(self.dim))   

                # Can't communicate with neighbor that is outside grid
                if np.any(neighbor_coord < 0):
                    continue

                neighbors.append(self.rankCoordToIndex(neighbor_coord))

        return neighbors

        
    def watershed_cf(self):
        neighbors = self.getRankNeighbors()

        print("Propagating critical points using water shed algorithm...")
        log2file("Propagating critical points using water shed algorithm...")
        cell_size = (self.dx[0]/self.cells[0])

        max_df = 0
        maxFlawRadius_rank = np.copy(self.df)
        for cp in self.crit_pts:
            if cp[2] > max_df:
                max_df = cp[2]
        max_n = int(max(math.ceil(max_df/cell_size), 1))

        cpIndex = np.empty(np.append(self.cells, self.dim))
        for i in range(0,self.cells[0]):
            for j in range(0,self.cells[1]):
                cpIndex[i][j] = np.array([i, j])

        i_range = np.arange(self.cell_range_x[0], self.cell_range_x[1], 1)
        j_range = np.arange(self.cell_range_y[0], self.cell_range_y[1], 1)

        for n in range(0,max_n+1):
            print("\tWatershed Iter: ", n, " of ", max_n)
            log2file("\tWatershed Iter: " + str(n) + " of " + str(max_n))
            np.random.shuffle(i_range)
            np.random.shuffle(j_range)

            for i0 in i_range:
                for j0 in j_range:
                    for i in range(i0-1, i0+2):
                        i = wrap(i, 0, self.cells[0]-1)
                        for j in range(j0-1, j0+2):
                            j = wrap(j, 0, self.cells[1]-1)
                            if maxFlawRadius_rank[i][j] > maxFlawRadius_rank[i0][j0]: 
                                cpi = cpIndex[i][j]
                                i00 = cpi[0]
                                j00 = cpi[1]

                                d = max(math.ceil(maxFlawRadius_rank[i][j]/cell_size), 1)
                                
                                di = i00-i0
                                dj = j00-j0

                                if self.periodic[0]:
                                  di = min([abs(di), abs(di-self.cells[0]), abs(di+self.cells[0])])
                                
                                if self.periodic[1]:
                                  dj = min([abs(dj), abs(dj-self.cells[1]), abs(dj+self.cells[1])])
                                
                                if di*di + dj*dj < d*d and maxFlawRadius_rank[i0][j0] > 0:
                                    maxFlawRadius_rank[i0][j0] = maxFlawRadius_rank[i][j] 
                                    cpIndex[i0][j0] = cpIndex[i][j]

            # Between iterations need to check cells on outer edge against neighbors
            for r in range(self.num_ranks_used):
                flaw_buffer = np.empty(self.cells, dtype=float) 
                cp_buffer = np.empty(np.append(self.cells, 2))
                if r == g_rank:
                    flaw_buffer = np.copy(maxFlawRadius_rank)
                    cp_buffer = np.copy(cpIndex)

                flaw_buffer = self.comm.bcast(flaw_buffer, root=r)
                cp_buffer = self.comm.bcast(cp_buffer, root=r)

                if r in neighbors:
                    neighbor_coord = self.rankIndexToCoord(r)
                    dr = neighbor_coord - self.rank_coord
                    for dd in range(self.dim):
                        # Bit of a hack for neighbors across periodic boundaries
                        # since neighbors must be within one cell distance
                        if dr[dd] > 1:
                            dr[dd] = -1

                        if dr[dd] < -1:
                            dr[dd] = 1
                    
                    # Merge field and cpIndex (should be a way to condense the following)
                    # X lower boundary
                    xx = wrap(self.cell_range_x[0]-1,0,self.cells[0]-1)
                    if np.all(dr == np.array([-1, 0])):
                        for u in j_range:
                            if flaw_buffer[xx][u] > maxFlawRadius_rank[xx][u]:
                                maxFlawRadius_rank[xx][u] = flaw_buffer[xx][u]
                                cpIndex[xx][u] = cp_buffer[xx][u]

                    # X upper boundary
                    xx = wrap(self.cell_range_x[1],0,self.cells[0]-1)
                    if np.all(dr == np.array([1, 0])):
                        for u in j_range:
                             if flaw_buffer[xx][u] > maxFlawRadius_rank[xx][u]:
                                maxFlawRadius_rank[xx][u] = flaw_buffer[xx][u]
                                cpIndex[xx][u] = cp_buffer[xx][u]

                    # Y lower boundary
                    yy = wrap(self.cell_range_y[0]-1,0,self.cells[1]-1)
                    if np.all(dr == np.array([0, -1])):
                        for v in i_range:
                             if flaw_buffer[v][yy] > maxFlawRadius_rank[v][yy]:
                                maxFlawRadius_rank[v][yy] = flaw_buffer[v][yy]
                                cpIndex[v][yy] = cp_buffer[v][yy]

                    # Y upper boundary
                    yy = wrap(self.cell_range_y[1],0,self.cells[1]-1)
                    if np.all(dr == np.array([0, 1])):
                        for v in i_range:
                             if flaw_buffer[v][yy] > maxFlawRadius_rank[v][yy]:
                                maxFlawRadius_rank[v][yy] = flaw_buffer[v][yy]
                                cpIndex[v][yy] = cp_buffer[v][yy]

                self.comm.Barrier()

        self.maxFlawRadius_complete = np.empty(self.cells, dtype=float)
        self.comm.Allreduce([maxFlawRadius_rank, MPI.DOUBLE], [self.maxFlawRadius_complete, MPI.DOUBLE], op=MPI.MAX)
    

    def syncField(self):
        self.maxFlawRadius = np.empty(self.cells, dtype=float) 
        if g_rank == 0:
            self.maxFlawRadius = self.maxFlawRadius_complete

        g_comm_world.bcast(self.maxFlawRadius,root=0)
        

    def write2File(self):
        # Save max flaw radius field to file
        with open(self.maxFlawRadius_file, 'w') as f:
            f.write(str(self.cells[0]) + " " + str(self.cells[1]) + "\n")
            for i in range(self.cells[0]):
                for j in range(self.cells[1]):
                        f.write(str(i) + " " + str(j) + " " + str(self.maxFlawRadius[i][j]) + "\n")

    def readFromFile(self):
        num_lines = countFileLines(self.maxFlawRadius_file)-1 # First line gives the number of cells
        line_range = np.ceil(num_lines*np.array([g_rank, g_rank+1])/g_num_ranks)
        file_pos = fileOffsetFromLineNum(self.maxFlawRadius_file, line_range[0]+1)

        count = 0
        with open(self.maxFlawRadius_file, 'r') as f:
            first_line = f.readline()
            self.cells = [int(t) for t in first_line.split()]            

            #Preallocating array sizes speeds things up
            self.pts = np.empty((np.prod(self.cells),self.dim), dtype=float)
            self.maxFlawRadius = np.empty(np.prod(self.cells), dtype=float)

            f.seek(file_pos)
            for i in range(int(line_range[0]), int(line_range[1])):
                line = f.readline()
                if i % 10000 == 0:
                    print('MaxFlaw: (' + str(i) + '/' + str(num_lines) + ')...')
                    log2file('MaxFlaw: (' + str(i) + '/' + str(num_lines) + ')...')

                terms = line.split()

                self.maxFlawRadius[i] = float(terms[self.dim]) 
                self.pts[i] = np.array([[int(t) for t in terms[:self.dim]]])

        # Sync max flaw radius field between processes
        for r in range(g_num_ranks):
            print("Syncing max flaw radius field on rank " + str(r) + "...")
            log2file("Syncing max flaw radius field on rank " + str(r) + "...")
            data_range = np.ceil(num_lines*np.array([r, r+1])/g_num_ranks)

            flaw_buffer = None
            pts_buffer = None
            if r == g_rank:
                flaw_buffer = np.copy(self.maxFlawRadius)
                pts_buffer = np.copy(self.pts)

            flaw_buffer = g_comm_world.bcast(flaw_buffer, root=r)
            pts_buffer = g_comm_world.bcast(pts_buffer, root=r)

            self.maxFlawRadius[int(data_range[0]):int(data_range[1])] = flaw_buffer[int(data_range[0]):int(data_range[1])]
            self.pts[int(data_range[0]):int(data_range[1]),:] = pts_buffer[int(data_range[0]):int(data_range[1]),:]

        keep = self.maxFlawRadius > 0 
        log2file(str(keep.shape) + " " + str(self.maxFlawRadius.shape) + " " + str(self.pts.shape))
        self.maxFlawRadius = self.maxFlawRadius[keep]
        self.pts = self.pts[keep,:]
        self.kdt = KDTree(self.pts, leaf_size=int(np.ceil(len(self.pts)/2)), metric='euclidean')


    def getFlawSize(self, pt):
        # Grid position of particle
        gx = math.floor(self.cells[0]*pt[0])
        gy = math.floor(self.cells[1]*pt[1])
        _, index = self.kdt.query(np.array([gx, gy]).reshape(1,-1),k=1)
        return self.maxFlawRadius[index[0][0]]


# Grid object used for efficiently performing random sequential absorption and poisson disc sampling
class Grid:
    def __init__(self,
    cell_size,
    **kwargs):
        self.cell_size = cell_size
        self.dim = kwargs.get('dim', 3)
        self.x0 = np.asarray(kwargs.get('x0', np.zeros(self.dim)))
        self.dx = np.asarray(kwargs.get('dx', np.ones(self.dim)))
        self.periodic = kwargs.get('periodic', [True for k in range(self.dim)])
        self.num_cells = np.ceil(self.dx / self.cell_size).astype(int)
        self.has_fractional_cells = np.multiply(self.num_cells, self.dx) > self.dx
        self.cells = -np.ones(tuple(self.num_cells)).astype(int)  # Cell value of -1 indicates empty cell

    def outside(self, x):
        return outside(x, self.x0, self.dx, self.periodic)

    def pos2cell(self, x):
        return np.floor((x - self.x0) / self.cell_size).astype(int)

    def place(self, cell, x):
        self.cells[tuple(cell)] = x
        return self

    def at(self, cell):
        return self.cells[tuple(cell)]

    def occupied(self, cell):
        return self.cells[tuple(cell)] >= 0

    def cells_in_range(self, cell, window):
        search_min = cell - window
        search_max = cell + window

        # If direction doesn't have periodic boundaries then cap search window to min and max cells
        # If grid has fractional cells and candidate cell is within range of them add 1 to search direction
        search_min = (search_min - np.ones(self.dim).astype(int) * np.logical_and(self.has_fractional_cells, cell < 1)) * self.periodic + np.logical_not(self.periodic) * np.maximum(
            np.zeros(self.dim).astype(int), search_min)

        search_max = (search_max + np.ones(self.dim).astype(int) * np.logical_and(self.has_fractional_cells, cell >= self.num_cells - 3)) * self.periodic + np.logical_not(self.periodic) * np.minimum(
            self.num_cells-1, search_max)

        search_range = search_max-search_min
        larger_than_cells = search_range > ( self.num_cells - 1 )

        search_min = np.zeros(self.dim).astype(int) * larger_than_cells + np.logical_not(larger_than_cells) * search_min
        search_max = (self.num_cells - 1) * larger_than_cells + np.logical_not(larger_than_cells) * search_max

        # Create list of all cells indices in search area
        cells_to_search_ranges = [np.array([]) for k in range(self.dim)]
        for k in range(self.dim):
            cells_to_search_ranges[k] = np.array(list(range(search_min[k], search_max[k] + 1))).astype(int)
        cells_to_search = combvec(cells_to_search_ranges).astype(int)


        cell_vals = []
        for c in cells_to_search:
            # Check if c index needs to be wrapped in periodic direction
            c = wrap(c, np.zeros(self.dim), self.num_cells - 1)

            val = self.cells[tuple(c)]
            if val >= 0:
                cell_vals.append(val)

        return cell_vals


# wraps integer index i to range of min and max (both inclusive)
# accepts nd arrays
def wrap(i, lower, upper):
    # Make sure all inputs are np arrays
    i = np.array(i).astype(int)
    lower = np.array(lower).astype(int)
    upper = np.array(upper).astype(int)

    r = upper - lower + 1
    i = i - lower

    i = i * (i < r) + (i >= r) * (lower + np.mod(i, r))
    i = i * (i >= 0) + (i < 0) * (upper - np.mod(np.absolute(i + 1), r))

    return i.astype(int)


# Generates all possible combinations of elements from list of offsets for each direction
def combvec(offsets):
    # dim = len(offsets)

    # combs = np.reshape(offsets[0], (len(offsets[0]), 1))
    # for i in range(1, dim):
    #     new_combs = np.empty((0, i+1))
    #     for j in offsets[i]:
    #         for k in combs:
    #             new_combs = np.append(new_combs, np.reshape(np.append(k, j), (1, len(k)+1)), axis=0)

    #     combs = new_combs

    # return combs

    # Appears to be a faster method
    return np.array(list(itertools.product(*offsets)))


# Randomly samples points on surface of unit n-sphere
def random_direction(dim: int=3):
    direction = np.random.normal(size=dim)
    return direction / np.linalg.norm(direction)


# Adds images of pores for nearest neighbors
def add_images(pores, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)
    x0 = np.asarray(kwargs.get('x0', np.zeros(dim)))
    dx = np.asarray(kwargs.get('dx', np.ones(dim)))
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    num_pores = len(pores)
    boundary_pores = np.empty((0, dim + 1))
    offsets = [np.array([0, -1, 1]) for i in range(dim)]
    offset_combs = combvec(offsets)
    for ff in offset_combs:
        for p in pores:
            boundary_pores = np.append(boundary_pores, np.reshape(p + np.append(np.multiply(ff, dx), [0]), (1, dim + 1)), axis=0)

    return boundary_pores


# Adds pores images of those that lie on the boundary of the RVE
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# pores = list of pores, format = x, y, z, r
def add_pore_images(pores, **kwargs):
  # Setup parameters
  dim = kwargs.get('dim', 3)
  x0 = np.asarray(kwargs.get('x0', np.zeros(dim)))
  dx = np.asarray(kwargs.get('dx', np.ones(dim)))
  periodic = kwargs.get('periodic', [True for k in range(dim)])

  num_pores = len(pores)

  indices = []
  boundary_pores = np.empty((0, dim + 1))
  for i in range(num_pores):
      offsets = [np.array([0]) for k in range(dim)]

      for k in range(dim):
          if periodic[k]:
              # Pore on - boundary
              if pores[i, k] - pores[i, dim] < x0[k]:
                  offsets[k] = np.append(offsets[k], 1.0)

              # Pore on + boundary
              if pores[i, k] + pores[i, dim] > x0[k] + dx[k]:
                  offsets[k] = np.append(offsets[k], -1.0)

      offset_combs = combvec(offsets)
      for ff in offset_combs:
          indices.append(i)
          boundary_pores = np.append(boundary_pores, np.reshape(pores[i, :] + np.append(np.multiply(ff, dx), [0]), (1, dim+1)), axis=0)

  return boundary_pores, indices


# Imports pores from a LAMMPS dump file
# pack_file = location of pack.dump file
# offset = translation of coordinates for simulation (default [0.0, 0.0, 0.0])
# add_boundary_images = boolean whether to include pores that lie on boundaries on other periodic surface 
# (e.g for generating the foam meshes in the particleFileWriter)
def import_pores(pack_file, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)

    # x0 = np.asarray(kwargs.get('x0', np.zeros(dim)))
    # dx = np.asarray(kwargs.get('dx', np.ones(dim)))
    offset = np.asarray(kwargs.get('offset', np.zeros(dim)))
    scale = kwargs.get('scale', 1.0)
    add_boundary_images = kwargs.get('add_boundary_images', True)
    radius_scale = kwargs.get('radius_scale', 1.0)

    print("Reading pores from " + pack_file)

    file = open(pack_file, "r")
    lines = file.readlines()
    file.close()

    # Stores current data block being processed and list of voids
    curr_block = ""
    pores = []
    box = []
    curr_iter = 0
    num_atoms = 0
    count = 0
    for line in lines:
        count += 1
        # print("Line{}: {}".format(count, line.strip()))
        if line == "" or line == "\n":
            continue

        # Split string by spaces
        line_terms = line.split()

        # Determined what to do based on first entry in line (ITEM: = some header info, number
        if line_terms[0] == "ITEM:":
            curr_block = line_terms[1]

            if curr_block == "BOX":
                # For box get periodicity from same line
                periodic = [line_terms[3] == 'pp', line_terms[4] == 'pp', line_terms[5] == 'pp']
        elif curr_block == "ITERATION":
            curr_iter = int(line)
        elif curr_block == "NUMBER":
            num_atoms = int(line)
        elif curr_block == "BOX":
            box.append([float(i) for i in line_terms])
        elif curr_block == "ATOMS":
            atom = [float(x) for x in line_terms]
            # Output format radius, x, y, z
            pores.append([scale*(atom[2] + offset[0]), scale*(atom[3] + offset[1]), scale*(atom[4] + offset[2]), scale * radius_scale * atom[5]])

    pores = np.asarray(pores)  # convert python array to numpy array

    x0 = [scale*(box[0][0] + offset[0]), scale*(box[1][0] + offset[1]), scale*(box[2][0] + offset[2])]
    dx = [scale*box[0][1], scale*box[1][1], scale*box[2][1]]

    print("Box after offset and scaling:")
    print("\t" + str(x0[0]) + " " + str(dx[0]))
    print("\t" + str(x0[1]) + " " + str(dx[1]))
    print("\t" + str(x0[2]) + " " + str(dx[2]))
    print("Periodicity: " + str(periodic[0]) + " " + str(periodic[1]) + " " + str(periodic[2]))
    print("Pores from file: " + str(len(pores)))
    if add_boundary_images:
        output_pores = add_pore_images(pores, x0=x0, dx=dx, dim=dim, periodic=periodic)
        print("Pore with images: " + str(len(output_pores)))
    else:
        output_pores = pores

    # Reformt pores so radius is first, then x, y, z
    radii = output_pores[:,dim].copy()
    output_pores[:,1:dim+1] = output_pores[:,0:dim] # shift x, y, z right one column
    output_pores[:,0] = radii

    return np.array(output_pores)


# # Generates densely packed random points in box domain with minimum separation
# # Generalized for 2 and 3 dimensions
# def poisson(spacing, **kwargs):
#     # Setup parameters
#     seed = kwargs.get('seed', np.random.randint(0, 1e8))
#     dim = kwargs.get('dim', 3)
#     x0 = np.asarray(kwargs.get('x0', np.zeros(dim)))
#     dx = np.asarray(kwargs.get('dx', np.ones(dim)))
#     trials = kwargs.get('trials', 30)
#     periodic = kwargs.get('periodic', [True for k in range(dim)])

#     x0 = x0[:dim]
#     dx = dx[:dim]

#     # Grid cell size if determined given longest cell size diagonal is radius of spacing
#     grid = Grid(spacing / math.sqrt(dim), dim=dim, x0=x0, dx=dx, periodic=periodic)

#     # First seed point is at center of domain
#     seed_pts = np.array([dx/2 + x0])  # Check that this makes the right size array
#     pores = np.empty((0, dim+1))

#     # Set seed for random number generation
#     np.random.seed(seed)

#     # Terminates when there are no more seed points left
#     while seed_pts.shape[0] != 0:
#         # Get the first point in the seed array
#         curr_seed_pt = seed_pts[0, :]
#         seed_pts = np.delete(seed_pts, 0, axis=0)

#         num_rejected = 0
#         while num_rejected < trials:
#             # Pick random distance and random angle from seed point
#             d = spacing*(1 + np.random.random() * 0.1)  # random spacing between seed and candidate point
#             direction = random_direction(dim)
#             candidate_pt = curr_seed_pt + d * direction

#             accept = True  # Accept candidate point unless found otherwise

#             # If periodic boundaries are not enabled test point must be inside domain
#             is_outside, candidate_pt = grid.outside(candidate_pt)
#             if is_outside:
#                 num_rejected = num_rejected + 1
#                 continue

#             candidate_cell = grid.pos2cell(candidate_pt)
#             if grid.occupied(candidate_cell):
#                 accept = False
#             else:
#                 neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))

#                 # print( "Neighbors", neighbor_pores)
#                 # for n in neighbor_pores:
#                 #     print(n, pores[n,:-1], grid.pos2cell(pores[n,:-1]))
#                 # print("Done Neighbors")

#                 # Check every cell for occupant, if any are within spacing reject candidate point
#                 for n in neighbor_pores:

#                     # Check distance between candidate point and point already there
#                     ddx = candidate_pt - pores[n, :-1]

#                     # Use shortest image distance if that point lies across periodic boundary
#                     ddxx = np.logical_not(periodic) * ddx + periodic * np.minimum(np.absolute(ddx + dx),
#                                                                                   np.minimum(np.absolute(ddx),
#                                                                                              np.absolute(ddx - dx)))
#                     # print(candidate_pt, pores[n, :-1], ddxx, np.sqrt(np.sum(ddxx**2)), np.sqrt(np.sum(ddxx**2)) < spacing)

#                     if np.sum(ddxx**2) < spacing**2:
#                         accept = False
#                         break  # Once a conflicting point is found no need to search the remaining cells

#             if accept:
#                 # Add pore to list of pores and seed points plus mark grid cell as occupied
#                 grid = grid.place(tuple(candidate_cell), pores.shape[0])
#                 new_pore = np.append(candidate_pt, spacing/2)
#                 pores = np.append(pores, np.reshape(new_pore, (1, dim+1)), axis=0)
#                 seed_pts = np.append(seed_pts, np.reshape(candidate_pt, (1, dim)), axis=0)
#             else:
#                 num_rejected = num_rejected + 1

#     return pores

# Generates densely packed random points in box domain with minimum separation
# Generalized for 2 and 3 dimensions
def poisson(spacing, **kwargs):
    seed = kwargs["seed"] if "seed" in kwargs else np.random.randint(0, 1e8)
    dim = int(kwargs.get("dim", 3))
    x0 = np.asarray(kwargs.get("x0", np.zeros(dim)), dtype=float)[:dim]
    dx = np.asarray(kwargs.get("dx", np.ones(dim)), dtype=float)[:dim]
    trials = int(kwargs.get("trials", 30))
    periodic = np.asarray(kwargs.get("periodic", [True for _ in range(dim)]), dtype=bool)[:dim]

    if dim not in (2, 3):
        raise ValueError("poisson() currently supports dim=2 or dim=3")

    if spacing <= 0.0:
        raise ValueError("spacing must be positive")

    if np.any(dx <= 0.0):
        raise ValueError("all dx entries must be positive")

    # Cell diagonal is about spacing.
    cell_size = spacing / math.sqrt(dim)
    inv_cell_size = 1.0 / cell_size
    num_cells = np.ceil(dx * inv_cell_size).astype(np.int64)

    cells = -np.ones(tuple(num_cells), dtype=np.int64)

    spacing2 = spacing * spacing
    pore_radius = 0.5 * spacing

    # Local RNG: reproducible without resetting global np.random state.
    rng = np.random.default_rng(seed)
    twopi = 2.0 * math.pi

    # Python lists avoid O(N^2) np.append / np.delete behavior.
    points = []
    active = [tuple((x0 + 0.5 * dx).tolist())]   # preserve old behavior: center is only an active seed
    active_head = 0

    # ------------------------------------------------------------------
    # 2D specialized path
    # ------------------------------------------------------------------
    if dim == 2:
        x00, y00 = float(x0[0]), float(x0[1])
        lx, ly = float(dx[0]), float(dx[1])
        xmax, ymax = x00 + lx, y00 + ly
        nx, ny = int(num_cells[0]), int(num_cells[1])
        perx, pery = bool(periodic[0]), bool(periodic[1])

        while active_head < len(active):
            sx, sy = active[active_head]
            active_head += 1

            num_rejected = 0
            while num_rejected < trials:
                d = spacing * (1.0 + 0.1 * rng.random())

                theta = twopi * rng.random()
                cx = sx + d * math.cos(theta)
                cy = sy + d * math.sin(theta)

                # Boundary handling.
                if perx:
                    cx = ((cx - x00) % lx) + x00
                elif cx < x00 or cx > xmax:
                    num_rejected += 1
                    continue

                if pery:
                    cy = ((cy - y00) % ly) + y00
                elif cy < y00 or cy > ymax:
                    num_rejected += 1
                    continue

                ic = int((cx - x00) * inv_cell_size)
                jc = int((cy - y00) * inv_cell_size)

                # Guard against exact upper-bound / floating point edge cases.
                if ic < 0:
                    ic = 0
                elif ic >= nx:
                    ic = nx - 1

                if jc < 0:
                    jc = 0
                elif jc >= ny:
                    jc = ny - 1

                accept = cells[ic, jc] < 0

                if accept:
                    # Need a slightly wider search near periodic boundaries when
                    # the domain length is not an integer multiple of cell_size.
                    if perx:
                        ilo = ic - 2 - (1 if ic < 2 else 0)
                        ihi = ic + 2 + (1 if ic >= nx - 3 else 0)
                    else:
                        ilo = max(0, ic - 2)
                        ihi = min(nx - 1, ic + 2)

                    if pery:
                        jlo = jc - 2 - (1 if jc < 2 else 0)
                        jhi = jc + 2 + (1 if jc >= ny - 3 else 0)
                    else:
                        jlo = max(0, jc - 2)
                        jhi = min(ny - 1, jc + 2)

                    for ii0 in range(ilo, ihi + 1):
                        ii = ii0 % nx if perx else ii0

                        for jj0 in range(jlo, jhi + 1):
                            jj = jj0 % ny if pery else jj0

                            n = cells[ii, jj]
                            if n >= 0:
                                px, py = points[n]

                                ddx = cx - px
                                ddy = cy - py

                                if perx:
                                    ddx = abs(ddx)
                                    ddx = min(ddx, lx - ddx)

                                if pery:
                                    ddy = abs(ddy)
                                    ddy = min(ddy, ly - ddy)

                                if ddx * ddx + ddy * ddy < spacing2:
                                    accept = False
                                    break

                        if not accept:
                            break

                if accept:
                    n = len(points)
                    cells[ic, jc] = n

                    p = (cx, cy)
                    points.append(p)
                    active.append(p)
                else:
                    num_rejected += 1

    # ------------------------------------------------------------------
    # 3D specialized path
    # ------------------------------------------------------------------
    else:
        x00, y00, z00 = float(x0[0]), float(x0[1]), float(x0[2])
        lx, ly, lz = float(dx[0]), float(dx[1]), float(dx[2])
        xmax, ymax, zmax = x00 + lx, y00 + ly, z00 + lz
        nx, ny, nz = int(num_cells[0]), int(num_cells[1]), int(num_cells[2])
        perx, pery, perz = bool(periodic[0]), bool(periodic[1]), bool(periodic[2])

        while active_head < len(active):
            sx, sy, sz = active[active_head]
            active_head += 1

            num_rejected = 0
            while num_rejected < trials:
                d = spacing * (1.0 + 0.1 * rng.random())

                # Uniform random direction on the sphere, faster than normalizing
                # three Gaussian samples for every trial.
                zdir = 2.0 * rng.random() - 1.0
                theta = twopi * rng.random()
                rxy = math.sqrt(max(0.0, 1.0 - zdir * zdir))

                cx = sx + d * rxy * math.cos(theta)
                cy = sy + d * rxy * math.sin(theta)
                cz = sz + d * zdir

                # Boundary handling.
                if perx:
                    cx = ((cx - x00) % lx) + x00
                elif cx < x00 or cx > xmax:
                    num_rejected += 1
                    continue

                if pery:
                    cy = ((cy - y00) % ly) + y00
                elif cy < y00 or cy > ymax:
                    num_rejected += 1
                    continue

                if perz:
                    cz = ((cz - z00) % lz) + z00
                elif cz < z00 or cz > zmax:
                    num_rejected += 1
                    continue

                ic = int((cx - x00) * inv_cell_size)
                jc = int((cy - y00) * inv_cell_size)
                kc = int((cz - z00) * inv_cell_size)

                # Guard against exact upper-bound / floating point edge cases.
                if ic < 0:
                    ic = 0
                elif ic >= nx:
                    ic = nx - 1

                if jc < 0:
                    jc = 0
                elif jc >= ny:
                    jc = ny - 1

                if kc < 0:
                    kc = 0
                elif kc >= nz:
                    kc = nz - 1

                accept = cells[ic, jc, kc] < 0

                if accept:
                    # Need a slightly wider search near periodic boundaries when
                    # the domain length is not an integer multiple of cell_size.
                    if perx:
                        ilo = ic - 2 - (1 if ic < 2 else 0)
                        ihi = ic + 2 + (1 if ic >= nx - 3 else 0)
                    else:
                        ilo = max(0, ic - 2)
                        ihi = min(nx - 1, ic + 2)

                    if pery:
                        jlo = jc - 2 - (1 if jc < 2 else 0)
                        jhi = jc + 2 + (1 if jc >= ny - 3 else 0)
                    else:
                        jlo = max(0, jc - 2)
                        jhi = min(ny - 1, jc + 2)

                    if perz:
                        klo = kc - 2 - (1 if kc < 2 else 0)
                        khi = kc + 2 + (1 if kc >= nz - 3 else 0)
                    else:
                        klo = max(0, kc - 2)
                        khi = min(nz - 1, kc + 2)

                    for ii0 in range(ilo, ihi + 1):
                        ii = ii0 % nx if perx else ii0

                        for jj0 in range(jlo, jhi + 1):
                            jj = jj0 % ny if pery else jj0

                            for kk0 in range(klo, khi + 1):
                                kk = kk0 % nz if perz else kk0

                                n = cells[ii, jj, kk]
                                if n >= 0:
                                    px, py, pz = points[n]

                                    ddx = cx - px
                                    ddy = cy - py
                                    ddz = cz - pz

                                    if perx:
                                        ddx = abs(ddx)
                                        ddx = min(ddx, lx - ddx)

                                    if pery:
                                        ddy = abs(ddy)
                                        ddy = min(ddy, ly - ddy)

                                    if perz:
                                        ddz = abs(ddz)
                                        ddz = min(ddz, lz - ddz)

                                    if ddx * ddx + ddy * ddy + ddz * ddz < spacing2:
                                        accept = False
                                        break

                            if not accept:
                                break

                        if not accept:
                            break

                if accept:
                    n = len(points)
                    cells[ic, jc, kc] = n

                    p = (cx, cy, cz)
                    points.append(p)
                    active.append(p)
                else:
                    num_rejected += 1

    pores = np.empty((len(points), dim + 1), dtype=float)
    if points:
        pores[:, :dim] = np.asarray(points, dtype=float)
        pores[:, dim] = pore_radius

    return pores

# Provides fast generation of nonuniform poisson disc sampling in 2 dimensions
# An arbitrary density profile can be defined via a function or lambda with bounds [0,1]
def fast_nonuniform_poisson_2(**kwargs):
    c = kwargs.get('c', [0.0, 0.0])
    dx = kwargs.get('dx', [1.0, 1.0])
    min_dist = kwargs.get('min_dist', 0.1)
    max_dist = kwargs.get('max_dist', min_dist)
    num_trials = kwargs.get('num_trials', 30)
    density = kwargs.get('density', lambda x, y: 1)
    seed = kwargs.get('seed', random.randint(0, 1000000))
    
    random.seed(seed)

    width = dx[0]
    height = dx[1]

    cell_size = min_dist / math.sqrt(2)
    search_size = int(math.ceil(max_dist / cell_size))

    ngx = int(width / cell_size) + 1  # of grid cells in x-direction
    ngy = int(height / cell_size) + 1  # of grid cells in y-direction
    grid = [[None for _ in range(ngx)] for _ in range(ngy)]

    # x0 = random.uniform(0, width)
    # y0 = random.uniform(0, height)
    x0 = random.uniform((width-min_dist)/2, (width+min_dist)/2)
    y0 = random.uniform((height-min_dist)/2, (height-min_dist)/2)

    samples = [[x0, y0]]
    seed_list = [[x0, y0]]
    sample_dists = [min_dist + (max_dist - min_dist) * (1 - density(x0 / width, y0 / height))]

    # Start with a random point
    gx = int(x0 / cell_size)
    gy = int(y0 / cell_size)
    grid[gx][gy] = 0

    while seed_list:
        idx = random.randint(0, len(seed_list) - 1)
        x, y = seed_list[idx]
        local_dist = sample_dists[idx]
        found = False
        for _ in range(num_trials):
            angle = random.uniform(0, 2 * math.pi)
            r = random.uniform(local_dist, 2 * local_dist)
            nx = x + r * math.cos(angle)
            ny = y + r * math.sin(angle)

            if not (0 <= nx < width and 0 <= ny < height):
                continue

            n_local_dist = min_dist + (max_dist - min_dist) * (1 - density(nx / width, ny / height))
            gx = int(nx / cell_size)
            gy = int(ny / cell_size)
            accept = True
            for i in range(max(0, gx - search_size), min(ngx, gx + search_size + 1)):
                for j in range(max(0, gy - search_size), min(ngy, gy + search_size + 1)):
                    n = grid[i][j]
                    if n is not None:
                        neighbor = samples[n]

                        dx = neighbor[0] - nx
                        dy = neighbor[1] - ny
                        distSqr = dx * dx + dy * dy

                        min_allowed = min(sample_dists[n], n_local_dist)
                        if distSqr < min_allowed * min_allowed:
                            accept = False
                            break
                if not accept:
                    break
            if accept:
                grid[gx][gy] = len(samples)
                samples.append([nx, ny])
                seed_list.append([nx, ny])
                sample_dists.append(n_local_dist)

                found = True
                break
        if not found:
            seed_list.pop(idx)

    # Offset points
    for s in samples:
        s[0] = s[0] + c[0]
        s[1] = s[1] + c[1]

    return samples

# Provides fast generation of nonuniform poisson disc sampling in 3 dimensions
# An arbitrary density profile can be defined via a function or lambda with bounds [0,1]
def fast_nonuniform_poisson_3(**kwargs):
    c = kwargs.get('c', [0.0, 0.0, 0.0])
    dx = kwargs.get('dx', [1.0, 1.0, 1.0])
    min_dist = kwargs.get('min_dist', 0.1)
    max_dist = kwargs.get('max_dist', min_dist)
    num_trials = kwargs.get('num_trials', 30)
    density = kwargs.get('density', lambda x, y, z: 1)
    seed = kwargs.get('seed', random.randint(0, 1000000))

    random.seed(seed)

    width = dx[0]
    height = dx[1]
    length = dx[2]

    cell_size = min_dist / math.sqrt(2)
    search_size = int(math.ceil(max_dist / cell_size))

    ngx = int(width / cell_size) + 1  # of grid cells in x-direction
    ngy = int(height / cell_size) + 1  # of grid cells in y-direction
    ngz = int(length / cell_size) + 1  # of grid cells in y-direction
    grid = [[[None for _ in range(ngx)] for _ in range(ngy)] for _ in range(ngz)]

    # x0 = random.uniform(0, width)
    # y0 = random.uniform(0, height)
    x0 = random.uniform((width-min_dist)/2, (width+min_dist)/2)
    y0 = random.uniform((height-min_dist)/2, (height-min_dist)/2)
    z0 = random.uniform((length - min_dist) / 2, (length - min_dist) / 2)

    samples = [[x0, y0, z0]]
    seed_list = [[x0, y0, z0]]
    sample_dists = [min_dist + (max_dist - min_dist) * (1 - density(x0 / width, y0 / height, z0 / length))]

    # Start with a random point
    gx = int(x0 / cell_size)
    gy = int(y0 / cell_size)
    gz = int(z0 / cell_size)
    grid[gx][gy][gz] = 0

    while seed_list:
        idx = random.randint(0, len(seed_list) - 1)
        x, y, z = seed_list[idx]
        local_dist = sample_dists[idx]
        found = False
        for _ in range(num_trials):
            d = random_direction(3)
            r = random.uniform(local_dist, 2 * local_dist)
            nx = x + r * d[0]
            ny = y + r * d[1]
            nz = z + r * d[2]

            if not (0 <= nx < width and 0 <= ny < height and 0 <= nz <= length):
                continue

            n_local_dist = min_dist + (max_dist - min_dist) * (1 - density(nx / width, ny / height, nz / length))
            gx = int(nx / cell_size)
            gy = int(ny / cell_size)
            gz = int(nz / cell_size)
            accept = True
            for i in range(max(0, gx - search_size), min(ngx, gx + search_size + 1)):
                for j in range(max(0, gy - search_size), min(ngy, gy + search_size + 1)):
                    for k in range(max(0, gz - search_size), min(ngz, gz + search_size + 1)):
                        n = grid[i][j][k]
                        if n is not None:
                            neighbor = samples[n]

                            dx = neighbor[0] - nx
                            dy = neighbor[1] - ny
                            dz = neighbor[2] - nz
                            distSqr = dx * dx + dy * dy + dz * dz

                            min_allowed = min(sample_dists[n], n_local_dist)
                            if distSqr < min_allowed * min_allowed:
                                accept = False
                                break
                    if not accept:
                        break
                if not accept:
                    break
            if accept:
                grid[gx][gy][gz] = len(samples)
                samples.append([nx, ny, nz])
                seed_list.append([nx, ny, nz])
                sample_dists.append(n_local_dist)

                found = True
                break
        if not found:
            seed_list.pop(idx)

    # Offset points
    for s in samples:
        s[0] = s[0] + c[0]
        s[1] = s[1] + c[1]
        s[2] = s[2] + c[2]

    return samples


# Creates pores which some minimum spacing via random sequential adsorption (RSA)
# N = number of desired particles
# spacing = minimum separation between spheres
# r = function alias for sampling radii (To be implemented, right now only monodisperse)
# seed = random seed for RSA
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# trials = number of rejected trials before terminating RSA
def rsa(N, spacing, **kwargs):
    # Setup parameters
    seed = kwargs.get('seed', np.random.randint(0, 1e8))
    dim = kwargs.get('dim', 3)
    x0 = np.asarray(kwargs.get('x0', np.zeros(dim)))
    dx = np.asarray(kwargs.get('dx', np.ones(dim)))
    trials = kwargs.get('trials', 30)
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    # Grid cell size if determined given longest cell size diagonal is radius of spacing
    grid = Grid(spacing / math.sqrt(dim), dim=dim)

    pores = np.empty((0, dim + 1))

    num_rejected = 0
    np.random.seed(seed)
    while pores.shape[0] < N and num_rejected <= trials:
        x = np.empty(dim + 1)
        x[dim] = spacing/2
        for k in range(dim):
            if periodic[k]:
                x[k] = np.random.uniform(0.0, dx[k])
            else:
                x[k] = np.random.uniform(0.0+spacing/2, dx[k]-spacing/2)

        candidate_cell = grid.pos2cell(x[:-1])

        accept = True
        if grid.occupied(candidate_cell):
            accept = False
        else:
            neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))
            for p in neighbor_pores:
                d = x[:dim] - pores[p, 0:-1]
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(d), np.abs(d + dx))
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(dd), np.abs(d - dx))

                if np.sum(np.square(dd)) < spacing ** 2:
                    accept = False
                    break

        if accept:
            grid = grid.place(candidate_cell, pores.shape[0])
            pores = np.append(pores, np.reshape(x, (1, dim + 1)), axis=0)
            num_rejected = 0
        else:
            num_rejected = num_rejected + 1

    for i in range(len(pores)):
        pores[i, 0:dim] = pores[i, 0:dim] + x0

    return pores


#############################################
# def voronoi_volumes(points, v0):
#   # computes volumes for a voronoi tesselation of given points
#   # or returns v0 for edge points without enough data to compute.
#     v = Voronoi(points)
#     vol = np.zeros(v.npoints)
#     for i, reg_num in enumerate(v.point_region):
#         indices = v.regions[reg_num]
#         if -1 in indices: # some regions can be opened
#             vol[i] = v0
#         else:
#             vol[i] = ConvexHull(v.vertices[indices]).volume
#     return vol


# ===========================================
# END UTILITY FUNCTIONS
# ===========================================

# ===========================================
# GEOMETRY OBJECTS
# ===========================================

#############################################
class Geometry:
  """
  Base class of PFW geometry objects
  """
  @abstractmethod
  def __init__(self,
               name,
               vel=_defaultVelocity,
               mat: int=_defaultMat,
               group: int=_defaultGroup,
               particleType=_defaultParticleType,
               damage=_defaultDamage,
               porosity=_defaultPorosity,
               temperature=_defaultTemperature,
               surfaceTraction=_defaultSurfaceTraction):
    # # Input checks
    # type_check_scalar(name, "Object name", str)
    # type_check_scalar(mat, "Material", int)
    # type_check_scalar(group, "Group", int)
    # type_check_scalar(particleType, "Particle type", int)
    # type_check_scalar(damage, "Damage", float)
    # type_check_scalar(porosity, "Porosity", float)
    # type_check_scalar(temperature, "Temperature", float)

    # type_check_array(surfaceTraction, "Surface traction", 3, floatType, canBeFunc=True)
    # type_check_array(vel, "Velocity", (3, 2), floatType, canBeFunc=True)

    self.name = name
    self.vel = vel if callable(vel) else np.asarray(vel) # Only need to check if callable if converting it array
    self.mat = mat
    self.group = group
    self.damage = damage
    self.particleType = particleType
    self.porosity = porosity
    self.temperature = temperature
    self.surfaceTraction = surfaceTraction if callable(surfaceTraction) else np.asarray(surfaceTraction)
  
  @abstractmethod
  def isInterior(self, pt):
    pass

  @abstractmethod
  def getSurfaceNormal(self, pt):
    return _defaultSurfaceNormal

  @abstractmethod
  def getSurfacePosition(self, pt):
    return _defaultSurfacePosition

  @abstractmethod
  def xMin(self):
    return -np.inf

  @abstractmethod
  def xMax(self):
    return np.inf

  @abstractmethod
  def getSubregions(self):
    return [(self.mat, self.particleType)]


#############################################
class shell(Geometry):
  """
  Geometry object for creating a spherical shell defined by center and inner and outer radii
  """
  def __init__(self,
               name,
               x0,
               ri,
               ro,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    

    self.x0 = np.asarray(x0)
    self.ri = ri
    self.ro = ro

  def isInterior(self, pt, skinDepth):
    x = np.asarray(pt) - self.x0
    #  Is point interior?
    if ( self.ri * self.ri <= np.dot(x,x) < self.ro * self.ro ):
      #  Is point a surface?
      if ( ( np.dot(x,x) > (self.ro - skinDepth) * (self.ro - skinDepth) ) or ( np.dot(x,x) < (self.ri + skinDepth) * (self.ri + skinDepth) ) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1 

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt) - self.x0
    xn = np.linalg.norm(x)

    n = x / xn 
    if np:
      n = -n
    return n

  def getSurfacePosition(self,pt):
    x = np.asarray(pt) - self.x0
    xn = np.linalg.norm(x)
    n = x / xn 
    di = xn - self.ri
    do = self.ro - xn

    if di < do:
      return di*-n
    else:
      return do*n

  def xMin(self):
    return self.x0[0] - self.ro

  def xMax(self):
    return self.x0[0] + self.ro


#############################################
class sphere(Geometry):
  """
  Geometry object for creating a sphere defined by center and radius
  """
  def __init__(self,
               name,
               x0,
               r,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)
    self.r = r
    self.rSqr = r*r

  def isInterior(self,pt,skinDepth):
    x = np.asarray(pt) - self.x0
    if np.dot(x,x) < self.rSqr:
      if np.dot(x,x) > (self.r - skinDepth)*(self.r - skinDepth):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1    
  
  def getSurfaceNormal(self,pt):
    x = np.asarray(pt) - self.x0
    return x / np.linalg.norm( x )
  
  def getSurfacePosition(self,pt):
    x = np.asarray(pt) - self.x0
    mag = np.linalg.norm(x)
    return (self.r/mag-1)*x

  def xMin(self):
    return self.x0[0] - self.r

  def xMax(self):
    return self.x0[0] + self.r


#############################################
class sphericalInclusion(Geometry):
  """
  Geometry object for creating a sphere defined by center and radius, but where surface is always false (useful for inclusions)
  """
  def __init__(self,
               name,
               x0,
               r,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)
    self.r = r
    self.rSqr = r*r

  def isInterior(self,pt):
    x = np.asarray(pt) - self.x0
    if npt.dot(x, x) < self.rSqr:
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self,pt):
    return np.array([0.0,0.0,0.0])

  def xMin(self):
    return self.x0[0] - self.r

  def xMax(self):
    return self.x0[0] + self.r


#############################################
class czSphericalPrill(Geometry):
    """
    Geometry object for creating a spherical prill with voronoi crystals bound by cohesive zones and defined by center, radius, and grain size

    Applies to 2D and 3D prills with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 center,
                 radius,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 vel=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 explicitGroups=False,
                 dim: int=3):
        super().__init__(name,
                         vel = vel,
                         mat = mat,
                         group = group,
                         particleType = particleType)
        self.dim = dim
        self.center = np.asarray(center)  # center
        self.center = self.center[:self.dim] # For 2D problem only need x and y
        self.radius = radius  # prill radius inluding shell
        self.radiusSquared = radius * radius
        self.x0 = self.center - self.radius
        self.x1 = self.center + self.radius
        self.grainDiameter = grainDiameter
        self.prill_porosity = porosity
        self.radialBias = radialBias
        self.seed = seed
        self.explicitGroups = explicitGroups

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.x0, dx=self.x1-self.x0, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        if self.explicitGroups:
            self.group = []
            for i in range(self.npts):
                self.group.append(i)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        
        self.matDirs = []
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)
            # Discourage porosity near the prill surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)
           
            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

    def isInterior(self, pt, skinDepth):
        # is the point within the object
        x = np.asarray(pt[:self.dim])
        xx = x - self.center

        # Check if point is inside sphere
        if np.dot(xx, xx) < self.radiusSquared:
            # Find voronoi cell closest to point
            _, index = self.kdt.query(x.reshape(1, -1), k=1)
            index = index[0,0]
            # If voronoi cell is not porosity
            if self.phase[index] != 0:
              surfaceFlag = 0 # Particle is interior unless otherwise determined

              # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
              minSurfaceDist = np.inf
              ridgePts = self.voronoi.ridge_points
              for i in range(len(ridgePts)):
                  p1 = ridgePts[i][0]
                  p2 = ridgePts[i][1]
                  if index == p1 or index == p2:
                      if index == p1:
                          p = p2
                      else:
                          p = p1

                      n = self.vpts[p, :] - self.vpts[index, :]
                      n = n / 2
                      d = np.linalg.norm(n)
                      n = n / d

                      dv = x - self.vpts[index, :]
                      dvc = np.dot(n, dv)  # component of points along voronoi face normal
                      if dvc > 0.0 and dvc > d - skinDepth:
                          surfaceFlag = 3
                          break

              # If within skinDepth of prill surface must be a surface
              xx = x - self.center
              if np.dot(xx,xx) > (self.radius - skinDepth)*(self.radius - skinDepth):
                if surfaceFlag != 3:
                  surfaceFlag = _defaultSurfaceFlag

              return surfaceFlag
        
        return -1
       

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfaceNormal = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        # Compute distance to surface of prill and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfaceNormal = xx / xx_mag

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, 0.0)

        return surfaceNormal

    def getSurfacePosition(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfacePosition = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfacePosition = surfaceDistance * n

        # Compute distance to surface of prill and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfacePosition = surfaceDistance * xx / xx_mag

        if self.dim == 2:
          surfacePosition = np.append(surfacePosition, 0.0)

        return surfacePosition

    def getGroup(self, pt):
        if self.explicitGroups:
          x = np.asarray(pt[:self.dim])
          _, index = self.kdt.query(x.reshape(1, -1), k=1)
          return self.group[index[0,0]]
        else:
          return  self.group

    def getMatDir(self, pt):
        x = np.asarray(pt[:self.dim])
        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        return self.matDirs[index[0,0]]

    def xMin(self):
        return self.x0[0]

    def xMax(self):
        return self.x1[0]


#############################################
class explicitBinderSphericalPrill(Geometry):
    """
    Geometry object for creating a spherical prill with voronoi crystals bound by cohesive zones and defined by center, radius, and grain size

    Applies to 2D prills with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 center,
                 radius,
                 grainDiameter,
                 binderThickness,
                 shellThickness,
                 porosity,
                 radialBias,
                 density,
                 seed,
                 vel=_defaultVelocity,
                 grainMat=_defaultMat,
                 binderMat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim: int=3):
        super().__init__(name,
                         vel = vel,
                         mat = grainMat,
                         group = group,
                         particleType = particleType)
        self.dim = dim
        self.center = np.asarray(center)  # center
        self.center = self.center[:self.dim] # For 2D problem only need x and y
        self.radius = radius  # prill radius inluding shell
        self.radiusSquared = radius * radius
        self.x0 = self.center - self.radius
        self.x1 = self.center + self.radius
        self.grainDiameter = grainDiameter
        self.binderThickness = binderThickness
        self.shellThickness = shellThickness
        self.grainMat = grainMat
        self.binderMat = binderMat
        self.prill_porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = np.array(fast_nonuniform_poisson_2(c=self.x0, dx=self.x1-self.x0, min_dist=0.2*grainDiameter, max_dist=self.grainDiameter, seed=self.seed, density=density))
        self.vpts = self.vpts[:, 0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        # self.group = []
        self.matDirs = []
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)
            # Discourage porosity near the prill surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)
            # self.group.append(i)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

    def isInterior(self, pt, skinDepth):
        # is the point within the object
        x = np.asarray(pt[:self.dim])
        xx = x - self.center

        # Check if point is inside sphere
        if np.dot(xx, xx) < self.radiusSquared:
          # Find voronoi cell closest to point
          _, index = self.kdt.query(x.reshape(1, -1), k=1)
          index = index[0,0]

          # If voronoi cell is not porosity
          if self.phase[index] != 0:
            surfaceFlag = 0 # Particle is interior unless otherwise determined

            # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
            ridgePts = self.voronoi.ridge_points
            for i in range(len(ridgePts)):
                p1 = ridgePts[i][0]
                p2 = ridgePts[i][1]
                if index == p1 or index == p2:
                    if index == p1:
                        p = p2
                    else:
                        p = p1
                    
                    # If neighbor is a void
                    if self.phase[p] == 0:
                      n = self.vpts[p, :] - self.vpts[index, :]
                      n = n / 2
                      d = np.linalg.norm(n)
                      n = n / d

                      dv = x - self.vpts[index, :]
                      dvc = np.dot(n, dv)  # component of points along voronoi face normal
                      if dvc > 0.0 and dvc > d - skinDepth:
                          surfaceFlag = _defaultSurfaceFlag
                          break

            # If within skinDepth of prill surface must be a surface
            xx = x - self.center
            if np.dot(xx,xx) > (self.radius - skinDepth)*(self.radius - skinDepth):
              surfaceFlag = _defaultSurfaceFlag

            return surfaceFlag
        
        return -1
       
    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfaceNormal = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                if self.phase[p] == 0:
                  n = self.vpts[p, :] - self.vpts[index, :]
                  n = n / 2
                  d = np.linalg.norm(n)
                  n = n / d

                  dv = (x - self.vpts[index, :])
                  surfaceDistance = (d - np.dot(dv, n))
                  if minSurfaceDist > surfaceDistance:
                      minSurfaceDist = surfaceDistance
                      surfaceNormal = n

        # Compute distance to surface of prill and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfaceNormal = xx / xx_mag

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, 0.0)

        return surfaceNormal

    def getSurfacePosition(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfacePosition = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                
                if self.phase[p] == 0:
                  n = self.vpts[p, :] - self.vpts[index, :]
                  n = n / 2
                  d = np.linalg.norm(n)
                  n = n / d

                  dv = (x - self.vpts[index, :])
                  surfaceDistance = (d - np.dot(dv, n))
                  if minSurfaceDist > surfaceDistance:
                      minSurfaceDist = surfaceDistance
                      surfacePosition = surfaceDistance * n

        # Compute distance to surface of prill and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfacePosition = surfaceDistance * xx / xx_mag

        if self.dim == 2:
          surfacePosition = np.append(surfacePosition, 0.0)

        return surfacePosition

    def getGroup(self, pt):
      return self.group

    def getMat(self, pt):
      # is the point within the object
      x = np.asarray(pt[:self.dim])
      xx = x - self.center

      # Find voronoi cell closest to point
      _, index = self.kdt.query(x.reshape(1, -1), k=1)
      index = index[0,0]

      # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
      ridgePts = self.voronoi.ridge_points
      for i in range(len(ridgePts)):
          p1 = ridgePts[i][0]
          p2 = ridgePts[i][1]
          if index == p1 or index == p2:
              if index == p1:
                  p = p2
              else:
                  p = p1

              n = (self.vpts[p, :] - self.vpts[index, :])/2
              d = np.linalg.norm(n)
              n = n / d

              dv = x - self.vpts[index, :]
              dvc = np.dot(n, dv)  # component of points along voronoi face normal
              # print(x, n, d, dv, dvc, dvc > 0.0, dvc > d - self.binderThickness/2, self.binderThickness)
              if dvc > 0.0 and dvc > d - self.binderThickness/2:
                  return self.binderMat
      
      # If within skinDepth of prill surface must be a surface
      if np.dot(xx,xx) > (self.radius - self.shellThickness) * (self.radius - self.shellThickness):
        return self.binderMat

      return self.grainMat

    def getMatDir(self, pt):
        x = np.asarray(pt[:self.dim])
        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        return self.matDirs[index[0,0]]

    def xMin(self):
        return self.x0[0]

    def xMax(self):
        return self.x1[0]

    def getSubregions(self):
      return [(self.grainMat, self.particleType),
              (self.binderMat, self.particleType)]


#############################################
class czCylindricalPrill(Geometry):
    """
    Geometry object for creating a cylindrical prill with voronoi crystals bound by cohesive zones and defined by center, radius, and grain size

    Applies to 2D and 3D prills with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 x1,
                 x2,
                 radius,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 vel=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim: int=3):
        super().__init__(name,
                         vel = vel,
                         mat = mat,
                         group = group,
                         particleType = particleType)
        self.dim = dim
        self.x1 = np.asarray(x0)
        self.x2 = np.asarray(x1)
        self.a = self.x1-self.x0
        self.a = self.a/np.linalg.norm(self.a)
        self.center = 0.5*(self.x0+self.x1)  # center
        self.radius = radius  # prill radius inluding shell
        self.radiusSquared = radius * radius

        self.grainDiameter = grainDiameter
        self.porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.x0, dx=self.x1-self.x0, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        self.matDirs = []
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)
            # Discourage porosity near the prill surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

    def insideCylinder(self, pt, l):
      # THIS doesn't seem quite right, FIXME
      # Seems like it's a check of surface flag not a check of interiorness?
      x = np.asarray(pt) - self.x1
      z = np.dot(x, self.axis)                                   # z-coordinate of test point
      r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

      if (z<l or z>(self.h-l)) or (r > self.r-l):
        return True
      
      return False

    def isInterior(self, pt, skinDepth):
      # is the point within the object
      x = np.asarray(pt[:self.dim])

      # Check if point is inside sphere
      if self.insideCylinder(pt, 0.0):
          # Find voronoi cell closest to point
          _, index = self.kdt.query(x.reshape(1, -1), k=1)
          index = index[0,0]

          # If voronoi cell is not porosity
          if self.phase[index] != 0:
            surfaceFlag = 0 # Particle is interior unless otherwise determined

            # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
            minSurfaceDist = np.inf
            ridgePts = self.voronoi.ridge_points
            for i in range(len(ridgePts)):
                p1 = ridgePts[i][0]
                p2 = ridgePts[i][1]
                if index == p1 or index == p2:
                    if index == p1:
                        p = p2
                    else:
                        p = p1

                    n = self.vpts[p, :] - self.vpts[index, :]
                    n = n / 2
                    d = np.linalg.norm(n)
                    n = n / d

                    dv = x - self.vpts[index, :]
                    dvc = np.dot(n, dv)  # component of points along voronoi face normal
                    if dvc > 0.0 and dvc > d - skinDepth:
                        surfaceFlag = 3
                        break

            # If within skinDepth of prill surface must be a surface
            if isInterior(pt, skinDepth):
              if surfaceFlag != 3:
                surfaceFlag = _defaultSurfaceFlag

            return surfaceFlag
      
      return -1      

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfaceNormal = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        # Compute distance to surface of prill and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfaceNormal = xx / xx_mag

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, 0.0)

        return surfaceNormal

    def getSurfacePosition(self, pt):
      # assumes the point is interior and a surface
      x = np.asarray(pt[:self.dim])

      _, index = self.kdt.query(x.reshape(1, -1), k=1)
      index = index[0,0]

      minSurfaceDist = np.inf
      surfacePosition = np.inf*np.ones((1,self.dim))
      ridgePts = self.voronoi.ridge_points
      for i in range(len(ridgePts)):
          p1 = ridgePts[i][0]
          p2 = ridgePts[i][1]
          if index == p1 or index == p2:
              if index == p1:
                  p = p2
              else:
                  p = p1
              n = self.vpts[p, :] - self.vpts[index, :]
              n = n / 2
              d = np.linalg.norm(n)
              n = n / d

              dv = (x - self.vpts[index, :])
              surfaceDistance = (d - np.dot(dv, n))
              if minSurfaceDist > surfaceDistance:
                  minSurfaceDist = surfaceDistance
                  surfacePosition = surfaceDistance * n

      x = np.asarray(pt)-self.x1
      z = np.dot(x,self.axis)  # z-coordinate of test point
      xr = x - z*self.axis
      r = np.linalg.norm( xr )  # r coordinate of test point

      dist_from_wall = self.r-r
      dist_from_top = self.h-z
      dist_from_bot = z

      min_surf_dist = min(dist_from_wall, min(dist_from_top, dist_from_bot))

      if min_surf_dist == dist_from_top:
        return dist_from_top * self.axis
      
      if min_surf_dist == dist_from_bot:
        return dist_from_bot * -self.axis

      if min_surf_dist == dist_from_wall:
        return (self.r/r - 1) * xr

      # Compute distance to surface of prill and check if closer
      xx = x - self.center
      xx_mag = np.linalg.norm(xx)
      surfaceDistance = self.radius - xx_mag
      if minSurfaceDist > surfaceDistance:
          surfacePosition = surfaceDistance * xx / xx_mag

      if self.dim == 2:
        surfacePosition = np.append(surfacePosition, 0.0)

      return surfacePosition

    def getGroup(self, pt):
        return  self.group

    def getMatDir(self, pt):
        x = np.asarray(pt[:self.dim])
        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        return self.matDirs[index[0,0]]


#############################################
class czPrill(Geometry):
    """
    Geometry object for creating a prill (box) with voronoi crystals bound by cohesive zones and defined by minimum corner, maximum corner, and grain size

    Applies to 2D and 3D prills with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 xmin,
                 xmax,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 vel=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 bondedSurfaceFraction=1.0,
                 neatSurfaceFraction=0.0,
                 dim: int=3):
        super().__init__(name,
                         vel = vel,
                         mat = mat,
                         group = group,
                         particleType = particleType)
        self.dim = dim
        self.xmin = np.asarray(xmin[:self.dim])
        self.xmax = np.asarray(xmax[:self.dim])
        self.dx = self.xmax - self.xmin
        self.center = 0.5 * ( self.xmin + self.xmax )
        self.grainDiameter = grainDiameter
        self.porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        self.bondedSurfaceFraction = bondedSurfaceFraction
        self.neatSurfaceFraction = neatSurfaceFraction

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.xmin, dx=self.dx, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        self.matDirs = []
        radius = np.min((self.xmax-self.xmin)*0.5)
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)

            # Discourage porosity near the prill surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

        # This should define a surface flag for each cell face based on the specified fraction of bonded surfaces.
        self.ridgePtFlags = []
        self.neatPtFlags = []
        for r in self.voronoi.ridge_points:
          p0 = r[0]
          p1 = r[1]
          
          rFlag = 2 
          if self.phase[p0] != 0 and self.phase[p1] != 0:
            rFlag = (3 if ( np.random.random() <= self.bondedSurfaceFraction ) else 2 )
          self.ridgePtFlags.append(rFlag)

          if rFlag == 2:
            self.neatPtFlags.append(_defaultCZTag)
          else:
            self.neatPtFlags.append( 1 if ( np.random.random() <= self.bondedSurfaceFraction * self.neatSurfaceFraction ) else 2 )
        
        # midpoints = []
        # for p in self.voronoi.ridge_points:
        #   midpoints.append(0.5*(self.vpts[p[0]] + self.vpts[p[1]]))
        # midpoints = np.vstack(midpoints) 
        # self.interface_kdt = KDTree(midpoints, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')

    def isInterior(self, pt, skinDepth):
        # is the point within the object
        x = np.asarray(pt[:self.dim])

        # Check if point is inside bounding
        if inside_box(x, self.xmin, self.dx, [False for d in range(self.dim)]):
          # Find voronoi cell closest to point
          _, index = self.kdt.query(x.reshape(1, -1), k=1)
          index = index[0,0]
          
          # If voronoi cell is not porosity
          if self.phase[index] != 0:
            # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
            minSurfaceDist = np.inf
            ridgePts = self.voronoi.ridge_points
            for i in range(len(ridgePts)):
              p1 = ridgePts[i][0]
              p2 = ridgePts[i][1]
              if index == p1 or index == p2:
                if index == p1:
                  p = p2
                else:
                  p = p1

                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = x - self.vpts[index, :]
                dvc = np.dot(n, dv)  # component of points along voronoi face normal
                if dvc > 0.0 and dvc > d - skinDepth:
                        return self.ridgePtFlags[i]
            return 0

        return -1
       

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfaceNormal = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, np.array([0.0]))

        return surfaceNormal

    def getSurfacePosition(self, pt):
        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        surfacePosition = np.inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfacePosition = surfaceDistance * n

        if self.dim == 2:
          surfacePosition = np.append(surfacePosition, np.array([0.0]))

        return surfacePosition

    def getGroup(self, pt):
        return  self.group

    def getMatDir(self, pt):
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        return self.matDirs[index]

    def xMin(self):
        return self.xmin[0]

    def xMax(self):
        return self.xmax[0]

    def getCZTag(self, pt):
        # x = np.asarray(pt[:self.dim])
        # _, index = self.interface_kdt.query(x.reshape(1, -1), k=1)
        # index = index[0,0]
        # return self.neatPtFlags[index]

        # assumes the point is interior and a surface
        x = np.asarray(pt[:self.dim])

        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.inf
        czTag = _defaultCZTag
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    czTag = self.neatPtFlags[i]

        return czTag


#############################################
class ellipsoid(Geometry):
  """
  Geometry object for creating a grid aligned ellipsoid defined by center and three lengths
  """
  def __init__(self,
               name,
               x0,
               a,
               b,
               c,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)
    self.a = a
    self.b = b
    self.c = c
    self.tolerance = 1e-6

  def isInterior(self,pt,skinDepth):
    x = pt[0] - self.x0[0]
    y = pt[1] - self.x0[1]
    z = pt[2] - self.x0[2]

    if (x/self.a)*(x/self.a) + (y/self.b)*(y/self.b) + (z/self.c)*(z/self.c) < 1:
      a = max( skinDepth/100, self.a - skinDepth)
      b = max( skinDepth/100, self.b - skinDepth)
      c = max( skinDepth/100, self.c - skinDepth)

      if ( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) > 1 ):
        return _defaultSurfaceFlag
      
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    p = self.getSurfacePosition(pt)
    n = np.divide(p, np.power(np.array([self.a,self.b,self.c]),2))
    return n / np.linalg.norm(n)

  def getSurfacePosition(self,pt):
    x = np.asarray(pt) - self.x0
    L = max(self.a,max(self.b,self.c))
    p1 = np.array([0.0,0.0,0.0])
    p2 = x
    n = p2 / np.linalg.norm(p2)
    p3 = p2 + L*n

    # Bisection method (numerical method to find surface point on ellipsoid)
    max_iterations = 100
    for i in range(max_iterations):
      g2 = np.sum(np.divide(np.power(p2,2),np.power(np.array([self.a,self.b,self.c]),2)))
      if abs(g2 - 1.0) < self.tolerance:
        return p2 - x

      if g2 < 1:
        p1 = p2
      else:
        p3 = p2

      p2 = (p1+p3)/2
    
    return None # Should through error    

  def xMin(self):
    return -self.a

  def xMax(self):
    return self.a
  

#############################################
class crystal(Geometry):
  """
  Generates faceted crystal analog in 3D defined by a height and min and max face offsets
  """
  def __init__(self,
               name,
               center,
               axis,
               height,
               rmin,
               rmax,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.center = center
    self.axis = axis
    self.height = height
    self.rmin = rmin
    self.rmax = rmax

    n0 = np.asarray(self.axis)
    n0 = n0 / np.linalg.norm(n0)

    n1 = np.array([3/(5*np.sqrt(2)), 2*np.sqrt(2)/5, 1/np.sqrt(2)]) 
    if np.all(n0 == n1):
      n1 = np.array([3/(5*np.sqrt(2)), 1/np.sqrt(2), 2*np.sqrt(2)/5])

    n1 = np.dot(np.identity(3)-np.tensordot(n0, n0, axes=0), n1)
    n1 = n1 / np.linalg.norm(n1)
    n2 = np.cross(n0, n1)
    self.matDir = np.vstack((n0, n1, n2))

    faceNormals = []

    # Face set 1
    for theta in np.linspace(0.0, 5*math.pi/3, 6): #math.pi/3):
      nn = np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 2
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): # 2*math.pi/3):
      nn = n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 3
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): #2*math.pi/3):
      nn = -n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    faceNormals.append(n0)
    faceNormals.append(-n0)

    self.faceNormals = faceNormals

    self.ci = np.vstack((np.zeros((6,3)), 
                         np.tile(self.height/2*n0,(3,1)),
                         np.tile(-self.height/2*n0,(3,1)),
                         np.zeros((2,3))))
    self.ri = (self.rmax-self.rmin)*np.random.rand(14,1)+self.rmin
    self.ri[-2:] = self.ri[-2:] + self.height

  def isInterior(self, pt, skinDepth):
    pc = np.asarray(pt) - self.center
    surfaceFlag = 0 # Surface flag defaults to interior unless otherwise determined
    for i in range(len(self.faceNormals)):
      dist2Surf = np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if dist2Surf > self.ri[i]:
        return -1
      
      if dist2Surf > self.ri[i] - skinDepth:
        surfaceFlag = _defaultSurfaceFlag

    return surfaceFlag

  def getSurfaceNormal(self,pt):
    pc = np.asarray(pt) - self.center
    min_distance_2_surf = np.inf
    surface_normal = np.empty((1,3))
    for i in range(len(self.faceNormals)):
      m = self.ri[i] - np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.faceNormals[i]

    return surface_normal

  def getSurfacePosition(self,pt):
    pc = np.asarray(pt) - self.center
    min_distance_2_surf = np.inf
    surface_normal = np.empty((1,3))
    for i in range(len(self.faceNormals)):
      m = self.ri[i] - np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.faceNormals[i]

    return min_distance_2_surf * surface_normal

  def xMin(self):
    return -rmax

  def xMax(self):   
    return rmax


#############################################
class indentor(Geometry):
  """
  Geometry object for creating a grid-aligned indentor defined by angle and number of facets'

  note: alpha=65.3 deg for berkovich with nFaces=3
  note: alpha=68 deg for vickers with nFaces=4
  """
  def __init__(self,
               name,
               x0,
               nFaces,
               alpha,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)         # tip coordinate, assumes identor is aligned with z-axis with tip in -z direction
    self.nFaces = nFaces # number of facets
    self.alpha = np.radians(alpha)   # facet angle, input in degrees, stored in radians, 
    
    # create set of unit vectors orthogonal to each facet.
    theta = 2.*np.pi/nFaces # angle between normals about indentor axis
    self.n = []
    k=0.
    for i in range(0,nFaces):
      k=k+1.
      ny = -np.sin(self.alpha)
      r = np.cos(self.alpha)
      nx = r*np.cos(theta*(k-0.5))
      nz = r*np.sin(theta*(k-0.5))

      # nz = -np.sin(self.alpha)
      # r = np.cos(self.alpha)
      # nx = r*np.cos(theta*(k-0.5))
      # ny = r*np.sin(theta*(k-0.5))

      (self.n).append(np.array([nx,ny,nz]))

  def isInterior(self, pt, skinDepth):  
    x = np.asarray(pt)-self.x0
    surfaceFlag = 0
    for i in range(0,self.nFaces):
      dist2Surf = x.dot(self.n[i]) 
      if dist2Surf > 0.0:
        return -1

      if dist2Surf > -skinDepth:
        surfaceFlag = _defaultSurfaceFlag
    
    return surfaceFlag

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)-self.x0
    min_distance_2_surf = np.inf
    surface_normal = np.empty((1,3))
    for i in range(0,self.nFaces):
      m =  -x.dot(self.n[i])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.n[i]

    return surface_normal

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-self.x0
    min_distance_2_surf = np.inf
    surface_normal = np.empty((1,3))
    for i in range(0,self.nFaces):
      m =  -x.dot(self.n[i])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.n[i]

    return min_distance_2_surf * surface_normal


#############################################
class spherical_indenter(Geometry):
  """
  Geometry object for creating a grid-aligned spherical tip indentor defined by radius and angle
  """
  def __init__(self,
               name,
               x0,
               radius,
               alpha,
               vel=_defaultVelocity,
               mat: int=_defaultMat,
               group: int=_defaultGroup,
               particleType=_defaultParticleType,
               dim: int=3):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)          # tip coordinate (bottom of spherical tip), assumes identor is aligned with y-axis with tip in -y direction
    self.radius = radius
    self.rSqr = radius * radius
    self.alpha = np.radians(alpha)   # tip half angle, input in degrees, stored in radians
    self.dim = dim

    self.g = self.x0[1] + self.radius * (1.0 - np.sin(self.alpha) )
    self.vtip = self.g - np.cos(self.alpha) * self.radius / np.tan(self.alpha) # height of virtual tip of cone shape
    
  def isInterior(self, pt, skinDepth):  
    x = np.asarray(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    surfaceFlag = 0
    if x[1] >= self.g:
      x[1] = x[1] - self.vtip
      rr = x[0] * x[0] + x[2] * x[2]
      if rr < (x[1] * np.tan(self.alpha))**2:
        if rr >= (x[1] * np.tan(self.alpha) - skinDepth * np.cos(self.alpha))**2:
          return _defaultSurfaceFlag
        else:
          return 0
    else:
      x[1] = x[1] - self.radius - self.x0[1]
      if np.dot(x,x) < self.rSqr:
        if np.dot(x,x) > (self.radius - skinDepth) * (self.radius - skinDepth):
          return _defaultSurfaceFlag
        else:
          return 0
    
    return -1    

  def getSurfaceNormal(self, pt):
    x = np.asarray(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      surfaceNormal = n / np.linalg.norm(n)
    else:
      x[1] = x[1] - self.radius - self.x0[1]
      surfaceNormal = x / np.linalg.norm(x)
    
    return surfaceNormal

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      n = n / np.linalg.norm(n)
      x[1] =  x[1] - self.vtip
      surfacePosition = -np.dot(n, x) * n
    else:
      x[1] = x[1] - ( self.radius + self.x0[1] )
      xmag = np.linalg.norm(x)
      n = x / xmag
      surfacePosition = (self.radius - xmag) * n 
    
    return surfacePosition


#############################################
class flatpunch_indenter(Geometry):
  """
  Geometry object for creating a grid-aligned flatpunch tip indentor defined by radius and angle
  """
  def __init__(self,
               name,
               x0,
               radius,
               alpha,
               vel=_defaultVelocity,
               mat: int=_defaultMat,
               group: int=_defaultGroup,
               particleType=_defaultParticleType,
               dim: int=3):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x0 = np.asarray(x0)          # tip coordinate (bottom of spherical tip), assumes identor is aligned with y-axis with tip in -y direction
    self.radius = radius
    self.rSqr = radius * radius
    self.alpha = np.radians(alpha)   # tip half angle, input in degrees, stored in radians
    self.dim = dim
    self.y0 = self.radius / np.tan(self.alpha)
    self.vtip = self.x0 - np.array([0.0, self.y0, 0.0])
    
  def isInterior(self, pt, skinDepth):  
    x = np.asarray(pt)-self.vtip
    if self.dim == 2:
      x[2] = 0.0

    surfaceFlag = 0
    if x[1] >= self.y0:
      rr = x[0] * x[0] + x[2] * x[2]
      if rr < (x[1] * np.tan(self.alpha))**2:
        if rr >= (x[1] * np.tan(self.alpha) - skinDepth * np.cos(self.alpha))**2 or x[1] < self.y0 + skinDepth:
          return _defaultSurfaceFlag
        else:
          return 0
    
    return -1    

  def getSurfaceNormal(self, pt):
    x = np.asarray(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      surfaceNormal = n / np.linalg.norm(n)
    else:
      x[1] = x[1] - self.radius - self.x0[1]
      surfaceNormal = x / np.linalg.norm(x)
    
    return surfaceNormal

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      n = n / np.linalg.norm(n)
      x[1] =  x[1] - self.vtip
      surfacePosition = -np.dot(n, x) * n
    else:
      x[1] = x[1] - ( self.radius + self.x0[1] )
      xmag = np.linalg.norm(x)
      n = x / xmag
      surfacePosition = (self.radius - xmag) * n 
    
    return surfacePosition


class simple_tensile_gripper(Geometry):
  """
  Geometry object representing one side of the gripper claw
  default assume microns and base unit of length is mm
  """
  def __init__(self,
               name,
               x0,
               inside_width=30.0e-3,
               hs=5.19615e-3,
               hg=8.5e-3,
               gap=12.0e-3,
               outer_width=55.0e-3,
               thickness=200.0e-3,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel=vel,
                     mat=mat,
                     group=group,
                     particleType=particleType)
    self.x0 = np.asarray(x0) # top of the gripper gap
    self.inside_width = inside_width
    self.hg=hg
    self.hs=hs
    self.gap=gap
    self.outer_width=outer_width
    self.thickness = thickness

  def isInterior(self, pt, skinDepth):
    x = pt-self.x0
    hh = x[1]
    xx = np.abs(x[0])

    # If below can't be inside
    if hh < -self.hg or xx > self.outer_width/2 or x[2] > self.thickness/2:
      return -1

    # Is particle oustide width surface particle or bottom of gripper
    if xx > self.outer_width/2 - skinDepth or xx < -self.hg + skinDepth:
      return 2

    # Is particle in bottom grip section?
    if hh < 0.0 and hh >= -self.hg and xx > self.gap/2:
      if xx < self.gap/2 + skinDepth or hh < -self.hg + skinDepth:
        return 2
      else:
        return 0

    # Is particle in slope of grips
    if hh >= 0 and hh < self.hs and xx > self.gap/2 + 0.5*(self.inside_width-self.gap)*hh/self.hs:
      angle = np.arctan(self.hs*2/(self.inside_width-self.gap))
      xxx = skinDepth/np.sin(angle)
      if xx < self.gap/2 + xxx + 0.5*(self.inside_width-self.gap)*hh/self.hs:
        return 2
      else:
        return 0
   
    # Is particle in straight section
    if hh >= self.hs and xx > self.inside_width/2:
      if xx < self.inside_width/2 + skinDepth:
        return 2
      else:
        return 0     
    
    return -1

  def getSurfaceNormal(self, pt):
    x = pt-self.x0
    hh = x[1]
    xx = np.abs(x[0])

    min_dist = np.inf
    surf_norm = _defaultSurfaceNormal

    if hh > -self.hg and hh < 0:
      # Gripper gap
      new_dist = xx - self.gap/2
      new_norm = np.array([-x[0]/xx, 0.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = np.copy(new_norm)

      # Gripper bottom
      new_dist = np.abs(hh + self.hg)
      new_norm = np.array([0.0, -1.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm

    # Gripper slope
    if hh >= 0 and hh < self.hs:
      xh = self.gap/2 + 0.5*(self.inside_width-self.gap)*hh/self.hs
      angle = np.arctan(self.hs*2/(self.inside_width-self.gap))
      dy = np.tan(np.pi/2- angle)
      m = np.sqrt(1+dy*dy)

      new_dist = np.sin(angle)*np.abs(xx-xh)
      new_norm = np.array([-x[0]/xx, dy, 0.0])/m
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm

    if hh >= self.hs:
      # Interior straight walls
      new_dist = xx - self.inside_width/2
      new_norm = np.array([-x[0]/xx, 0.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm
    
    new_dist = self.outer_width/2 - xx
    new_norm = np.array([x[0]/xx, 0.0, 0.0])
    if new_dist < min_dist:
      min_dist = new_dist
      surf_norm = new_norm

    # If none of the above should be outside walls
    return surf_norm

  def getSurfacePosition(self, pt):
    x = pt-self.x0
    hh = x[1]
    xx = np.abs(x[0])

    min_dist = np.inf
    surf_norm = _defaultSurfaceNormal

    if hh > -self.hg and hh < 0:
      # Gripper gap
      new_dist = xx - self.gap/2
      new_norm = np.array([-x[0]/xx, 0.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = np.copy(new_norm)

      # Gripper bottom
      new_dist = np.abs(hh + self.hg)
      new_norm = np.array([0.0, -1.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm

    # Gripper slope
    if hh >= 0 and hh < self.hs:
      xh = self.gap/2 + 0.5*(self.inside_width-self.gap)*hh/self.hs
      angle = np.arctan(self.hs*2/(self.inside_width-self.gap))
      dy = np.tan(np.pi/2- angle)
      m = np.sqrt(1+dy*dy)

      new_dist = np.sin(angle)*np.abs(xx-xh)
      new_norm = np.array([-x[0]/xx, dy, 0.0])/m
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm

    if hh >= self.hs:
      # Interior straight walls
      new_dist = xx - self.inside_width/2
      new_norm = np.array([-x[0]/xx, 0.0, 0.0])
      if new_dist < min_dist:
        min_dist = new_dist
        surf_norm = new_norm
    
    new_dist = self.outer_width/2 - xx
    new_norm = np.array([x[0]/xx, 0.0, 0.0])
    if new_dist < min_dist:
      min_dist = new_dist
      surf_norm = new_norm

    # If none of the above should be outside walls
    return min_dist*surf_norm

  def xMin(self):
    return -self.outer_width/2 + self.x0[0]

  def xMax(self):
    return self.outer_width/2 + self.x0[0] 


#############################################
class tensile_gripper(Geometry):
  """
  Tensile silicon gripper for micromechanical experiments (currently uses dimensions of Alemnis medium tensile gripper)
  """
  def __init__(self,
               name,
               x0,
               r=15.0,
               h=20.0,
               h2=5.19615,
               h3=8.5,
               gap=12.0,
               width=55.0,
               thickness=200.0,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel=vel,
                     mat=mat,
                     group=group,
                     particleType=particleType)
    self.x0 = np.asarray(x0)
    self.r=r
    self.h=h
    self.h2=h2
    self.h3=h3
    self.gap=gap
    self.width=width
    self.thickness=thickness

  def isInterior(self, pt, skinDepth):
    inside = False
    surface = False

    x = pt-self.x0

    # Is within thickness?
    if np.abs(x[2]) < self.thickness/2 and np.abs(x[0]) < self.width/2:
      hh = x[1]
      # Is it above all interior features?
      if hh > self.r + self.h + self.h2 + skinDepth:
        inside=True

      # Is particle in bottom grip section?
      if hh < 0.0 and hh > -self.h3 and np.abs(x[0]) > self.gap/2:
        inside = True
        if hh < -self.h3 + skinDepth or np.abs(x[0]) < self.gap/2 + skinDepth:
          surface = True

      # Is particle in slope of grips
      if hh >= 0 and hh < self.h2 and np.abs(x[0]) > self.gap/2 + (self.r-self.gap/2) * hh/self.h2:
        inside = True
        angle = np.arctan(self.h2/(self.r-self.gap*0.5))
        xx = skinDepth/np.sin(angle)
        if np.abs(x[0]) < self.gap/2 + xx + (self.r-self.gap/2) * hh/self.h2:
          surface = True

      # Is particle in straight section
      if hh >= self.h2 and hh < self.h2 + self.h and np.abs(x[0]) > self.r:
        inside = True
        if np.abs(x[0]) < self.r + skinDepth:
          surface = True
      
      # Is particle in top arch?
      if hh >= self.h + self.h2:
        rad = np.array([x[0], x[1]-self.h-self.h2, 0.0])
        rad = np.inner(rad, rad)
        if rad > self.r  * self.r:
          inside = True
          if rad < (self.r + skinDepth)*(self.r + skinDepth):
            surface = True          
      
      if inside:
        # Is within thickness surface depth?
        if np.abs(x[2]) > self.thickness/2 - skinDepth or np.abs(x[0]) > self.width/2 - skinDepth:
          surface = 2

    if inside:
      if surface:
        return 2
      else:
        return 0
    
    return -1

  def getSurfaceNormal(self, pt):
    return _defaultSurfaceNormal

  def getSurfacePosition(self, pt):
    return _defaultSurfacePosition

  def xMin(self):
    return -self.width/2

  def xMax(self):
    return self.width/2


#############################################
class box(Geometry):
    """
    Geometry object for creating a grid-aligned box defined by two corners.
    """
    def __init__(self, 
                 name, 
                 x0, 
                 x1, 
                 vel=_defaultVelocity, 
                 mat=_defaultMat, 
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      super().__init__(name,
                       vel = vel,
                       mat = mat,
                       group = group,
                       particleType = particleType)
      type_check_scalar(dim, "Dimension", int)
      type_check_array(x0, "x0", dim, floatType)
      type_check_array(x1, "x1", dim, floatType)
      type_check_array(flaggedSurfaces, "Flagged surfacecs", 2*dim, bool)

      self.dim = dim
      self.x0 = np.asarray(x0[:self.dim])
      self.x1 = np.asarray(x1[:self.dim])
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.asarray(flaggedSurfaces[:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.asarray(pt[:self.dim])
      if np.all( np.logical_and( x >= self.xmin, x < self.xmax) ):
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfaceNormal = np.array([0.0, 0.0, 0.0])
      surfaceNormal[d] = s
      return surfaceNormal

    def getSurfacePosition(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfacePosition = np.array([0.0, 0.0, 0.0])
      surfacePosition[d] = dx[minI]
      return surfacePosition

      # x = np.asarray(pt[:self.dim])
      # dx = np.vstack((self.x0 - x, self.x1 - x))

      # dxI = np.argmin(np.absolute(dx), axis=0)
      # dxMin = []
      # for d in range(self.dim):
      #   dxMin.append(dx[dxI[d]][d])
      # dxMin = np.array(dxMin)

      # minI = np.argmin(np.absolute(dxMin))

      # surfacePos = np.array([0.0, 0.0, 0.0])
      # surfacePos[minI] = dx[dxI[minI]][minI]

      # return surfacePos

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]

#############################################
class notchedBar(Geometry):
    """
    Grid-aligned box defined by two corners with an edge notch in the +y face, having a 45degree
    angle and specified depth

    NOTE, TODO: Surface normals around notch have not been defined.
    """
    def __init__(self, 
                 name, 
                 x0, # corner 1 
                 x1, # corner 2
                 h, # notch depth
                 vel=_defaultVelocity, 
                 mat: int=_defaultMat, 
                 group: int=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim: int=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      super().__init__(name,
                       vel = vel,
                       mat = mat,
                       group = group,
                       particleType = particleType)
      self.dim = dim
      self.x0 = np.asarray(x0[:self.dim])
      self.x1 = np.asarray(x1[:self.dim])
      self.h = h
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.asarray(flaggedSurfaces[:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.asarray(pt[:self.dim])
      if ( np.all( np.logical_and( x >= self.xmin, x < self.xmax) ) and ( x[1] < self.xmax[1] - self.h + np.abs(x[0] - 0.5*( self.xmin[0] + self.xmax[0] ) ) ) ):    
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfaceNormal = np.array([0.0, 0.0, 0.0])
      surfaceNormal[d] = s
      return surfaceNormal

    def getSurfacePosition(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfacePosition = np.array([0.0, 0.0, 0.0])
      surfacePosition[d] = dx[minI]
      return surfacePosition

      # x = np.asarray(pt[:self.dim])
      # dx = np.vstack((self.x0 - x, self.x1 - x))

      # dxI = np.argmin(np.absolute(dx), axis=0)
      # dxMin = []
      # for d in range(self.dim):
      #   dxMin.append(dx[dxI[d]][d])
      # dxMin = np.array(dxMin)

      # minI = np.argmin(np.absolute(dxMin))

      # surfacePos = np.array([0.0, 0.0, 0.0])
      # surfacePos[minI] = dx[dxI[minI]][minI]

      # return surfacePos

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]

#############################################
class box2(Geometry):
    """
    Geometry object for creating a grid-aligned box defined by two corners, mixed normals at corners for contact testing
    """
    def __init__(self, 
                 name, 
                 x0, 
                 x1, 
                 vel=_defaultVelocity, 
                 mat=_defaultMat, 
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim: int=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      super().__init__(name,
                       vel = vel,
                       mat = mat,
                       group = group,
                       particleType = particleType)
      self.dim = dim
      self.x0 = np.asarray(x0[:self.dim])
      self.x1 = np.asarray(x1[:self.dim])
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.asarray(flaggedSurfaces[:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.asarray(pt[:self.dim])
      if np.all( np.logical_and( x >= self.xmin, x < self.xmax) ):
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      dx_min = np.min(np.absolute(dx))
      nearest = ( np.absolute( dx ) - dx_min )< 1e-6
      

      surfaceNormal = np.array( [dx[0] if nearest[0] else 0.0 + dx[2] if nearest[2] else 0.0, dx[1] if nearest[1] else 0.0 + dx[3] if nearest[3] else 0.0, 0.0])
      # print(dx, dx_min, nearest, surfaceNormal, np.linalg.norm(surfaceNormal))
      # surfaceNormal = np.array( [nearest[0] * -dx[0] + nearest[3] * dx[3], nearest[1] * -dx[1] + nearest[4] * dx[4], nearest[2] * -dx[2] + nearest[5] * dx[5]])
      # if self.dim == 2:
      #   surfaceNormal[2] = 0.0
      return surfaceNormal / np.linalg.norm(surfaceNormal)

    def getSurfacePosition(self, pt):
      x = np.asarray(pt[:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      dx_min = np.min(np.absolute(dx))
      nearest = ( np.absolute( dx ) - dx_min )< 1e-2

      surfaceNormal = np.array( [dx[0] if nearest[0] else 0.0 + dx[2] if nearest[2] else 0.0, dx[1] if nearest[1] else 0.0 + dx[3] if nearest[3] else 0.0, 0.0])
      # if self.dim == 2:
      #   surfaceNormal[2] = 0.0
      return surfaceNormal

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]


#############################################
class polygon(Geometry):
  """
  Geometry object for creating a polygon described by ordered vertices.'
  """
  def __init__(self,
               name,
               plist,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=None):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.plist = np.asarray(plist)
    self.flaggedSurfaces = flaggedSurfaces
    if self.flaggedSurfaces is None:
      self.flaggedSurfaces = [True for i in range(self.plist.shape[0])]

  def ccw(self,x1,x2,x3):
    return (x3[1]-x1[1])*(x2[0]-x1[0]) - (x2[1]-x1[1])*(x3[0]-x1[0]) > -10**-16

  def intersect(self,A,B,C,D):
    return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

  def isInside(self,vertices,point):
    v_arr = np.asarray(vertices)
    xmin = min(v_arr[:,0])
    ymin = min(v_arr[:,1])
    xmax = max(v_arr[:,0])
    ymax = max(v_arr[:,1])
    dx = xmax - xmin
    dy = ymax - ymin
    outside = [xmin,ymin]
    outside[0] -= dx
    outside[1] -= dy
    nv = v_arr.shape[0]
    nIntersections = 0
    for i in range(nv):
      p1 = vertices[i]
      p2 = vertices[np.mod(i+1,nv)]
      if(self.intersect(p1,p2,outside,point)):
        nIntersections += 1
    if np.mod(nIntersections,2)==0:
      return False
    else:
      return True

  def isInterior(self, point, skinDepth):
    x = np.asarray(point)
    if(self.isInside(self.plist,x)):
      vertices = self.plist
      nv = vertices.shape[0]

      nearestI = -1
      nearestEdgeD = np.inf
      for i in range(nv):
        if not self.flaggedSurfaces[i]:
          continue
        
        j = np.mod(i+1,nv)
        v = x - vertices[i,:]
        w = vertices[j,:] - vertices[i,:]

        wNorm = np.linalg.norm(w)
        w = w / wNorm
        dw = np.dot(v,w)
        if dw >= 0.0 and dw <= wNorm:
          v = v - dw * w
          d = np.linalg.norm(v)
          if d < nearestEdgeD:
            nearestEdgeD = d
            nearestI = i

      if nearestEdgeD < skinDepth:
        return _defaultSurfaceFlag
      
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)
    # Assumes point is already internal
    # Find the nearest edge and use it's normal
    vertices = self.plist
    nv = vertices.shape[0]

    nearestEdgeD = np.inf
    surfaceNormal = np.empty((1,3))
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfaceNormal = -v / d

    return surfaceNormal
  
  def getSurfacePosition(self,pt):
    x = np.asarray(pt)
    # Assumes point is already internal
    # Find the nearest edge and vector to surface
    vertices = self.plist
    nv = vertices.shape[0]
    
    nearestEdgeD = np.inf
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfacePosition = -v

    return surfacePosition

  def signedDist(self, pt):
    # Computes the distance to the nearest surface
    x = np.asarray(pt)
    x = x[:2]

    vertices = self.plist
    nv = vertices.shape[0]
    
    nearestEdgeD = np.inf
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]
      w = w / np.linalg.norm(w)

      n = np.cross(np.array([0.0, 0.0, 1.0]), np.array([w[0], w[1], 0.0]))

      # print(vertices[i,:], vertices[j,:], v, w, n, np.dot(n[:2], v) )
      nearestEdgeD = min(nearestEdgeD, np.dot(n[:2], v) )

    return nearestEdgeD

  def xMin(self):
    arr = np.asarray(self.plist)
    xmin = min(arr[:,0])
    return xmin

  def xMax(self):
    arr = np.asarray(self.plist)
    xmax = max(arr[:,0])
    return xmax


#############################################
class foam(Geometry):
  """
  Geometry object for creating a grid-aligned box defined by two corners with spherical pores.
  
  This fills a box with spherical pores defined by the array pores
  assumes this is a numpy array [[r,x,y,z],[r,x,y,z],...]
  The surface flags should be correctly set and the searching of pores
  for isInterior uses a KD-tree. which should be faster than
  using particle file writer with each pore as its own object.
  """
  def __init__(self,
               name,
               x0,
               x1,
               pores,
               vel=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.pores = pores
    self.k = min(5, len(pores))

    pts = []
    for p in pores:
      pts.append([ p[1], p[2], p[3] ])
    
    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=int(np.ceil(len(pts)/2)), metric='euclidean')

  def isInterior(self,pt,skinDepth):   
    x = np.asarray(pt)
    
    if np.all(np.logical_and(x >= self.x0, x < self.x1)):
      _, index = self.kdt.query(x.reshape(1,-1), k=self.k)
      surfaceFlag = 0 # Interior unless otherwise determined
      for i in index[0]:
        p=self.pores[i]
        rSqr = ( x[0] - p[1] )**2 + ( x[1] - p[2] )**2 + ( x[2] - p[3] )**2
        if rSqr < p[0]**2:
          return -1 # Inside a pore
        
        if rSqr < ( p[0] + skinDepth )**2:
          surfaceFlag = _defaultSurfaceFlag # near internal pore surface
          
      return surfaceFlag
    
    return -1

  def xMin(self):
    return self.x0[0]

  def xMax(self):
    return self.x1[0]


#############################################
class poissonDiskFoam(Geometry):
  """
  Geometry object for a grid-aligned foam block with Poisson-disk pores.

  In three dimensions this object removes monodisperse spherical pores from the
  retained box.  With ``dim=2`` it removes circular plane-strain pores from the
  active x-y plane; the pores are implicitly extruded through the z thickness.
  The pore centers are generated with a minimum center spacing and optional
  periodicity, so pore cuts match across periodic faces.

  Parameters
  ----------
  x0, x1 : array-like
      Opposite corners of the foam block.  The first ``dim`` entries are used.
  poreRadius or poreDiameter : float
      Monodisperse pore size.  If poreDiameter is supplied, poreRadius is taken
      as half of it.
  porosity : float
      Target void volume fraction for ``dim=3`` or area fraction for ``dim=2``.
      The realized value is stored as ``realizedPorosity`` because an integer
      number of pores is generated.
  periodic : array-like bool
      Periodicity of the pore structure in the first ``dim`` directions.  For
      nonperiodic directions, pore centers are kept at least one pore radius away
      from the block faces.
  minLigament : float
      Extra clearance added to the minimum pore-center spacing.  The minimum
      center spacing is ``2*poreRadius + minLigament``.
  flaggedSurfaces : array-like bool
      Exterior foam-block surfaces to flag, ordered consistently with ``box`` as
      [x-, y-, z-, x+, y+, z+] for ``dim=3`` and [x-, y-, x+, y+] for ``dim=2``.
      Pore surfaces are always flagged.
  """
  def __init__(self,
               name,
               x0,
               x1,
               poreRadius=None,
               poreDiameter=None,
               porosity=0.1,
               seed=1,
               periodic=[False, False, False],
               minLigament=0.0,
               maxTrialsPerPore=10000,
               flaggedSurfaces=None,
               dim=3,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)

    self.dim = int(dim)
    if self.dim not in (2, 3):
      raise ValueError("poissonDiskFoam supports dim=2 or dim=3")

    x0 = np.asarray(x0, dtype=float)
    x1 = np.asarray(x1, dtype=float)
    if x0.size < self.dim or x1.size < self.dim:
      raise ValueError("poissonDiskFoam requires x0 and x1 with length at least dim")
    self.x0 = np.minimum(x0[:self.dim], x1[:self.dim])
    self.x1 = np.maximum(x0[:self.dim], x1[:self.dim])
    self.dx = self.x1 - self.x0

    if np.any(self.dx <= 0.0):
      raise ValueError("poissonDiskFoam requires x1 to be greater than x0 in all active directions")

    if poreDiameter is not None:
      poreDiameter = float(poreDiameter)
      if poreRadius is None:
        poreRadius = 0.5 * poreDiameter
      elif not math.isclose(float(poreDiameter), 2.0 * float(poreRadius), rel_tol=1.0e-10, abs_tol=1.0e-14):
        raise ValueError("poreDiameter must equal 2*poreRadius when both are supplied")

    if poreRadius is None:
      raise ValueError("poissonDiskFoam requires poreRadius or poreDiameter")

    self.poreRadius = float(poreRadius)
    self.poreDiameter = 2.0 * self.poreRadius
    self.porosity = float(porosity)
    self.seed = int(seed)
    periodic = np.asarray(periodic, dtype=bool)
    if periodic.size < self.dim:
      raise ValueError("poissonDiskFoam requires periodic with length at least dim")
    self.periodic = periodic[:self.dim]
    self.minLigament = float(minLigament)
    self.minCenterSpacing = 2.0 * self.poreRadius + self.minLigament
    self.maxTrialsPerPore = int(maxTrialsPerPore)

    if self.poreRadius <= 0.0:
      raise ValueError("poreRadius must be positive")
    if self.porosity < 0.0 or self.porosity >= 1.0:
      raise ValueError("porosity must satisfy 0 <= porosity < 1")
    if self.minLigament < 0.0:
      raise ValueError("minLigament must be non-negative")
    if self.maxTrialsPerPore <= 0:
      raise ValueError("maxTrialsPerPore must be positive")

    for d in range(self.dim):
      if self.poreDiameter >= self.dx[d]:
        raise ValueError("poreDiameter must be smaller than the domain length in every active direction")

    if flaggedSurfaces is None:
      flaggedSurfaces = []
      for d in range(self.dim):
        flaggedSurfaces.append(not self.periodic[d])
      for d in range(self.dim):
        flaggedSurfaces.append(not self.periodic[d])
    flaggedSurfaces = np.asarray(flaggedSurfaces, dtype=bool)
    if flaggedSurfaces.size < 2*self.dim:
      raise ValueError("poissonDiskFoam requires flaggedSurfaces with length at least 2*dim")
    self.flaggedSurfaces = flaggedSurfaces[:2*self.dim]

    self.measure = float(np.prod(self.dx))
    if self.dim == 3:
      self.poreMeasure = 4.0 * math.pi * self.poreRadius**3 / 3.0
      self.measureName = "volume fraction"
    else:
      self.poreMeasure = math.pi * self.poreRadius**2
      self.measureName = "area fraction"

    self.targetPoreCount = int(round(self.porosity * self.measure / self.poreMeasure)) if self.porosity > 0.0 else 0
    if self.porosity > 0.0:
      self.targetPoreCount = max(1, self.targetPoreCount)

    self.centers = self._samplePoreCenters(self.targetPoreCount)
    self.realizedPoreCount = int(self.centers.shape[0])
    self.realizedPorosity = self.realizedPoreCount * self.poreMeasure / self.measure

    if self.realizedPoreCount > 0:
      self.pores = np.column_stack((np.full(self.realizedPoreCount, self.poreRadius), self.centers))
    else:
      self.pores = np.empty((0, self.dim + 1), dtype=float)

    self._buildPoreSearchTree()

    if g_rank == 0:
      print(f"poissonDiskFoam '{self.name}': target {self.measureName} = {self.porosity:.6g}, "
            f"realized {self.measureName} = {self.realizedPorosity:.6g}, pores = {self.realizedPoreCount}, "
            f"pore diameter = {self.poreDiameter:.6g}, dim = {self.dim}")

  def _samplePoreCenters(self, targetCount):
    if targetCount == 0:
      return np.empty((0, self.dim), dtype=float)

    centerX0 = self.x0.copy()
    centerDx = self.dx.copy()
    for d in range(self.dim):
      if not self.periodic[d]:
        centerX0[d] += self.poreRadius
        centerDx[d] -= 2.0 * self.poreRadius

    if np.any(centerDx <= 0.0):
      raise ValueError("No valid pore-center domain remains after applying nonperiodic pore-radius margins")

    rng = np.random.default_rng(self.seed)
    centers = []
    rejectedSinceLastAccept = 0
    minSpacing2 = self.minCenterSpacing * self.minCenterSpacing

    while len(centers) < targetCount:
      if rejectedSinceLastAccept >= self.maxTrialsPerPore:
        raise RuntimeError(
          "poissonDiskFoam could only place " + str(len(centers)) + " of " + str(targetCount) +
          " requested pores. Reduce porosity, reduce pore size/minLigament, enlarge the domain, " +
          "or increase maxTrialsPerPore."
        )

      candidate = centerX0 + rng.random(self.dim) * centerDx

      if len(centers) == 0:
        centers.append(candidate)
        rejectedSinceLastAccept = 0
        continue

      centerArray = np.asarray(centers)
      delta = np.abs(centerArray - candidate)
      for d in range(self.dim):
        if self.periodic[d]:
          delta[:, d] = np.minimum(delta[:, d], self.dx[d] - delta[:, d])

      if np.all(np.einsum('ij,ij->i', delta, delta) >= minSpacing2):
        centers.append(candidate)
        rejectedSinceLastAccept = 0
      else:
        rejectedSinceLastAccept += 1

    return np.asarray(centers, dtype=float)

  def _buildPoreSearchTree(self):
    if self.realizedPoreCount == 0:
      self.queryCenters = np.empty((0, self.dim), dtype=float)
      self.queryRadii = np.empty((0,), dtype=float)
      self.kdt = None
      return

    offsets = []
    for d in range(self.dim):
      if self.periodic[d]:
        offsets.append([-self.dx[d], 0.0, self.dx[d]])
      else:
        offsets.append([0.0])

    centers = []
    radii = []
    for offset in itertools.product(*offsets):
      offsetArray = np.asarray(offset, dtype=float)
      centers.append(self.centers + offsetArray)
      radii.append(np.full(self.realizedPoreCount, self.poreRadius))

    self.queryCenters = np.vstack(centers)
    self.queryRadii = np.concatenate(radii)
    leafSize = max(1, int(np.ceil(self.queryCenters.shape[0] / 20)))
    self.kdt = KDTree(self.queryCenters, leaf_size=leafSize, metric='euclidean')

  def _nearestPoreWithin(self, pt, searchDistance):
    if self.kdt is None:
      return np.inf, None, 0.0, np.inf

    x = np.asarray(pt[:self.dim], dtype=float)
    idx = self.kdt.query_radius(x.reshape(1, -1), r=float(searchDistance))[0]
    if len(idx) == 0:
      return np.inf, None, 0.0, np.inf

    delta = self.queryCenters[idx] - x
    d2 = np.einsum('ij,ij->i', delta, delta)
    j = int(np.argmin(d2))
    qidx = int(idx[j])
    centerDistance = math.sqrt(float(d2[j]))
    radius = float(self.queryRadii[qidx])
    return centerDistance - radius, self.queryCenters[qidx], radius, centerDistance

  def _nearestPore(self, pt):
    if self.kdt is None:
      return np.inf, None, 0.0, np.inf

    x = np.asarray(pt[:self.dim], dtype=float)
    dist, idx = self.kdt.query(x.reshape(1, -1), k=1)
    qidx = int(idx[0][0])
    centerDistance = float(dist[0][0])
    radius = float(self.queryRadii[qidx])
    return centerDistance - radius, self.queryCenters[qidx], radius, centerDistance

  def _nearestOuterSurface(self, pt):
    x = np.asarray(pt[:self.dim], dtype=float)
    bestDistance = np.inf
    bestNormal = _defaultSurfaceNormal.copy()
    bestPosition = _defaultSurfacePosition.copy()

    for d in range(self.dim):
      if self.flaggedSurfaces[d]:
        distance = x[d] - self.x0[d]
        if distance < bestDistance:
          bestDistance = distance
          bestNormal = np.zeros(3)
          bestNormal[d] = -1.0
          bestPosition = np.zeros(3)
          bestPosition[d] = -distance

      if self.flaggedSurfaces[d + self.dim]:
        distance = self.x1[d] - x[d]
        if distance < bestDistance:
          bestDistance = distance
          bestNormal = np.zeros(3)
          bestNormal[d] = 1.0
          bestPosition = np.zeros(3)
          bestPosition[d] = distance

    return bestDistance, bestNormal, bestPosition

  def isInterior(self, pt, skinDepth):
    x = np.asarray(pt[:self.dim], dtype=float)
    if not np.all(np.logical_and(x >= self.x0, x < self.x1)):
      return -1

    poreSurfaceDistance, _, _, _ = self._nearestPoreWithin(x, self.poreRadius + skinDepth)
    if poreSurfaceDistance < 0.0:
      return -1

    surfaceFlag = 0
    if poreSurfaceDistance <= skinDepth:
      surfaceFlag = _defaultSurfaceFlag

    outerDistance, _, _ = self._nearestOuterSurface(x)
    if outerDistance <= skinDepth:
      surfaceFlag = _defaultSurfaceFlag

    return surfaceFlag

  def getSurfaceNormal(self, pt):
    x = np.asarray(pt[:self.dim], dtype=float)
    poreSurfaceDistance, center, radius, centerDistance = self._nearestPore(x)
    outerDistance, outerNormal, _ = self._nearestOuterSurface(x)

    if center is not None and poreSurfaceDistance <= outerDistance:
      if centerDistance <= 0.0:
        return _defaultSurfaceNormal
      normal = np.zeros(3)
      normal[:self.dim] = (center - x) / centerDistance
      return normal

    return outerNormal

  def getSurfacePosition(self, pt):
    x = np.asarray(pt[:self.dim], dtype=float)
    poreSurfaceDistance, center, radius, centerDistance = self._nearestPore(x)
    outerDistance, _, outerPosition = self._nearestOuterSurface(x)

    if center is not None and poreSurfaceDistance <= outerDistance:
      if centerDistance <= 0.0:
        return _defaultSurfacePosition
      pos = np.zeros(3)
      pos[:self.dim] = center + (x - center) * (radius / centerDistance) - x
      return pos

    return outerPosition

  def xMin(self):
    return self.x0[0]

  def xMax(self):
    return self.x1[0]

  def yMin(self):
    return self.x0[1]

  def yMax(self):
    return self.x1[1]

  def zMin(self):
    return self.x0[2] if self.dim == 3 else 0.0

  def zMax(self):
    return self.x1[2] if self.dim == 3 else 0.0


#############################################
class twoMaterialBox(Geometry):
  """
  Geometry object for creating a grid-aligned box defined by two corners with random assignment of group 1 or group 2.
  """
  def __init__(self,
               name,
               x0,
               x1,
               vel=_defaultVelocity,
               mat1=0,
               mat2=1,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     group = group,
                     particleType = particleType)
    assert len(x0) == 3, "x0 must have length of 3"
    for i in x0:
       assert isinstance(i, floatType), "x0 elements must be floats"
    assert len(x1) == 3, "x1 must have length of 3"
    for i in x1:
       assert isinstance(i, floatType), "x1 elements must be floats"

    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    
    self.mats = [mat1,mat2]
    
  def getMat(self,pt):
    mat = random.choice(self.mats)
    print(' assigning mat ', mat)
    return mat # will error if map doesn't contain index as key

  def isInterior(self, pt, skinDepth):
    x = np.asarray(pt)
    # is the point within the object?
    if np.all( np.logical_and(x > self.x0, x < self.x1) ):
      # is point on the surface?
      if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def xMin(self):
    return min(self.x0[0], self.x1[0])

  def xMax(self):
    return max(self.x0[0], self.x1[0])


#############################################
class twoFieldBox(Geometry):
  """
  Geometry object for creating a grid-aligned box defined by two corners with random assignment of group 1 or group 2.
  """
  def __init__(self,
               name,
               x0,
               x1,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group1=0,
               group2=1,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    assert len(x0) == 3, "x0 must have length of 3"
    for i in x2:
       assert isinstance(i, floatType), "x0 elements must be floats"
    assert len(x1) == 3, "x1 must have length of 3"
    for i in x1:
       assert isinstance(i, floatType), "x1 elements must be floats"

    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self, pt, skinDepth):
    x = np.asarray(pt)
    # is the point within the object?
    if np.all( np.logical_and(x > self.x0, x < self.x1) ):
      # is point on the surface?
      if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def xMin(self):
    return min(self.x0[0], self.X1[0])

  def xMax(self):
    return max(self.x0[0], self.X1[0])


#############################################
class cone(Geometry):
  """
  Geometry object for creating a cone defined by two points and base radius
  """
  def __init__(self,
               name,
               x1,
               x2,
               r,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    type_check_array(x1, "x1", 3, floatType)
    type_check_array(x2, "x2", 3, floatType)
    type_check_scalar(r, "r", floatType)

    self.x1 = np.asarray(x1)
    self.x2 = np.asarray(x2)
    self.r = r

    self.h = np.linalg.norm(self.x2-self.x1) # height of cone axis
    self.axis = (self.x2-self.x1)/self.h    

  def isInterior(self, pt, skinDepth):    
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

    if (z >= 0 and z < self.h):
      r_h = self.r * (1-z / self.h)
      if (r < r_h):
          theta = np.arccos(self.r/np.sqrt(self.h * self.h + self.r * self.r))
          s = skinDepth / np.sin(theta) 
          if r >= r_h-s or z < skinDepth:
            return _defaultSurfaceFlag
          else:
            return 0    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x-z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    theta = np.arctan(self.h/self.r)

    rh = self.r * (1-z / self.h)
    min_dist = (rh - r) * np.sin(theta)

    surf_norm = xr + self.axis*r/np.tan(theta)
    surf_norm = surf_norm/np.linalg.norm(surf_norm)  

    new_dist = z
    new_norm = -self.axis

    if new_dist < min_dist:
      surf_norm = new_norm

    return surf_norm
   
  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x-z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    theta = np.arctan(self.h/self.r)

    rh = self.r * (1-z / self.h)
    min_dist = (rh - r) * np.sin(theta)

    surf_norm = xr + self.axis*r/np.tan(theta)
    surf_norm = surf_norm/np.linalg.norm(surf_norm)  

    new_dist = z
    new_norm = -self.axis

    if new_dist < min_dist:
      surf_norm = new_norm
      min_dist = new_dist

    return surf_norm*min_dist
 
  # def xMin(self):
  #   return min(self.x1[0], self.x2[0])-self.r

  # def xMax(self):
  #   return max(self.x1[0], self.x2[0])+self.r


#############################################
class stl(Geometry):
  """
  Geometry object for a closed triangulated STL surface.

  This object follows the same interface used by analytic geometry objects such
  as sphere and cylinder:

    * isInterior(pt, skinDepth) returns -1 outside the STL volume, 0 for
      interior particles, and _defaultSurfaceFlag for interior particles within
      skinDepth of the STL surface.
    * getSurfaceNormal(pt) returns the outward unit normal of the nearest STL
      triangle.  Pass flipNormal=True if the STL winding is inward.
    * getSurfacePosition(pt) returns the vector from pt to the nearest point on
      the STL surface.  This matches the relative-vector convention used by
      cylinder.getSurfacePosition().

  The inside/outside test assumes the STL is a reasonably closed manifold.  It
  casts rays along rayAxis and uses a 2D binning acceleration structure over the
  transverse coordinates.  Surface projection uses a cKDTree of triangle
  centroids followed by exact closest-point checks on the nearest candidates.

  Parameters
  ----------
  fileName / stlFile / filename : str
      STL file path.  Binary and ASCII STL files are supported.  Relative paths
      are interpreted relative to the current PFW run directory, where PFW also
      copies files listed with #[pfw_dependency].
  x0 : array-like length 3
      Translation applied to STL coordinates after scaling.
  scale : float or array-like length 3
      Coordinate scale factor.  The MPM examples use mm, microseconds, and
      milligrams, so STL files authored in mm usually use scale=1.0.
  center : array-like length 3 or None
      Optional target bounding-box center after scale/x0.
  rayAxis : int
      Ray-casting axis; 0 means +x rays through yz bins.
  binCounts : int, tuple, or None
      Number of bins in the two transverse ray-casting directions.  None chooses
      a size based on triangle count.
  kNearest : int
      Number of centroid-nearest triangles to check for surface projection.
  flipNormal : bool
      Flip all STL triangle normals after loading.
  """
  def __init__(self,
               name,
               fileName=None,
               stlFile=None,
               filename=None,
               x0=_defaultSurfacePosition,
               scale=1.0,
               center=None,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               surfaceFlag=_defaultSurfaceFlag,
               flaggedSurfaces=True,
               rayAxis=0,
               binCounts=None,
               kNearest=64,
               flipNormal=False,
               eps=1.0e-10):
    super().__init__(name,
                     vel=vel,
                     mat=mat,
                     group=group,
                     particleType=particleType)

    if fileName is None:
      fileName = stlFile
    if fileName is None:
      fileName = filename
    if fileName is None:
      raise ValueError("stl geometry requires fileName, stlFile, or filename")

    self.fileName = str(fileName)
    self.surfaceFlag = surfaceFlag
    self.flaggedSurfaces = flaggedSurfaces
    self.rayAxis = int(rayAxis)
    if self.rayAxis not in (0, 1, 2):
      raise ValueError("rayAxis must be 0, 1, or 2")
    self.projAxes = [axis for axis in (0, 1, 2) if axis != self.rayAxis]
    self.kNearest = max(1, int(kNearest))
    self.eps = float(eps)

    triangles, normals = self._read_stl(self.fileName)

    scale = np.asarray(scale, dtype=float)
    if scale.shape == ():
      scale = np.array([float(scale), float(scale), float(scale)])
    if scale.shape != (3,):
      raise ValueError("scale must be a scalar or length-3 array")

    shift = np.asarray(x0, dtype=float)
    if shift.shape != (3,):
      raise ValueError("x0 must be a length-3 array")

    triangles = triangles * scale.reshape((1, 1, 3)) + shift.reshape((1, 1, 3))
    if center is not None:
      center = np.asarray(center, dtype=float)
      if center.shape != (3,):
        raise ValueError("center must be None or a length-3 array")
      bb_center = 0.5 * (np.min(triangles, axis=(0, 1)) + np.max(triangles, axis=(0, 1)))
      triangles = triangles + (center - bb_center).reshape((1, 1, 3))

    # Recompute normals after any anisotropic scale so surface normals are
    # consistent with the transformed geometry.
    normals = np.cross(triangles[:, 1, :] - triangles[:, 0, :],
                       triangles[:, 2, :] - triangles[:, 0, :])
    nmag = np.linalg.norm(normals, axis=1)
    keep = nmag > self.eps
    if not np.all(keep):
      triangles = triangles[keep]
      normals = normals[keep]
      nmag = nmag[keep]
    if len(triangles) == 0:
      raise ValueError("STL mesh contains no non-degenerate triangles: " + self.fileName)
    normals = normals / nmag[:, None]
    if flipNormal:
      normals = -normals

    self.triangles = triangles
    self.normals = normals
    self.triMins = np.min(triangles, axis=1)
    self.triMaxs = np.max(triangles, axis=1)
    self.boundsMin = np.min(triangles, axis=(0, 1))
    self.boundsMax = np.max(triangles, axis=(0, 1))
    self.centroids = np.mean(triangles, axis=1)
    self._surfaceCache = {}

    self._build_ray_bins(binCounts)
    self._centroidTree = scipy.spatial.cKDTree(self.centroids)

  @staticmethod
  def _read_stl(path):
    import os
    import struct

    path = os.path.expandvars(os.path.expanduser(str(path)))
    if not os.path.exists(path):
      raise FileNotFoundError("STL file not found: " + str(path))

    size = os.path.getsize(path)
    with open(path, "rb") as f:
      header = f.read(80)
      count_bytes = f.read(4)
      if len(count_bytes) == 4:
        n_tri = struct.unpack("<I", count_bytes)[0]
        expected = 84 + 50 * n_tri
        if n_tri > 0 and expected == size:
          raw = f.read(50 * n_tri)
          dtype = np.dtype([
              ("normal", "<f4", (3,)),
              ("vertices", "<f4", (3, 3)),
              ("attr", "<u2"),
          ])
          data = np.frombuffer(raw, dtype=dtype, count=n_tri)
          return data["vertices"].astype(float), data["normal"].astype(float)

    # Fall back to ASCII STL.  This is intentionally simple and permissive.
    vertices = []
    normals = []
    current_normal = np.array([0.0, 0.0, 0.0])
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
      for line in f:
        parts = line.strip().split()
        if not parts:
          continue
        key = parts[0].lower()
        if key == "facet" and len(parts) >= 5 and parts[1].lower() == "normal":
          current_normal = np.array([float(parts[2]), float(parts[3]), float(parts[4])], dtype=float)
        elif key == "vertex" and len(parts) >= 4:
          vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
          if len(vertices) % 3 == 0:
            normals.append(current_normal.copy())

    if len(vertices) % 3 != 0 or len(vertices) == 0:
      raise ValueError("Could not parse STL file as binary or ASCII: " + str(path))
    triangles = np.asarray(vertices, dtype=float).reshape((-1, 3, 3))
    normals = np.asarray(normals, dtype=float)
    return triangles, normals

  def _build_ray_bins(self, binCounts):
    tri_count = len(self.triangles)
    if binCounts is None:
      n = int(np.clip(np.sqrt(tri_count) / 3.0, 32, 192))
      binCounts = (n, n)
    elif isinstance(binCounts, int):
      binCounts = (binCounts, binCounts)
    else:
      binCounts = tuple(binCounts)

    if len(binCounts) != 2:
      raise ValueError("binCounts must be None, an integer, or a length-2 tuple")

    self.binCounts = (max(1, int(binCounts[0])), max(1, int(binCounts[1])))
    self.projMin = self.boundsMin[self.projAxes].astype(float)
    self.projMax = self.boundsMax[self.projAxes].astype(float)
    span = self.projMax - self.projMin
    span[span <= self.eps] = 1.0
    self.projSpan = span

    self._rayBins = [[[] for _ in range(self.binCounts[1])] for _ in range(self.binCounts[0])]

    for idx in range(tri_count):
      mn = self.triMins[idx, self.projAxes]
      mx = self.triMaxs[idx, self.projAxes]
      lo = np.floor((mn - self.projMin) / self.projSpan * np.asarray(self.binCounts)).astype(int)
      hi = np.floor((mx - self.projMin) / self.projSpan * np.asarray(self.binCounts)).astype(int)
      lo = np.maximum(lo, 0)
      hi = np.minimum(hi, np.asarray(self.binCounts) - 1)
      for i in range(lo[0], hi[0] + 1):
        for j in range(lo[1], hi[1] + 1):
          self._rayBins[i][j].append(idx)

  def _surface_is_flagged(self):
    if isinstance(self.flaggedSurfaces, (list, tuple, np.ndarray)):
      return bool(self.flaggedSurfaces[0]) if len(self.flaggedSurfaces) else False
    return bool(self.flaggedSurfaces)

  def _ray_bin_indices(self, pt):
    q = np.asarray(pt, dtype=float)[self.projAxes]
    if np.any(q < self.projMin - self.eps) or np.any(q > self.projMax + self.eps):
      return None
    ij = np.floor((q - self.projMin) / self.projSpan * np.asarray(self.binCounts)).astype(int)
    ij = np.maximum(ij, 0)
    ij = np.minimum(ij, np.asarray(self.binCounts) - 1)
    return int(ij[0]), int(ij[1])

  def _ray_intersection_coordinates(self, pt):
    bin_ij = self._ray_bin_indices(pt)
    if bin_ij is None:
      return []

    p = np.asarray(pt, dtype=float)
    q = p[self.projAxes]
    a = self.rayAxis
    b, c = self.projAxes
    xs = []

    for idx in self._rayBins[bin_ij[0]][bin_ij[1]]:
      tri = self.triangles[idx]
      # Fast transverse-bounding-box rejection.
      if q[0] < self.triMins[idx, b] - self.eps or q[0] > self.triMaxs[idx, b] + self.eps:
        continue
      if q[1] < self.triMins[idx, c] - self.eps or q[1] > self.triMaxs[idx, c] + self.eps:
        continue

      v0 = tri[0, self.projAxes]
      v1 = tri[1, self.projAxes]
      v2 = tri[2, self.projAxes]
      e1 = v1 - v0
      e2 = v2 - v0
      den = e1[0] * e2[1] - e1[1] * e2[0]
      if abs(den) < self.eps:
        continue
      rhs = q - v0
      u = (rhs[0] * e2[1] - rhs[1] * e2[0]) / den
      v = (e1[0] * rhs[1] - e1[1] * rhs[0]) / den
      w = 1.0 - u - v
      if u >= -self.eps and v >= -self.eps and w >= -self.eps:
        xint = w * tri[0, a] + u * tri[1, a] + v * tri[2, a]
        if xint > p[a] + self.eps:
          xs.append(float(xint))

    if not xs:
      return []
    xs.sort()

    # Shared triangle edges can produce duplicate intersections at the same
    # coordinate.  Count each cluster once for robust odd-even parity.
    unique = [xs[0]]
    tol = 100.0 * self.eps * max(1.0, abs(xs[0]), np.linalg.norm(self.boundsMax - self.boundsMin))
    for x in xs[1:]:
      if abs(x - unique[-1]) > tol:
        unique.append(x)
    return unique

  @staticmethod
  def _closest_point_on_triangle(p, a, b, c):
    # Christer Ericson, Real-Time Collision Detection, closest point on triangle.
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = np.dot(ab, ap)
    d2 = np.dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
      return a

    bp = p - b
    d3 = np.dot(ab, bp)
    d4 = np.dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
      return b

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
      v = d1 / (d1 - d3)
      return a + v * ab

    cp = p - c
    d5 = np.dot(ab, cp)
    d6 = np.dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
      return c

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
      w = d2 / (d2 - d6)
      return a + w * ac

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
      w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
      return b + w * (c - b)

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    return a + ab * v + ac * w

  def _nearest_surface(self, pt):
    p = np.asarray(pt, dtype=float)
    key = tuple(np.round(p, 12))
    cached = self._surfaceCache.get(key)
    if cached is not None:
      return cached

    k = max(1, min(self.kNearest, len(self.triangles)))
    _, idxs = self._centroidTree.query(p, k=k)
    idxs = np.atleast_1d(idxs)

    best_idx = int(idxs[0])
    best_point = self.triangles[best_idx, 0, :]
    best_dist2 = np.inf
    for idx in idxs:
      idx = int(idx)
      tri = self.triangles[idx]
      cp = self._closest_point_on_triangle(p, tri[0], tri[1], tri[2])
      dist2 = float(np.dot(cp - p, cp - p))
      if dist2 < best_dist2:
        best_dist2 = dist2
        best_point = cp
        best_idx = idx

    result = (best_point, self.normals[best_idx], best_dist2, best_idx)
    # Keep the cache bounded; PFW often calls normal and position for the same
    # surface-flagged particle back-to-back.
    if len(self._surfaceCache) > 8192:
      self._surfaceCache.clear()
    self._surfaceCache[key] = result
    return result

  def isInterior(self, pt, skinDepth=0.0):
    p = np.asarray(pt, dtype=float)
    if np.any(p < self.boundsMin - self.eps) or np.any(p > self.boundsMax + self.eps):
      return -1

    intersections = self._ray_intersection_coordinates(p)
    if (len(intersections) % 2) == 0:
      return -1

    if self._surface_is_flagged() and skinDepth is not None and skinDepth > 0.0:
      _, _, dist2, _ = self._nearest_surface(p)
      if dist2 <= skinDepth * skinDepth:
        return self.surfaceFlag
    return 0

  def getSurfaceNormal(self, pt):
    p = np.asarray(pt, dtype=float)
    cp, n, _, _ = self._nearest_surface(p)
    n = np.asarray(n, dtype=float)
    nmag = np.linalg.norm(n)
    if nmag < self.eps:
      return _defaultSurfaceNormal
    n = n / nmag
    # STL normals should be outward.  If an interior near-surface point lies on
    # the positive side of the nearest triangle normal, the local winding is
    # inward or inconsistent, so flip the returned normal for this query.
    if np.dot(n, p - cp) > 0.0:
      n = -n
    return n

  def getSurfacePosition(self, pt):
    p = np.asarray(pt, dtype=float)
    cp, _, _, _ = self._nearest_surface(p)
    return cp - p

  def xMin(self):
    return float(self.boundsMin[0])

  def xMax(self):
    return float(self.boundsMax[0])

  def getSubregions(self):
    return [(self.mat, self.particleType)]


# Upper-case alias for users who prefer geometry class names to mirror the file
# format name.  Existing PFW objects use lower-case class names, so geom.stl is
# the canonical spelling used in the examples.
STL = stl

#############################################
class cylinder(Geometry):
  """
  Geometry object for creating a cylinder defined by two points and radius

  flaggedSurfaces = [bottom, top, outside walls, inside walls]
  """
  def __init__(self,
               name,
               x1,
               x2,
               r,
               ri=0.0,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=[True, True, True, True],
               n2 = None):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    type_check_array(x1, "x1", 3, floatType)
    type_check_array(x2, "x2", 3, floatType)
    type_check_scalar(r, "r", floatType)
    type_check_scalar(r, "ri", floatType)

    self.x1 = np.asarray(x1)
    self.x2 = np.asarray(x2)
    self.r = r
    self.ri = ri #inner radius, defaults to 0 for solid cylinder
    self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks
    self.flaggedSurfaces = flaggedSurfaces

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.axis = (self.x2-self.x1)/self.h    

    # option to define a slide plane passing through x2 with outward normal n2.
    self.n2 = np.asarray(n2)
    if (self.n2.any() == None) or (n2 is None):
      self.n2 = self.axis

    # for now we will default to cylinder having axis-aligned fiber direction,
    # may make more general later, but this allows for backwards compatibility

    if abs(np.dot(self.axis, np.array([0.0,0.0,1.0]))-1) < 1e-12:
      m2 = np.cross(np.array([0.0,1.0,0.0]),self.axis)
    else:
      m2 = np.cross(np.array([0.0,0.0,1.0]),self.axis)
    m2 = m2 / np.linalg.norm(m2)
    m3 = np.cross(self.axis,m2) # TODO double check this

    self.matDir = np.vstack((self.axis, m2, m3))

  def isInterior(self, pt, skinDepth):    
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

    dist_from_top = np.dot(self.n2, pt - self.x2 )

    if (z >= 0 and dist_from_top <= 0 ) and (r > self.ri and r < self.r):
      if ( (z<skinDepth and self.flaggedSurfaces[0]) or ( dist_from_top > -skinDepth and self.flaggedSurfaces[1]) or ( r > self.r-skinDepth and self.flaggedSurfaces[2]) or (self.ri > 0.0 and r < self.ri+skinDepth and self.flaggedSurfaces[3]) ): 
        return _defaultSurfaceFlag
      else:
        return 0    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point

    dist_from_inner_wall = r - self.ri
    dist_from_outer_wall = self.r - r
    dist_from_top = -np.dot( self.n2, pt - self.x2 )
    dist_from_bot = z

    if (self.ri > 0):
        min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )
    else:
        min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return self.n2

    if min_surf_dist == dist_from_bot:
      return -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return xr / r

    if min_surf_dist == dist_from_inner_wall:
      return -xr / r

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    dist_from_top = -np.dot( self.n2, pt - self.x2 )

    dist_from_inner_wall = r - self.ri
    dist_from_outer_wall = self.r - r
    dist_from_top = self.h-z
    dist_from_bot = z

    if (self.ri > 0):
        min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )
    else:
        min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return dist_from_top * self.n2
    
    if min_surf_dist == dist_from_bot:
      return dist_from_bot * -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return (self.r/r - 1) * xr

    if min_surf_dist == dist_from_inner_wall:
      return (self.ri/r - 1) * xr

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r

#############################################
class pillar(Geometry):
  """
  Geometry object for creating a micro pillar compression sample with 
  a cylinder gage length and a cylinder base defined by two points and radii
  where an elliptical radius profile is constructed to have dr/dz=0 at z2
  and a wider base at z1.   (r1>r2)
  
  material has a preffered direction, and defined strength scale layers

  flaggedSurfaces = [bottom, top, outside walls]
  
  WARNING: surface normals and position are radially outward, not
           correct for curved surface
           
  TODO: add the cut n2 capability and inner radius option from cylinder
  """
  def __init__(self,
               name,
               x1,
               x2,
               gage,
               r1,
               r2,
               matDir,           
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=[True, False, True],
               n2 = None):
    Geometry.__init__(self,
                      name,
                      vel = vel,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    type_check_array(x1, "x1", 3, floatType)
    type_check_array(x2, "x2", 3, floatType)
    type_check_array(gage, "gage", floatType)
    type_check_scalar(r1, "r1", floatType)
    type_check_scalar(r2, "r2", floatType)
    type_check_scalar(matDir, "matDir", 3, floatType)
    type_check_scalar(weibullLayerThickness, "weibullLayerThickness", floatType)
    type_check_scalar(weibullReferenceVolume, "weibullReferenceVolume", floatType)
    type_check_scalar(weibullModulus, "weibullModulus", floatType)
    
    self.x1 = np.asarray(x1) # base point
    self.x2 = np.asarray(x2) # top point
        
    self.r1 = r1 # base radius
    self.r2 = r2 # top radius.
    
    self.weibullLayerThickness = weibullLayerThickness
    self.weibullReferenceVolume = weibullReferenceVolume
    self.weibullModulus = weibullModulus
        
    self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks      
    self.flaggedSurfaces = flaggedSurfaces

    # gage length
    self.gage = gage
    # total height
    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.axis = (self.x2-self.x1)/self.h    
    if(self.gage > self.h):
      print( 'pillar defined with gage > h')
    
    # height of base 
    self.base = self.h - self.gage

    # for now we will default to cylinder having axis-aligned fiber direction,
    # may make more general later, but this allows for backwards compatibility
    matDir = np.asarray(matDir)
    if abs(np.dot(matDir, np.array([0.0,0.0,1.0]))-1) < 1e-12:
      m2 = np.cross(np.array([0.0,1.0,0.0]),matDir)
    else:
      m2 = np.cross(np.array([0.0,0.0,1.0]),matDir)
    m2 = m2 / np.linalg.norm(m2)
    m3 = np.cross(matDir,m2)

    self.matDir = np.vstack((matDir, m2, m3))
    
  ## elliptical profile, z is distance from the base
  def R(self, z):
    if(z <= 0):
      return self.r1
    elif (z>= self.base):
      return self.r2
    else:
      return 0.5*( 2.*self.r1*self.base + 2.*( self.r2 - self.r1 )*np.sqrt( z*( 2.*self.base - z ) ) ) / self.base
  
  def isInterior(self, pt, skinDepth: float) -> int:    
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point
    rs = self.R(z) # surface radius at z
       
    if (z >= 0 and z < self.h) and r < rs :
      if ( ( z < skinDepth and self.flaggedSurfaces[0] ) or ( z > self.h - skinDepth and self.flaggedSurfaces[1] ) or ( r > rs-skinDepth and self.flaggedSurfaces[2] ) ): 
        return _defaultSurfaceFlag
      else:
        return 0    
    return -1
  
  

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    
    rs = self.R(z) # surface radius at z

    dist_from_outer_wall = rs- r
    dist_from_top = self.h-z
    dist_from_bot = z

    min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return self.axis

    if min_surf_dist == dist_from_bot:
      return -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return xr / r

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    rs = self.R(z) # surface radius at z

    dist_from_outer_wall = rs - r
    dist_from_top = self.h-z
    dist_from_bot = z

    min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return dist_from_top * self.axis
    
    if min_surf_dist == dist_from_bot:
      return dist_from_bot * -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return (rs/r - 1) * xr

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r1

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r1
  
#############################################
class pillarBase(Geometry):
  """
  Geometry object for creating a cylinder base defined by two points and radii
  where an elliptical radius profile is constructed to have dr/dz=0 at z2
  and a wider base at z1.   (r1>r2)

  flaggedSurfaces = [bottom, top, outside walls]
  
  WARNING: surface normals and position are radially outward, not
           correct for curved surface
  """
  def __init__(self,
               name,
               x1,
               x2,
               r1,
               r2,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=[True, False, True]):
    Geometry.__init__(self,
                      name,
                      vel = vel,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    type_check_array(x1, "x1", 3, floatType)
    type_check_array(x2, "x2", 3, floatType)
    type_check_scalar(r1, "r1", floatType)
    type_check_scalar(r2, "r2", floatType)

    self.x1 = np.asarray(x1) # base point
    self.x2 = np.asarray(x2) # top point
    self.r1 = r1 # base radius
    self.r2 = r2 # top radius.
    self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks
    self.flaggedSurfaces = flaggedSurfaces

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.axis = (self.x2-self.x1)/self.h    

    # for now we will default to cylinder having axis-aligned fiber direction,
    # may make more general later, but this allows for backwards compatibility
    if abs(np.dot(self.axis, np.array([0.0,0.0,1.0]))-1) < 1e-12:
      m2 = np.cross(np.array([0.0,1.0,0.0]),self.axis)
    else:
      m2 = np.cross(np.array([0.0,0.0,1.0]),self.axis)
    m2 = m2 / np.linalg.norm(m2)
    m3 = np.cross(self.axis,m2) # TODO double check this

    self.matDir = np.vstack((self.axis, m2, m3))

  ## elliptical profile:
  def R(self, z):
    if(z <= 0):
      return self.r1
    elif (z>= self.h):
      return self.r2
    else:
      return 0.5*( 2.*self.r1*self.h + 2.*( self.r2 - self.r1 )*np.sqrt( z*( 2.*self.h - z ) ) ) / self.h
  
  def isInterior(self, pt, skinDepth: float) -> int:    
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point
    rs = self.R(z) # surface radius at z
       
    if (z >= 0 and z < self.h) and r < rs :
      if ( (z<skinDepth and self.flaggedSurfaces[0]) or (z>self.h-skinDepth and self.flaggedSurfaces[1]) or ( r > rs-skinDepth and self.flaggedSurfaces[2])  ): 
        return _defaultSurfaceFlag
      else:
        return 0    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    
    rs = self.R(z) # surface radius at z

    dist_from_outer_wall = rs- r
    dist_from_top = self.h-z
    dist_from_bot = z

    min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return self.axis

    if min_surf_dist == dist_from_bot:
      return -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return xr / r

  def getSurfacePosition(self,pt):
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    rs = self.R(z) # surface radius at z

    dist_from_outer_wall = rs - r
    dist_from_top = self.h-z
    dist_from_bot = z

    min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return dist_from_top * self.axis
    
    if min_surf_dist == dist_from_bot:
      return dist_from_bot * -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return (rs/r - 1) * xr

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r1

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r1


# #############################################
# class expandingRing(Geometry):
#   """
#   Geometry object for creating a cylinder defined by two points and radius
#   """
#   def __init__(self,
#                name,
#                x1,
#                x2,
#                innerRadius,
#                outerRadius,
#                radialVelocity,
#                vel=_defaultVelocity,
#                mat=_defaultMat,
#                group=_defaultGroup,
#                particleType=_defaultParticleType):
#     super().__init__(name,
#                      vel = vel,
#                      mat = mat,
#                      group = group,
#                      particleType = particleType)
#     self.x1 = np.asarray(x1)
#     self.x2 = np.asarray(x2)
#     self.innerRadius = innerRadius
#     self.outerRadius = outerRadius

#     self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks
#     self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
#     self.axis = (self.x2-self.x1)/self.h    

#     # for now we will default to cylinder having axis-aligned fiber direction,
#     # may make more general later, but this allows for backwards compatibility

#     if abs(np.dot(self.axis, np.array([0.0,0.0,1.0]))-1) < 1e-12:
#       m2 = np.cross(np.array([0.0,1.0,0.0]),self.axis)
#     else:
#       m2 = np.cross(np.array([0.0,0.0,1.0]),self.axis)
#     m2 = m2 / np.linalg.norm(m2)
#     m3 = np.cross(self.axis,m2) # TODO double check this

#     self.matDir = np.vstack((self.axis, m2, m3))

#   def isInterior(self, pt, skinDepth):    
#     x = np.asarray(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

#     if (z >= 0 and z < self.h) and (r < self.outerRadius) and (r > self.innerRadius):
#       if ( (z<skinDepth) or (z>self.h-skinDepth) or ( r > self.outerRadius-skinDepth) or (r < self.innerRadius+skinDepth) ):
#         return _defaultSurfaceFlag
#       else:
#         return 0
    
#     return -1

#   def getSurfaceNormal(self,pt):
#     x = np.asarray(pt)-self.x1
#     z = np.dot(x, self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     dist_from_inner_wall = r - self.innerRadius
#     dist_from_outer_wall = self.outerRadius - r
#     dist_from_top = self.h-z
#     dist_from_bot = z

#     min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )

#     if min_surf_dist == dist_from_top:
#       return self.axis

#     if min_surf_dist == dist_from_bot:
#       return -self.axis

#     if min_surf_dist == dist_from_inner_wall:
#       return -xr / r

#     if min_surf_dist == dist_from_outer_wall:
#       return xr / r

#   def getSurfacePosition(self,pt):
#     x = np.asarray(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     dist_from_inner_wall = r - self.innerRadius
#     dist_from_outer_wall = self.outerRadius - r
#     dist_from_top = self.h-z
#     dist_from_bot = z

#     min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )

#     if min_surf_dist == dist_from_top:
#       return dist_from_top * self.axis
    
#     if min_surf_dist == dist_from_bot:
#       return dist_from_bot * -self.axis

#     if min_surf_dist == dist_from_inner_wall:
#       return (self.innerRadius/r - 1) * xr

#     if min_surf_dist == dist_from_outer_wall:
#       return (self.outerRadius/r - 1) * xr

#   def getVelocity(self,pt):
#     x = np.asarray(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     if ( r > 1.e-12):
#       er = xr/r
#       v = self.v + self.radialVelocity*er
#     else:
#       v = self.v

#     return v

#   def xMin(self):
#     return min(self.x1[0], self.x2[0]) - max(self.innerRadius,self.outerRadius) 

#   def xMax(self):
#     return max(self.x1[0], self.x2[0]) + max(self.innerRadius,self.outerRadius) 

#############################################
class whiskers(Geometry):
  """
  Geometry object for creating a cylinder defined by two points and radius
  """
  def __init__(self,
               name,
               x1,
               x2,
               r,
               numWhiskers,
               volFracWhiskers,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.x1 = np.asarray(x1)
    self.x2 = np.asarray(x2)
    self.r = r
    self.nw = numWhiskers
    self.VF = volFracWhiskers

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.e1 = (self.x2-self.x1)/self.h       # unit vector along cylinder axis

    # random vectors orthogonal to e1 and each other
    self.e2 = np.random.rand(3)
    self.e2 = self.e2 - np.dot(self.e2,self.e1)*self.e1
    self.e2 = self.e2/np.linalg.norm(self.e2)
    self.e3 = np.cross(self.e1,self.e2)

  def isInterior(self, pt, skinDepth):    
    x = np.asarray(pt)
    z = np.dot(x-self.x1,self.e1)                                   # z-coordinate of test point
    pMinusX1InPlane = (x-self.x1) - np.dot(x-self.x1,self.e1)*self.e1
    r = np.linalg.norm( pMinusX1InPlane )  # r coordinate of test point
  
    ex = np.dot( self.e2,pMinusX1InPlane )
    ey = np.dot( self.e3,pMinusX1InPlane )

    theta = np.arctan2(ex,ey)
    dtheta = 2.0*np.pi/self.nw

    if ( ( (z>=0 and z<self.h) and (r < self.r) ) and ( (theta%dtheta)/dtheta < self.VF) ):
      if ( (z<skinDepth or z>(self.h-skinDepth)) or (r > self.r-skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r


#############################################
class toroid(Geometry):
  """
  Geometry object for creating a circular toroid defined by a point, direction ring radius and revolved circle radius
  """
  def __init__(self,
               name,
               x0,
               n,
               r1,
               r2,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    assert len(x0) == 3, "x0 must have length of 3"
    for i in x0:
       assert isinstance(i, floatType), "x0 elements must be floats"
    assert len(n) == 3, "n must have length of 3"
    for i in n:
       assert isinstance(i, floatType), "n elements must be floats"
    assert len(r1) == 3, "r1 must have length of 3"
    for i in r1:
       assert isinstance(i, floatType), "r1 elements must be floats"
    assert len(r2) == 3, "r2 must have length of 3"
    for i in r2:
       assert isinstance(i, floatType), "r2 elements must be floats"

    self.x0 = np.array(x0) # center point
    self.n = np.array(n) # axes unit vector
    self.n = self.n / np.linalg.norm(self.n)
    self.r1 = r1         # ring radius
    self.r2 = r2         # revolved circle radius

  def isInterior(self,pt,skinDepth):    
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    r = np.linalg.norm(x - z*self.n)

    if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) < self.r2):
      if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) > self.r2 - skinDepth):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    xr = x - z*self.n
    r = np.linalg.norm(xr)
    xrr = x - self.r1 * xr / r
    xrr_norm = np.linalg.norm(xrr)

    return xrr / xrr_norm

  def getSurfacePosition(self,pt):
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    xr = x - z*self.n
    r = np.linalg.norm(xr)
    xrr = x - self.r1 * xr / r 
    xrr_norm = np.linalg.norm( xrr )

    return ( self.r2 / xrr_norm - 1) * xrr

  def xMin(self):
    return self.x0[0]-self.r1-self.r2

  def xMax(self):
    return self.x0[0]+self.r1+self.r2


#############################################
class spinodal(Geometry):
  """
  Geometry object for generating arbitrary tailorable periodic spinodals microstructures using wavelet method

  Because N must be large it is too expensive to generate spinodal directly. Instead distance field is generated and interpolated to 
  generate spinodal microstructure

  # rho = target relative density
  # seed  = random seed for generating wave vectors and phases
  # num_waves = number of waves should be large (1000-100000)
  # a = cell size specified as a scalar or vector [ax, ay, az]
  # shell = boolean determines if geometry is a shell or solid spinodal
  """
  def __init__(self, 
               name, 
               rho, 
               seed, 
               shell=False, 
               grid_size=[100,100,100], 
               a=[0.33333,0.33333,0.33333], 
               num_waves=1000,
               vel=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    
    self.rho = rho
    self.seed = seed
    self.a = np.array(a)
    self.cellSize = np.array([a,a,a])
    self.shell = shell
    self.num_waves = num_waves
    self.offset = np.array([0.0,0.0,0.0])
    if self.shell:
      rho_s = self.rho / 2 + 0.5
    else:
      rho_s = self.rho
    self.level_set = np.sqrt( 2.0 ) * scipy.special.erfinv(2.0*rho_s-1.0)
    
    print(self.rho, rho_s, self.level_set)

    self.A = np.sqrt( 2.0 ) / self.num_waves #np.sqrt( 2.0 / self.num_waves )  # Prefactor for level set function (avoid computing for every point)

    self.phase = []
    self.q = []
    for n in range(self.num_waves):
      self.q.append(np.divide(2*np.pi*random_direction(), self.a).tolist())
      self.phase.append(np.random.uniform(0.0, 2*np.pi))
    
    self.phase  = np.array(self.phase)
    self.q = np.array(self.q)
    print(self.q.shape, self.phase.shape)

  def fphi(self, u):
    c = self.A*np.sum(np.cos(np.inner(self.q, u[:,np.newaxis].T)+self.phase))
    # print(c, self.A, np.sum(np.cos(np.inner(self.q, u[:,np.newaxis].T)+self.phase)), self.level_set, np.inner(self.q, u[:,np.newaxis].T).shape)
    
    return c

  def fdphi(self, u):
    c = 0
    for n in range(self.num_waves):
      c = c - self.kappa * np.sin(self.kappa * np.dot(self.k[n], u) + self.phase[n])
    
    return (1.41421356237 * c / self.num_waves) 

  # def pt2coords(self, pt):
  #   return 2*np.pi*np.divide(np.array(pt) - self.offset, self.cellSize) 

  def vec2Surf(self, pt):
    # Currently seem to work with gyroid may need to be adjusted later
    tolerance = 0.01
    maxIter = 1000
    alpha = 0.1
    
    c = self.level_set
    p = pt
    phi = self.fphi(p)
    dphi = self.fdphi(p)
    # print("i:", 0, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
    for i in range(maxIter):
      s = 1
      if phi <= 0:
        c = -abs(c)
        s = 0
      dpt = (-1)**(1+s)*alpha*abs(phi - c)*dphi
      phi = self.fphi(p)
      dphi = self.fdphi(p)
      p = p + dpt
      # print("i:", i, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
      if abs(phi - c) < tolerance:
        # print("Converged!")
        return p - pt
      
    print('Did not converge!')

  def isInterior(self, pt, skinDepth):
    # u = self.pt2coords(pt)
    u = np.array(pt)

    # Use level set method
    if self.shell:
      if abs(self.fphi(u)) < self.level_set:
        return True
    else:
      # print(u, self.fphi(u), self.level_set)
      if self.fphi(u) < self.level_set:
        return True

    return False

  def getSurfaceNormal(self, pt):
    # u = self.pt2coords(pt)
    # v = self.vec2Surf(u)
    # s = self.fdphi(u+v)
    # if self.fphi(u) <= 0:
    #   s = -s
    # return s / np.linalg.norm(s)
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self, pt):
    # u = self.pt2coords(pt)
    # v = self.vec2Surf(u)
    # return np.multiply(v, self.cellSize)/(2*np.pi)
    return np.array([0.0,0.0,0.0])


#############################################
class tpms(Geometry):
  """
  Geometry object for creating grid aligned triply periodic minimal surface

  # tpms_type = type of triply periodic minimal surface (TPMS) or spinodal
  # rho = target relative density
  # cellSize = array of size 3 that determines cellSize in x, y, z optionally also a single scalar
  # offset = array of size 3 to offset the cell of the TPMS
  """
  def __init__(self,
               name,
               tpms_type,
               rho=0.5,
               cellSize=[1.0,1.0,1.0],
               offset=[0.0,0.0,0.0],
               vel=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.tpms_type = tpms_type
    self.rho = rho
    self.cellSize= np.array(cellSize)
    self.offset = np.array(offset)

    # When python3 is updated to 3.10, we can use match instead of ifs
    if self.tpms_type == "gyroid":
      self.level_set = np.interp(self.rho, _gyroid_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.sin(u[0])*np.cos(u[1])+np.sin(u[1])*np.cos(u[2])+np.sin(u[2])*np.cos(u[0])
      self.dphi = lambda u : np.array([np.cos(u[0])*np.cos(u[1])-np.sin(u[2])*np.sin(u[0]),
                                                -np.sin(u[0])*np.sin(u[1])+np.cos(u[1])*np.cos(u[2]),
                                                -np.sin(u[1])*np.sin(u[2])+np.cos(u[2])*np.cos(u[0])])
      return

    if self.tpms_type == "shwarz diamond":
      self.level_set = np.interp(self.rho, _schwarzDiamond_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.cos(u[0])*np.cos(u[1])*np.cos(u[2]) - np.sin(u[0])*np.sin(u[1])*np.sin(u[2])
      self.dphi = lambda u : np.array([-np.sin(u[0])*np.cos(u[1])*np.cos(u[2]) - np.cos(u[0])*np.sin(u[1])*np.sin(u[2]),
                                        -np.cos(u[0])*np.sin(u[1])*np.cos(u[2]) - np.sin(u[0])*np.cos(u[1])*np.sin(u[2]),
                                        -np.cos(u[0])*np.cos(u[1])*np.sin(u[2]) - np.sin(u[0])*np.sin(u[1])*np.cos(u[2])])
      return

    if self.tpms_type == "schwarz primitive":
      self.level_set = np.interp(self.rho, _schwarzPrimitive_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.cos(u[0]) + np.cos(u[1]) + np.cos(u[2])
      self.dphi = lambda u : np.array([-np.sin(u[0]),
                                        -np.sin(u[1]),
                                        -np.sin(u[2])])
      return
    
    if self.tpms_type == "schwarz hexagonal":
      print("schwarz hexagonal tpms not currently implemented!")
      return
    
    if self.tpms_type == "schwarz clp":
      print("schwarz hexagonal tpms not currently implemented!")
      return
        
    print('No TPMS of type specified!')

  def pt2coords(self, pt):
    return 2*np.pi*np.divide(np.array(pt) - self.offset, self.cellSize) 

  def vec2Surf(self, pt):
    # Currently seem to work with gyroid may need to be adjusted later
    tolerance = 0.01
    maxIter = 1000
    alpha = 0.1
    
    c = self.level_set
    p = pt
    phi = self.phi(p)
    dphi = self.dphi(p)
    # print("i:", 0, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
    for i in range(maxIter):
      s = 1
      if phi <= 0:
        c = -abs(c)
        s = 0
      dpt = (-1)**(1+s)*alpha*abs(phi - c)*dphi
      phi = self.phi(p)
      dphi = self.dphi(p)
      p = p + dpt
      # print("i:", i, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
      if abs(phi - c) < tolerance:
        # print("Converged!")
        return p - pt
      
    print('Did not converge!')

  def isInterior(self, pt, skinDepth):
    u = self.pt2coords(pt)

    # Use level set method
    if abs(self.phi(u)) < self.level_set:
      v = self.vec2Surf(u)
      vMag = np.linalg.norm(np.multiply(v, self.cellSize))/(2*np.pi)

      if vMag < skinDepth:
        return _defaultSurfaceFlag

      return 0

    return -1

  def getSurfaceNormal(self, pt):
    u = self.pt2coords(pt)
    v = self.vec2Surf(u)
    s = self.dphi(u+v)
    if self.phi(u) <= 0:
      s = -s
    return s / np.linalg.norm(s)

  def getSurfacePosition(self, pt):
    u = self.pt2coords(pt)
    v = self.vec2Surf(u)
    return np.multiply(v, self.cellSize)/(2*np.pi)


class polygonInclusions(Geometry):
  """
  Generates polygon inclusions with surface flags
  points in polygons should be arrange in CCW order
  """
  def __init__(self,
               name,
               inclusions,
               vel=_defaultVelocity,
               mat=_defaultMat,
               fillGroup=_defaultGroup,
               inclusionGroup=_defaultGroup):
    super().__init__(name,
                     mat=mat,
                     group=fillGroup)
    self.fillGroup = fillGroup
    self.inclusionGroup = inclusionGroup

    self.polygons = []
    self.num_polygons = len(inclusions)
    for p in range(self.num_polygons):
      self.polygons.append(polygon(str(p), inclusions[p]))
  
  def isInterior(self, pt, skinDepth):
    for p in range(self.num_polygons):
      d = self.polygons[p].signedDist(pt)
      if np.abs(d) < skinDepth:
        return _defaultSurfaceFlag

    return 0

  def getGroup(self, pt):
    for p in range(self.num_polygons):
      d = self.polygons[p].signedDist(pt)
      if d >= 0.0:
        return self.inclusionGroup
    
    return self.fillGroup


#############################################
class fill(Geometry):
  """
  Geometry object for creating particles to fill the background.
  """
  def __init__(self,
               name,
               mat=_defaultMat,
               group=_defaultGroup):
    super().__init__(name,
                     mat=mat,
                     group=group)
    print('Warning: Fill should only be used with untransformed objects.')

  def isInterior(self, pt, skinDepth):
    return 0

  def getSurfaceNormal(self, pt):
    return np.array([0.0, 0.0, 0.0])

  def getSurfacePosition(self, pt):
    return np.array([0.0, 0.0, 0.0])


#############################################
class VCCTL(Geometry):
  """
  Geometry object for creating a 3D object from a VCCTL voxelized dataset

  data        text file with header stripped, containing one line per voxel and integer
              values indicating phase
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group

  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               nk,
               x0,
               x1,
               map,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):    
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0, x1)
    self.x1 = np.maximum(x0, x1)

  def isInterior(self,pt,skinDepth):
    x = np.array(pt)
    if np.all( np.logical_and( x > self.x0, x < self.x1 ) ):
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >-1):
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
        
        # check if point is next to porosity      
        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )
        kmin = int( max(0,k-1) )
        kmax = int( min(self.nk-1,k+1) )
        for i in range(imin,imax):
          for j in range(jmin,jmax):
            for k in range(kmin,kmax):
              n = i*self.nj*self.nk + j*self.nk + k
              index = self.data[n]
              mat = self.map[index] # will error if map doesn't contain index as key

              if (mat < 0):
                return _defaultSurfaceFlag

        return 0 # interior to the domain and not porosity
      
      return -1  # internal porosity
    
    return -1 # Out of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = self.data[n]
    mat = self.map[index]
    return mat # will error if map doesn't contain index as key

#############################################
class GrainStack(Geometry):
  def __init__(self, 
               name, 
               ct_array, 
               x0, 
               x1, 
               vel=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup, 
               particleType=_defaultParticleType, 
               flaggedSurfaces=None): 
    super().__init__(name, 
                     vel = vel, 
                     mat = mat, 
                     group = group, 
                     particleType = particleType) 
    # flaggedSurfaces = [bottom, top, outside walls, inside walls] 
    # ct_array is a numpy array with dimensions equal to num voxels in x, y, and z 
    # values indicate is void or material 
    # x0, x1 coordinates of the -,+ corners of the object in the domain [x-, y-, z-] 
    # mat default material type 
    # group contact group 
    # vel velocity vector 
    # damage scalar damage value for material 
    # ni, nj, nk number of voxels in the x,y, and z directions. 

    self.flaggedSurfaces = flaggedSurfaces if flaggedSurfaces is not None else [True, True, True, False] 
    type_check_array(x0, "x0", 3, floatType) 
    type_check_array(x1, "x1", 3, floatType) 
    self.name = name 
    self.type = 'GrainStack' 
    self.ct_array = ct_array 
    self.ni, self.nj, self.nk = ct_array.shape 
    self.x0 = np.array(x0) 
    self.x1 = np.array(x1) 
    self.vel = vel if callable(vel) else np.array(vel) 
    self.mat = mat 
    self.group = group 

    self.dx = (self.x1[0]-self.x0[0])/(self.ni - 1.0 ) if self.ni > 1 else 1.0 
    self.dy = (self.x1[1]-self.x0[1])/(self.nj - 1.0 ) if self.nj > 1 else 1.0 
    self.dz = (self.x1[2]-self.x0[2])/(self.nk - 1.0 ) if self.nk > 1 else 1.0

  def inBoundingBox(self, pt, eps=1e-8):
    pt = np.array(pt)
    return np.all(pt >= self.x0 - eps) and np.all(pt <= self.x1 + eps)

  def castToIJK(self, pt):
    if not self.inBoundingBox(pt):
        raise ValueError(f"Point {pt} is outside the voxel bounding box.")
    i = int(np.floor((pt[0] - self.x0[0]) / self.dx))
    j = int(np.floor((pt[1] - self.x0[1]) / self.dy))
    k = int(np.floor((pt[2] - self.x0[2]) / self.dz))
    return [i, j, k]

  def isInterior(self, pt, skinDepth):
    if not self.inBoundingBox(pt):
        return -1

    i, j, k = self.castToIJK(pt)
    
    # Check if current voxel is void
    if self.ct_array[i, j, k] == 0:
        return -1

    material_id = self.ct_array[i, j, k]
    
    # Compute number of voxels to check around the point in each direction
    di = max(int(np.ceil(skinDepth / self.dx)), 1)
    dj = max(int(np.ceil(skinDepth / self.dy)), 1)
    dk = max(int(np.ceil(skinDepth / self.dz)), 1)

    # Get bounds of local neighborhood, clamped to valid index range
    i_min = max(i - di, 0)
    i_max = min(i + di + 1, self.ni)
    j_min = max(j - dj, 0)
    j_max = min(j + dj + 1, self.nj)
    k_min = max(k - dk, 0)
    k_max = min(k + dk + 1, self.nk)

    # Check neighbors
    for ii in range(i_min, i_max):
        for jj in range(j_min, j_max):
            for kk in range(k_min, k_max):
                neighbor_id = self.ct_array[ii, jj, kk]
                if neighbor_id != material_id:
                    return _defaultSurfaceFlag

    return 0

#############################################
class CT(Geometry):
  """
  Geometry object for creating a 3D object from a CT voxelized dataset
  
  data are generated by a mathematica script, flattening a 3D array with
  ordering [z,y,x]

  This is similar (possibly identical) to the VCCTL object, but kept separate to allow
  for specialization without compromising backwards compatibility.
  data        text file with header stripped, containing one line per voxel and integer
              values indicating phase
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group
    
  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               nk,
               x0,
               x1,
               map,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self,pt, skinDepth):
    x = np.array(pt)
    
    if np.all( np.logical_and( x > self.x0, x < self.x1 ) ):
      i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )

      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >= 0): # interior to the domain and not porosity
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
    
        # check if point is next to porosity
        i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
        j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
        k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
        i=max(min(self.ni-1,i),0)
        j=max(min(self.nj-1,j),0)
        k=max(min(self.nk-1,k),0)
        
        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )
        kmin = int( max(0,k-1) )
        kmax = int( min(self.nk-1,k+1) )
        poreSurface = False
        for i in range(imin,imax):
          for j in range(jmin,jmax):
            for k in range(kmin,kmax):
              n = i*self.nj*self.nk + j*self.nk + k
              index = int(self.data[n])
              mat = self.map.get(index) # will error if map doesn't contain index as key
              if (mat <0):
                return _defaultSurfaceFlag
        return 0   
      
      return -1  # internal porosity
    
    return -1    # outside of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])

    i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = int(self.data[n])
    mat = self.map.get(index)
    return mat # will error if map doesn't contain index as key


#############################################
class bitmap(Geometry):
  """
  Geometry object for creating a 2D object from an image file

  data        binary array
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group

  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               x0,
               x1,
               map,
               vel=_defaultVelocity,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self, pt, skinDepth):
    x = np.array(pt)
    if np.all( np.logical_and(x >= self.x0 , x < self.x1) ):
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      index = int(self.data[j,i])

      if (index >0):
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
        
        # check if point is next to porosity
        i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
        j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )

        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )

        for i in range(imin,imax):
          for j in range(jmin,jmax):
            index = int(self.data[j,i])
            if (self.data[n]==0):
              return _defaultSurfaceFlag
        
        return 0   # interior to the domain and not porosity
    
    return -1   

  def getMat(self,pt):
    x = np.array(pt)
    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    #n = i*self.nj + j
    index = int(self.data[j,i])

    #index = self.data[n]
    return self.map[index] # will error if map doesn't contain index as key


# ===========================================
# END GEOMETRY OBJECTS
# ===========================================


# ===========================================
# GEOMETRY WRAPPERS 
# (Modify other geometry objects as inputs)
# ===========================================

#############################################
class BaseWrapper:
  @abstractmethod
  def __init__(self,
               name,
               subObject):
    self.name = name
    self.subObject = subObject
    self.particleType = subObject.particleType

  def isInterior(self, pt, skinDepth):
    return self.subObject.isInterior(pt, skinDepth)

  def getSurfaceNormal(self, pt):
    if hasattr(self, "surfaceNormal"):
        return self.surfaceNormal(self, pt) if callable(self.surfaceNormal) else self.surfaceNormal

    if hasattr(self.subObject, 'getSurfaceNormal'):
        return self.subObject.getSurfaceNormal(pt)

    if hasattr(self.subObject, 'surfaceNormal'):
        return self.subObject.surfaceNormal

    return _defaultSurfaceNormal
    
  def getSurfacePosition(self,pt):
    if hasattr(self, "surfacePosition"):
        return self.surfacePosition(self, pt) if callable(self.surfacePosition) else self.surfacePosition

    if hasattr(self.subObject, 'getSurfacePosition'):
        return self.subObject.getSurfacePosition(pt)

    if hasattr(self.subObject, 'surfacePosition'):
        return self.subObject.surfacePosition
    
    return _defaultSurfacePosition

  def getStrengthScale(self,pt):
    if hasattr(self, "strengthScale"):
        return self.strengthScale(self, pt) if callable(self.strengthScale) else self.strengthScale
    
    if hasattr(self.subObject, 'getStrengthScale'):
        return self.subObject.getStrengthScale(pt)

    if hasattr(self.subObject, 'strengthScale'):
        return self.subObject.strengthScale
    
    return _defaultStrengthScale

  def getVelocity(self,pt):
    if hasattr(self, "vel"):
        return self.vel(self, pt) if callable(self.vel) else self.vel

    if hasattr(self.subObject, 'getVelocity'):
        return self.subObject.getVelocity( pt )
    else:
        return self.subObject.vel( self.subObject, pt ) if callable( self.subObject.vel ) else self.subObject.vel

    return _defaultVelocity

  def getGroup(self, pt):
    if hasattr(self, "group"):
        return self.group(self, pt) if callable(self.group) else self.group
    
    if hasattr(self.subObject, 'getGroup'):
        return self.subObject.getGroup(pt)

    if hasattr(self.subObject, 'group'):
        return self.subObject.group

    return _defaultGroup

  def getMat(self, pt):
    if hasattr(self, "mat"):
        return self.mat(self, pt) if callable(self.mat) else self.mat
    
    if hasattr(self.subObject, 'getMat'):
        return self.subObject.getMat(pt)

    if hasattr(self.subObject, 'mat'):
        return self.subObject.mat

    return _defaultMat

  def getMatDir(self, pt):
    if hasattr(self, "matDir"):
        return self.matDir(self, pt) if callable(self.matDir) else self.matDir
    
    if hasattr(self.subObject, 'getMatDir'):
        return self.subObject.getMatDir(pt)

    if hasattr(self.subObject, 'matDir'):
        return self.subObject.matDir

    return _defaultMatDir

  def getDamage(self, pt):
    if hasattr(self, "damage"):
        return self.damage(self, pt) if callable(self.damage) else self.damage
    
    if hasattr(self.subObject, 'getDamage'):
        return self.subObject.getDamage(pt)

    if hasattr(self.subObject, 'damage'):
        return self.subObject.damage

    return _defaultDamage

  def getPorosity(self, pt):
    if hasattr(self, "porosity"):
        return self.porosity(self, pt) if callable(self.porosity) else self.porosity
    
    if hasattr(self.subObject, 'getPorosity'):
        return self.subObject.getPorosity(pt)

    if hasattr(self.subObject, 'porosity'):
        return self.subObject.porosity

    return _defaultPorosity

  def getTemperature(self, pt):
    if hasattr(self, "temperature"):
        return self.temperature(self, pt) if callable(self.temperature) else self.temperature
    
    if hasattr(self.subObject, 'getTemperature'):
        return self.subObject.getTemperature(pt)

    if hasattr(self.subObject, 'temperature'):
        return self.subObject.temperature

    return _defaultTemperature

  def getSurfaceTraction(self, pt):
    if hasattr(self, "surfaceTraction"):
        return self.surfaceTraction(self, pt) if callable(self.surfaceTraction) else self.surfaceTraction
    
    if hasattr(self.subObject, 'getSurfaceTraction'):
        return self.subObject.getSurfaceTraction(pt)

    if hasattr(self.subObject, 'surfaceTraction'):
        return self.subObject.surfaceTraction

    return _defaultSurfaceTraction

  def getCZTag(self, pt):
    if hasattr(self, "czTag"):
        return self.czTag(self, pt) if callable(self.czTag) else self.czTag
    
    if hasattr(self.subObject, 'getCZTag'):
        return self.subObject.getCZTag(pt)

    if hasattr(self.subObject, 'czTag'):
        return self.subObject.czTag

    return _defaultCZTag

  def xMin(self):
    return self.subObject.xMin()

  def xMax(self):
    return self.subObject.xMax()
  
  def getSubregions(self):
    return self.subObject.getSubregions()
  

#############################################
class materialDirectionWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               matDir):
    super().__init__(name,
                     subObject)

    if np.shape(matDir) == (3,):
      # matDir will be a row-wise 3x3 where the first row is u1=matDir_3x1, and rows 2 and 3 are orthogonal.
      u1, u2, u3 = generate_orthonormal_axes(matDir)
      self.matDir = np.vstack((u1, u2, u3))
    elif np.shape(matDir) == (3, 3):
      self.matDir = np.asarray(matDir)
    else:
      self.matDir = _defaultMatDir
      log2file("Unsupported material direction size...")


############################################
class surfaceFlagWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               flagType):
    super().__init__(name,
                     subObject)
    self.flagType = flagType

  def isInterior(self, pt, skinDepth):
    flag = self.subObject.isInterior(pt, skinDepth)
        
    if flag > 0:
      return self.flagType

    return flag


############################################
class czTagWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               czTag):
    super().__init__(name,
                     subObject)
    self.czTag = czTag
  
  def getCZTag( self, pt ):
    return self.czTag( pt) if callable(self.czTag) else self.czTag
    

class functionalDeletionWrapper(BaseWrapper):
  def __init__(self, name, subObject, deleteFunction):
    super().__init__(name, subObject)
    self.deleteFunction = deleteFunction  # A function that takes pt and returns -1 or 0

  def isInterior(self, pt, skinDepth):
    # Let subObject decide initial flag
    flag = self.subObject.isInterior(pt, skinDepth)

    # Only consider interior or surface points
    if flag >= 0:
      # Use function to determine if we override flag
      delete_flag = self.deleteFunction(pt)
      if delete_flag == -1:
        return -1  # Delete the point
    return flag


############################################
class surfaceNormalWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               surfaceNormal):
    super().__init__(name,
                     subObject)
    self.surfaceNormal = np.asarray(surfaceNormal)


############################################
class surfacePositionWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               surfacePosition):
    super().__init__(name,
                     subObject)
    self.surfacePosition = np.asarray(surfacePosition)


############################################
class shrinkageFlagWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               flag):
    super().__init__(name,
                     subObject)
    self.flag = flag


############################################
class strengthScaleWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               strengthScale):
    super().__init__(name,
                     subObject)
    self.strengthScale = strengthScale


############################################
class damageWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               damage):
    super().__init__(name,
                     subObject)
    self.damage=damage


############################################
class porosityWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               porosity):
    super().__init__(name,
                     subObject)
    self.porosity=porosity


############################################
class temperatureWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               temperature):
    super().__init__(name,
                     subObject)
    self.temperature = temperature


class voronoiWrapperWithLayeredStrengthScale(BaseWrapper):
  """
  Voronoi tessellation + layered strengthScale inside each cell.

  strengthScale(pt) = cellSliceStrengthScale[cellId(pt)][layerId(pt)]
  """

  def __init__(self,
               name,
               subObject,
               x0,
               x1,
               flawSize,
               weibullVolume,
               weibullMinVolume,
               weibullReferenceVolume,
               weibullModulus,
               weibullSeed: int,
               weibullLayerThickness,
               nIntegrationPoints=10,
               matDir=np.array([1,0,0]),
               bondedSurfaceFraction=1.0,
               vpts=None,
               dim: int = 3,
               randomMatDir: bool = False):

    super().__init__(name, subObject)
    self.object = subObject
    self.dim = dim

    self.bondedSurfaceFraction = bondedSurfaceFraction

    # ---- store inputs you use later ----
    self.seed = int(weibullSeed)
    self.weibullLayerThickness = float(weibullLayerThickness)
    self.nIntegrationPoints = int(nIntegrationPoints)

    # Keep these around (used in slice-strength computation)
    self.weibullMinVolume = float(weibullMinVolume)
    self.weibullReferenceVolume = float(weibullReferenceVolume)
    self.weibullModulus = float(weibullModulus)

    self.randomMatDir = bool(randomMatDir)

    # ---- bounding box used for generating Voronoi seeds ----
    self.x0_full = np.array(x0, dtype=float)
    self.x1_full = np.array(x1, dtype=float)
    self.x0 = self.x0_full[:self.dim]
    self.x1 = self.x1_full[:self.dim]
    self.dx = (self.x1_full - self.x0_full)

    # ---- global layering direction (fallback) ----
    if matDir is None and (not self.randomMatDir):
      matDir = np.array([0.0, 0.0, 1.0])
    self.globalMatDir = np.array(matDir, dtype=float)
    n = np.linalg.norm(self.globalMatDir)
    if n <= 0.0:
      raise ValueError("matDir must be nonzero")
    self.globalMatDir /= n

    # ---- build Voronoi seeds ----
    if vpts is None:
      vpts_raw = poisson(flawSize, x0=self.x0, dx=(self.x1 - self.x0), seed=self.seed, dim=self.dim)
    else:
      vpts_raw = np.array(vpts, dtype=float)

    self.vpts = np.asarray(vpts_raw, dtype=float)[:, :self.dim]
    self.npts = self.vpts.shape[0] #number of voronoi points and soon to be cells


    # ---- nearest-site lookup ----
    self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(self.npts / 2)), metric='euclidean')

    # ---- Voronoi diagram ----
    self.voronoi = Voronoi(self.vpts)
    vor = self.voronoi
 
    # ---- cell volume + centroid + radius ----
    v0 = np.prod(self.x1 - self.x0) / self.npts # default volume for unbounded and edge case cells
    # set up empty arrays for volume, centroid, radius attributes
    ################# self.cellVolume  = np.zeros(vor.npoints, dtype=float)
    self.cellCentroid = np.zeros((vor.npoints, self.dim), dtype=float)   
    self.cellRadius  = np.zeros(vor.npoints, dtype=float)


    # iterate through all voronoi cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
  
      # for the unbounded cells 
      if (len(indices) == 0) or (-1 in indices):
          # volume fallback
          ################# self.cellVolume[i] = v0
  
          # assign cellCentroid to be the current seed point
          self.cellCentroid[i, :self.dim] = vor.points[i, :self.dim]
          
          #original did not have a way to get cell radius for unbounded cells
          k_neigh = min(6, self.npts)  #ask kd tree for the 6 closest neighboring points. can adjust 6 as needed
          dists, _ = self.kdt.query(vor.points[i, :self.dim].reshape(1, -1), k=k_neigh) # #ask kdt for  distances to those points - automatically sorts from closest to furthest
          R = float(dists[0, -1]) #take the furthest of the neighpors and make that distance a floate
          self.cellRadius[i] = max(R, 1e-12) # this is a safety in case there is only one seed
  
          continue
  
      #for bounded cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      verts = vor.vertices[indices]  # safe because no -1
  
      # If we have an edge case, like a cell in a corner:
      if (self.dim == 3) and (verts.shape[0] < 4):
          ################# self.cellVolume[i] = v0
          
          # "centroid of the indices" interpretted as mean of the available vertices:
          self.cellCentroid[i, :self.dim] = np.mean(verts[:, :self.dim], axis=0)
          
          # get radius based on mean of available verts (or 0 if somehow empty, but it shouldn't be here)
          # subtract the centroid from each vertex
          # then we get vectors from the centroid
          # get distances  of each centroid-vertex vector.
          # take the largest and assign to the cell 'radius'
          # and if no vertices (perhaps if there is one cell for the entire object), make radius 0 to not crash program
          self.cellRadius[i] = float(np.max(np.linalg.norm(
              verts[:, :self.dim] - self.cellCentroid[i, :self.dim], axis=1
          ))) if verts.shape[0] > 0 else 0.0
          
          # 
          self.cellRadius[i] = max(self.cellRadius[i], 1e-12) # this is a safety in case there is only one seed
          continue
  
      # Normal bounded cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # in the original, we have volume for only  the unbounded cells. 
      # adding in volume calc for bounded cells too
      # work from sizeEffectVoronoiWeibullBox
      
      ##############self.cellVolume[i] = max(ConvexHull(verts).volume, vMin)
      
      #compute centroid from the vertices
      c = self.voronoi_cell_centroid(verts)  # should return (2,) or (3,)
      #assign to object attribute
      self.cellCentroid[i, :self.dim] = c[:self.dim]
  
      # get distances  of each centroid-vertex vector.
      # take the largest and assign to the cell 'radius'
      self.cellRadius[i] = float(np.max(np.linalg.norm(
          verts[:, :self.dim] - self.cellCentroid[i, :self.dim], axis=1
      )))
      self.cellRadius[i] = max(self.cellRadius[i], 1e-12)
    
    self.nCentroidsInside = int(np.sum([subObject.isInterior(c, 0.0) >= 0
                                   for c in self.cellCentroid]))

    # This should define a surface flag for each cell face based on the specified fraction of bonded surfaces.
    self.ridgePtFlags = [ (3 if ( np.random.random() <= self.bondedSurfaceFraction ) else 2 ) for p in self.voronoi.ridge_points ]

    print("Centroids inside cylinder:", self.nCentroidsInside, "out of", self.npts)
               
    # ---- per-cell material direction ----
    self.cellMatDir = None
    self.cellMatDirMat = []
    if self.randomMatDir:
      cellMatDir = []
      for i in range(self.npts):
        d = random_direction(dim=self.dim)
        if self.dim == 2:
          d = np.append(d, 0.0)
        d = np.asarray(d, dtype=float)
        d /= max(np.linalg.norm(d), 1e-15)
        u1, u2, u3 = generate_orthonormal_axes(d)
        self.cellMatDirMat.append(np.vstack((u1, u2, u3)))
        cellMatDir.append(d)
      self.cellMatDir = np.asarray(cellMatDir, dtype=float)
      self.cellMatDirMat = np.asarray(self.cellMatDirMat, dtype=float)
    else:
      # If you don't want random per-cell directions, use globalMatDir for all cells
      self.cellMatDir = np.tile(self.globalMatDir, (self.npts, 1))

    # ---- number of slices per cell ----
    self.cellNumSlices = np.zeros(self.npts, dtype=int)
    for i in range(self.npts):
       self.cellNumSlices[i] = max(
        1,
        int(np.ceil(2.0 * float(self.cellRadius[i]) / float(self.weibullLayerThickness)))
      )


    # ---- compute slice strength scales for each cell (ragged) ----
    self.cellSliceStrengthScale = [None] * self.npts 
    rng = np.random.default_rng(self.seed)

    cellVolsFromSlices = []
    for i in range(self.npts): 
      m1 = np.asarray(self.cellMatDir[i], dtype=float)

      # construct m2, m3 spanning plane normal to m1
      if abs(np.dot(m1, np.array([0.0, 0.0, 1.0])) - 1.0) < 1e-12:
        m2 = np.cross(np.array([0.0, 1.0, 0.0]), m1)
      else:
        m2 = np.cross(np.array([0.0, 0.0, 1.0]), m1)
      m2 = m2 / max(np.linalg.norm(m2), 1e-15)
      m3 = np.cross(m1, m2)

      nSlices_i = int(self.cellNumSlices[i]) 
      Rcell = float(self.cellRadius[i])

      # If cellRadius is 0 (edge/unbounded fallback), still define 1 slice value
      if Rcell <= 0.0:
        R = rng.uniform(1e-20, 1.0)
        s = ((self.weibullReferenceVolume / self.weibullMinVolume) *
             (np.log(R) / np.log(0.5))) ** (1.0 / self.weibullModulus)
        self.cellSliceStrengthScale[i] = np.array([s], dtype=float)
        continue

      slicez = np.linspace(-Rcell, Rcell, nSlices_i)
      points = np.linspace(-Rcell, Rcell, self.nIntegrationPoints)

      sliceStrengthScale = np.zeros(nSlices_i, dtype=float) 

 
      layerVolsList=[]
      for j in range(nSlices_i):        
        p0 = self.cellCentroid[i] + slicez[j] * m1 #get center origin of slice

        count = 0
        for r2 in points:
          for r3 in points:
            pt = p0 + r2 * m2 + r3 * m3 #test point on the plane
            # if point is interior to subobject and is nearest to the voronoi cell i
            if (subObject.isInterior(pt, 0.0) >= 0) and (self._cellId(pt) == i):
              count += 1
        # below, we have a grid of nintegration points, hence squaring of nIntegration points 
        # the area of the slice is approximated as (2*Rcell)^2
        layerVol = max(
          self.weibullMinVolume,
          self.weibullLayerThickness *
          (count / float(self.nIntegrationPoints**2)) *
          (2.0 * Rcell)**2
        )
        
        layerVolsList.append(layerVol)

        R = rng.uniform(1e-20, 1.0)
        s = ((self.weibullReferenceVolume / layerVol) *
             (np.log(R) / np.log(0.5))) ** (1.0 / self.weibullModulus)

        sliceStrengthScale[j] = s  

      cellVolsFromSlices.append(np.sum(layerVolsList))
      self.cellSliceStrengthScale[i] = sliceStrengthScale
      
      

    ## after all cells are processed:
    self.cellVolumeFromSlices = np.asarray(cellVolsFromSlices, dtype=float)
    #
    ## per-cell "reference radius" = radius of equal-volume sphere
    self.cellRefRadius = (3.0 * self.cellVolumeFromSlices / (4.0 * np.pi)) ** (1.0/3.0)
    
    ## reference radius of all cells
    ## average reference radius
    #self.avgRefRadius = float(np.mean(self.cellRefRadius))
    #print("avg reference radius:", self.avgRefRadius, "| n cells:", self.cellRefRadius.size)
    
    # reference radius for only interior cells
    inside = np.array([subObject.isInterior(c, 0.0) >= 0 for c in self.cellCentroid], dtype=bool)
    
    # avoid empty slice
    if np.any(inside):
        self.avgRefRadiusInside = float(np.mean(self.cellRefRadius[inside]))
    else:
        self.avgRefRadiusInside = float("nan")
    
    print("avg reference radius (inside):", self.avgRefRadiusInside)
  
  def isInterior(self, pt, skinDepth):
      # is the point within the object
      x = np.asarray(pt[:self.dim])

      surf_flag = self.subObject.isInterior(pt, skinDepth)
      if surf_flag == 0:
        # Find voronoi cell closest to point
        _, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]
        
        # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
        minSurfaceDist = np.inf
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
          p1 = ridgePts[i][0]
          p2 = ridgePts[i][1]
          if index == p1 or index == p2:
            if index == p1:
              p = p2
            else:
              p = p1

            n = self.vpts[p, :] - self.vpts[index, :]
            n = n / 2
            d = np.linalg.norm(n)
            n = n / d

            dv = x - self.vpts[index, :]
            dvc = np.dot(n, dv)  # component of points along voronoi face normal
            if dvc > 0.0 and dvc > d - skinDepth:
              return self.ridgePtFlags[i]
      return surf_flag

  # ---------------------------
  # point -> nearest Voronoi site id
  def _cellId(self, pt) -> int:
      x = np.asarray(pt[:self.dim], dtype=float)
      dist, idx = self.kdt.query(x.reshape(1, -1), k=5)  # show top-5 nearest
  
      # print only a few times
      if not hasattr(self, "_dbg_cellid_calls"):
          self._dbg_cellid_calls = 0
      if self._dbg_cellid_calls < 20:
          self._dbg_cellid_calls += 1
  
      return int(idx[0][0])   # closest seed = Voronoi cell id

  # ---------------------------
  # centroid from vertices (ConvexHull-based; only use when verts form a real polytope)
  def voronoi_cell_centroid(self, vertices):
    vertices = np.asarray(vertices)
    dim = vertices.shape[1]
    if dim not in (2, 3):
      raise ValueError("Only 2D or 3D vertices are supported")
    hull = ConvexHull(vertices)
    if dim == 2:
      return self.centroid_2d(vertices, hull)
    else:
      return self.centroid_3d(vertices, hull)

  def centroid_2d(self, vertices, hull):
    pts = vertices[hull.vertices]
    x = pts[:, 0]; y = pts[:, 1]
    x_next = np.roll(x, -1); y_next = np.roll(y, -1)
    cross = x * y_next - x_next * y
    area = 0.5 * np.sum(cross)
    if np.isclose(area, 0):
      return np.mean(pts, axis=0)
    cx = np.sum((x + x_next) * cross) / (6 * area)
    cy = np.sum((y + y_next) * cross) / (6 * area)
    return np.array([cx, cy])

  def centroid_3d(self, vertices, hull):
    centroid = np.zeros(3)
    volume = 0.0
    origin = np.zeros(3)
    for simplex in hull.simplices:
      a, b, c = vertices[simplex]
      v = np.dot(a, np.cross(b, c)) / 6.0
      tetra_centroid = (origin + a + b + c) / 4.0
      centroid += v * tetra_centroid
      volume += v
    if np.isclose(volume, 0):
      return np.mean(vertices, axis=0)
    return centroid / volume

  # ---------------------------
  # PFW hooks
  def getMatDir(self, pt):
    cellId = self._cellId(pt)
    return self.cellMatDirMat[cellId]

  def getStrengthScale(self, pt):
    i = self._cellId(pt)
    m1 = self.cellMatDir[i]
    z = float(self.cellRadius[i] + np.dot(m1, pt - self.cellCentroid[i]))
    k = int(np.floor(z / self.weibullLayerThickness))
    k = max(0, min(int(self.cellNumSlices[i]) - 1, k))
    return float(self.cellSliceStrengthScale[i][k])

  def getDamage(self, pt):
    if hasattr(self.object, "getDamage"):
      return self.object.getDamage(pt)
    if hasattr(self.object, "damage"):
      dmg = self.object.damage(self.object, pt) if callable(self.object.damage) else self.object.damage
      return dmg
    return 0.0

############################################
class porousPolygon(Geometry):
  """
  Geometry object for creating a polygon described by ordered vertices with 2d pores defined by radius and porosity.'
  """
  def __init__(self,
               name,
               plist,
               porosity,
               poreSize,
               seed = 1,
               vel=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=None):
    super().__init__(name,
                     vel = vel,
                     mat = mat,
                     group = group,
                     particleType = particleType)
    self.plist = np.asarray(plist)
    self.porosity = porosity
    self.poreSize = poreSize
    self.seed = seed # random seed for porosity
    self.flaggedSurfaces = flaggedSurfaces
    if self.flaggedSurfaces is None:
      self.flaggedSurfaces = [True for i in range(self.plist.shape[0])]
      
    self.dim = 2
    
    # microPores:
    self.porosity = porosity
    self.poreSize = poreSize
    
    # percent area occupied by RCP uniform circles in a plane is .82-.84
    packingFraction = 0.82 
        
    self.poreRadiusSquared = 0.25*poreSize*poreSize
    # compute the effective disk-packing diameter to give the desired porosity and pore size.
    self.packingDiameter = poreSize*np.sqrt(packingFraction/porosity)
    
    self.dim = 2
    self.x0 = np.array([ min([ float(v[0]) for v in plist]),min([ float(v[1]) for v in plist]),min([ float(v[2]) for v in plist])])
    self.x1 = np.array([ max([ float(v[0]) for v in plist]),max([ float(v[1]) for v in plist]),max([ float(v[2]) for v in plist])])
    self.dx = self.x1-self.x0
        
    self.pores = poisson(self.packingDiameter, x0=self.x0, dx=self.dx[:self.dim], seed=self.seed, dim=self.dim)
    self.k = min(5, len(self.pores))
        
    self.kdt = KDTree(self.pores, leaf_size=int(np.ceil(len(self.pores) / 2)), metric='euclidean')

  def ccw(self,x1,x2,x3):
    return (x3[1]-x1[1])*(x2[0]-x1[0]) - (x2[1]-x1[1])*(x3[0]-x1[0]) > -10**-16

  def intersect(self,A,B,C,D):
    return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

  def isInside(self,vertices,point):
    v_arr = np.asarray(vertices)
    xmin = min(v_arr[:,0])
    ymin = min(v_arr[:,1])
    xmax = max(v_arr[:,0])
    ymax = max(v_arr[:,1])
    dx = xmax - xmin
    dy = ymax - ymin
    outside = [xmin,ymin]
    outside[0] -= dx
    outside[1] -= dy
    nv = v_arr.shape[0]
    nIntersections = 0
    for i in range(nv):
      p1 = vertices[i]
      p2 = vertices[np.mod(i+1,nv)]
      if(self.intersect(p1,p2,outside,point)):
        nIntersections += 1
    if np.mod(nIntersections,2)==0:
      return False
    else:
      return True

  def isInterior(self, point, skinDepth):
    x = np.asarray(point)
    
    # Is inside polygon
    if(self.isInside(self.plist,x)):
      
      # check if in a pore:
      _, index = self.kdt.query(x.reshape(1,-1), k=self.k)
      surfaceFlag = 0 # Interior unless otherwise determined
      for i in index[0]:
        p=self.pores[i]
        rSqr = ( x[0] - p[0] )**2 + ( x[1] - p[1] )**2  #+ ( x[2] - p[2] )**2
        if rSqr < self.poreRadiusSquared:
          return -1 # Inside a pore
        if rSqr < ( 0.5*self.poreSize + skinDepth )**2:
          surfaceFlag = _defaultSurfaceFlag # near internal pore surface
          
      # not a pore surface, check if polygon surface
      if surfaceFlag != _defaultSurfaceFlag:   
        vertices = self.plist
        nv = vertices.shape[0]

        nearestI = -1
        nearestEdgeD = np.inf
        for i in range(nv):
          if not self.flaggedSurfaces[i]:
            continue
          
          j = np.mod(i+1,nv)
          v = x - vertices[i,:]
          w = vertices[j,:] - vertices[i,:]

          wNorm = np.linalg.norm(w)
          w = w / wNorm
          dw = np.dot(v,w)
          if dw >= 0.0 and dw <= wNorm:
            v = v - dw * w
            d = np.linalg.norm(v)
            if d < nearestEdgeD:
              nearestEdgeD = d
              nearestI = i

          if nearestEdgeD < skinDepth:
            surfaceFlag = _defaultSurfaceFlag
      
      return surfaceFlag
    return -1

  def getSurfaceNormal(self,pt):
    x = np.asarray(pt)
    # Assumes point is already internal
    # Find the nearest edge and use it's normal
    vertices = self.plist
    nv = vertices.shape[0]
    nearestEdgeD = np.inf
    surfaceNormal = np.empty((1,3))
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfaceEdgeNormal = -v / d
    
    # check if in a pore:
    _, index = self.kdt.query(x.reshape(1,-1), k=self.k)
    nearestPoreD = np.inf # distance to nearest pore
    for i in index[0]:
      p=self.pores[i]
      rSqr = ( x[0] - p[0] )**2 + ( x[1] - p[1] )**2 #+ ( x[2] - p[2] )**2
      if rSqr < nearestPoreD*nearestPoreD:
        nearestPoreD = np.sqrt(rSqr)
        nearestPoreSurfaceD = nearestPoreD - 0.5*self.poreSize # distance to surface of nearest pore.
        surfacePoreNormal = np.array(p-x) / nearestPoreD
    
    # select closest surface, edge or pore
    if (nearestEdgeD <= nearestPoreSurfaceD):
      return surfaceEdgeNormal
    else:
      return surfacePoreNormal
  
  def getSurfacePosition(self,pt):
    x = np.asarray(pt)
    # Assumes point is already internal
    # Find the nearest edge and vector to surface
    vertices = self.plist
    nv = vertices.shape[0]
    
    nearestEdgeD = np.inf
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfaceEdgePosition = -v
          
    # check if in a pore:
    _, index = self.kdt.query(x.reshape(1,-1), k=self.k)
    nearestPoreD = np.inf # distance to nearest pore
    for i in index[0]:
      p=np.array(self.pores[i])
      rSqr = ( x[0] - p[0] )**2 + ( x[1] - p[1] )**2 #+ ( x[2] - p[2] )**2
      if rSqr < nearestPoreD*nearestPoreD:
        nearestPoreD = np.sqrt(rSqr) 
        surfacePoreNormal = np.array((p-x))/nearestPoreD
        nearestPoreSurfaceD = nearestPoreD - 0.5*self.poreSize # distance to surface of nearest pore.
        # relative position.
        surfacePorePosition = (p - np.array(surfacePoreNormal*0.5*self.poreSize)) - x
    
    # select closest surface, edge or pore
    if (nearestEdgeD <= nearestPoreSurfaceD):
      return surfaceEdgePosition
    else:
      return surfacePorePosition

  def xMin(self):
    arr = np.asarray(self.plist)
    xmin = min(arr[:,0])
    return xmin

  def xMax(self):
    arr = np.asarray(self.plist)
    xmax = max(arr[:,0])
    return xmax
  
############################################
class layeredVoronoiWeibullWrapper(BaseWrapper):
  """
  Voronoi tessellation + layered strengthScale inside each cell.

  strengthScale(pt) = cellSliceStrengthScale[cellId(pt)][layerId(pt)]
  """

  def __init__(self,
               name,
               subObject,
               x0,
               x1,
               flawSize,
               weibullReferenceVolume,
               weibullModulus,
               weibullSeed,
               minWeibullVolume,
               weibullLayerThickness,
               nIntegrationPoints=10,
               matDirs=None,
               vpts=None,
               dim: int = 3,
               randomMatDir=False):

    super().__init__(name, subObject)
    self.object = subObject
    self.dim = dim

    # ---- store inputs you use later ----
    self.seed = int(weibullSeed)
    self.weibullLayerThickness = float(weibullLayerThickness)
    self.nIntegrationPoints = int(nIntegrationPoints)

    # Keep these around (used in slice-strength computation)
    self.minWeibullVolume = float(minWeibullVolume)
    self.weibullReferenceVolume = float(weibullReferenceVolume)
    self.weibullModulus = float(weibullModulus)

    self.matDirs = matDirs # If None use subobject matDirs as passthrough
    self.randomMatDir = bool(randomMatDir)

    # ---- bounding box used for generating Voronoi seeds ----
    self.x0 = np.asarray(x0, dtype=float)
    self.x1 = np.asarray(x1, dtype=float)
    self.dx = (self.x1 - self.x0)
    self.x0 = self.x0[:self.dim]
    self.x1 = self.x1[:self.dim]
    
    # ---- build Voronoi seeds ----
    if vpts is None:
      vpts_raw = poisson(flawSize, x0=self.x0, dx=(self.x1 - self.x0), seed=self.seed, dim=self.dim)
    else:
      vpts_raw = np.array(vpts, dtype=float)

    self.vpts = np.asarray(vpts_raw, dtype=float)[:, :self.dim]
    self.npts = self.vpts.shape[0] #number of voronoi points and soon to be cells

    # ---- nearest-site lookup ----
    self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(self.npts / 2)), metric='euclidean')

    # ---- Voronoi diagram ----
    self.voronoi = Voronoi(self.vpts)
    vor = self.voronoi
 
    # ---- cell volume + centroid + radius ----
    v0 = np.prod(self.x1 - self.x0) / self.npts # default volume for unbounded and edge case cells
    # set up empty arrays for volume, centroid, radius attributes
    ################# self.cellVolume  = np.zeros(vor.npoints, dtype=float)
    self.cellCentroid = np.zeros((vor.npoints, self.dim), dtype=float)   
    self.cellRadius  = np.zeros(vor.npoints, dtype=float)

    # iterate through all voronoi cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
  
      # for the unbounded cells 
      if (len(indices) == 0) or (-1 in indices):
          # volume fallback
          ################# self.cellVolume[i] = v0
  
          # assign cellCentroid to be the current seed point
          self.cellCentroid[i, :self.dim] = vor.points[i, :self.dim]
          
          #original did not have a way to get cell radius for unbounded cells
          k_neigh = min(6, self.npts)  #ask kd tree for the 6 closest neighboring points. can adjust 6 as needed
          dists, _ = self.kdt.query(vor.points[i, :self.dim].reshape(1, -1), k=k_neigh) # #ask kdt for  distances to those points - automatically sorts from closest to furthest
          R = float(dists[0, -1]) #take the furthest of the neighpors and make that distance a floate
          self.cellRadius[i] = max(R, 1e-12) # this is a safety in case there is only one seed
  
          continue
  
      #for bounded cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      verts = vor.vertices[indices]  # safe because no -1
  
      # If we have an edge case, like a cell in a corner:
      if (self.dim == 3) and (verts.shape[0] < 4):
          ################# self.cellVolume[i] = v0
          
          # "centroid of the indices" interpretted as mean of the available vertices:
          self.cellCentroid[i, :self.dim] = np.mean(verts[:, :self.dim], axis=0)
          
          # get radius based on mean of available verts (or 0 if somehow empty, but it shouldn't be here)
          # subtract the centroid from each vertex
          # then we get vectors from the centroid
          # get distances  of each centroid-vertex vector.
          # take the largest and assign to the cell 'radius'
          # and if no vertices (perhaps if there is one cell for the entire object), make radius 0 to not crash program
          self.cellRadius[i] = float(np.max(np.linalg.norm(
              verts[:, :self.dim] - self.cellCentroid[i, :self.dim], axis=1
          ))) if verts.shape[0] > 0 else 0.0
          
          # 
          self.cellRadius[i] = max(self.cellRadius[i], 1e-12) # this is a safety in case there is only one seed
          continue
  
      # Normal bounded cells: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # in the original, we have volume for only  the unbounded cells. 
      # adding in volume calc for bounded cells too
      # work from sizeEffectVoronoiWeibullBox
      
      ##############self.cellVolume[i] = max(ConvexHull(verts).volume, vMin)
      
      #compute centroid from the vertices
      c = self.voronoi_cell_centroid(verts)  # should return (2,) or (3,)
  
      # get distances  of each centroid-vertex vector.
      # take the largest and assign to the cell 'radius'
      self.cellRadius[i] = float(np.max(np.linalg.norm(
          verts[:, :self.dim] - self.cellCentroid[i, :self.dim], axis=1
      )))
      self.cellRadius[i] = max(self.cellRadius[i], 1e-12)
    
    if self.dim == 2:
      self.cellCentroid = np.concatenate((self.cellCentroid, np.zeros((self.cellCentroid.shape[0], 1))), axis=1)

    self.nCentroidsInside = int(np.sum([subObject.isInterior( c , 0.0) >= 0
                                   for c in self.cellCentroid]))
               
    # ---- global layering direction (fallback) ----
    if matDirs is None:
      if not self.randomMatDir:
        matDirs = _defaultMatDir
        rowVec = False
      else:
        # ---- per-cell material direction ----
        cellMatDir = []
        for i in range(self.npts):
          d = random_direction(dim=self.dim)
          if self.dim == 2:
            d = np.append(d, 0.0)
          d = np.asarray(d, dtype=float)
          d /= max(np.linalg.norm(d), 1e-15)
          cellMatDir.append(d)
        self.cellMatDir = np.asarray(cellMatDir, dtype=float)
        rowVec = True       
    else:
      # First figure out shape of material directions
      matDirs = np.asarray(matDirs)
      shape = matDirs.shape
      match len(shape):
        case 1:
          self.cellMatDir=matDirs.reshape(1, 3, 1)
          rowVec = True
        case 2:
          self.cellMatDir=np.expand_dims(matDirs, axis=2)
          rowVec = True
        case 3:
          self.cellMatDir = matDirs
          rowVec = False
        case _:
          raise "Unknown material direction shape" 

      numMatDirs = self.cellMatDir.shape[0]
      
      # If only one material direction is specified treat as global, otherwise there should be the same number of material
      # directions as voronoi cell centers
      if numMatDirs == 1:
        # If you don't want random per-cell directions, use globalMatDir for all cells
        self.cellMatDir = np.tile(self.cellMatDir, (self.npts, 1))
      else:
        assert numMatDirs == self.npts, "When specifying more than one material direction, the number of directions must match the number of voronoi cells"

    # Make matrices of material directions
    if rowVec:
      self.cellMatDirMat = []
      for d in self.cellMatDir:
        print(d)
        u1, u2, u3 = generate_orthonormal_axes(d)
        self.cellMatDirMat.append(np.vstack((u1, u2, u3)))
      self.cellMatDirMat = np.asarray(self.cellMatDirMat, dtype=float)
    else:
      self.cellMatDirMat = self.cellMatDir

    # ---- number of slices per cell ----
    self.cellNumSlices = np.zeros(self.npts, dtype=int)
    for i in range(self.npts):
       self.cellNumSlices[i] = max(
        1,
        int(np.ceil(2.0 * float(self.cellRadius[i]) / float(self.weibullLayerThickness)))
      )

    # ---- compute slice strength scales for each cell (ragged) ----
    self.cellSliceStrengthScale = [None] * self.npts 
    rng = np.random.default_rng(self.seed)

    cellVolsFromSlices = []
    for i in range(self.npts): 
      m1 = np.asarray(self.cellMatDir[i][0], dtype=float)
      m2 = np.asarray(self.cellMatDir[i][1], dtype=float)
      m3 = np.asarray(self.cellMatDir[i][2], dtype=float)

      # # construct m2, m3 spanning plane normal to m1
      # if abs(np.dot(m1, np.array([0.0, 0.0, 1.0])) - 1.0) < 1e-12:
      #   m2 = np.cross(np.array([0.0, 1.0, 0.0]), m1)
      # else:
      #   m2 = np.cross(np.array([0.0, 0.0, 1.0]), m1)
      # m2 = m2 / max(np.linalg.norm(m2), 1e-15)
      # m3 = np.cross(m1, m2)

      nSlices_i = int(self.cellNumSlices[i]) 
      Rcell = float(self.cellRadius[i])

      # If cellRadius is 0 (edge/unbounded fallback), still define 1 slice value
      if Rcell <= 0.0:
        R = rng.uniform(1e-20, 1.0)
        s = ((self.weibullReferenceVolume / self.weibullMinVolume) *
             (np.log(R) / np.log(0.5))) ** (1.0 / self.weibullModulus)
        self.cellSliceStrengthScale[i] = np.array([s], dtype=float)
        continue

      slicez = np.linspace(-Rcell, Rcell, nSlices_i)
      points = np.linspace(-Rcell, Rcell, self.nIntegrationPoints)

      sliceStrengthScale = np.zeros(nSlices_i, dtype=float) 

      layerVolsList=[]
      for j in range(nSlices_i):  
        p0 = self.cellCentroid[i] + slicez[j] * m1 #get center origin of slice

        count = 0
        for r2 in points:
          for r3 in points:
            pt = p0 + r2 * m2 + r3 * m3 #test point on the plane

            # if point is interior to subobject and is nearest to the voronoi cell i
            if (subObject.isInterior(pt, 0.0) >= 0) and (self._cellId(pt) == i):
              count += 1
        # below, we have a grid of nintegration points, hence squaring of nIntegration points 
        # the area of the slice is approximated as (2*Rcell)^2
        layerVol = max(
          self.minWeibullVolume,
          self.weibullLayerThickness *
          (count / float(self.nIntegrationPoints**2)) *
          (2.0 * Rcell)**2
        )
        
        layerVolsList.append(layerVol)

        R = rng.uniform(1e-20, 1.0)
        s = ((self.weibullReferenceVolume / layerVol) *
             (np.log(R) / np.log(0.5))) ** (1.0 / self.weibullModulus)

        sliceStrengthScale[j] = s  

      cellVolsFromSlices.append(np.sum(layerVolsList))
      self.cellSliceStrengthScale[i] = sliceStrengthScale      

    ## after all cells are processed:
    self.cellVolumeFromSlices = np.asarray(cellVolsFromSlices, dtype=float)
    
    ## per-cell "reference radius" = radius of equal-volume sphere
    self.cellRefRadius = (3.0 * self.cellVolumeFromSlices / (4.0 * np.pi)) ** (1.0/3.0)
    
    # reference radius for only interior cells
    inside = np.array([subObject.isInterior(c, 0.0) >= 0 for c in self.cellCentroid], dtype=bool)
    
    # avoid empty slice
    if np.any(inside):
        self.avgRefRadiusInside = float(np.mean(self.cellRefRadius[inside]))
    else:
        self.avgRefRadiusInside = float("nan")
    
    print("avg reference radius (inside):", self.avgRefRadiusInside)
  
  def isInterior(self, pt, skinDepth):
      return self.object.isInterior(pt, skinDepth)

  # ---------------------------
  # point -> nearest Voronoi site id
  def _cellId(self, pt) -> int:
      x = np.asarray(pt[:self.dim], dtype=float)
      dist, idx = self.kdt.query(x.reshape(1, -1), k=5)  # show top-5 nearest
  
      # print only a few times
      if not hasattr(self, "_dbg_cellid_calls"):
          self._dbg_cellid_calls = 0
      if self._dbg_cellid_calls < 20:
          self._dbg_cellid_calls += 1
  
      return int(idx[0][0])   # closest seed = Voronoi cell id

  # ---------------------------
  # centroid from vertices (ConvexHull-based; only use when verts form a real polytope)
  def voronoi_cell_centroid(self, vertices):
    vertices = np.asarray(vertices)
    dim = vertices.shape[1]
    if dim not in (2, 3):
      raise ValueError("Only 2D or 3D vertices are supported")
    hull = ConvexHull(vertices)
    if dim == 2:
      return self.centroid_2d(vertices, hull)
    else:
      return self.centroid_3d(vertices, hull)

  def centroid_2d(self, vertices, hull):
    pts = vertices[hull.vertices]
    x = pts[:, 0]; y = pts[:, 1]
    x_next = np.roll(x, -1); y_next = np.roll(y, -1)
    cross = x * y_next - x_next * y
    area = 0.5 * np.sum(cross)
    if np.isclose(area, 0):
      return np.mean(pts, axis=0)
    cx = np.sum((x + x_next) * cross) / (6 * area)
    cy = np.sum((y + y_next) * cross) / (6 * area)
    return np.array([cx, cy])

  def centroid_3d(self, vertices, hull):
    centroid = np.zeros(3)
    volume = 0.0
    origin = np.zeros(3)
    for simplex in hull.simplices:
      a, b, c = vertices[simplex]
      v = np.dot(a, np.cross(b, c)) / 6.0
      tetra_centroid = (origin + a + b + c) / 4.0
      centroid += v * tetra_centroid
      volume += v
    if np.isclose(volume, 0):
      return np.mean(vertices, axis=0)
    return centroid / volume

  def getMatDir(self, pt):
    cellId = self._cellId(pt)
    return self.cellMatDirMat[cellId]

  def getStrengthScale(self, pt):
    i = self._cellId(pt)
    m1 = self.cellMatDir[i][0]
    z = self.cellRadius[i] + np.dot(m1, pt - self.cellCentroid[i])
    k = int(np.floor(z / self.weibullLayerThickness))
    k = max(0, min(int(self.cellNumSlices[i]) - 1, k))
    return self.cellSliceStrengthScale[i][k]


############################################
class voronoiWeibullBoxWrapper(BaseWrapper):
  """
  Box wrapper for another object that will be used to assign voronoi-cell weibull distribution of strength scale.
  Works for 2D or 3D cases.
  Box should be bigger than the subobject.
  """
  def __init__(self,
               name,
               subObject,
               x0,
               x1,
               flawSize,
               weibullVolume,
               weibullModulus,
               weibullSeed: int,
               vMin,
               vpts=None,
               dim: int=3,
               randomMatDir: bool=False):
    super().__init__(name,
                     subObject)
    self.object = subObject
    self.dim = dim
    self.x0 = np.asarray(x0)
    self.x1 = np.asarray(x1)
    self.dx = self.x1-self.x0
    self.x0 = self.x0[:self.dim]
    self.x1 = self.x1[:self.dim]
    self.randomMatDir = randomMatDir
    self.seed = weibullSeed

    if vpts is None:
      self.vpts = poisson(flawSize, x0=self.x0, dx=self.dx[:self.dim], seed=self.seed, dim=self.dim)
    else:
      self.vpts = vpts

    self.vpts = self.vpts[:,:self.dim] # Remove spacing from points
    self.npts = self.vpts.shape[0]
    self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')
    self.voronoi = Voronoi(self.vpts)

    # Estimate average volume per Voronoi cell
    v0 = np.prod(self.x1 - self.x0) / self.npts
    vor = self.voronoi
    vol = np.zeros(vor.npoints)

    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
      if (-1 in indices) or (vor.vertices[indices].shape[0] < 1):
        vol[i] = v0
      else:
        vol[i] = ConvexHull(vor.vertices[indices]).volume

        # Adjust volume by intersection with interior
        numInteriorVertices = 0
        numVertices = len(indices)
        for v in vor.vertices[indices]:
          if (dim==2):
            v = np.append(v, np.array([0.0])) # TODO There could be situation where the user unintentionaly sets the domain depth so that it does not straddle the origin
          if subObject.isInterior(v,0.0) >= 0:
            numInteriorVertices += 1

        vol[i] = vol[i]*numInteriorVertices/numVertices*(1.0 if self.dim==2 else self.dx[self.dim-1])

      vol[i] = max(vol[i], vMin)

    # Compute strength scale for each cell
    cellStrengthScale = []
    for i in range(self.npts):
      s = ((weibullVolume / vol[i]) * 
          (np.log(np.random.uniform(1e-20, 1.0)) / np.log(0.5))) ** (1.0 / weibullModulus)
      cellStrengthScale.append(s)
    self.cellStrengthScale = np.array(cellStrengthScale)

    # cells can have random orientation or can inheret from subObject
    cellMatDir=[]
    if ( self.randomMatDir ): 
      for i in range(0,self.npts):
        d = random_direction(dim=self.dim)
        if self.dim == 2:
          d = np.append(d, 0.0)
        u1, u2, u3 = generate_orthonormal_axes(d)
        d = np.vstack((u1, u2, u3))
        cellMatDir.append(d)
    self.cellMatDir = np.array(cellMatDir)

  def getMatDir(self,pt):
    if ( self.randomMatDir ):  
      x = np.asarray(pt[:self.dim])
      _, index = self.kdt.query(x.reshape(1,-1), k=1)
      matDir = self.cellMatDir[index[0][0]]
    else:
      if hasattr( self.object, 'getMatDir' ):
        matDir = self.object.getMatDir( pt )
      elif hasattr( self.object, 'matDir' ):
        matDir = self.object.matDir( self.object, matDir ) if callable( self.object.matDir ) else self.object.matDir
      else:
        matDir = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
    return matDir

  def getStrengthScale(self, pt):
    x = np.asarray(pt[:self.dim])
    _, index = self.kdt.query(x.reshape(1, -1), k=1)
    return self.cellStrengthScale[index[0][0]]

  def getDamage(self, pt):
    return self.object.getDamage(pt)


############################################
class pointwisePorosityWrapper(BaseWrapper):
  """
  Box wrapper for another object that will randomly define isInterior=False for a fraction of points equal
  to the specified porosity.  Not deterministic and doesn't create surface flags.  Used to introduce
  a simple explicit representation of under-resolved porosity.

  Works for 2D or 3D cases
  """
  def __init__(self,
               name,
               subObject,
               porosity):
    super().__init__(name,
                     subObject)
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 
    self.porosity = porosity

  def getMat(self, pt):
    if np.random.random() < self.porosity:
        return -1
    else:
        return self.object.getMat(pt)


############################################
class voronoiMatDirBoxWrapper(BaseWrapper):
  """
  Box wrapper for another object that will be used to assign voronoi-cell random distribution of material direction
  Works for 2D or 3D cases
  Box should be bigger than the subobject
  """
  def __init__(self,
               name,
               subObject,
               x0,
               x1,
               flawSize,
               seed=1,
               dim: int=3):
    super().__init__(name,
                     subObject)
    self.dim = dim
    self.x0 = np.asarray(x0)
    self.x1 = np.asarray(x1)
    self.dx = self.x1-self.x0
    self.x0 = self.x0[:self.dim]
    self.x1 = self.x1[:self.dim]
    self.seed = seed
    
    self.vpts = poisson(flawSize, x0=self.x0, dx=self.dx[:self.dim], seed=self.seed, dim=self.dim)
    self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
    self.npts = self.vpts.shape[0]
    self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(len(self.vpts) / 2)), metric='euclidean')

    # define the value of strength scale that will be assigned to each cell's particles.
    cellMatDir=[]
    for i in range(0,self.npts):
      d = random_direction()
      cellMatDir.append(d)

    self.cellMatDir = np.array(cellMatDir)

  def getMatDir(self,pt):
    x = np.asarray(pt[:self.dim])
    _, index = self.kdt.query(x.reshape(1,-1), k=1)
    matDir = self.cellMatDir[index[0][0]]
    return matDir

#############################################
class sizeEffectVoronoiWeibullBoxWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               flawfield_file,
               grid,
               x0,
               x1,
               pores,
               gamma,
               max_strength_scale,
               voronoi_file,
               weibullVolume,
               weibullModulus,
               weibullSeed,
               vMin,
               periodic = [False, False, False]):
    super().__init__(name,
                     subObject)
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 

    self.x0 = np.asarray(x0) # [x0, y0] should define box bigger than subObject
    self.x1 = np.asarray(x1) # [x1, y1] should define box bigger than subObject
    self.dx = self.x1-self.x0
    self.vel = subObject.vel
    self.mat = subObject.mat
    self.group = subObject.group
    self.strengthScale = 'map'  
    self.periodic = periodic
    self.gamma = gamma
    self.max_strength_scale = max_strength_scale
    self.weibullModulus = weibullModulus

    log2file("Creating sizeEffectVoronoiWeibullBox geometry object...")

    log2file("Creating max flaw radius field...")
    self.maxFlawRadiusField = maxFlawRadius(flawfield_file, 
                                            grid, 
                                            grid, 
                                            np.array([0.0,0.0,0.0]), 
                                            np.array([1.0,1.0,1.0]), 
                                            pores, 
                                            dim = 3,
                                            periodic = self.periodic, 
                                            readFromFile=True)
    log2file("Finished creating max flaw radius field")
    
    num_lines = countFileLines(voronoi_file)
    line_range = np.ceil(num_lines*np.array([g_rank, g_rank+1])/g_num_ranks)
    file_pos = fileOffsetFromLineNum(voronoi_file, line_range[0])

    print("Reading voronoi weibull center lines " + str(int(line_range[0])) + " to " + str(int(line_range[1])) + "on rank " + str(g_rank) + "...")
    log2file("Reading voronoi weibull center lines " + str(int(line_range[0])) + " to " + str(int(line_range[1])) + "...")

    pts = np.empty((num_lines, 4), dtype=float)
    vcFlawSize = np.empty(num_lines, dtype=float)
    with open(voronoi_file, 'r') as f:
        f.seek(file_pos)
        for i in range(int(line_range[0]), int(line_range[1])):
            line = f.readline()
            if i % 1000 == 0:
                print('Voronoi: (' + str(i) + '/' + str(num_lines) + ')...')
                log2file('Voronoi: (' + str(i) + '/' + str(num_lines) + ')...')
            terms = line.split()
            c = np.array([float(j) for j in terms])
            c = np.multiply(c, self.dx) + self.x0

            a = self.maxFlawRadiusField.getFlawSize(np.divide(np.array(c)-self.x0, self.dx))    
            vcFlawSize[i] = self.gamma/(math.sqrt(a)*self.dx[0])

            pts[i,:] = np.append(c, 0.0)

    print("Finished reading voronoi centers from file")
    log2file("Finished reading voronoi centers from file")

    for r in range(g_num_ranks):
        print("Syncing voronoi cell centers on rank " + str(r) + "...")
        log2file("Syncing voronoi cell centers on rank " + str(r) + "...")
        data_range = np.ceil(num_lines*np.array([r, r+1])/g_num_ranks)

        vcflaw_buffer = None
        pts_buffer = None
        if r == g_rank:
            vcflaw_buffer = np.copy(vcFlawSize)
            pts_buffer = np.copy(pts)

        vcflaw_buffer = g_comm_world.bcast(vcflaw_buffer, root=r)
        pts_buffer = g_comm_world.bcast(pts_buffer, root=r)

        vcFlawSize[int(data_range[0]):int(data_range[1])] = vcflaw_buffer[int(data_range[0]):int(data_range[1])]
        pts[int(data_range[0]):int(data_range[1]),:] = pts_buffer[int(data_range[0]):int(data_range[1]),:]
        
    self.npts = len(pts)
    v0=np.prod(self.dx)/self.npts # average volume, will be assigned to edge cells.

    log2file("Adding images of voronoi cell centers...")
    pts = add_images(pts, dim=3, x0=self.x0, dx=self.dx, periodic=self.periodic)
    nipts = len(pts) # Update number of points to include those of images too
    pts = pts[:,0:3] # Removes spacing from end of array  
    log2file("Finished adding images of voronoi cell centers")

    #voronoi tesselation of points
    log2file("Creating voronoi cells using scipy...")
    # vor=Voronoi(pts)
    pv = ParallelVoronoi(self.x0, self.x1, pts, self.periodic)
    log2file("Finished creating voronoi cells using scipy")

    # neighbor list for points:
    log2file("Creating KDTree for voronoi cell centers...")
    cell_centers_w_images = add_images(np.hstack((pv.cell_centers, np.zeros((len(pv.cell_centers),1)))), dim=3, x0=self.x0, dx=self.dx, periodic=self.periodic)
    cell_centers_w_images = cell_centers_w_images[:,0:3]
    self.kdt = KDTree(cell_centers_w_images, leaf_size=int(np.ceil(len(pts)/2)), metric='euclidean')
    log2file("Finished creating KDTree for voronoi cell centers")

    # self.vor = vor #CC: Temporary fix to grab from particle file writer for surface detection of voronoi cells

    self.rr = [np.random.uniform(1e-20,1.0) for i in range(self.npts)] # Store random weibull perturbations for strength query

    # compute volume of each voronoi cell
    # vol = np.zeros(vor.npoints)
    vol = np.zeros(self.npts)
    self.cellStrengthScale=np.empty(self.npts, dtype=float)
    
    pts_range = np.ceil(self.npts*np.array([g_rank, g_rank+1])/g_num_ranks)

    # for i in range(self.npts):
    for i in range(int(pts_range[0]), int(pts_range[1])):
      if i % 10000 == 0:
        if g_rank == 0:
            print("Computing volume of voronoi cell " + str(i+1) + "/" + str(self.npts))
        log2file("Computing volume of voronoi cell " + str(i+1) + "/" + str(self.npts))

      # reg_num = vor.point_region[i]
      indices = pv.cells[i] #vor.regions[reg_num]
      vertices = pv.vertices[indices] #vor.vertices[indices]
      numVertices = vertices.shape[0]

      if ( (-1 in indices) or ( numVertices < 1 ) ): # some regions can be opened
        vol[i] = v0
      else:
        vol[i] = ConvexHull(vertices).volume
        numInteriorVertices = 0
                
        for v in vertices:
          v3d = mapToRange(v, self.x0, self.x1)

          if subObject.isInterior(v3d, fromFile=True):
            numInteriorVertices += 1
        
        vol[i] = vol[i]*numInteriorVertices/numVertices

      vol[i] = max( vol[i], vMin )

      # define the value of strength scale that will be assigned to each cell's particles.
      # s = ( weibullVolume/vol[i] )**(1.0/weibullModulus)
      s = min(self.max_strength_scale, max(( weibullVolume/vol[i] )**(1.0/weibullModulus), vcFlawSize[i]))*(np.log( self.rr[i] )/np.log(0.5) )**(1.0/self.weibullModulus)
      self.cellStrengthScale[i] = s

    # Sync cellStrengths
    for r in range(g_num_ranks):
        print("Syncing cell strength scales on rank " + str(r) + "...")
        log2file("Syncing cell strength scales on rank " + str(r) + "...")
        buffer_range = np.ceil(self.npts*np.array([r, r+1])/g_num_ranks)

        strength_buffer = None
        if r == g_rank:
            strength_buffer = np.copy(self.cellStrengthScale)

        strength_buffer = g_comm_world.bcast(strength_buffer, root=r)

        self.cellStrengthScale[int(buffer_range[0]):int(buffer_range[1])] = strength_buffer[int(buffer_range[0]):int(buffer_range[1])]

    log2file("Finished constructing sizeEffectVoronoiWeibullBox")
    print("Finished constructing sizeEffectVoronoiWeibullBox")

    log2file("Finished creating sizeEffectVoronoiWeibullBox geometry object")


  def getMat(self, pt):
    return self.object.getMat(pt)

  def isInterior(self, pt, skinDepth):
    return self.object.isInterior(pt, skinDepth)

  def getStrengthScale(self, pt) -> float:
    # Define the value of strength scale that will be assigned to each cell's particles.

    # # Precompute inverse square of flaw field to speed up job setup?
    # a = self.maxFlawRadiusField.getFlawSize(np.divide(np.array(pt)-self.x0, self.dx))    
    # s_size = self.gamma/(math.sqrt(a)*self.dx[0]) # maxFlawRadius uses square unit domain, need to new domain size

    # Voronoi Weibull Strength Scale
    x = np.asarray(pt)
    _, index = self.kdt.query(x.reshape(1,-1), k=1)
    # #KDT Tree includes image centers of cells so need to find cell volume 
    index[0][0] = index[0][0] % self.npts # Ensure image indices are wrapped into interior cell indices
    s_weibull = self.cellStrengthScale[index[0][0]]
    # return min(self.max_strength_scale, max(s_weibull, s_size))*(np.log( self.rr[index[0][0]] )/np.log(0.5) )**(1.0/self.weibullModulus)
    
    return s_weibull # CC debug trying strengh scale with size effect determined only by center of voronoi cell
  
  def xMin(self):
    # This might be used for binning objects
    return self.object.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.object.xMax()


#############################################
class surfaceFlagBoxWrapper(BaseWrapper):
  def __init__(self,
               name,
               x0,
               x1,
               surfaceFlag: int,
               subObject):
    super().__init__(name,
                     subObject)
    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.surfaceFlag = surfaceFlag
  
  def isInterior(self, pt, skinDepth):
    x = np.asarray(pt)
    flag = self.subObject.isInterior(pt, skinDepth)
    if flag > 0 and np.all( np.logical_and(x > self.x0, x < self.x1)):
      return self.surfaceFlag
    
    return flag

#############################################
class damageBoxWrapper(BaseWrapper):
  def __init__(self,
               name,
               subObject,
               x0,
               x1,
               damage):
    super().__init__(name,
                     subObject)
    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.boxDamage = damage

  def getDamage(self, pt) -> float:
    x = np.asarray(pt)
    if np.all( np.logical_and(x > self.x0, x < self.x1)):
      return self.boxDamage
    
    if hasattr(self.subObject, 'damage'):
      return self.subObject.damage

    if hasattr(self.subObject, 'getDamage'):
      return self.subObject.getDamage(pt)
    
    return 0.0
  

class pennyShapedCrackWrapper(BaseWrapper):
  def __init__(self,
               name: str,
               subObject,
               x1,
               x2,
               r,
               damage: float):
    super().__init__(name,
                     subObject)
    self.x1 = np.asarray(x1)
    self.x2 = np.asarray(x2)
    self.r = r
    self.crackDamage = damage
    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.axis = (self.x2-self.x1)/self.h    

  def getDamage(self, pt) -> float:
    x = np.asarray(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point 
    
    if (z>0 and z<self.h and r<self.r ):
      return self.crackDamage
    
    if hasattr(self.subObject, 'damage'):
      return self.subObject.damage
    
    if hasattr(self.subObject, 'getDamage'):
      return self.subObject.getDamage(pt)
  
    return 0.0
  
class layeredStrengthScaleBoxWrapper(BaseWrapper):
  # defines strength scale for a box based on layers relative to a direction
  # subObject isInterior is used to estimate true layer volume
  def __init__(self,
               name: str,
               subObject,
               x1,
               x2,
               matDir = None,
               weibullReferenceVolume = 1,
               weibullModulus = 1,
               weibullSeed = 1,
               weibullLayerThickness = 1,
               weibullMinVolume = 0,
               nIntegrationPoints = 10,
               layerNormal = None,
               dim = 3,
               domainThickness = None
               ):
    super().__init__(name,
                     subObject)
    
    x1 = np.asarray(x1) # bounding box corner
    x2 = np.asarray(x2) # bounding box corner
    self.x0 = 0.5*(x1+x2) # box center.
    
    # if specified define self.matDir, otherwise it will be inherited from subobject
    if (matDir is not None):
      self.matDir = matDir
   
    # corner to corner vector for bounding box, used to define slices.
    h = np.linalg.norm(x2[:dim]-x1[:dim]) 
    
    if ( dim == 2 and domainThickness is None):
      print('ERROR you need to specify domain Thickness if dim=2')
    else:
      self.domainThickness = domainThickness

    # layerNormal is single direction used to define slices.  but material may have
    # multiple directions, matDirs that will be stored in the matDir variable 
    # read by pfw.  By default we let layer normal be the first direction.    
    if layerNormal is None:    
      if (np.size(matDir) == 9):
        self.layerNormal = np.asarray(matDir[0])
      elif (np.size(matDir) == 3):
        self.layerNormal = np.asarray(matDir)
    else:
      self.layerNormal = layerNormal
    norm = np.linalg.norm(self.layerNormal)
    self.layerNormal = self.layerNormal / norm
       

    # number of slices for strength scale layers based on prescribed thickness
    # adjusted to have an integer number of slices.
    nslices = int( h  / weibullLayerThickness )

    # relative position on layer normal axis relative to x0
    self.slicez = np.linspace( -0.5*h, 0.5*h, nslices)

    # for each slice defined by a point along the box axis,
    # create an array of test points on a plane normal to matDir
    # estimate volume based on fraction of those points interior to subobject.
    m1 = self.layerNormal
    if (dim==3):
      if abs(np.dot( m1, np.array([0.0,0.0,1.0]))-1) < 1e-12:
        m2 = np.cross(np.array([0.0,1.0,0.0]),m1)
      else:
        m2 = np.cross(np.array([0.0,0.0,1.0]),m1)
        m2 = m2 / np.linalg.norm(m2)
      m3 = np.cross(m1,m2)
    elif (dim == 2 ):
      m1[2] = 0.0
      m3 = np.array([0.0,0.0,1.0])
      m2 = np.cross(m3,m1)
      m2[2]=0.0
    else:
      print( " Bad dim input to weibull layer wrapper.")        
       

    
    # define the value of strength1 scale that will be assigned to each cell's particles.
    sliceStrengthScale=[]
    
    # test points in transverse plane.
    points = np.linspace( -0.5*h, 0.5*h, nIntegrationPoints)
    
    np.random.seed(weibullSeed)
    for i in range(0,nslices):
      p0 = self.x0 + self.layerNormal*self.slicez[i]  # Position on layer normal axis of layer-center
      count = 0     
      if (dim == 3):
        for r2 in points:
          for r3 in points:
            if ( subObject.isInterior( p0 + r2*m2 + r3*m3 , 0.0 ) >= 0 ):           
              count += 1
        vol = max( weibullMinVolume, weibullLayerThickness * (count / nIntegrationPoints**2 ) * h**2 )
      if (dim == 2):
        for r2 in points:
          if ( subObject.isInterior( p0 + r2*m2, 0.0 ) >= 0 ):           
            count += 1
        vol = max( weibullMinVolume, weibullLayerThickness * (count / nIntegrationPoints ) * domainThickness * h )
                
      R = np.random.uniform()
      s = ( ( weibullReferenceVolume / vol )*( np.log( R ) / np.log(0.5) ) )**(1.0/weibullModulus)
      sliceStrengthScale.append(s)
      print("slice ",i," has vol = ",vol,"minVol = ",weibullMinVolume," and strengthScale = ",s)
    self.sliceStrengthScale = np.array(sliceStrengthScale)
    
  def getStrengthScale(self, pt) -> float:
    # z-coordinate of test point relative to box-corner axis
    z = np.dot( np.asarray(pt) - self.x0, self.layerNormal )
    
    # index of the slice containing the point.
    idx = np.argmin(np.abs(self.slicez - z ))
    
    return self.sliceStrengthScale[idx]

##################################################################

# class voronoiWrapperWithLayeredStrengthScale(BaseWrapper):
#   """
#   Voronoi tessellation + layered strengthScale inside each cell.

#   **This may be deprecated.  The same can be achieved combining voronoi tesselation
#   wrappers and layered strength scale wrappers, by using the same voronoi points for the 
#   two cases.**

#   strengthScale(pt) = cellSliceStrengthScale[cellId(pt), layerId(pt)]
#   For each particle:
#   Determine which Voronoi cell particle is in               | _cellId.
#   Determine the cell’s material direction                   | _matDirVec_at from cellMatDir, wrapped object, or globalMatDir)
#   Turn material direction into a 3×3 basis                  | getMatDir.
#   Determine where the particle is along that direction      |  getStrengthScale (compute k)
#   Determine which layer that particle is in                 | getStrengthScale (compute k)
#   Determine Weibull strength for the cell and layer         | cellSliceStrengthScale[cellId, k]
  
#   """

#   def __init__(self,
#                name,
#                subObject,                   # object subjected to strenght scaling
#                x0,                          # x0 and x1 defin dimentions for voronoi polygons
#                x1,                          # x0 and x1 affect voronoi cell volumes and where intra-cell layer boundaries fall
#                flawSize,                    # controls how far away seeds are form one another
#                weibullVolume,               # webiull volume, moudulus and seed affect s, the strengthscale
#                weibullModulus,
#                weibullSeed: int,
#                vMin,                        # vmin is the lower bound on volume used in strength scaling
#                weibullLayerThickness,       # thickness of layers within each cell
#                nIntegrationPoints = 10,     # points in each direction to estimate layer volume
#                matDir=np.array([1,0,0]),    # "None": no layers, "random": each cell has random orientation, [nx,yn,zn]: prescribed direction
#                vpts=None,                   # optional voronoi points
#                dim: int = 3,                # used for cell material directions
#                randomMatDir: bool = False): # allows random distribution of strengthscale layers in the cells


#     super().__init__(name, subObject)       # initialize base wrapper.  
#     self.object = subObject                 # wrapper “acts like” subObject, but override things like getStrengthScale() and getMatDir().
#     self.dim = dim

#     # Define the bounding box dimentions of the Voronoi Cells.
#     #The Voronoi cells will be created inside this box, 
#     # their sizes will vary depending on the site distribution (vpts) and flawSize.
#     self.x0_full = np.array(x0, dtype=float)     #store input box corners, x0 and x1 in arrays
#     self.x1_full = np.array(x1, dtype=float)
#     self.x0 = self.x0_full[:self.dim]            #adjust the arrays depending on the dimensions of the problem
#     self.x1 = self.x1_full[:self.dim]
#     self.dx = (self.x1_full - self.x0_full)   #box size vector for the cells | dx=x1−x0=[Lx​, Ly​, Lz​]


#     #assign other traits to object
#     self.randomMatDir = bool(randomMatDir)    # assignes the flag for random layering directions in each cell
#     self.seed = int(weibullSeed)              #converts the input seed to an integer | initialize the random number generator(s), 
#                                               #so you get reproducible Voronoi patterns and Weibull draws.
#     self.weibullLayerThickness = float(weibullLayerThickness) #assigns thickness of layers within each cell
#                                               # later will be in k=⌊layerThickness(pt−layerOrigin)⋅m1​​⌋


#     # Global layer direction (used if randomMatDir=False and subObject doesn't provide getMatDir)
#     # the default direction is z
#     if matDir is None and randomMatDir == False:
#       matDir = np.array([0.0, 0.0, 1.0]) #make horizontal layers
#     self.globalMatDir = np.array(matDir, dtype=float) # convert what ever matDir is, to an array to use later
#     self.globalMatDir /= np.linalg.norm(self.globalMatDir)
#                               #normalize the direction vector to unit length 
#                               # ensure that z measures distance in physical units along the direction
#                               # will be needed for dot product with m1 later. 




#     # Build Voronoi points~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # the seeds points are not provided, created them. 
#     # if they are provided, ensure they are in an array format
#     if vpts is None:
#       #us a posson distribution to define the seed points.
#       vpts_raw = poisson(flawSize, x0=self.x0, dx=(self.x1 - self.x0), seed=self.seed, dim=self.dim)
#     else:
#       vpts_raw = np.array(vpts, dtype=float)
    
#     # helps with formatting if vpts varies in dimensions from one problem input to another
#     self.vpts = np.asarray(vpts_raw, dtype=float)[:, :self.dim]
#     self.npts = self.vpts.shape[0]               #number of voronoi points
    
#     # Create our look up structure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
#     # Even though the Voronoi diagram defines the cell boundaries, 
#     # it’s expensive to test “which polyhedron contains this point” directly.
#     # So instead we use the key property of Voronoi tessellations:
#     # A point belongs to the cell of the seed it is closest to.
#     # KDTree makes that closest-seed query fast for millions of particles.
#     # good reference for KDTree tutorial: https://www.youtube.com/watch?v=TLxWtXEbtFE
#     self.kdt = KDTree(self.vpts, leaf_size=int(np.ceil(self.npts / 2)), metric='euclidean')
    
    
#     # Build the Voronoi tessellation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # refer here for more information. https://docs.scipy.org/doc//scipy-1.6.3/reference/generated/scipy.spatial.Voronoi.html
#     self.voronoi = Voronoi(self.vpts)

#     # Estimate average volume | take product of box lx, ly and possibly lz / number of points
#     v0 = np.prod(self.x1 - self.x0) / self.npts
#     vor = self.voronoi

#     # vol = np.zeros(vor.npoints)
#     centroid = np.zeros(shape=(vor.npoints, self.dim))
#     cellRadius = np.zeros(vor.npoints)
        
#     for i, reg_num in enumerate(vor.point_region):    # get the region id for each voronoi cell
#       indices = vor.regions[reg_num]      
#       if (-1 in indices) or (len(indices) < 1):        # for unbounded cells, set their volume to the average volume
#         vol[i] = v0
#       else:
#         verts = vor.vertices[indices]                     # get the voronoi vertices to pass to ConvexHull(verts) to estimate area/volume.
#       # compute the centroid for voronoi cell with corner indices from "indices"
#       centroid[i] = self.voronoi_cell_centroid(verts)    # TODO: check this is real.
#       # Find radius of bounding sphere used to define slices and integration points:
#       cellRadius[i] = 0.0
#       for v in verts:
#         dist = np.linalg.norm( v - centroid[i,:self.dim] )
#         if dist > cellRadius[i]:
#           cellRadius[i] = dist
      
#     self.cellCentroid = centroid # (npts, dim)
#     self.cellRadius = cellRadius
    
#     # Cell material directions (only used if randomMatDir=True)
#     # Cell material directions (only used if randomMatDir=True)
#     # --------------------------------------------------------
#     # If randomMatDir is enabled, give each Voronoi cell its own random
#     # unit direction vector. This direction is later used to:
#     #   - define the layering direction inside that cell (getStrengthScale)
#     #   - define the material direction basis returned to the particleFileWriter (getMatDir)
#     #
#     # If randomMatDir is False, we do NOT create per-cell directions here; instead
#     # direction will come from:
#     #   - either the wrapped subObject.getMatDir(pt) or
#     #   - the wrapper's globalMatDir fallback.
#     # ** will need to update if we want to assign differnt directions to cells 
#     # (like set 50% of directions to x normal and 50% to y normal)
#     self.cellMatDir = None
#     if self.randomMatDir:
#       cellMatDir = []
#       for i in range(self.npts):            # for each voronoi cell  
#         d = random_direction(dim=self.dim)  # generate a random direction vector
#         if self.dim == 2:                   # if 2d, append 0 to make it a 3d vector
#           d = np.append(d, 0.0)
#         d = np.array(d, dtype=float)        # ensure it is an array
#         d /= np.linalg.norm(d)              # normalize to unit length
#         cellMatDir.append(d)                # add to list of cell directions
#       self.cellMatDir = np.array(cellMatDir, dtype=float)    # Convert Python list -> NumPy array of shape (npts, 3)

#     # Decide number of layers in each voronoi cell ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # remember that we want layers of weibullLayerThickness.
#     cellNumSlices = []
#     for i in range(self.npts):            # for each voronoi cell  
#       cellNumSlices[i] = np.ceil( 2.*self.cellRadius[i] / self.weibullLayerThickness ) 
#     self.cellNumSlices = cellNumSlices

#     # for each point, create an array of test points on a plane normal to matDir
#     # estimate volume based on fraction of those points interior to subobject.
#     for i in range(self.npts):            # for each voronoi cell  
#       m1 = self.matDir[i]
#       if abs(np.dot( m1, np.array([0.0,0.0,1.0]))-1) < 1e-12:
#         m2 = np.cross(np.array([0.0,1.0,0.0]),m1)
#       else:
#         m2 = np.cross(np.array([0.0,0.0,1.0]),m1)
#         m2 = m2 / np.linalg.norm(m2)
#       m3 = np.cross(m1,m2)
      
#       # integration points
#       slicez = np.linspace( -self.cellRadius[i] , self.cellRadius[i], self.cellNumSlices[i] )
#       points = np.linspace( -self.cellRadius[i] , self.cellRadius[i], self.nIntegrationPoints )

#       # define the value of strength scale that will be assigned to each cell's particles.
#       sliceStrengthScale=[]
#       np.random.seed(weibullSeed)
#       for j in range(0,self.nslices):
#         count = 0
#         p0 = self.cellCentroid[i] + slicez[j]*self.cellMatdir[i]  
#         for r2 in points:
#           for r3 in points: 
#             # if point is interior to subobject and is nearest to the voronoi cell i
#             pt = p0 + r2*m2 + r3*m3 # test point
#             if ( ( subObject.isInterior( pt , 0.0 ) >= 0 ) and (self._cellID( pt ) == i ) ):           
#               count += 1
      
#         vol = max( weibullMinVolume, weibullLayerThickness * (count / nIntegrationPoints**2 ) * (2.*self.cellRadius[i])**2 )
#         R = np.random.uniform()
#         s = ( ( weibullReferenceVolume / vol )*( np.log( R ) / np.log(0.5) ) )**(1.0/weibullModulus)
#         sliceStrengthScale[j].append(s)
#         print("cell", i,", slice ",j," has vol = ",vol,"minVol = ",weibullMinVolume," and strengthScale = ",s)
        
#       # slice strength scale is array of arrays
#       self.cellSliceStrengthScale[i] = sliceStrengthScale
        

  # # Establish which cell a point belongs to ~~~~~~~~~~~~~~~~~~~~~~~~
  # def _cellId(self, pt) -> int:
  #   x = np.array(pt[:self.dim], dtype=float)  # ensure this will work in 2 or 3 dimentions
  #   _, idx = self.kdt.query(x.reshape(1, -1), k=1) # query the KDTree to find closest voronoi seed point
  #   # query() returns (distance(s), index/indices). We ignore the distances with "_".
  #   # x.reshape(1, -1) makes it a 2D array of shape (1, self.dim), as KDTree expects.
  #   return int(idx[0][0])

  # def voronoi_cell_centroid(vertices):
  #     """
  #     AI SLOP, check this.
      
  #     Compute the centroid of a convex Voronoi cell in 2D or 3D
  #     given only its vertices.

  #     Parameters
  #     ----------
  #     vertices : (N,2) or (N,3) array-like
  #         Voronoi cell vertices (unordered)

  #     Returns
  #     -------
  #     centroid : ndarray
  #         Centroid of the cell (area centroid in 2D, volume centroid in 3D)
  #     """
  #     vertices = np.asarray(vertices)
  #     dim = vertices.shape[1]

  #     if dim not in (2, 3):
  #         raise ValueError("Only 2D or 3D vertices are supported")

  #     hull = ConvexHull(vertices)

  #     if dim == 2:
  #         return _centroid_2d(vertices, hull)
  #     else:
  #         return _centroid_3d(vertices, hull)


  # def _centroid_2d(vertices, hull):
  #     """
  #     AI SLOP, check this.
      
  #     Area centroid of a convex polygon from its hull
  #     """
  #     pts = vertices[hull.vertices]
  #     x = pts[:, 0]
  #     y = pts[:, 1]

  #     # Shoelace formula
  #     x_next = np.roll(x, -1)
  #     y_next = np.roll(y, -1)

  #     cross = x * y_next - x_next * y
  #     area = 0.5 * np.sum(cross)

  #     if np.isclose(area, 0):
  #         return np.mean(pts, axis=0)

  #     cx = np.sum((x + x_next) * cross) / (6 * area)
  #     cy = np.sum((y + y_next) * cross) / (6 * area)

  #     return np.array([cx, cy])


  # def _centroid_3d(vertices, hull):
  #     """
  #     AI SLOP, check this.
      
  #     Volume centroid of a convex polyhedron via tetrahedralization
  #     """
  #     centroid = np.zeros(3)
  #     volume = 0.0
  #     origin = np.zeros(3)

  #     for simplex in hull.simplices:
  #         a, b, c = vertices[simplex]

  #         # Tetrahedron (origin, a, b, c)
  #         v = np.dot(a, np.cross(b, c)) / 6.0
  #         tetra_centroid = (origin + a + b + c) / 4.0

  #         centroid += v * tetra_centroid
  #         volume += v

  #     if np.isclose(volume, 0):
  #         return np.mean(vertices, axis=0)

  #     return centroid / volume


  # # # Helper function that gets the direction for the individual cell
  # # def _matDirVec_at(self, pt, cellId: int):
  
  # #   # Priority rule (what to use if multiple options exist):
  # #   #   1) If doing random orientations per Voronoi cell, use the
  # #   #      precomputed random direction for this specific cell.
  # #   if self.randomMatDir and self.cellMatDir is not None:
  # #     # return the random direction gibven to the cell
  # #     return self.cellMatDir[cellId]

  # #   #   2) Otherwise, inherit a direction from the wrapped subObject
  # #   #      (e.g., a materialDirectionWrapper or something with getMatDir()).
  # #   if hasattr(self.object, "getMatDir"):
  # #     md = np.asarray(self.object.getMatDir(pt))
      
      
  # #     #correct for the different shape of MatDir that are in the input files
  # #     # Case 1: subObject returns a full 3x3 orthonormal basis matrix.: ex. (3,3)
  # #     if md.ndim == 2 and md.shape == (3,3):
  # #       v = md[0]
  # #     # case 2: subObject returns a single direction vector (length 3): ex ((3,))
  # #     elif md.ndim == 1 and md.size == 3:
  # #       v = md
  # #     #Case 3: subObject returned some other shape/type: ex. (1,3), (3,1)
  # #     else:
  # #       v = self.globalMatDir
        
  # #     v = np.asarray(v, dtype=float)
  # #     # Normalize the direction to unit length 
  # #     # then dot-products measure “distance along direction” consistent with length units.
  # #     return v / np.linalg.norm(v)
    
    
  # #   # made to handle matDir input
  # #   if hasattr(self.object, "matDir"):
  # #       # `matDir` might be:
  # #       #   (a) a function of (object, pt) that returns a direction at this point, OR
  # #       #   (b) a stored constant direction (vector or basis).
  # #       # This line handles both cases:
  # #     md = self.object.matDir(self.object, pt) if callable(self.object.matDir) else self.object.matDir
  # #     md = np.asarray(md)
      
  # #     # If in 3x3 basis matrix, first row is primary direction vector (a1). use for layering / orientation.
  # #     if md.ndim == 2 and md.shape == (3,3):
  # #       v = md[0]
  # #     # If is a single 3-vector, use it directly as the direction.
  # #     elif md.ndim == 1 and md.size == 3:
  # #       v = md
  # #     # otherwise use single global direction.
  # #     else:
  # #       v = self.globalMatDir
      
  # #     v = np.asarray(v, dtype=float)
  # #     return v / np.linalg.norm(v)
 
  # #   #   3) Priority Rule 3: use the wrapper’s globalMatDir.
  # #   return self.globalMatDir
  

    
  # # finally turn “one direction / cell” into a 3×3 matrix of row directions that pfw can writes out~~~~~~~~~~~~~~~
  # #earlier,   In __init__, we estiblished the direction per Voronoi cell
  # # In _matDirVec_at(pt, cellId) we established which vector to use r a direction at each point. 
  # # either subObject.getMatDir | subObject.matDir | globalMatDir
  # # now we assign the directions to the rows
    
  # def getMatDir(self, pt):
  #   cellId = self._cellId(pt)                           #establish  to which cell we belong
  #   return self.cellMatDir[cellId]

  #   # TODO: add support for prescribed list of matDir corresponding to prescribed voronoi points list.
  #   # TODO: add support for 3 material directions.
  #   # a1 = self._matDirVec_at(pt, cellId)                 # get direction vector. Either random, from bObject.getMatDIr , subObject.matDir, or global dir
    
  #   # # now need to build: [[a1x a2x a3x], [a1y a2y a3y], [a1z a2z a3z]]
  #   # # where a1 is the primary direction vector, and a2, a3 are orthogonal to it.
  #   # # If a1 is nearly aligned with z-axis, use y-axis instead for better numerical stability.
  #   #   tmp = np.array([0.0, 0.0, 1.0])                    # make a helper vector alined with z axis
  #   # if abs(np.dot(a1, tmp)) > 0.99:                     # check if a1 is parallel to z axis
  #   #   tmp = np.array([0.0, 1.0, 0.0])                   # if so, use y axis instead for helper vector
  #   # a2 = np.cross(tmp, a1); a2 /= np.linalg.norm(a2)    # make a2 the cross product of helper and a1, normalize to unit length
  #   # a3 = np.cross(a1, a2)                               # make a3 perpendiculr to a1 and a2 (it's already unit length)
  #   # return np.vstack((a1, a2, a3))                      # make 3x3 matrix with rows a1, a2, a3


  # # for a Voronoi cell, for a layer in that cell, determine the Weibull strength scale ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # def getStrengthScale(self, pt):
    
  #   i = self._cellId(pt)                         # get cellID from kdtree search
  #   m1 = self.cellMatDir[i]                      # get the material direction, either random, subObject.getMatDIr , subObject.matDir, or global dir
  #   # poistion along slice-normal axis
  #   z = float( self.cellRadius[i] + np.dot( m1, pt - self.cellCentroid[i] ) )

  #   # slice index (elsewhere called j in initialization)
  #   k = int(np.floor(z / self.weibullLayerThickness)) # tells us which lyaer the point is in
  #   k = max( 0, min( self.cellNumSlices[i] - 1, k ) )     # for edge cases where z may be slightly outside the furthest layer
  #   return float(self.cellSliceStrengthScale[i, k])  # return strength scale for cell and layer


  # def getDamage(self, pt):                        #damage getter
  #   return self.object.getDamage(pt)

# ===========================================
# END GEOMETRY WRAPPERS
# ===========================================

# ===========================================
# TRANSFORMS
#
# ===========================================

#############################################
# Homogeneous Transform Matrices
def translate(dx):
  return np.array([[1.0, 0.0, 0.0, dx[0]],
                   [0.0, 1.0, 0.0, dx[1]],
                   [0.0, 0.0, 1.0, dx[2]],
                   [0.0, 0.0, 0.0, 1.0]])


def scale(ds):
  if len(ds) == 1:
    ds = np.array([ds,ds,ds])
  
  return np.array([[ds[0], 0.0, 0.0, 0.0],
                   [0.0, ds[1], 0.0, 0.0],
                   [0.0, 0.0, ds[2], 0.0],
                   [0.0, 0.0,   0.0, 1.0]])


# a0 = normal of reflection plane
# x0 = center of reflection
def reflect(a0, x0=np.array([0.0,0.0,0.0])):
  x0 = np.asarray(x0)
  a0 = np.asarray(a0)
  a0 = a0 / np.linalg.norm(a0)
  reflection = np.array([[1-2*a0[0]**2,     -2*a0[0]*a0[1],  -2*a0[0]*a0[2], 0.0],
                          [-2*a0[0]*a0[1], 1-2*a0[1]**2,     -2*a0[1]*a0[2], 0.0],
                          [-2*a0[0]*a0[2],  -2*a0[1]*a0[2], 1-2*a0[2]**2,    0.0],
                          [0.0,                        0.0,          0.0,    1.0]])
  return np.matmul(translate(x0), np.matmul(reflection, translate(-x0)))


# a0 = axis of rotation
# alpha = angle of rotation (radians)
# x0 = center of rotation
def rotate(alpha,a0=np.array([0.0,0.0,1.0]),x0=np.array([0.0,0.0,0.0])):
  x0 = np.asarray(x0)
  a0 = np.asarray(a0)
  aa = np.outer(a0,a0)
  A = np.array([ [ 0,     -a0[2],      a0[1]], 
                 [ a0[2],      0,     -a0[0]], 
                 [-a0[1],  a0[0],          0] ])
  R = np.identity(4)
  R[:3,0:3] = np.cos(alpha)*np.identity(3) + (1-np.cos(alpha))*aa + np.sin(alpha)*A
  return np.matmul(translate(x0), np.matmul(R, translate(-x0)))


#############################################
class transform(BaseWrapper):
  @abstractmethod
  def __init__(self,
               name,
               subObject,
               transform):
    super().__init__(name,
                     subObject)
    self.transform = np.asarray(transform)
    self.inverse = np.linalg.inv(self.transform[:3,:3])

  def transformPoint(self, pt):
    pt = np.asarray(pt)
    pt = np.append(pt, 1.0)
    pt = np.matmul(self.transform, pt)
    return pt[:3]
  
  def transformVector(self, vec):
    return np.matmul(self.inverse, np.asarray(vec))

  def isInterior(self, pt, skinDepth):
    return super().isInterior(self.transformPoint(pt), skinDepth)

  def getSurfaceNormal(self, pt):
    surfaceNormal = super().getSurfaceNormal(self.transformPoint(pt))
    return self.transformVector(surfaceNormal)
    
  def getSurfacePosition(self,pt):
    surfacePosition = super().getSurfacePosition(self.transformPoint(pt))
    return self.transformVector(surfacePosition)

  def getGroup(self, pt):
    return super().getGroup(self.transformPoint(pt))

  def getMatDir(self, pt):
    matDir = super().getMatDir(self.transformPoint(pt))
    return self.transformVector(matDir)

  def getDamage(self, pt) -> float:
    return super().getDamage(self.transformPoint(pt))

  def getPorosity(self, pt) -> float:
    return super().getPorosity(self.transformPoint(pt))

  def getTemperature(self, pt) -> float:
    return super().getPorosity(self.transformPoint(pt))

  def getSurfaceTraction(self, pt):
    surfaceTraction = super().getSurfaceTraction(self.transformPoint(pt))
    return self.transformVector(surfaceTraction)

  def xMin(self):
    return -np.inf

  def xMax(self):
    return np.inf


# ===========================================
# END TRANSFORMS
# ===========================================

# ===========================================
# SET OPERATIONS
# 
# ===========================================
class SetOperation(Geometry):
  @abstractmethod
  def __init__(self,
               name,
               subObjA,
               subObjB,
               defaultToA: bool=True):
    self.subObjA = subObjA
    self.subObjB = subObjB
    self.defaultToA = defaultToA

    if self.defaultToA:
      super().__init__(name,
                       vel=self.subObjA.vel,
                       mat=self.subObjA.mat,
                       group=self.subObjA.group,
                       particleType=self.subObjA.particleType)
    else:
      super().__init__(name,
                       vel=self.subObjB.vel,
                       mat=self.subObjB.mat,
                       group=self.subObjB.group,
                       particleType=self.subObjA.particleType)

  @abstractmethod
  def isInterior(self,pt,skinDepth):
    pass

  @abstractmethod
  def getSurfaceNormal(self, pt):
    if self.defaultToA:
      return self.subObjA.getSurfaceNormal(pt)
    else:
      return self.subObjB.getSurfaceNormal(pt)

  @abstractmethod
  def getSurfacePosition(self,pt):
    if self.defaultToA:
      return self.subObjA.getSurfacePosition(pt)
    else:
      return self.subObjB.getSurfacePosition(pt)

  @abstractmethod
  def getGroup(self, pt):
    obj = self.subObjA
    if not self.defaultToA:
      obj = self.subObjB

    if hasattr(obj,'getGroup'):
        group = obj.getGroup( pt )
    else:
        if hasattr(obj, 'group'):
          group = obj.group
        else:
          group = 0
    
    return group

  @abstractmethod
  def getDamage(self, pt):
    obj = self.subObjA
    if not self.defaultToA:
      obj = self.subObjB

    if hasattr(obj,'getDamage'):
        damage = obj.getDamage( pt )
    else:
        if hasattr(obj, 'damage'):
          damage = obj.damage
        else:
          damage = 0
    
    return damage

  @abstractmethod
  def getStrengthScale(self, pt):
    obj = self.subObjA
    if not self.defaultToA:
      obj = self.subObjB

    if hasattr(obj,'getStrengthScale'):
        strengthScale = obj.getStrengthScale( pt )
    else:
        if hasattr(obj, 'strengthScale'):
          strengthScale = obj.strengthScale
        else:
          strengthScale = _defaultStrengthScale
    
    return strengthScale

  @abstractmethod
  def getMatDir(self, pt):
    if self.defaultToA:
      return self.subObjA.getMatDir(pt)
    else:
      return self.subObjB.getMatDir(pt)

  @abstractmethod
  def getSubregions(self):
    subregions = self.subObjA.getSubregions()
    subregions.extend(self.subObjB.getSubregions())
    return subregions
  
  @abstractmethod
  def getVelocity(self, pt):
    obj = self.subObjA
    if not self.defaultToA:
      obj = self.subObjB

    if hasattr(obj,'getVelocity'):
        vel = obj.getVelocity( pt )
    else:
        if hasattr(obj, 'vel'):
          vel = obj.vel(pt) if callable(obj.vel) else obj.vel
        else:
          vel = _defaultVelocity
    return vel

  # @abstractmethod
  # def getDamage(self, pt):
  #   pass

  # @abstractmethod
  # def getPorosity(self, pt):
  #   pass

  # @abstractmethod
  # def getTemperature(self, pt):
  #   pass

  # @abstractmethod
  # def getSurfaceTraction(self, pt):
  #   pass


#############################################
class union(SetOperation):
  """
  Geometry object for creating a union of two objects
  """
  def __init__(self,
               name,
               A,
               B):
    SetOperation.__init__(self,
                          name,
                          A,
                          B)

  def isInterior(self, pt, skinDepth):
    self.intA = self.subObjA.isInterior(pt, skinDepth)
    self.intB = self.subObjB.isInterior(pt, skinDepth)

    intA = self.intA
    intB = self.intB

    if intA==0 or intB==0:
      return 0
    
    if intA==2 or intB==2:
      return 2

    return -1 # Shouldn't reach here, if it did something went wrong

  def getSurfaceNormal(self, pt):
    intA = self.intA
    intB = self.intB
    
    surf_norm = np.array([0.0, 0.0, 0.0])

    if intA == 2 and intB == 2:
      surf_norm = (self.subObjA.getSurfaceNormal(pt) + self.subObjB.getSurfaceNormal(pt))/2.0

    if intA == 2 and intB < 2:
      surf_norm = self.subObjA.getSurfaceNormal(pt)

    if intB == 2 and intA < 2:
      surf_norm = self.subObjB.getSurfaceNormal(pt)

    return surf_norm

  def getSurfacePosition(self, pt):
    intA = self.intA
    intB = self.intB
    
    surf_pos = np.array([0.0, 0.0, 0.0])

    if intA == 2 and intB == 2:
      surf_pos = (self.subObjA.getSurfacePosition(pt) + self.subObjB.getSurfacePosition(pt))

    if intA == 2 and intB < 2:
      surf_pos = self.subObjA.getSurfacePosition(pt)

    if intB == 2 and intA < 2:
      surf_pos = self.subObjB.getSurfacePosition(pt)

    return surf_pos


#############################################
class intersection(SetOperation):
  """
  Geometry object for creating an intersection of two objects
  """
  def __init__(self,
               name,
               A,
               B):
    SetOperation.__init__(self,
                          name,
                          A,
                          B)

  def isInterior(self, pt, skinDepth):
    sA = self.subObjA.isInterior(pt, skinDepth)
    sB = self.subObjB.isInterior(pt, skinDepth)
    if sA >= 0 and sB >= 0:
      return max(sA, sB) # Should work to always return correct surface flag
    
    return -1

  def getSurfaceNormal(self, pt):
    # return _defaultSurfaceNormal
    # # Cases
    # # Both A and B define surface normals
    # # Only one defines surface normals
    # # Both define a surface normal for the same voxel (use surface normal of )
    
    sA = self.subObjA.getSurfacePosition(pt)
    sB = self.subObjB.getSurfacePosition(pt)

    if np.dot(sA, sA) <= np.dot(sB, sB):
      return self.subObjA.getSurfaceNormal(pt)
    else:
      return self.subObjB.getSurfaceNormal(pt)

  def getSurfacePosition(self, pt):
    sA = self.subObjA.getSurfacePosition(pt)
    sB = self.subObjB.getSurfacePosition(pt)

    if np.dot(sA, sA) <= np.dot(sB, sB):
      return sA
    else:
      return sB
  
  def getCZTag(self, pt):
    return self.subObjA.getCZTag(pt) # Temporarily hardcoded


#############################################
class difference(SetOperation):
  """
  Geometry object for the set difference A - B.

  The retained material points are those that are inside object A and outside
  object B.  The SetOperation base stores the operands as subObjA/subObjB; this
  class therefore intentionally uses those names rather than the legacy
  self.A/self.B attributes.

  Surface information is approximated from the closest retained boundary.  On
  the outer boundary of A, normals/positions come from A.  On the cut boundary
  produced by subtracting B, the surface position comes from B and the normal is
  reversed because the outward normal of A - B points into the removed region.
  """
  def __init__(self,
               name,
               A,
               B):
    SetOperation.__init__(self,
                          name,
                          A,
                          B)
    self.intA = -1
    self.intB = -1

  def isInterior(self, pt, skinDepth):
    self.intA = self.subObjA.isInterior(pt, skinDepth)
    self.intB = self.subObjB.isInterior(pt, skinDepth)

    # True set difference: keep A only where B is absent.  The legacy XOR logic
    # incorrectly retained points that were inside B but outside A.
    if self.intA < 0:
      return -1

    if self.intB >= 0:
      return -1

    # Preserve surface flags from A and also flag points close to the cut
    # surface of B when the sub-object can provide a surface-position vector.
    try:
      sB = self.subObjB.getSurfacePosition(pt)
      if np.linalg.norm(sB) <= skinDepth:
        return _defaultSurfaceFlag
    except Exception:
      pass

    if self.intA > 0:
      return self.intA

    return 0

  def _surfacePositionFromA(self, pt):
    try:
      return self.subObjA.getSurfacePosition(pt)
    except Exception:
      return _defaultSurfacePosition

  def _surfacePositionFromB(self, pt):
    try:
      return self.subObjB.getSurfacePosition(pt)
    except Exception:
      return _defaultSurfacePosition

  def _closestSurfaceIsB(self, pt):
    sA = self._surfacePositionFromA(pt)
    sB = self._surfacePositionFromB(pt)
    return np.dot(sB, sB) < np.dot(sA, sA)

  def getSurfaceNormal(self, pt):
    if self._closestSurfaceIsB(pt):
      try:
        return -self.subObjB.getSurfaceNormal(pt)
      except Exception:
        return _defaultSurfaceNormal

    try:
      return self.subObjA.getSurfaceNormal(pt)
    except Exception:
      return _defaultSurfaceNormal

  def getSurfacePosition(self, pt):
    if self._closestSurfaceIsB(pt):
      return self._surfacePositionFromB(pt)

    return self._surfacePositionFromA(pt)

  def xMin(self):
    return self.subObjA.xMin()

  def xMax(self):
    return self.subObjA.xMax()
# ===========================================
# END SET OPERATIONS
# ===========================================
