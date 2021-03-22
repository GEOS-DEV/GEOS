import segyio 
import sys
import numpy as np



def export_to_segy(pressure, rcvCoord, ishot, dt_cycle):
    """Export the pressure value calculated by GEOSX to a segy file
    
    Parameters
    ----------
    pressure :
        Numpy array with pressure values at all receivers coordinates
    
    rcvCoord :
        Receiver coordinates
    
    ishot :
        Number of current shot
   
    dt_cycle :
        Frequency of value export
    """ 
    
    spec = segyio.spec()
    ilines  = [x[0] for x in rcvCoord]   
    xlines  = [x[1] for x in rcvCoord] 
    spec.ilines = list(set(ilines))
    spec.xlines = list(set(xlines))
    spec.samples = list(range(int(pressure[:,0].size/dt_cycle)))   
    spec.sorting = 2
    spec.format  = 1

    with segyio.create("/home/m3d/Desktop/pygeosx/sismoTrace/sismoTraceShot"+str(ishot)+".sgy", spec) as f:
        for i in range(len(rcvCoord)):
            f.header[i] = {segyio.su.offset : 1, segyio.su.iline  :int(rcvCoord[i][0]*1000), segyio.su.xline  : int(rcvCoord[i][1]*1000)}
            f.trace[i]  = pressure[::dt_cycle,i]      #valeur de pression
        f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING)
       
