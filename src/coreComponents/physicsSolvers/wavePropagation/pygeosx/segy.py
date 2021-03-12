import segyio 
import sys
import numpy as np



def export_to_segy(pressure, rcvCoord, ishot, dt_cycle):
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
       
       
'''
def export_to_segy(pressure, rcvCoord, ishot):
   spec = segyio.spec()
   spec.ilines  = [1, 2, 3, 4]   #il * xl = nb receptors 
   spec.xlines  = [11, 12, 13]
   spec.samples = list(range(50))   #nombre de dt saved
   spec.sorting = 2
   spec.format  = 1

   with segyio.create("seismic_trace_shot_i="+str(ishot)+".sgy", spec) as f:
       start = 0.0
       step  = 0.00001
		# fill a trace with predictable values: left-of-comma is the inline
		# number. Immediately right of comma is the crossline number
		# the rightmost digits is the index of the sample in that trace meaning
		# looking up an inline's i's jth crosslines' k should be roughly equal
		# to i.j0k
       trace = np.arange(start = start, 
                         stop  = start + step * len(spec.samples),
                         step  = step,
                         dtype = np.single)
       tr=0
       for il in spec.ilines:
           for xl in spec.xlines:
               f.header[tr] = {segyio.su.offset : 1, segyio.su.iline  : il, segyio.su.xline  : xl}
               f.trace[tr] = trace + (xl / 100.0) + il      #valeur de pression
               tr += 1
       print('Amplitude Inline range: ' + str(np.amin(f.ilines)) + ' - ' +str(np.amax(f.ilines))) 
       print('Amplitude Crossline range: ' + str(np.amin(f.xlines)) + ' - ' +str(np.amax(f.xlines)))
       f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING)
'''
