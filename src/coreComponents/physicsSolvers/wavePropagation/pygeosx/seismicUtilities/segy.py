import numpy as np
import os
import segyio
from time import sleep

def exportToSegy(table, shot, filename, directory, rank=0):
    segyfile = os.path.join(directory, filename+".sgy")
    nsamples = table.shape[0]
    rcvCoord = [receiver.coords for receiver in shot.receivers.receivers_list]

    if rank == 0:
        if os.path.exists(directory):
            pass
        else:
            os.mkdir(directory)

        srcCoord = shot.sources.source_list[0].coords
        spec = segyio.spec()
        ilines  = [x[0] for x in rcvCoord]
        xlines  = [x[1] for x in rcvCoord]
        spec.ilines = list(set(ilines))
        spec.xlines = list(set(xlines))
        spec.samples = np.arange(nsamples)
        spec.sorting = 2
        spec.format  = 1

        with segyio.create(segyfile, spec) as f:
            for i in range(len(rcvCoord)):
                f.header[i] = {segyio.su.scalco : -100,
                               segyio.su.scalel : -100,
                               segyio.su.sx : int(srcCoord[0]*100),
                               segyio.su.sy : int(srcCoord[1]*100),
                               segyio.su.sdepth : int(srcCoord[2]*100),
                               segyio.su.gx : int(rcvCoord[i][0]*100),
                               segyio.su.gy : int(rcvCoord[i][1]*100),
                               segyio.su.gelev : int(rcvCoord[i][2]*100)}
                f.trace[i] = np.zeros(nsamples, dtype=np.float32)
            f.bin.update(tsort = segyio.TraceSortingFormat.INLINE_SORTING)
            f.close()
    else:
        sleep(1)

    with segyio.open(segyfile, 'r+', ignore_geometry=True) as f:
        for i in range(len(rcvCoord)):
            if any(table[1:,i])==True:
                f.trace[i] = np.ascontiguousarray(table[:, i], dtype = np.float32)
        f.close()
