import segyio
import os

def create_segy(shot_list, nsamples, tracePath):
   
    if os.path.exists(tracePath):
        pass
    else:
        os.mkdir(tracePath)

    for ishot in range(len(shot_list)):
        rcvCoord = shot_list[ishot].getReceiverSet().getSetCoord()
        srcCoord = shot_list[ishot].getSource().getCoord()
        spec = segyio.spec()
        ilines  = [x[0] for x in rcvCoord]
        xlines  = [x[1] for x in rcvCoord]
        spec.ilines = list(set(ilines))
        spec.xlines = list(set(xlines))
        spec.samples = list(range(nsamples))
        spec.sorting = 2
        spec.format  = 1

        with segyio.create(tracePath + "/sismoTraceShot"+str(ishot)+".sgy", spec) as f:
            for i in range(len(rcvCoord)):
                f.header[i] = {segyio.su.scalco : -100,
                               segyio.su.scalel : -100,
                               segyio.su.sx : int(srcCoord[0]*100),
                               segyio.su.sy : int(srcCoord[1]*100),
                               segyio.su.sdepth : int(srcCoord[2]*100),
                               segyio.su.gx : int(rcvCoord[i][0]*100),
                               segyio.su.gy : int(rcvCoord[i][1]*100),
                               segyio.su.gelev : int(rcvCoord[i][2]*100)}
                f.trace[i] = [0.0] * nsamples
            f.bin.update(tsort = segyio.TraceSortingFormat.INLINE_SORTING)


def export_to_segy(pressure, rcvCoord, ishot, tracePath):
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
    
    with segyio.open(tracePath + "/sismoTraceShot"+str(ishot)+".sgy", 'r+', ignore_geometry=True) as f:
        for i in range(len(rcvCoord)):
            if any(pressure[:,i])==True:
                f.trace[i] = pressure[:, i]


def export_for_acquisition(shot_list, acq_path, acq_name):
    if os.path.exists(acq_path):
        pass
    else:
        os.mkdir(acq_path)

    if os.path.exists(acq_path + acq_name):
        pass
    else:
        os.mkdir(acq_path + acq_name)

    for j in range(len(shot_list)):
        srcCoord = shot_list[j].getSource().getCoord()
        rcvCoord = shot_list[j].getReceiverSet().getSetCoord()

        spec = segyio.spec()
        ilines  = [x[0] for x in rcvCoord]
        xlines  = [x[1] for x in rcvCoord]
        spec.ilines = list(set(ilines))
        spec.xlines = list(set(xlines))
        spec.samples = [0]
        spec.sorting = 2
        spec.format  = 1

        with segyio.create(acq_path + acq_name + "/acquisition_" + str(j) + ".sgy", spec) as f:
            for i in range(len(rcvCoord)):
                f.header[i] = {segyio.su.scalco : -100,
                               segyio.su.scalel : -100,
                               segyio.su.sx : int(srcCoord[0]*100),
                               segyio.su.sy : int(srcCoord[1]*100),
                               segyio.su.sdepth : int(srcCoord[2]*100),
                               segyio.su.gx : int(rcvCoord[i][0]*100),
                               segyio.su.gy : int(rcvCoord[i][1]*100),
                               segyio.su.gelev : int(rcvCoord[i][2]*100)}
                f.bin.update(tsort = segyio.TraceSortingFormat.INLINE_SORTING)
