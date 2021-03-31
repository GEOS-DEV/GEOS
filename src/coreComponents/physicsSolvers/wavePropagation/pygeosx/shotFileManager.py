import ast

from shot import *
from source import *
from receiver import *


def exportShotList(proc, shot_list):
    shot_file = "/home/m3d/Desktop/pygeosx/shot_list/shot_list_" + str(proc) + ".txt"
    f = open(shot_file, 'w+')
    wavelet = shot_list[0].getSource().getFunction()
    dt = shot_list[0].getSource().getTimeStep()
    f.write(str(dt) + '\n')
    f.write(str(wavelet.tolist()) +'\n')
    rcvlist = []
    for shot in (shot_list):
        src = shot.getSource().getCoord().tolist()
        rcvs = shot.getReceiverSet().getSetCoord()
        for rcv in rcvs:
            rcvlist.append(rcv.tolist())

        f.write("{}\v{}\n".format(str(src), str(rcvlist)))

    f.close()
    return shot_file

def readShotList(shot_file):
    print(shot_file)
    f = open(shot_file, 'r')
    dt = float(f.readline())
    wavelet = ast.literal_eval(f.readline())
    shot_list = []
    for line in f.readlines():
        src = ast.literal_eval(line.splitlines()[0])
        rcvs = ast.literal_eval(line.splitlines()[1])
        shot_list.append(Shot(Source(src, wavelet, dt), ReceiverSet([Receiver(rcv) for rcv in(rcvs)])))

    f.close()
    return shot_list
