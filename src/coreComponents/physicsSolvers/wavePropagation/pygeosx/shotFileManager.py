import ast

from shot import *
from source import *
from receiver import *

import os

def exportShotList(proc, shot_list):
    path = os.path.abspath(os.getcwd()) + "/shots_lists/"
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    shot_file = path + "shot_list" + str(proc) + ".txt"
    f = open(shot_file, 'w+')
    dt = shot_list[0].getSource().getTimeStep()
    wavelet = shot_list[0].getSource().getFunction()
    f.write(str(dt) + '\n')
    f.write(str(wavelet.tolist()) +'\n')
    rcvlist = []
    for shot in (shot_list):
        src = shot.getSource().getCoord().tolist()
        rcvs = shot.getReceiverSet().getSetCoord()
        for rcv in rcvs:
            rcvlist.append(rcv.tolist())

        f.write("{}\v{}\n".format(str(src), str(rcvlist)))
        rcvlist = []

    f.close()
    return shot_file

def readShotList(shot_file):
    path = os.path.abspath(os.getcwd()) + "/shots_lists/"
    f = open(shot_file, 'r')
    dt = float(f.readline())
    wavelet = ast.literal_eval(f.readline())
    shot_list = []
    for line in f.readlines():
        src = ast.literal_eval(line.splitlines()[0])
        rcvs = ast.literal_eval(line.splitlines()[1])
        shot_list.append(Shot(Source(src, wavelet, dt), ReceiverSet([Receiver(rcv) for rcv in(rcvs)])))

    f.close()
    os.remove(shot_file)
    os.rmdir(path)
    return shot_list
