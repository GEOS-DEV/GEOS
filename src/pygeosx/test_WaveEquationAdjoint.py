import pytest
import os
import sys
import numpy as np
from mpi4py import MPI
from pygeosx import initialize,apply_initial_conditions,_finalize,pylvarray

comm = MPI.COMM_WORLD
nRank= comm.Get_size()
rank = comm.Get_rank()


def oneStepForward(solver, time, dt, save_snapshot, shot_index = 0):
    oneStep(solver, time, dt, True, save_snapshot, shot_index)

def oneStepBackward(solver, time, dt, save_snapshot, shot_index = 0):
    oneStep(solver, time, dt, False, save_snapshot, shot_index)

def oneStep(solver, time, dt, isForward, save_snapshot, shot_index = 0):
    solver.get_wrapper("forward").value()[0] = 1 if isForward else 0
    solver.get_wrapper("saveFields").value()[0] = 1 if save_snapshot else 0
    solver.get_wrapper("shotIndex").value()[0] = shot_index
    solver.execute(time, dt)

def saveGradient(geosx):
    geosx.get_group("Tasks/partialGradientCollection").collect(0.0, 0)
    geosx.get_group("Outputs/partialGradientOutput").output(0.0, 0)
    geosx.get_group("Tasks/partialGradientCollection1").collect(0.0, 0)
    geosx.get_group("Outputs/partialGradientOutput1").output(0.0, 0)


def getPressureAtReceivers(solver):

    """Get the pressures values at receivers coordinates

    Return
    ------
    pressureAtReceivers : Numpy Array
        Array containing the pressure values at all time step at all receivers coordinates
    """

    pressureAtReceivers = solver.get_wrapper("pressureNp1AtReceivers").value()
    return pressureAtReceivers.to_numpy()

def getPressure():
    return solver.get_wrapper("/domain/MeshBodies/mesh/meshLevels/Level0/nodeManager/pressure_n").value().to_numpy()

def updateSourceAndReceivers(solver, sources_list=[], receivers_list=[]):

    """Update Sources and Receivers positions in GEOSX.

       Parameters
       ----------
       sources_list : list
           List of coordinates for the sources

       receivers_list : list
           List of coordinates for the receivers
       """

    src_pos_geosx = solver.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pylvarray.RESIZEABLE)

    rcv_pos_geosx = solver.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pylvarray.RESIZEABLE)


    src_pos_geosx.resize(len(sources_list))
    if len(sources_list) == 0:
        src_pos_geosx.to_numpy()[:] = np.zeros((0,3))
    else:
        src_pos_geosx.to_numpy()[:] = sources_list[:]

    rcv_pos_geosx.resize(len(receivers_list))
    if len(receivers_list) == 0:
        rcv_pos_geosx.to_numpy()[:] = np.zeros((0,3))
    else:
        rcv_pos_geosx.to_numpy()[:] = receivers_list[:]

    solver.reinit()

def resetWaveField(solver,geosx):

    """Reinitialize all pressure values on the Wavefield to 0.
    """

    meshName="mesh"
    solver.get_wrapper("indexSeismoTrace").value()[0] = 0

    nodeManagerPath = "domain/MeshBodies/"+meshName+"/meshLevels/Level0/nodeManager/"

    pressure_nm1 = geosx.get_wrapper(nodeManagerPath + "pressure_nm1").value()
    pressure_nm1.set_access_level(pylvarray.MODIFIABLE)

    pressure_n = geosx.get_wrapper(nodeManagerPath + "pressure_n").value()
    pressure_n.set_access_level(pylvarray.MODIFIABLE)

    pressure_np1 = geosx.get_wrapper(nodeManagerPath + "pressure_np1").value()
    pressure_np1.set_access_level(pylvarray.MODIFIABLE)

    pressure_nm1.to_numpy()[:] = 0.0
    pressure_n.to_numpy()[:]   = 0.0
    pressure_np1.to_numpy()[:] = 0.0

def getSourceValuesSize(solver):
    """ Retrieve the number of time steps in the source """
    src_value = solver.get_wrapper("sourceValue").value()
    return src_value.to_numpy().shape[0]
def getSource(solver):
    """ Retrieve the number of time steps in the source """
    src_value = solver.get_wrapper("sourceValue").value()
    return src_value.to_numpy()[:]

def updateSourceValue(solver, value):

    """Update the value of the source in GEOSX

    Parameters
    ----------
    value : array/list
        List/array containing the value of the source at each time step
    """

    src_value = solver.get_wrapper("sourceValue").value()
    src_value.set_access_level(pylvarray.MODIFIABLE)

    src_value.to_numpy()[:] = value[:]

def gatherSeismos(solver):
    comm.Barrier()
    # get sismos.
    sismos = getPressureAtReceivers(solver)
    n1=sismos.shape[0]
    n2=sismos.shape[1]
    # send sismos to rank 0
    tmp=np.zeros((n1,n2))
    tmp[:,:]=sismos[:,:]
    if rank !=0:
        comm.Send(tmp,0,1)
    else:
        tmp1=np.zeros((n1,n2))
        for i in range(1,nRank):
            comm.Recv(tmp1,i,1)
            for j in range(0,n2):
                if np.abs(np.min(tmp[:,j]))<0.0000001 and np.abs(np.max(tmp[:,j]))<0.00000001: 
                    tmp[:,j]+=tmp1[:,j]

    comm.Barrier()
    return tmp

def print_group(group, indent=0):
    print("{}{}".format(" " * indent, group))

    indent += 4
    print("{}wrappers:".format(" " * indent))

    for wrapper in group.wrappers():
        print("{}{}".format(" " * (indent + 4), wrapper))
        print_with_indent(str(wrapper.value(False)), indent + 8)

    print("{}groups:".format(" " * indent))

    for subgroup in group.groups():
        print_group(subgroup, indent + 4)

def print_with_indent(msg, indent):
    indent_str = " " * indent
    print(indent_str + msg.replace("\n", "\n" + indent_str))

def seismo_to_txt(seismo, filename):
    with open(filename,  "w") as f:
        for x in np.nditer(seismo):
            f.write(f"{x}\n")

def test_WaveEquationAcousticSolverAdjoint():
    argv = [ "test_WaveEquationAcousticSolverAdjoint.py", "-i", "test_WaveEquationAcousticSolverAdjoint.xml",
             "-z", str(nRank) ]
    geosx=initialize(rank,argv)
    #print_group(geosx)
    #exit(1)
    # get parameters form XML file
    solver=geosx.get_group("/Solvers/acousticSolver")
    srcPos=solver.get_wrapper("sourceCoordinates").value()
    recvPos=solver.get_wrapper("receiverCoordinates").value()
    event=geosx.get_group("/Events")
    maxTime=event.get_wrapper("maxTime").value()
    eventSolver=geosx.get_group("/Events/solverApplications")
    dt=eventSolver.get_wrapper("forceDt").value()

    vel_info=False
    if vel_info: 
        # get velocity infornation: min,max.
        objectPath="Region" 
        meshName="modelbicouche"
        velocity =solver.get_wrapper("/domain/MeshBodies/"+meshName+"/meshLevels/Level0/ElementRegions/elementRegionsGroup/"+objectPath+"/elementSubRegions/hexahedra/mediumVelocity").value()
        minVal=np.zeros(1)
        maxVal=np.zeros(1)
        minVal[0]=np.min(velocity.to_numpy())
        maxVal[0]=np.max(velocity.to_numpy())
        comm.Barrier()
        comm.Allreduce(MPI.IN_PLACE, maxVal, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, minVal, op=MPI.MIN)
        comm.Barrier()

    if rank ==0 :
       print()
       print("===============================================================")
       print(" Running pygeosx with solver:",solver)
       print(" Source coordinatesi : ",srcPos.to_numpy())
       print(" Receiver coordinates: ",recvPos.to_numpy())
       print(" simulation duration : ",maxTime)
       print(" Time Step           : ",dt)        
       if vel_info:
           print(" Min Vel             : ",minVal)
           print(" Max Vel             : ",maxVal)
       print("===============================================================")

    apply_initial_conditions()
    time=0
    cycle=0
    nsteps = maxTime / dt
    src_pos = [1000.01, 1000.01, 1000.01]
    rcv_pos = [1060.01, 1060.01, 1060.01 ]
    src_pos2 = [1000.01, 1000.01, 1000.01]
    rcv_pos2 = [1060.01, 1060.01, 1060.01 ]
    updateSourceAndReceivers(solver, sources_list = [src_pos], receivers_list = [rcv_pos])
    source1 = []
    source2 = []
    source = getSource(solver)
#    sourceSize = getSourceValuesSize(solver)
#    while sourceSize > 0:
#        source1.append([10.0])
#        source2.append([10.0])
#        sourceSize-=1
    source1 = source
    source2 = source*2
    print(source)
    #source1[1] = [10.0]
    #source2[1] = [30.0]
    updateSourceValue(solver, source1)

    while time < maxTime:
        if rank == 0 and cycle%10 == 0:
            print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", cycle+1)
        oneStepForward(solver, time,dt, False, 0)
        time+=dt
        cycle+=1
    solver.cleanup(maxTime)
    seismos1 = gatherSeismos(solver)
    seismo_to_txt(seismos1, "seismo1.txt")
    seismo_to_txt(source1, "source1.txt")
    seismo_to_txt(source2, "source2.txt")

    time=0
    cycle=0
    resetWaveField(solver,geosx)
    updateSourceAndReceivers(solver, sources_list = [src_pos], receivers_list = [rcv_pos])
    updateSourceValue(solver, source1)
    while time < maxTime:
        if rank == 0 and cycle%10 == 0:
            print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", cycle+1)
        oneStepForward(solver, time,dt, False, 0)
        time+=dt
        cycle+=1
    solver.cleanup(maxTime)
    seismos11 = gatherSeismos(solver)
    seismo_to_txt(seismos11, "seismo1bis.txt")

#    exit(1)
    resetWaveField(solver,geosx)
    updateSourceAndReceivers(solver, sources_list = [src_pos2], receivers_list = [rcv_pos2])
    updateSourceValue(solver, source2)
    dir(solver)

    while time>=0:
        if rank == 0 and cycle%10 == 0:
            print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", cycle)
        oneStepBackward(solver, time,dt, False, 0)
        time-=dt
        cycle-=1
    solver.cleanup(maxTime)

    seismos2 = gatherSeismos(solver)
    seismo_to_txt(seismos2, "seismo2_rev.txt")
    seismos2_rev = np.copy(seismos2)
    for i in range(seismos2.shape[0]):
        seismos2_rev[-i] = seismos2[i]

    seismo_to_txt(seismos2_rev, "seismo2_fixed.txt")
    print(seismos2[0:10])
    print(seismos2[-10:-1])
    print(seismos2_rev[0:10])
    print(seismos2_rev[-10:-1])
    prod21 = np.dot(np.array(seismos2_rev, float).transpose(), np.array(source1, float))
    prod12 = np.dot(np.array(seismos1, float).transpose(), np.array(source2, float))
    print(prod12)
    print(prod21)
    assert prod21 == pytest.approx(prod12)

    _finalize()

if __name__ == "__main__":
    test_WaveEquationAcousticSolverAdjoint()
