import pygeosx
import numpy as np


def recomputeSourceAndReceivers(solver, sources, receivers):
    updateSourceAndReceivers(solver, sources, receivers)

    solver.reinit()


def updateSourceAndReceivers( solver, sources_list=[], receivers_list=[] ):
    src_pos_geosx = solver.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)

    rcv_pos_geosx = solver.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)


    src_pos_geosx.resize(len(sources_list))
    if len(sources_list) == 0:
        src_pos_geosx.to_numpy()[:] = np.zeros((0,3))
    else:
        src_pos = [source.coords for source in sources_list]
        src_pos_geosx.to_numpy()[:] = src_pos[:]

    rcv_pos_geosx.resize(len(receivers_list))
    if len(receivers_list) == 0:
        rcv_pos_geosx.to_numpy()[:] = np.zeros((0,3))
    else:
        rcv_pos = [receiver.coords for receiver in receivers_list]
        rcv_pos_geosx.to_numpy()[:] = rcv_pos[:]

    solver.reinit()

def updateSourceValue( solver, value ):
    src_value = solver.get_wrapper("sourceValue").value()
    src_value.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    src_value.to_numpy()[:] = value[:]


def residualLinearInterpolation(rtemp, maxTime, dt, dtSeismoTrace):
    r = np.zeros((int(maxTime/dt)+1, np.size(rtemp,1)))
    for i in range(np.size(rtemp,1)):
        r[:,i] = np.interp(np.linspace(0, maxTime, int(maxTime/dt)+1), np.linspace(0, maxTime, int(maxTime/dtSeismoTrace)+1), rtemp[:,i])

    return r


def resetWaveField(group):
    group.get_wrapper("Solvers/acousticSolver/indexSeismoTrace").value()[0] = 0
    nodeManagerPath = "domain/MeshBodies/mesh/Level0/nodeManager/"

    pressure_nm1 = group.get_wrapper(nodeManagerPath + "pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_n = group.get_wrapper(nodeManagerPath + "pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_np1 = group.get_wrapper(nodeManagerPath + "pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_geosx = group.get_wrapper("Solvers/acousticSolver/pressureNp1AtReceivers").value()
    pressure_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_nm1.to_numpy()[:] = 0.0
    pressure_n.to_numpy()[:]   = 0.0
    pressure_np1.to_numpy()[:] = 0.0
    pressure_geosx.to_numpy()[:] = 0.0



def setTimeVariables(problem, maxTime, dt, dtSeismoTrace):
    problem.get_wrapper("Events/maxTime").value()[0] = maxTime
    problem.get_wrapper("Events/solverApplications/forceDt").value()[0] = dt
    problem.get_wrapper("/Solvers/acousticSolver/dtSeismoTrace").value()[0] = dtSeismoTrace


def print_pressure(pressure, ishot):
    print("\n" + "Pressure value at receivers for configuration " + str(ishot) + " : \n")
    print(pressure)
    print("\n")

def print_shot_config(shot_list, ishot):
    print("\n \n" + "Shot configuration number " + str(ishot) + " : \n")
    print(shot_list[ishot])

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


def print_flag(shot_list):
    i = 0
    for shot in shot_list:
        print("Shot " + str(i) + " status : " + shot.flag)
        i += 1
    print("\n")
