import pygeosx
import numpy as np
import os
import h5py
import segyio


def forward(solver, shot, outputWaveField, residual_flag=False, datafile=None, rank=0):

    """Solves forward problem, and compute the residual and costFunction

    Parameters
    ----------
    solver : AcousticSolver
        AcousticSolver object

    shot : Shot
        Contains all informations on current shot

    outputWaveField : int
        Output interval for the Wavefield

    rank : int
        Process rank


    Return
    ------
    residual : array
        Contains the residual between the simulation result and the real data
    """

    if residual_flag:
        if datafile is None:
            raise ValueError("If 'residual_flag' is set to True, you must sepecify 'datafile' for residual computation.")

    if rank == 0 :
        print("\nForward")

    costDir = "partialCostFunction"

    #Time loop
    time = 0.0
    dt = solver.dt
    cycle = 0
    shot.flag = "In Progress"

    while time < solver.maxTime:
        if rank == 0:
            print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", cycle+1)

        #Execute one time step
        solver.execute(time)

        time += dt
        cycle += 1

        #Collect/export waveField values to hdf5
        if not cycle % outputWaveField:
            if not rank:
                print("Outputing Wavefield...")
            solver.outputWaveField(time)

    residual = 0
    #Residual computation
    if residual_flag:
        pressure = solver.getPressureAtReceivers()
        residualTemp = computeResidual(datafile, pressure)
        residual = residualLinearInterpolation(residualTemp, solver.maxTime, dt, solver.dtSeismo)

        if rank == 0:
            computePartialCostFunction(costDir, residualTemp, shot)

    #Reset waveField
    solver.resetWaveField()

    return residual


def backward(solver, shot, outputWaveField, rank=0):

    """Solves backward problem

    Parameters
    ----------
    solver : AcousticSolver
        AcousticSolver object

    shot : Shot
        Contains all informations on current shot

    outputWaveField : int
        Output interval for the Wavefield

    rank : int
        Process rank
    """

    if rank == 0:
        print("\nbackward")

    #Reverse time loop
    time = solver.maxTime
    dt = solver.dt
    cycle = int(solver.maxTime / dt)
    while time > 0:
        if rank == 0:
            print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", cycle)

        #Execute one time step backward
        solver.execute(time)

        #Collect/export waveField values to hdf5
        if not cycle % outputWaveField:
            if rank == 0:
                print("Outputing Wavefield...")
            solver.outputWaveField(time)

        time -= dt
        cycle -= 1


def computeResidual(filename, data_obs):

    """Compute the difference between the result of the forward problem and the real data

    Parameters
    ----------
    filename : str
        Path to the observed data file

    data_obs : array
        Simulated data

    Return
    ------
    residual : array
        Contains the residual between the simulated data and the observed data
    """

    residual = np.zeros(data_obs.shape)
    with segyio.open(filename, 'r+', ignore_geometry=True) as f:
        for i in range(data_obs.shape[1]):
            residual[:, i] = f.trace[i] - data_obs[:, i]
        f.close()

    return residual


def residualLinearInterpolation(rtemp, maxTime, dt, dtSeismoTrace):

    """Interpolate the residual for each simulation timestep

    Parameters
    ----------
    rtemp : array
        Residual at dtSeismoTrace timestep

    maxTime : fload
        End time of simulation

    dt : float
        Time step of simulation

    dtSeismoTrace : float
        Time step of output seismoTrace

    Return
    ------
    residual : array
        Contains the residual for each simulation time step
    """

    r = np.zeros((int(maxTime/dt)+1, np.size(rtemp,1)))
    for i in range(np.size(rtemp,1)):
        r[:,i] = np.interp(np.linspace(0, maxTime, int(maxTime/dt)+1), np.linspace(0, maxTime, int(maxTime/dtSeismoTrace)+1), rtemp[:,i])

    return r


def computePartialCostFunction(directory_in_str, residual, shot):

    """Compute cost function for the current shot.

    Parameters
    ----------
    directory_in_str : str
        Path to the directory where to ouput result

    residual : array
        Contains value of residual

    shot : Shot
        Current shot
    """

    if os.path.exists(directory_in_str):
        pass
    else:
        os.mkdir(directory_in_str)

    partialCost = 0
    for i in range(residual.shape[1]):
        partialCost += np.linalg.norm(residual[:, i])

    partialCost = partialCost/2

    with open(os.path.join(directory_in_str, "partialCost_"+shot.id+".txt"), 'w') as f:
        f.write(str(partialCost))

    f.close()


def computeFullCostFunction(directory_in_str, acquisition):

    """Compute the full cost function by summing all partial cost function

    Parameters
    ----------
    directory_in_str : str
        Path to the directory containing all partial cost functions value

    acquisition : Acquisition
        Contains informations on the full acquisition

    Return
    ------
    fullCost : float
        Sum of all partial cost functions values
    """

    directory = os.fsencode(directory_in_str)
    nfiles = len(acquisition.shots)

    n=0
    fullCost=0
    while n < nfiles:
        file_list = os.listdir(directory)
        if len(file_list) != 0:
            filename = os.fsdecode(file_list[0])
            costFile = open(os.path.join(directory_in_str,filename), "r")
            fullCost += float(costFile.readline())
            costFile.close()
            os.remove(os.path.join(directory_in_str,filename))
            n += 1
        else:
            continue

    os.rmdir(directory_in_str)

    return fullCost


def computePartialGradient(directory_in_str, shot):

    """Compute partial gradient for the current shot.

    Parameters
    ----------
    directory_in_str : str
        Path to the directory where to ouput result

    shot : Shot
        Current shot
    """

    if os.path.exists(directory_in_str):
        pass
    else:
        os.mkdir(directory_in_str)

    h5fNp1 = h5py.File(os.path.join(directory_in_str,"forwardWaveFieldNp1_"+shot.id+".hdf5"), "r+")
    h5bNp1 = h5py.File(os.path.join(directory_in_str,"backwardWaveFieldNp1_"+shot.id+".hdf5"), "r+")

    h5fN = h5py.File(os.path.join(directory_in_str,"forwardWaveFieldN_"+shot.id+".hdf5"), "r+")
    h5bN = h5py.File(os.path.join(directory_in_str,"backwardWaveFieldN_"+shot.id+".hdf5"), "r+")

    h5fNm1 = h5py.File(os.path.join(directory_in_str,"forwardWaveFieldNm1_"+shot.id+".hdf5"), "r+")
    h5bNm1 = h5py.File(os.path.join(directory_in_str,"backwardWaveFieldNm1_"+shot.id+".hdf5"), "r+")

    keysNp1 = list(h5fNp1.keys())
    keysN = list(h5fN.keys())
    keysNm1 = list(h5fNm1.keys())

    for i in range(len(keysNp1)-2):
        if ("ReferencePosition" in keysNp1[i]) or ("Time" in keysNp1[i]):
            del h5fNp1[keysNp1[i]]
            del h5bNp1[keysNp1[i]]

            del h5fN[keysNp1[i]]
            del h5bN[keysNp1[i]]

            del h5fNm1[keysNp1[i]]
            del h5bNm1[keysNp1[i]]

    keysNp1 = list(h5fNp1.keys())

    h5p = h5py.File(os.path.join(directory_in_str,"partialGradient_"+shot.id+".hdf5"), "w")

    h5p.create_dataset("ReferencePosition", data = h5fNp1[keysNp1[-2]][0], dtype='d')
    h5p.create_dataset("Time", data = h5fNp1[keysNp1[-1]], dtype='d')
    h5p.create_dataset("partialGradient", data = np.zeros(h5bNp1[keysNp1[0]][0].shape[0]), dtype='d')


    keyp = list(h5p.keys())

    n = len( list( h5fNp1[keysNp1[0]] ) )
    for i in range(n):
        h5p["partialGradient"][:] -= (h5fNp1["pressure_np1"][i, :] - 2*h5fN["pressure_n"][i, :] + h5fNm1["pressure_nm1"][i, :]) * h5bN["pressure_n"][-i-1, :] / (shot.dt*shot.dt)


    h5fNp1.close()
    h5bNp1.close()
    h5fN.close()
    h5bN.close()
    h5fNm1.close()
    h5bNm1.close()
    h5p.close()

    os.remove(os.path.join(directory_in_str,"forwardWaveFieldNp1_"+shot.id+".hdf5"))
    os.remove(os.path.join(directory_in_str,"backwardWaveFieldNp1_"+shot.id+".hdf5"))
    os.remove(os.path.join(directory_in_str,"forwardWaveFieldN_"+shot.id+".hdf5"))
    os.remove(os.path.join(directory_in_str,"backwardWaveFieldN_"+shot.id+".hdf5"))
    os.remove(os.path.join(directory_in_str,"forwardWaveFieldNm1_"+shot.id+".hdf5"))
    os.remove(os.path.join(directory_in_str,"backwardWaveFieldNm1_"+shot.id+".hdf5"))


def computeFullGradient(directory_in_str, acquisition):

    """Compute the full gradient by summing all partial gradients

    Parameters
    ----------
    directory_in_str : str
        Path to the directory containing all partial gradients

    acquisition : Acquisition
        Contains informations on the full acquisition
    """

    directory = os.fsencode(directory_in_str)

    limited_aperture_flag = acquisition.limited_aperture
    nfiles = len(acquisition.shots)

    totalNodes = (acquisition.nx+1)*(acquisition.ny+1)*(acquisition.nz+1)

    h5F = h5py.File("fullGradient.hdf5", "w")
    h5F.create_dataset("fullGradient", data = np.zeros(totalNodes), dtype='d', chunks=True, maxshape=(totalNodes,))
    h5F.create_dataset("ReferencePosition", data = np.zeros((totalNodes, 3)), chunks=True, maxshape=(totalNodes, 3))
    h5F.create_dataset("Time", data = np.zeros(1), chunks=True, maxshape=(None,))

    n=0
    while True:
        file_list = os.listdir(directory)
        if n == nfiles:
            break
        elif len(file_list) > 0:
            for file in file_list:
                filename = os.fsdecode(file)
                if filename.startswith("partialGradient"):
                    h5p = h5py.File(os.path.join(directory_in_str, filename), 'r')
                    keysp = list(h5p.keys())

                    if limited_aperture_flag:
                        shotInd = int(os.path.join(directory_in_str, filename).rsplit('.', 1)[0][-3::]) - 1
                        shot = acquisition.shots[shotInd]

                        cx = (shot.boundary[0][1] - shot.boundary[0][0])/shot.nx
                        cy = (shot.boundary[1][1] - shot.boundary[1][0])/shot.ny
                        cz = (shot.boundary[2][1] - shot.boundary[2][0])/shot.nz


                        for i in range(h5p["ReferencePosition"].shape[0]):
                            coord = h5p["ReferencePosition"][i]

                            xi = int(coord[0]/cx)
                            yi = int(coord[1]/cy)
                            zi = int(coord[2]/cz)

                            ind = zi + yi*(acquisition.nz+1) + xi*(acquisition.nz+1)*(acquisition.ny+1)

                            h5F["fullGradient"][ind] += h5p["partialGradient"][i]
                            h5F["ReferencePosition"][ind] = coord

                    else:
                        h5F["fullGradient"][:] += h5p["partialGradient"][:]

                    if n == 0:
                       h5F["Time"].resize(h5p["Time"].shape[0], axis=0)
                    h5F["Time"][:] = h5p["Time"][:,0]

                    h5p.close()
                    os.remove(os.path.join(directory_in_str,filename))
                    n+=1
                else:
                    continue

        else:
            continue

    h5F.close()
    os.rmdir(directory_in_str)

    return h5F



def armijoCondition(J, Jm, d, gradJ, alpha, k1=0.8):

    """Armijo condition checking

    Parameters
    ----------
    J : float
        Value of the cost function

    Jm : float
        Value of the cost function for perturbed model

    d : array
        Descent direction

    gradJ : array
        Gradient of the cost function

    alpha : float
        Coefficient

    k1 : float
        Coefficient

    Return
    ------
    success : bool
        0 if test fails, 1 if it passes
    """

    success = (Jm <= k1 * alpha * np.dot(d, gradJ))
    return success

def goldsteinCondition(J, Jm, d, gradJ, alpha, k2=0.2):

    """Armijo condition checking

    Parameters
    ----------
    J : float
        Value of the cost function

    Jm : float
        Value of the cost function for perturbed model

    d : array
        Descent direction

    gradJ : array
        Gradient of the cost function

    alpha : float
        Coefficient

    k1 : float
        Coefficient

    Return
    ------
    success : bool
        0 if test fails, 1 if it passes
    """

    success = (Jm >= k2 * alpha * np.dot(d, gradJ))
    return success


def descentDirection(gradJ, method="steepestDescent"):

    """Compute the descent direction

    Parameters
    ----------
    gradJ : array
        Gradient of the cost function

    method : str
        Method to use for descent direction calculation
    """

    if method == "steepestDescent":
        d = -gradJ
    else:
        print("to be implemented")

    return d


def linesearchDirection(J, gradJ, alpha, d, m, p1=1.5, p2=0.5):
    while True:
        m1 = m + alpha * d
        J1 = computeCostFunction(acqs, m1, client)
        succesArmijo = armijoCondition(J, J1, gradJ, alpha)
        if succesArmijo:
            succesGoldstein = goldsteinCondition(J, J1, gradJ, alpha)
            if succesGoldstein:
                return m1
            else:
                alpha*=p1
        else:
            alpha*=p2


def updateModel(acqs, m):
    print("To be implemented")


#==============================================================================================#


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
