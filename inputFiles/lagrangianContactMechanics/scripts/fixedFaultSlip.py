import numpy as np
import os
import sys
import xml.etree.ElementTree as ElementTree
import matplotlib
import matplotlib.pyplot as plt
import argparse
import geos.hdf5_wrapper as hdf5_wrapper
    
class SingularCrackSlip:

    def __init__(self, mechanicalParameters, length ):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]
        poisson_ratio= (3 * K - 2 * G) / (2 * (3 * K + G))
        
        mu_star = G /( 1 - poisson_ratio)
        self.tau_0 = 0.0
        self.tau_r = -1.0

        self.scaling = 2*(self.tau_0 - self.tau_r)/mu_star
        self.halfLength = length

    def computeSlip(self, x):
        return self.scaling * np.sqrt(self.halfLength**2 - x**2)

    def computeTraction(self, x):
        if x < -self.halfLength or x > self.halfLength:
            return self.tau_0 + (self.tau_0-self.tau_r) * ( np.abs(x)/np.sqrt(x**2 - self.halfLength**2) - 1 )
        else:
            return self.tau_r
class GaussianSlip:

    def __init__(self, peakSlip, length ):
        self.scaling =  peakSlip
        self.halfLength = length

    def computeSlip(self, x):
        denom = 1 / (self.halfLength/2)
        return self.scaling*np.exp(-0.5*((x)/denom)**2)

class GaussianSlip2D:
    def __init__(self, peakSlip, lengthX, lengthY):
        self.scaling = peakSlip
        self.sigma_x = lengthX / 2
        self.sigma_y = lengthY / 2

    def computeSlip(self, x, y):
        exponent = -0.5 * ((x / self.sigma_x) ** 2 + (y / self.sigma_y) ** 2)
        return self.scaling * np.exp(exponent)  

class RadialSlip:
    def __init__(self, stressDrop, mechanicalParameters, crackRadius):
        self.deltaTau = stressDrop
        self.K = mechanicalParameters["bulkModulus"]
        self.G = mechanicalParameters["shearModulus"]
        self.poisson_ratio= (3 * self.K - 2 * self.G) / (2 * (3 * self.K + self.G))
        print(f'nu is {self.poisson_ratio}')
        self.R = crackRadius
        self.prefactor = (24 * self.R * self.deltaTau) / (7 * np.pi * self.G)

    def computeSlip(self, x, y):
        r = np.sqrt(x**2 + y**2)
        slip = np.zeros_like(r)
        inside = r <= self.R
        slip[inside] = self.prefactor * np.sqrt(1 - (r[inside] / self.R)**2)
        return slip

class SinusoidalSlip:
    def __init__(self, amplitude, wavenumber, mechanicalParameters, direction="x"):
        """
        Initialize sinusoidal slip pattern.

        Parameters:
        - amplitude (float): Maximum slip amplitude A
        - wavenumber (float): Wavenumber k (1/m)
        - mechanicalParameters (dict): Dictionary with 'shearModulus' and optionally 'bulkModulus'
        - direction (str): 'x' or 'y', direction of slip variation
        """
        self.A = amplitude
        self.k = wavenumber
        assert direction in ("x", "y"), "direction must be 'x' or 'y'"
        self.direction = direction

        self.G = mechanicalParameters["shearModulus"]
        self.K = mechanicalParameters.get("bulkModulus", None)

        # Compute Poisson's ratio if bulk modulus is provided
        self.nu = None
        if self.K is not None:
            self.nu = (3 * self.K - 2 * self.G) / (2 * (3 * self.K + self.G))
            print(f"Poisson's ratio (nu) = {self.nu:.4f}")

    def computeSlip(self, x, y):
        """
        Compute slip field over a grid defined by x, y.

        Returns:
        - slip (np.ndarray): slip magnitude at each (x, y)
        """
        coord = x if self.direction == "x" else y
        slip = self.A * np.cos( self.k * coord )
        return slip

    def computeStressChange(self, x, y, mode="antiplane"):
        """
        Compute shear stress change from sinusoidal slip.

        Parameters:
        - mode (str): "antiplane" or "planeStrain"

        Returns:
        - stressChange (np.ndarray): Shear stress change at each point
        """
        slip = self.computeSlip(x, y)

        if mode == "antiplane":
            prefactor = -0.5 * self.G * self.k
        elif mode == "planeStrain":
            if self.nu is None:
                raise ValueError("Plane strain mode requires bulkModulus to compute Poisson's ratio.")
            prefactor = -self.G / (2 * (1 - self.nu)) * self.k
        else:
            raise ValueError("Mode must be 'antiplane' or 'planeStrain'")

        return prefactor * slip
    


def getMechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.find('Constitutive/ElasticIsotropic')

    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
    return mechanicalParameters

def getFractureLengthFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    rectangle = tree.find('Geometry/Box')
    xmin = rectangle.get("xMin")
    xmax = rectangle.get("xMax")
    xmin = [float(i) for i in xmin[1:-1].split(",")]
    xmax = [float(i) for i in xmax[1:-1].split(",")]
    length = ( xmax[0] - xmin[0] ) / 2
    origin = 0.0
   
    return length, origin

def curve_check_solution(**kwargs):
    #-------- Extract info from XML
    xmlFilePath = f'./LagrangeContactBubbleStab_FixedSlip_base.xml'

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)

    # Get length of the fracture
    xmlFilePath = f'./LagrangeContactBubbleStab_FixedSlip_smoke.xml'
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)
    halfLength = 2.0

    x = kwargs['traction elementCenter']
    x_geos = x[0, :, 0]

    return analytical_solution(x, mechanicalParameters, totalHalfLength, halfLength)

def analytical_solution(x, mechanicalParameters, totalHalfLength, halfLength):

    singularCrackSlipSolution = SingularCrackSlip(mechanicalParameters, halfLength)
    x = np.linspace(-totalHalfLength, totalHalfLength, 10000, endpoint=True)
    traction_analytical = np.zeros(len(x))
    i = 0
    for xCell in x:
        traction_analytical[i] = singularCrackSlipSolution.computeTraction(xCell)
        i += 1  
    return traction_analytical

def plot_traction_solution_1d(inputFileDirectory, outputDirectory):
    # Read HDF5
    import hdf5_wrapper
    hdf5File1Path = f'outputDirectory/traction.hdf5'

    # Read HDF5
    data = hdf5_wrapper.hdf5_wrapper(hdf5File1Path).get_copy()
    traction = data['traction']
    traction = np.asarray(traction)
    traction_geos = traction[0, :, 1]
    x = data['traction elementCenter']
    x_geos = x[0, :, 0]

     #-------- Extract info from XML
    xmlFilePath = f'{inputFileDirectory}/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_base.xml'

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)

    # Get length of the fracture
    xmlFilePath = f'{inputFileDirectory}lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_smoke.xml'
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)
    halfLength = 2.0

    traction_analytical = analytical_solution(x, mechanicalParameters, totalHalfLength, halfLength)    

    fsize = 30
    msize = 15
    lw = 2
    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")
    
    # Plot analytical (continuous line) and numerical (markers) aperture solution
    ax.plot(x, traction_analytical, color='r', label='Traction analytical', lw=lw)
    ax.plot(x_geos, traction_geos, color='k', label='geos', marker="o", lw=lw)

    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('Shear traction', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.75, 0.9), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("traction.png")

def plot_traction_solution_2d(inputFileDirectory, outputDirectory):
    # Read GEOS traction file
    hdf5FilePath = f'{outputDirectory}/traction.hdf5'
    data = hdf5_wrapper.hdf5_wrapper(hdf5FilePath).get_copy()

    # Extract traction and coordinates
    traction = np.asarray(data['traction'])[0]              # Shape: (N, 3)
    centers = np.asarray(data['traction elementCenter'])[0] # Shape: (N, 3)

    # Extract traction and coordinates
    traction = np.asarray(data['traction'])[0]              # Shape: (N, 3)
    centers = np.asarray(data['traction elementCenter'])[0] # Shape: (N, 3)

    # Extract x, y, τx, τy, τz
    x = centers[:, 0]
    y = centers[:, 1]
    tau_x = traction[:, 1]
    tau_y = traction[:, 2]
    tau_z = traction[:, 0]

    # Assume structured grid: try to reshape automatically
    try:
        N = int(np.sqrt(len(x)))  # assumes square grid
        X = x.reshape((N, N))
        Y = y.reshape((N, N))
        TX = tau_x.reshape((N, N))
        TY = tau_y.reshape((N, N))
        TZ = tau_z.reshape((N, N))
    except ValueError:
        print("Could not reshape arrays to square grid. Check data.")
        return

    print(f'max traction change {np.max(TX)/1.e6}')
    print(f'min traction change {np.min(TX)/1.e6}')        

    # Plot contour plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    titles = ["τₓ [Pa]", "τᵧ [Pa]", "τ_z [Pa]"]
    fields = [TX/1.e6, TY/1.e6, TZ/1.e6]

    for ax, field, title in zip(axes, fields, titles):
        c = ax.contourf(X, Y, field, cmap='viridis')
        fig.colorbar(c, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_aspect('equal')

    plt.savefig("traction_2d_contours.png")
    plt.show()

def output_tables_1d(x, slip, name):
    # Save x to x.csv with one value per row
    np.savetxt('x.csv', x, fmt='%f')

    # Save aperture_analytical to jump.csv with one value per row
    np.savetxt(f'{name}.csv', slip, fmt='%f') 

def output_tables_2d(x, y, slip, name):
    # Save 2D coordinates and slip as CSVs
    np.savetxt('x_2d.csv', x, fmt='%f')
    np.savetxt('y_2d.csv', y, fmt='%f')
    np.savetxt(f'{name}.csv', slip, fmt='%f')         

def generate_tables_1d(inputFileDirectory):
    #-------- Extract info from XML
    xmlFilePath = f'{inputFileDirectory}/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_base.xml'

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    appliedPressure = 1.0

    # Get length of the fracture
    xmlFilePath = f'{inputFileDirectory}/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_smoke.xml'
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)
    halfLength = 2.0

    # Initialize analytical solution
    singularCrackSlipSolution = SingularCrackSlip(mechanicalParameters, halfLength )
    peakStrength = 3.0
    gaussianSlipSolution = GaussianSlip( peakStrength, halfLength)

    # Plot analytical (continuous line) and numerical (markers) aperture solution
    x = np.linspace(-totalHalfLength, totalHalfLength, 10000, endpoint=True)
    singularCrackSlip = np.zeros(len(x))
    gaussianSlip = np.zeros(N*N)
    i = 0
    for xCell in x:
        if xCell > -halfLength and xCell < halfLength:
            singularCrackSlip[i] = singularCrackSlipSolution.computeSlip(xCell)
        gaussianSlip[i] = gaussianSlipSolution.computeSlip(xCell)
        i += 1

    fsize = 24
    msize = 15
    lw = 6
    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")

    ax.plot(x, singularCrackSlip , color='k', label='Analytical Solution', lw=lw)
    ax.grid()
    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('slip [m]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.7, 1), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("singularCrackSlip.png")

    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")

    ax.plot(x, gaussianSlip , color='k', label='Analytical Solution', lw=lw)
    ax.grid()
    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('slip [m]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.75, 0.9), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("gaussianSlip.png")

    output_tables_1d(x, singularCrackSlip, "singularCrackSlip")
    output_tables_1d(x, gaussianSlip, "gaussianSlip")

def generate_tables_2d(inputFileDirectory):
    #-------- Extract info from XML
    xmlFilePath = f'{inputFileDirectory}/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_3d_base.xml'
    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    crackRadius = 25
    stressDrop = 1.0e6

    # Get length of the fracture
    xmlFilePath = f'{inputFileDirectory}/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_3d_smoke.xml'
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)

    # Initialize analytical solution
    peakSlip = 1.0e-1
    gaussianSlipSolution = GaussianSlip2D( 1.0e-1, crackRadius, crackRadius)
    radialSlipSolution = RadialSlip( stressDrop, mechanicalParameters, crackRadius )
    sinusoidalSlipSolution = SinusoidalSlip( amplitude=1.0e-1, wavenumber=0.068, mechanicalParameters=mechanicalParameters, direction="x"  )

    # 2D domain
    N = 1000  # grid resolution
    x = np.linspace(-totalHalfLength, totalHalfLength, N)
    y = np.linspace(-totalHalfLength, totalHalfLength, N)
    X, Y = np.meshgrid(x, y)

    gaussianSlip = gaussianSlipSolution.computeSlip(X, Y)
    radialSlip = radialSlipSolution.computeSlip(X, Y)
    sinusoidalSlip = sinusoidalSlipSolution.computeSlip(X, Y)

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.contourf(X, Y, gaussianSlip, levels=50, cmap='viridis')
    fig.colorbar(c, ax=ax, label="Slip [m]")

    ax.set_title("2D Gaussian Slip Distribution")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig("gaussianSlip_2d.png")
    plt.show()


    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.contourf(X, Y, radialSlip, levels=50, cmap='viridis')
    fig.colorbar(c, ax=ax, label="Slip [m]")

    ax.set_title("2D Radial Slip Distribution")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig("radialSlip.png")
    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(8, 6))
    c0 = ax[0].contourf(X, Y, sinusoidalSlip, levels=50, cmap='viridis')
    c1 = ax[1].contourf(X, Y, sinusoidalSlipSolution.computeStressChange(X, Y, "planeStrain"), levels=50, cmap='viridis')
    fig.colorbar(c0, ax=ax[0], label="Slip [m]")
    fig.colorbar(c1, ax=ax[1], label="Stress [MPa]")

    ax[0].set_title("Sinusoidal Slip Distribution")
    ax[0].set_xlabel("X [m]")
    ax[0].set_ylabel("Y [m]")
    ax[0].set_aspect('equal')
    ax[1].set_title("Tangential stress Distribution")
    ax[1].set_xlabel("X [m]")
    ax[1].set_ylabel("Y [m]")
    ax[1].set_aspect('equal')
    plt.tight_layout()
    plt.savefig("sinusoidalSlip.png")
    plt.show()
    
    # Compute 2D slip
    gaussianSlip_flat = np.zeros(N*N)
    radialSlip_flat = np.zeros(N*N)
    sinusoidalSlip_flat = np.zeros(N*N)
    i=0
    for xCell in x:
        for yCell in y:
            gaussianSlip_flat[i] = gaussianSlipSolution.computeSlip(xCell, yCell)
            radialSlip_flat[i] = radialSlipSolution.computeSlip(xCell, yCell)
            sinusoidalSlip_flat[i] = sinusoidalSlipSolution.computeSlip(xCell, yCell)
            i += 1

    output_tables_2d(x, y, gaussianSlip_flat, "gaussianSlip")
    output_tables_2d(x, y, radialSlip_flat, "radialSlip")
    output_tables_2d(x, y, sinusoidalSlip_flat, "sinusoidalSlip")         

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--action', type=str, choices=['generate_tables', 'plotTractions'], required=True, help='Action to perform: generate_tables or plotTractions')
    parser.add_argument('-i', '--input-files-path', type=str, required=True, help='Path to the inputFilesFolder')
    parser.add_argument('-o', '--output-dir', type=str, help='Directory containing the output files')
    parser.add_argument('-d', '--dimension', type=str, choices=['1d', '2d'], required=True, help='Type of table to generate: 1d or 2d')

    args = parser.parse_args()
    
    if args.action == 'generate_tables':
        print("Generating tables...")
        if args.dimension == '1d':
            generate_tables_1d(os.path.normpath(args.input_files_path))
        elif args.dimension == '2d':
            generate_tables_2d(os.path.normpath(args.input_files_path))
        else:
            print("Invalid dimension. Use '1d' or '2d'.")
            sys.exit(1)
    elif args.action == 'plotTractions':
        print("Plotting tractions...")
        if args.dimension == '1d':
            plot_traction_solution_1d(os.path.normpath(args.input_files_path), os.path.normpath(args.output_dir))
        elif args.dimension == '2d':
            plot_traction_solution_2d( os.path.normpath(args.input_files_path), os.path.normpath(args.output_dir) )
        else:
            print("Invalid dimension. Use '1d' or '2d'.")
            sys.exit(1)
