
import argparse
import numpy as np
import xml.etree.ElementTree as ElementTree
import matplotlib.pyplot as plt

class analyticalSolution():
    def __init__( self, parameters ):
        self.cTau = 1e-3
        self.aSigma = 1e6
        self.t_a = 1e6/parameters['backgroundStressingRate']
        self.initial_normal_traction     = parameters['initialFaultNormalTraction']
        self.initial_fault_shear_traction = parameters['initialFaultShearTraction']

    def compute_shear_stress( self, time_n, dt ):
        return self.aSigma * np.log( self.cTau*( time_n + dt ) + 1 ) + self.initial_fault_shear_traction

    def compute_seismic_rate(self, time_n, dt, shear_stress):
        K = np.exp( shear_stress / self.aSigma - shear_stress / self.aSigma )
        denom = ( self.cTau*(time_n + dt - self.t_a ) + 1) * np.exp( (time_n+dt)  / self.t_a ) + self.cTau*self.t_a
        return K/denom

def getParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse( xmlFilePath )

    solver = tree.find('Solvers/DieterichSeismicityRate')
    
    parameters = dict.fromkeys(["initialFaultNormalTraction", "initialFaultShearTraction", "backgroundStressingRate"])
    parameters["initialFaultNormalTraction"] = float(solver.get("initialFaultNormalTraction"))
    parameters["initialFaultShearTraction"]  = float(solver.get("initialFaultShearTraction"))
    parameters["backgroundStressingRate"]  = float(solver.get("backgroundStressingRate"))

    return parameters

def curve_check_solution(**kwargs):
     
    xmlFilePath = "./SeismicityRate_verificationProblem_smoke.xml"
    parameters = getParametersFromXML( xmlFilePath )

    analytical_solution = analyticalSolution( parameters )

    time = 0.0
    dt = 3600.0
    final_time = 360000.0
    times = []
    analytical_rates = []
    while time < final_time:
        tau = analytical_solution.compute_shear_stress(time, dt)
        analytical_rate = analytical_solution.compute_seismic_rate(time, dt, tau)
        time += dt
        times.append(time)
        analytical_rates.append(analytical_rate)

    return analytical_rates    


def debug( xmlFilePath ):
    #-------- Extract info from XML
    parameters = getParametersFromXML( xmlFilePath )

    analytical_solution = analyticalSolution( parameters )

    time = 0.0
    dt = 3600.0
    final_time = 360000.0
    times = []
    analytical_rates = []
    while time < final_time:
        tau = analytical_solution.compute_shear_stress(time, dt)
        analytical_rate = analytical_solution.compute_seismic_rate(time, dt, tau)
        time += dt
        times.append(time)
        analytical_rates.append(analytical_rate)
        

    # Plot analytical (continuous line) and numerical (markers) aperture solution 
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.plot(times, analytical_rates)
    ax.set_xlabel('time [s]', weight="bold")
    ax.set_ylabel('seismic rate', weight="bold")
    plt.savefig("seismicRate.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--xml-file', type=str, help='Path to XML file')
    args = parser.parse_args()
    debug( args.xml_file )
