import argparse
import numpy as np
import csv
import matplotlib.pyplot as plt
import os
import math 
from itertools import product

class RateAndStateParameters:
    def __init__(self):
        self.amax = 0.016
        self.a0 = 0.004
        self.b = 0.01
        self.Vstar = 1.0e-6
        self.f = 0.6
        self.Rvw = 200
        self.Drs = 5.0e-3
        self.shearModulus = 32.04e9 
        self.cs = 3.464e3

    def a(self, r):
        if (r > self.Rvw ):
            return self.amax
        else:
            return self.a0    

    def shearImpedance(self):
        return self.shearModulus / (2*self.cs)            

class BP7:
    def __init__(self):
        self.x = [-500, 500]
        self.y = [-500, 500]
        self.z = [-500, 500]
        self.Rn = 150
        self.Vinit = 1.0e-9
        self.Vzero = 1.0e-20
        self.normalTraction = 25e6 
        self.T = 1
        self.deltaTau = 1.75e6
        self.rsParameters = RateAndStateParameters()

    def reference_shearTraction( self, r ): 
        inner_expression = self.Vinit / (2* self.rsParameters.Vstar) * math.exp( self.rsParameters.f + self.rsParameters.b * math.log( self.rsParameters.Vstar / self.Vinit ) / self.rsParameters.a( r ) )
        t0 = self.normalTraction * self.rsParameters.a(r) *  math.asinh( inner_expression )
        ty = t0 * self.Vinit / np.linalg.norm([self.Vinit, self.Vzero])
        tz = t0 * self.Vzero / np.linalg.norm([self.Vinit, self.Vzero])
        return [ty, tz]

    def nucleation_forcing( self, r, t ):
        return self.deltaTau * self.G1( r ) * self.G2( t )

    def G1(self, r):
        if ( r > self.Rn ):
            return 0
        else:
            return math.exp( r**2 / (r**2 - self.Rn**2) )

    def G2(self, t):
        if (t > self.T):
            return 1.0
        else:
            return math.exp((t-self.T)**2 / (t*(t-2*self.T)))        

# Function to save a numpy array to a .geos file
def save_to_file(filename, data):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        for value in data:
            writer.writerow([value])

def writeBP7Tables( data_to_save, dir ):    
    # Loop through the dictionary and save each file
    for filename, data in data_to_save.items():
        save_to_file(f'{dir}/{filename}.geos', data)

def calcualteBP7Parameters( dir ):
    bp7 = BP7()
    ycoords = np.linspace(bp7.y[0], bp7.y[1], 100)
    zcoords = np.linspace(bp7.z[0], bp7.z[1], 100)

    # Prepare lists to store data
    data = { 'x': [-1000, 1000], 'y': [], 'z': [], 'backgroungShearTraction_y': [], 'backgroungShearTraction_z': [] }
    for y, z in product(ycoords, zcoords):
        data['y'].append(y)
        data['z'].append(z)
        r = np.sqrt(y**2 + z**2)
        tau0 = bp7.reference_shearTraction( r )
        data['backgroungShearTraction_y'].append( tau0[0] )
        data['backgroungShearTraction_z'].append( tau0[1] )
    
    writeBP7Tables( data, dir )  

    forcing_data = { 'x': [-1000, 1000], 'forcing_y': [], 'forcing_z': [], 'forcing_t': [], 'forcing_value': [] }
    for t in np.linspace(0.1, 1.0, 10):
      for y, z in product(ycoords, zcoords):
        forcing_data['forcing_y'].append(y)
        forcing_data['forcing_z'].append(z)
        forcing_data['forcing_t'].append(t)
        r = np.sqrt(y**2 + z**2)
        forcing_data['forcing_value'].append( bp7.nucleation_forcing( r, t ) )

    writeBP7Tables( forcing_data, dir )  
   
     
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, help='output dir', default='.')
    args = parser.parse_args()
    normalized_dir = os.path.abspath( args.dir ) 
    calcualteBP7Parameters( normalized_dir )


    
  