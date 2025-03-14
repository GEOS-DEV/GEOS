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
        self.Drs = 5.0e-4
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
        inner_expression = self.Vinit / (2* self.rsParameters.Vstar) * math.exp( ( self.rsParameters.f + self.rsParameters.b * math.log( self.rsParameters.Vstar / self.Vinit ) ) / self.rsParameters.a( r ) )
        t0 = self.normalTraction * self.rsParameters.a(r) *  math.asinh( inner_expression ) + self.rsParameters.shearImpedance() * self.Vinit
        tx = t0 * self.Vinit / np.linalg.norm([self.Vinit, self.Vzero])
        ty = t0 * self.Vzero / np.linalg.norm([self.Vinit, self.Vzero])
        return [tx, ty]

    def nucleation_forcing( self, r, t ):
        return self.deltaTau * self.G1( r ) * self.G2( t )

    def G1(self, r):
        if ( r > self.Rn ):
            return 0.
        else:
            return math.exp( r**2 / (r**2 - self.Rn**2) )

    def G2(self, t):
        if (t > self.T):
            return 1.0
        else:
            return math.exp((t-self.T)**2 / (t*(t-2*self.T)))

    def initialState( self ):
        theta = self.rsParameters.Drs / self.Vinit
        print( f'theta = {theta}' )
        Psi = self.rsParameters.f + self.rsParameters.b * math.log( self.rsParameters.Vstar * theta / self.rsParameters.Drs );
        return Psi        

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

def calcualteBP7Parameters( dir, printTables ):
    bp7 = BP7()
    if (printTables):
        xcoords = np.linspace(bp7.x[0], bp7.x[1], 100).tolist()
        ycoords = np.linspace(bp7.y[0], bp7.y[1], 100).tolist()
        tcoords = np.linspace(0.1, 1.0, 1000).tolist()
        tcoords.append( 3.156e+8 )
        # Prepare lists to store data
        data = { 'z': [0.0], 'x': xcoords, 'y': ycoords, 'backgroungShearTraction_x': [], 'backgroungShearTraction_y': [], 'a': [] }
        for x, y in product(xcoords, ycoords):
            r = np.sqrt(x**2 + y**2)
            tau0 = bp7.reference_shearTraction( r )
            data['backgroungShearTraction_x'].append( tau0[0] )
            data['backgroungShearTraction_y'].append( tau0[1] )
            data['a'].append( bp7.rsParameters.a(r) )
    
        writeBP7Tables( data, dir )  
        forcing_data = { 'z': [0.0], 'x': xcoords, 'y': ycoords, 't': tcoords, 'backgroungShearTractionWithForcing_x': [] }

        xnc = -50.0
        ync = -50.0
        
        for t in tcoords:
            for x, y in product(xcoords, ycoords):
                rn = np.sqrt((x - xnc)**2 + (y - ync)**2)
                r = np.sqrt(x**2 + y**2)
                tau0 = bp7.reference_shearTraction( r )
                forcing_data['backgroungShearTractionWithForcing_x'].append( tau0[0] + bp7.nucleation_forcing( rn, t ) )

        writeBP7Tables( forcing_data, dir ) 

    print( f'initial state variable = {bp7.initialState()}' )
    print( f'shearImpedance = {bp7.rsParameters.shearImpedance() }' )
   
     
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-pt', '--print-tables', type=bool, help='whether to print tables or not', default=False)
    parser.add_argument('-d', '--dir', type=str, help='output dir', default='.')
    args = parser.parse_args()
    normalized_dir = os.path.abspath( args.dir ) 
    calcualteBP7Parameters( normalized_dir, args.print_tables )


    
  