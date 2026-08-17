import math
import numpy as np

class DPAnalyticalSolution:
    def __init__(self):
        self.bulkMod = 0.5e9 # Pa
        self.shearMod = 0.3e9 # Pa
        self.lambdaCoeff = self.bulkMod - 2*self.shearMod/3. # Pa 
        self.cohesion = 5e3 # Pa
        self.frictionAngle = 15.27 # deg
        self.thermalExpansionCoeff = 3e-7
        self.phi = 6 * np.sin(np.deg2rad(self.frictionAngle)) / ( 3 - np.sin(np.deg2rad(self.frictionAngle)) )
        self.C = 6 * self.cohesion * np.cos(np.deg2rad(self.frictionAngle)) / ( 3 - np.sin(np.deg2rad(self.frictionAngle)) )
    
    def compute_stress(self, deltaTemp):
        sigma_x = 0
        sigma_z = 0

        epsMech_y = - self.thermalExpansionCoeff * deltaTemp

        # We first assume all elastic deformation
        epsMech_x = -epsMech_y * self.lambdaCoeff / 2 / (self.lambdaCoeff + self.shearMod)
        epsMech_z = epsMech_x

        sigma_y = (self.lambdaCoeff + 2*self.shearMod) * epsMech_y + 2 * self.lambdaCoeff * epsMech_x

        P = sigma_y / 3. 
        Q = sigma_y

        yieldFunc = Q + self.phi * P - self.C

        if yieldFunc <= 0.0:
            eps_x = epsMech_x + self.thermalExpansionCoeff * deltaTemp

            return sigma_x, sigma_y, eps_x
        else: 
            sigma_y = self.C / (1 + self.phi/3.0)
            epsMech_x = (sigma_y - (3*self.lambdaCoeff + 2*self.shearMod)*epsMech_y )/ (6*self.lambdaCoeff + 4*self.shearMod)
            multipler = -(sigma_y - (self.lambdaCoeff + 2*self.shearMod)*epsMech_y - 2*self.lambdaCoeff*epsMech_x)/ (4/3 * self.shearMod) 
            eps_x = epsMech_x + self.thermalExpansionCoeff * deltaTemp

            return sigma_x, sigma_y, eps_x
        
    def compute_disp(self, deltaTemp):
        _, _, esp_x = self.compute_stress(deltaTemp)

        disp_y = 0
        disp_x = esp_x*1

        return disp_x, disp_y