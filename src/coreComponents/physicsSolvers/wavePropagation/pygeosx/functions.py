'''
Created on 8/02/2021

@author: macpro
'''

import numpy as np
import math as m
    

def ricker(maxT, dt, f0):
    
    T0 = 1.0/f0;
    fi = np.zeros(int(maxT/dt))
    
    for t in range(int(maxT/dt)):
        t0 = dt*t
        if t0 <= -0.9*T0 or t0 >= 2.9*T0:
            fi[t] = 0.0;
        else:
            tmp      = f0*t0-1.0
            f0tm1_2  = 2*(tmp*np.pi)*(tmp*np.pi)
            gaussian = m.exp( -f0tm1_2 )
            fi[t]    = -(t0-1)*gaussian
               
    return fi



#Calculate dt using order of space discretization method, Wave velocity, and the radius of the included sphere in element 
def calculDt(mesh):
   
    if mesh.getOrd()==1:
        nx = 1
        ny = 2
        nz = 4
    elif mesh.getOrd()==3:
        nx = 3
        ny = 12
        nz = 48
    elif mesh.getOrd()==5:
        nx = 5
        ny = 30
        nz = 180
    

    h    = np.linalg.norm(mesh.getElem_List()[0].getNode_List()[0].getCoords() - mesh.getElem_List()[0].getNode_List()[nx].getCoords())/2
    Vmax = mesh.getElem_List()[0].getSpeed()
    V    = mesh.getElem_List()[0].getVolume()
    
    for elem in (mesh.getElem_List()):
        
        if elem.getVolume()<V: 
        
            xradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nx].getCoords())/2
            yradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[ny].getCoords())/2
            zradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nz].getCoords())/2
            
            if xradius < h:
                h = xradius
            if yradius < h:
                h = yradius
            if zradius < h:
                h = zradius
            
        if Vmax < elem.getSpeed():
            Vmax = elem.getSpeed()
            
        
    dt = h/(Vmax*mesh.ord)           
            
    return dt
    

            
