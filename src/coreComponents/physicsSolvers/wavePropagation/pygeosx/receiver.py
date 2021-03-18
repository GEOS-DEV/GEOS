'''
Created on 10/02/2021

@author: macpro
'''
import numpy as np

class Receiver:
    def __init__(self,pos):
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
        self.coords = np.array([self.x, self.y, self.z])
        
    def __repr__(self):
        return '('+str(self.x)+','+str(self.y)+','+str(self.z)+') \n'
    
    def getCoord(self):
        return self.coords 
    
    def setReceiverPos(self, coord, value):
        if coord==0:
            self.x = value
        elif coord==1:
            self.y = value
        elif coord==2:
            self.z = value
    
    #Transpose the receiver position    
    def linearMove(self, coord, value):
        if coord==0:
            self.x += value
        elif coord==1:
            self.y += value
        elif coord==2:
            self.z += value
            
    
class ReceiverSet:
        
        
    def __init__(self,receiver_list=None):
        if receiver_list is None:
            self.receiver_list = []
            self.n = 0
        else:
            self.receiver_list = receiver_list
            self.n = len(receiver_list)
            
    def __repr__(self): 
        return str(self.receiver_list)
            
    def getInsideDomain(self,meshBoundaries):
        list=[]
        
        for r in self.receiver_list:     
            if meshBoundaries[0][0]<=r.x and r.x<=meshBoundaries[0][1] \
            and meshBoundaries[1][0]<=r.y and r.y<=meshBoundaries[1][1] \
            and meshBoundaries[2][0]<=r.z and r.z<=meshBoundaries[2][1]:
                list.append(r)
        
        final_receiver_set = ReceiverSet(list)
        return final_receiver_set
    
    def getList(self):
        return self.receiver_list
    
    def getNumberOfReceivers(self):
        return self.n
    
    def getSetCoord(self):
        listRcv = []
        for i in range(self.n):
            listRcv.append(np.array([self.receiver_list[i].x, self.receiver_list[i].y, self.receiver_list[i].z]))
        return listRcv
        
        
    #Transpose all receivers position    
    def linearSetMove(self,coord,value):
        for i in range(self.n):
            if coord==0:
                self.receiver_list[i].x += value
            elif coord==1:
                self.receiver_list[i].y += value
            elif coord==2:
                self.receiver_list[i].z += value
                
        
        
