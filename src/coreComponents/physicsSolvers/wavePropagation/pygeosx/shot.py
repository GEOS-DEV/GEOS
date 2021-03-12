'''
Created on 10/02/2021

@author: macpro
'''
from source import Source
from receiver import ReceiverSet

'''Create Shot object composed of one Source and a ReceiverSet'''
class Shot:
    def __init__(self,Source, ReceiverSet):
        self.source = Source
        self.receivers = ReceiverSet
        self.flag = "Undone"
        
    def __repr__(self):
        return 'Source position : \n'+str(self.source) +' \n\n' + 'Receivers positions : \n' + str(self.receivers) + '\n\n'
    
    def getSource(self):
        return self.source
    
    def getReceiverSet(self):
        return self.receivers
        
    def getFlag(self):
        return self.flag
        
    def flagUpdate(self, string):
        self.flag = string
        
        

    
