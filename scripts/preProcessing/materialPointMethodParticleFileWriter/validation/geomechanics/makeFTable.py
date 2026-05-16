import re 
import numpy as np
import matplotlib.pyplot as plt

#initialize plots
##plot strain versus time data from the experiment
timeData=[]
strainData=[]

with open ('/g/g19/homel1/pfwx/particleFileWriter/validation/geomechanics/strainvstime_IC1.csv') as csvfile:
	for line in csvfile:
		line = line.strip();
		values = line.split(',');
		
		#for some reason, there are other characters in the first column. Strip them.
		values[0]=re.sub(r'[^0-9.]','',values[0]);

		values[0]=float(values[0]);
		values[1]=float(values[1]);
		timeData.append(values[0]);
		strainData.append(values[1]);
  
timeDataNorm=[i/timeData[-1] for i in timeData];

stopTime = 2000
FTable = []
for i,t in enumerate(timeDataNorm):
    if (  i==0 or timeDataNorm[i] > timeDataNorm[ max(0,i-1) ] ):
        ev = strainData[i]
        J = np.exp(-ev)
        Fii = J**(1/3)
        FTable.append([t*stopTime,Fii,Fii,Fii])
        print([t*stopTime,Fii,Fii,Fii])

FTable = np.delete(FTable,list(range(25,150)),axis=0)
  
fig = plt.plot(FTable[:,1])
plt.savefig("temp.png", bbox_inches="tight")