def SaveTimeHistory( name, array ):
    with open(name, 'w') as filehandle:
        filehandle.write("       time   cellNumber\n")
        for row in array:
            for entry in row:
                filehandle.write("%10.4g" % entry)
            filehandle.write("\n")


def GetSuspendedProppantCells(ghostRanks, offsetDistance, proppantConcThresh, packVfs):
    SetActivePlots(2)    
    doms = GetDomains()
    for d in doms:
      TurnDomainsOff(d)

    dnum=1
    lengthMax=0.0
    cellNum=0
    for d in doms:
      TurnDomainsOn(d) 
      num = len(ghostRanks[dnum-1])
      for dd in range(0, num):
        v=PickByZone(domain=dnum, element=dd)
        if(int(ghostRanks[dnum-1][dd]) < 0) :
          coord=v['point']
          yc = float(coord[1])      
          zc = float(coord[2])
          proppantConcentration=float(v['Fracture_ElementFields/proppantConcentration'])
          if(yc > offsetDistance and proppantConcentration > proppantConcThresh and packVfs[dnum-1][dd] < 1.0) :
            cellNum += 1
      TurnDomainsOff(d)
      dnum += 1
    return cellNum


def GetGhostRank(ghostRanks):
    AddPlot("Pseudocolor", "Fracture_ElementFields/ghostRank")
    DrawPlots()
    SuppressMessages(2)    
    doms = GetDomains()
    for d in doms:
      TurnDomainsOff(d)

    dnum=1
    for d in doms:
      TurnDomainsOn(d) 
      Query("NumZones")
      a = GetQueryOutputValue()
      lnum = int(a[0])
      lgnum = int(a[1])  
      num = lnum + lgnum
      b=[]
      for dd in range(0, num):
        v=PickByZone(domain=dnum, element=dd)
        u=v['Fracture_ElementFields/ghostRank']
        b.append(u)
      ghostRanks.append(b)  
      TurnDomainsOff(d)
      dnum += 1


def GetPackVf(packVfs):
    SetActivePlots(1)
    doms = GetDomains()
    for d in doms:
      TurnDomainsOff(d)

    dnum=1
    for d in doms:
      TurnDomainsOn(d) 
      Query("NumZones")
      a = GetQueryOutputValue()
      lnum = int(a[0])
      lgnum = int(a[1])  
      num = lnum + lgnum
      b=[]
      for dd in range(0, num):
        v=PickByZone(domain=dnum, element=dd)
        u=float(v['Fracture_ElementFields/proppantPackVolumeFraction'])
        b.append(u)
      packVfs.append(b)  
      TurnDomainsOff(d)
      dnum += 1

      

if len(sys.argv) < 8:
    sys.exit('Usage: %s path/To/database database startTime dCycleNumber offsetDistance proppantConcThresh outputFile' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: path %s was not found!' % sys.argv[1])

database = "localhost:" + sys.argv[1] + "/" + sys.argv[2]
startTime = float(sys.argv[3])
dCycleNumber = int(sys.argv[4])
offsetDistance = float(sys.argv[5])
proppantConcThresh = float(sys.argv[6])
outputFile = sys.argv[7]

print "using " + database

OpenDatabase(database)

ghostRanks = []
SetTimeSliderState(0)
GetGhostRank(ghostRanks)

AddPlot("Pseudocolor", "Fracture_ElementFields/proppantPackVolumeFraction")
AddPlot("Pseudocolor", "Fracture_ElementFields/proppantConcentration")
DrawPlots()
SuppressMessages(2)

timehist = []

for state in range(TimeSliderGetNStates()):
    SetTimeSliderState(state)
    SetQueryFloatFormat("%g")
    time = float(Query("Time")[:-1].split(' ')[-1])
    if time == startTime:
        start = state
    if time >= startTime:    
        if (state - start) % dCycleNumber == 0:
            print "processing state, time ", state, time
            packVfs = []
            GetPackVf(packVfs)            
            cellNum = GetSuspendedProppantCells(ghostRanks, offsetDistance, proppantConcThresh, packVfs)
            timehist.append([float(time), float(cellNum)])
            
SaveTimeHistory(outputFile, timehist)

quit()    
