def SaveTimeHistory( name, array ):
    with open(name, 'w') as filehandle:
        filehandle.write("       time   lengthMax   cellNumber\n")
        for row in array:
            for entry in row:
                filehandle.write("%10.4g" % entry)
            filehandle.write("\n")


def GetPackLengthArea(ghostRanks):
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
          packVf=float(v['Fracture_ElementFields/proppantPackVolumeFraction'])
          if(packVf >= 1.0) : 
            cellNum += 1
            if(zc < 0.01) :
              if(yc > lengthMax) : lengthMax = yc
      TurnDomainsOff(d)
      dnum += 1
    return lengthMax, cellNum


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


if len(sys.argv) < 4:
    sys.exit('Usage: %s path/To/database database startTime dCycleNumber outputFile' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: path %s was not found!' % sys.argv[1])

database = "localhost:" + sys.argv[1] + "/" + sys.argv[2]
startTime = float(sys.argv[3])
dCycleNumber = int(sys.argv[4])
outputFile = sys.argv[5]

print "using " + database

OpenDatabase(database)

ghostRanks = []
SetTimeSliderState(0)
GetGhostRank(ghostRanks)

AddPlot("Pseudocolor", "Fracture_ElementFields/proppantPackVolumeFraction")
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
            lengthMax, cellNum = GetPackLengthArea(ghostRanks)
            timehist.append([float(time), float(lengthMax), float(cellNum)])

SaveTimeHistory(outputFile, timehist)

quit()    
