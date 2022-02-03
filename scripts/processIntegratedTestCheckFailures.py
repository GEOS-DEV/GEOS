#!/usr/bin/env python3
import sys
import os

def findFiles(folder):
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if( '.data' in filename ):
                yield os.path.join(root, filename)



matchStrings = [ 'Error:' ]
#exclusionStrings = [ 'sizedFromParent', 'different shapes' ]
exclusionStrings = [ 'sizedFromParent', 'different shapes', 'but not the' ]

numTrailingLines = 5


for fileName in findFiles(sys.argv[1]):
    #fileName = 'integratedTests/compositionalMultiphaseFlow/4comp_2ph_1d_01/4comp_2ph_1d_01.data'

    filteredErrors=''

    with open(fileName) as f:
        lines = f.readlines()
        
        for i in range(0,len(lines)):
            line = lines[i]
            if all(matchString in line for matchString in matchStrings):
                matchBlock = lines[i-1]
                matchBlock += line

                
                for j in range(1,numTrailingLines+1):
                    if not ('0: ********************************************************************************' in lines[i+j]):
                        matchBlock += lines[i+j]
                    else:
                        break
                i += j
#                print( j )


                if not any( excludeString in matchBlock for excludeString in exclusionStrings):
                    filteredErrors += matchBlock

        
    if( len( filteredErrors ) ):
        print( "IN ", fileName )
        print( filteredErrors, flush=True )

#for i in range(1,1+1):
#    print( i )
