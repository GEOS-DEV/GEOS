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

numTrailingLines = 1


for fileName in findFiles(sys.argv[1]):#'integratedTests'):
    #fileName = 'integratedTests/compositionalMultiphaseFlow/4comp_2ph_1d_01/4comp_2ph_1d_01.data'
    print( 'Processing ',fileName )

    filteredErrors=''

    with open(fileName) as f:
        lines = f.readlines()
        
        for i in range(0,len(lines)-numTrailingLines):
            line = lines[i]
            if all(matchString in line for matchString in matchStrings):
                matchBlock = line

                
                for j in range(1,numTrailingLines+1):
#                    print( len(lines),i,j, line)
                    matchBlock += lines[i+j]


                if not any( excludeString in matchBlock for excludeString in exclusionStrings):
                    filteredErrors += matchBlock

        
    if( len( filteredErrors ) ):
        print( filteredErrors, flush=True )

#for i in range(1,1+1):
#    print( i )
