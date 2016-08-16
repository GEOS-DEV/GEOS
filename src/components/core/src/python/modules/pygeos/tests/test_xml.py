
import sys
sys.path.append('../../')
import pygeos
import os
import glob

# Find the targets
os.chdir('xml')
targets = glob.glob('*.xml')

# Parse the xml files
for target in targets:
  # try:
  target_parsed = pygeos.PreprocessGEOSXML(target)
  # except:
  #   print('Error parsing file: %s' % target)

# Move the results to a tmp directory 
os.system('mkdir -p results')
os.system('rm -f results/*')
os.system('mv *int.xml results/')
os.system('mv *prep.xml results/')