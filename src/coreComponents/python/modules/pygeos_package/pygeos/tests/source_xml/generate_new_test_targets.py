#!/usr/gapps/GEOS/python/anaconda2/bin/python
import os
import pygeos


targets = ['basic.xml',
           'includes.xml',
           'symbolic.xml']


for t in targets:
  raw_file = pygeos.preprocessGEOSXML(t, verbose=0)
  pygeos.format_xml_file(raw_file, '../target_xml/target_%s' % (t))
  os.system('mv %s ../target_xml/raw_%s' % (raw_file, t))

