#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
from gdal2mbtiles import MBTilesBuilder
# .bashrc
# export PYTHONPATH=/usr/lib/mapmbtiles:$PYTHONPATH

logging.basicConfig(level=logging.DEBUG)
input_directory="source"
output_directory="output"
s_name="";
# zoom-levels 0-3
Source_filepath="%s/lidarplots.mbtiles" % input_directory
for i_loop in range(0,1):
 if i_loop == 0:
  # this db was to big to include in GitHub
  # remote: warning: File samples/source/lidarplots.mbtiles is 90.62 MB; this is larger than GitHub's recommended maximum file size of 50 MB
  # wget http://www.mj10777.de/public/download/mbtiles/test_mbtiles/lidarplots.mbtiles -O source/lidarplots.mbtiles
  Output_filepath="%s/lidarplots.mbtiles" % output_directory
 print  "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
 mb = MBTilesBuilder(mbtiles_input=Source_filepath, mbtiles_output=Output_filepath)
 # lidarplots.mbtiles db used 'tiles' as a table and not a view
 # - the result will be created with 'tiles' as a view with 'map' and 'images' as a table
 # the metadata table used mixed lettering for 'name' : 'minZoom'
 # - these will be recognized and changed to lower case
 # all of the 'metadata' of the 'input' will be placed in the 'output' mbtiles.db
 # for i_loop == 0 : the file will be rewritten and all metadata saved
 mb.run()



