#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging,os
from gdal2mbtiles import MBTilesBuilder
# .bashrc
# export PYTHONPATH=/usr/lib/mapmbtiles:$PYTHONPATH

logging.basicConfig(level=logging.DEBUG)
input_directory="source"
output_directory="output"
# i_parm=0: run programm ; = 1 show info only and return
i_parm=0
#----------------------------------------------------------------------------------
# Step 1 - we will subdivide an existing mbtiles
#----------------------------------------------------------------------------------
i_min_level=2
i_max_level=3
# zoom-levels 0-3
Source_filepath="%s/1861_Mercator_World.mbtiles" % input_directory
if os.path.exists(Source_filepath):
 for i_loop in range(i_min_level,i_max_level+1):
  if i_loop == 0:
   Output_filepath="%s/1861_Mercator_World.mbtiles" % output_directory
  if i_loop == 1:
   Output_filepath="%s/1861_Mercator_Europe.mbtiles" % output_directory
  if i_loop == 2:
   Output_filepath="%s/1861_Mercator_Europe_Africa.mbtiles" % output_directory
  print  "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
  mb_step1 = MBTilesBuilder(mbtiles_input=Source_filepath, mbtiles_output=Output_filepath)
  # all of the 'metadata' of the 'input' will be placed in the 'output' mbtiles.db
  # for i_loop == 0 : the file will be rewritten and all metadata saved
  # for i_loop == 1 : the area of Europe will be extracted
  # for i_loop == 2 : the area of Europe and Africa will be extracted
  # when finished with the inserts, the  the min/max zoom_levels, bounds and center will be calculated and saved
  if i_loop == 1:
   position_west=-8.0
   position_south=36.0
   position_east=80.0
   position_north=77.0
   # Add area and do zoom-level 3 only
   mb_step1.add_coverage(bbox=(position_west,position_south,position_east,position_north), zoomlevels=[3])
   # - the original 'name' and calulated center will be changed ; here 'Europe' and the Brandenburg Gate, Berlin as center
   # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
   metadata_list = [
          (u'name',u'1861 Mercator Europe'),
          (u'center',u'13.3777065575123,52.5162690144797,3'),
         ]
   mb_step1.add_metadata(metadata_list)
  if i_loop == 2:
   position_west=-8.0
   position_south=36.0
   position_east=80.0
   position_north=77.0
   # Add area and do zoom-levels 2 to 3 (inclusive)
   mb_step1.add_coverage(bbox=(position_west,position_south,position_east,position_north),zoomlevels=[2,3])
   position_west=-32.0
   position_south=-37.0
   position_east=53.0
   position_north=36.0
   # Add area and do zoom-levels 2 to 3 (inclusive)
   mb_step1.add_coverage(bbox=(position_west,position_south,position_east,position_north),zoomlevels=[2,3])
   # - the original 'name' and calulated center will be changed ; here 'Europe/Africa' and Piazza della Signoria, Florence as center
   # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
   metadata_list = [
          (u'name',u'1861 Mercator Europe/Africa'),
          (u'center',u'11.255867505557262,43.76961577178412,3'),
         ]
   mb_step1.add_metadata(metadata_list)
  # with 'False' we are telling MBTilesBuilder to exit if the file exists
  # with 'True' we are telling MBTilesBuilder to use an existing mbtiles file
  # i_parm=0: run programm ; = 1 show info only and return
  # - if no parmameters are given: 'False,0' is assumed
  # i_parm=1
  mb_step1.run(False,i_parm)
  del mb_step1
#----------------------------------------------------------------------------------
# Step 2 - we will add to an existing mbtiles a portion from another mbtiles file
#----------------------------------------------------------------------------------
# zoom-levels 0-8 [this file is the complete collection - it includes the same tiles as 1861_Mercator_World.mbtiles]
# 'wget http://www.mj10777.de/public/download/mbtiles/1861_Mercator_Countries.mbtiles -O source/1861_Mercator_Countries.mbtiles
Source_filepath="%s/1861_Mercator_Countries.mbtiles" % input_directory
if os.path.exists(Source_filepath):
 # we will be adding an area from another mbtiles file to the existing one
 Output_filepath="%s/1861_Mercator_Europe.mbtiles" % output_directory
 if os.path.exists(Output_filepath):
  print  "Step 2\nSource: ",Source_filepath," \nOutput: ",Output_filepath
  mb_step2 = MBTilesBuilder(mbtiles_input=Source_filepath,mbtiles_output=Output_filepath)
  # 12.656250,51.618017,14.062500,53.330873
  position_west=12.656250
  position_south=51.618017
  position_east=14.062500
  position_north=53.330873
  # Add area and do zoom-levels 4 to 8 (inclusive)
  mb_step2.add_coverage(bbox=(position_west,position_south,position_east,position_north), zoomlevels=[4,5,6,7,8])
  # with 'False' we are telling MBTilesBuilder to exit if the file exists
  # with 'True' we are telling MBTilesBuilder to use an existing mbtiles file
  # i_parm=0: run programm ; = 1 show info only and return
  # - if no parmameters are given: 'False,0' is assumed
  # i_parm=1
  mb_step2.run(True,i_parm)
  del mb_step2
#----------------------------------------------------------------------------------
# Step 3 - checking for logical errors
# - setting i_parm=1 will show some information about what is going to be done
# -- amount of tiles to be created
# -- if the output file will be overridden
# -- if given bounds are correct
#----------------------------------------------------------------------------------
# zoom-levels 0-8 [this file is the complete collection - it includes the same tiles as 1861_Mercator_World.mbtiles]
# 'wget http://www.mj10777.de/public/download/mbtiles/1861_Mercator_Countries.mbtiles -O source/1861_Mercator_Countries.mbtiles
Source_filepath="%s/1861_Mercator_Countries.mbtiles" % input_directory
if os.path.exists(Source_filepath):
 # we will be adding an area from another mbtiles file to the existing one
 Output_filepath="%s/1861_Mercator_Europe.mbtiles" % output_directory
 if os.path.exists(Output_filepath):
  print  "Step 3 - checking for logical errors\nSource: ",Source_filepath," \nOutput: ",Output_filepath
  mb_step3 = MBTilesBuilder(mbtiles_input=Source_filepath,mbtiles_output=Output_filepath)
  # 12.656250,51.618017,14.062500,53.330873
  position_west=12.656250
  position_south=51.618017
  position_east=14.062500
  position_north=-53.330873
  # Add area and do zoom-levels 4 to 8 (inclusive)
  mb_step3.add_coverage(bbox=(position_west,position_south,position_east,position_north), zoomlevels=[4,5,6,7,8])
  # with 'False' we are telling MBTilesBuilder to exit if the file exists
  # with 'True' we are telling MBTilesBuilder to use an existing mbtiles file
  # i_parm=0: run programm ; = 1 show info only and return
  # - if no parmameters are given: 'False,0' is assumed
  i_parm=1
  # mb_step3.run(True,i_parm)
  del mb_step3
#----------------------------------------------------------------------------------

