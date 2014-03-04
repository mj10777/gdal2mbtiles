#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging,os
from gdal2mbtiles import GDAL2MbTiles

logging.basicConfig(level=logging.DEBUG)
input_directory="source"
output_directory="output"
# Image is in WGS 84 / Pseudo-Mercator ; epsg:3857
# Europe
position_west= -939258.204 # -8.0
position_south=4226661.916 # 36.0
position_east=8922952.934 # 80.0
position_north=13932330.020 # 77.0
te_bounds="%f,%f,%f,%f"% (position_west,position_south,position_east,position_north)
s_name="";
Source_filepath="%s/1861_Mercator_Europe.tif" % output_directory
i_min_level=1
i_max_level=1
# zoom-levels 0-3
for i_loop in range(i_min_level,i_max_level+1):
 if i_loop == 0:
  if os.path.exists(Source_filepath):
   # this should tile the whole area of the Image
   Output_filepath="%s/1861_Mercator_Europe.mbtiles" % output_directory
   s_title=u'1861 Mercator Europe from Countries map'
   s_copyright='Johnston, Alexander Keith, 1804-1871'
   min_zoom=3
   max_zoom=4
   gdal2mbtiles_parms = ['--profile','mercator',
    '--tile-format','jpeg'
    ]
   gdal2mbtiles_parms.append('--mbtiles')
   # gdal2mbtiles_parms.append('--verbose')
   # gdal2mbtiles_parms.append('--resume')
   gdal2mbtiles_parms.extend(['--zoom',"%i-%i" % (min_zoom,max_zoom)])
   gdal2mbtiles_parms.extend(['--title',s_title])
   gdal2mbtiles_parms.extend(['--copyright',s_copyright])
   # gdal2mbtiles.extend(['--s_srs',"%i" % (3068)])
   # after the other parameters: the input file
   gdal2mbtiles_parms.append(Source_filepath.encode('utf-8'))
   # last parameter the output file
   # there is no output parameter
   print  "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
   gdal_tiles = GDAL2MbTiles(gdal2mbtiles_parms)
   gdal_tiles.process()
 if i_loop == 1:
  # zoom-levels 0-8 [this file is the complete collection - it includes the same tiles as 1861_Mercator_World.mbtiles]
  # 'wget http://www.mj10777.de/public/download/mbtiles/1861_Mercator_Countries.mbtiles -O source/1861_Mercator_Countries.mbtiles
  if os.path.exists(Source_filepath):
   output_directory="blombo"
   Output_filepath="%s/1861_Mercator_Berlin.mbtiles" % output_directory
   # this should tile only a portion of the Image
   # Berlin Area
   position_west=1252344.271 # 12.0
   position_south=6731350.459 # 52.2969487
   position_east=1565430.339 # 14.0
   position_north=7044436.527 # 52.7002333
   te_bounds="%f,%f,%f,%f"% (position_west,position_south,position_east,position_north)
   s_title=u'1861 Mercator Germany from Countries map'
   s_copyright='Johnston, Alexander Keith, 1804-1871'
   min_zoom=4
   max_zoom=8
   gdal2mbtiles_parms = ['--profile','mercator',
    '--tile-format','jpeg'
    ]
   gdal2mbtiles_parms.append('--mbtiles')
   gdal2mbtiles_parms.append('--verbose')
   # gdal2mbtiles_parms.append('--resume')
   gdal2mbtiles_parms.extend(['--zoom',"%i-%i" % (min_zoom,max_zoom)])
   gdal2mbtiles_parms.extend(['--title',s_title])
   gdal2mbtiles_parms.extend(['--copyright',s_copyright])
   gdal2mbtiles_parms.extend(['--te',te_bounds])
   # gdal2mbtiles.extend(['--s_srs',"%i" % (3068)])
   # after the other parameters: the input file
   gdal2mbtiles_parms.append(Source_filepath.encode('utf-8'))
   # last parameter the output file
   gdal2mbtiles_parms.append(Output_filepath.encode('utf-8'))
   print  "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
   gdal_tiles = GDAL2MbTiles(gdal2mbtiles_parms)
   gdal_tiles.process()



