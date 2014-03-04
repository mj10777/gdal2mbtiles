#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging,os
from gdal2mbtiles import ImageExporter
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
i_min_level=0
i_max_level=2
#----------------------------------------------------------------------------------
# Compression settings - only used when tif files are used, otherwise ignored
# - for tif files a GTif (Geofeferenced Tif) will be created
# -- will only be used if 'Output_filepath' ends with ".tif" or ".tiff"
#----------------------------------------------------------------------------------
# possible values: 'NONE' ; 'PACKBITS' ; 'DEFLATE' ; 'LZW'
s_tiff_compression="LZW"
# possible values: 1 or 2 for 'DEFLATE' or 'LZW', otherwise ignored
i_tiff_predictor=2
# possible values: 1 to 9 for 'DEFLATE', otherwise ignored
i_tiff_zlevel=9
#----------------------------------------------------------------------------------
# Compression settings - only used when jpg files are used, otherwise ignored
# -- will only be used if 'Output_filepath' ends with ".jpg" or ".jpeg"
#----------------------------------------------------------------------------------
# possible values: 1 to 95 [default=75]
# 1861_Mercator_Europe.jpg was not possible:
# IOError: encoder error -2 when writing image file
i_jpg_quality=75
#----------------------------------------------------------------------------------
# zoom-levels 0-3
Source_filepath="%s/1861_Mercator_World.mbtiles" % input_directory
if os.path.exists(Source_filepath):
 for i_loop in range(i_min_level,i_max_level+1):
#----------------------------------------------------------------------------------
  if i_loop == 0:
   Output_filepath="%s/1861_Mercator_World.tif" % output_directory
   print "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
   ie_step_1 = ImageExporter(mbtiles_input=Source_filepath,jpg_quality=i_jpg_quality,tiff_compression=s_tiff_compression,tiff_predictor=i_tiff_predictor,tiff_zlevel=i_tiff_zlevel)
   # 12.656250,51.618017,14.062500,53.330873
   position_west=-8.0
   position_south=36.0
   position_east=80.0
   position_north=77.0
   # mbtiles_metadata 'name'        will be set as 'TIFFTAG_DOCUMENTNAME'
   # mbtiles_metadata 'description' will be set as 'TIFFTAG_IMAGEDESCRIPTION'
   # - this values can be overridden here - otherwise the values of 'Source_filepath' will be used
   # if 'TIFFTAG_SOFTWARE' is not set: the bounds in wsg84 will be used
   # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
   metadata_list = [
   #        (u'name',''),
   #        (u'TIFFTAG_DOCUMENTNAME',''),
   #        (u'description',''),
   #        (u'TIFFTAG_IMAGEDESCRIPTION',''),           
   #        (u'TIFFTAG_SOFTWARE','mapmbtiles.ImageExporter'),
           (u'TIFFTAG_DATETIME','1861'),
           (u'TIFFTAG_ARTIST','Johnston, Alexander Keith, 1804-1871'),
           (u'TIFFTAG_HOSTCOMPUTER','http://www.davidrumsey.com/luna/servlet/detail/RUMSEY~8~1~21292~610098:Chart-of-the-world-on-Mercators-pro?sort=Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No&qvq=w4s:/what/Atlas%20Map/World%20Atlas/where/World;sort:Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No;lc:RUMSEY~8~1&mi=42&trs=583'),
           (u'TIFFTAG_COPYRIGHT','David Rumsey Map Collection'),
          ]
   ie_step_1.add_metadata(metadata_list)
   ie_step_1.export_image(bbox=(position_west,position_south,position_east,position_north),zoomlevel=3,imagepath=Output_filepath)
   del ie_step_1
#----------------------------------------------------------------------------------
  if i_loop == 1:
   Source_filepath="%s/1861_Mercator_Countries.mbtiles" % input_directory
   if os.path.exists(Source_filepath):
    Output_filepath="%s/1861_Mercator_Europe.tif" % output_directory
    print "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
    ie_step_2 = ImageExporter(mbtiles_input=Source_filepath,jpg_quality=i_jpg_quality,tiff_compression=s_tiff_compression,tiff_predictor=i_tiff_predictor,tiff_zlevel=i_tiff_zlevel)
    # -8.0,36.0,80.0,77.0
    position_west=-8.0
    position_south=36.0
    position_east=80.0
    position_north=77.0
    # mbtiles_metadata 'name'        will be set as 'TIFFTAG_DOCUMENTNAME'
    # mbtiles_metadata 'description' will be set as 'TIFFTAG_IMAGEDESCRIPTION'
    # - these values can be overridden here - otherwise the values of 'Source_filepath' will be used
    # if 'TIFFTAG_SOFTWARE' is not set: the bounds in wsg84 will be used
    # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
    metadata_list = [
           (u'name',u'1861 Mercator Europe'),
    #       (u'TIFFTAG_SOFTWARE','mapmbtiles.ImageExporter'),
           (u'TIFFTAG_DATETIME','1861'),
           (u'TIFFTAG_ARTIST','Johnston, Alexander Keith, 1804-1871'),
           (u'TIFFTAG_HOSTCOMPUTER','http://www.davidrumsey.com/luna/servlet/detail/RUMSEY~8~1~21292~610098:Chart-of-the-world-on-Mercators-pro?sort=Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No&qvq=w4s:/what/Atlas%20Map/World%20Atlas/where/World;sort:Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No;lc:RUMSEY~8~1&mi=42&trs=583'),
           (u'TIFFTAG_COPYRIGHT','David Rumsey Map Collection'),
          ]
    ie_step_2.add_metadata(metadata_list)
    ie_step_2.export_image(bbox=(position_west,position_south,position_east,position_north),zoomlevel=8,imagepath=Output_filepath)
    # -8.4375, 35.46066995149529, 80.15625, 77.157162522661
    del ie_step_2
#----------------------------------------------------------------------------------
  if i_loop == 2:
   Source_filepath="%s/1861_Mercator_Countries.mbtiles" % input_directory
   if os.path.exists(Source_filepath):
    Output_filepath="%s/1861_Mercator_Berlin.tif" % output_directory
    print "Loop: ",i_loop,"\nSource: ",Source_filepath," \nOutput: ",Output_filepath
    ie_step_3 = ImageExporter(mbtiles_input=Source_filepath,jpg_quality=i_jpg_quality,tiff_compression=s_tiff_compression,tiff_predictor=i_tiff_predictor,tiff_zlevel=i_tiff_zlevel)
    position_west=12.0
    position_south=52.2969487
    position_east=14.0
    position_north=52.7002333
    # mbtiles_metadata 'name'        will be set as 'TIFFTAG_DOCUMENTNAME'
    # mbtiles_metadata 'description' will be set as 'TIFFTAG_IMAGEDESCRIPTION'
    # - these values can be overridden here - otherwise the values of 'Source_filepath' will be used
    # if 'TIFFTAG_SOFTWARE' is not set: the bounds in wsg84 will be used
    # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
    metadata_list = [
           (u'name',u'1861 Mercator Berlin'),
    #       (u'TIFFTAG_SOFTWARE','mapmbtiles.ImageExporter'),
           (u'TIFFTAG_DATETIME','1861'),
           (u'TIFFTAG_ARTIST','Johnston, Alexander Keith, 1804-1871'),
           (u'TIFFTAG_HOSTCOMPUTER','http://www.davidrumsey.com/luna/servlet/detail/RUMSEY~8~1~21292~610098:Chart-of-the-world-on-Mercators-pro?sort=Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No&qvq=w4s:/what/Atlas%20Map/World%20Atlas/where/World;sort:Pub_List_No_InitialSort%2CPub_Date%2CPub_List_No%2CSeries_No;lc:RUMSEY~8~1&mi=42&trs=583'),
           (u'TIFFTAG_COPYRIGHT','David Rumsey Map Collection'),
          ]
    ie_step_3.add_metadata(metadata_list)
    ie_step_3.export_image(bbox=(position_west,position_south,position_east,position_north),zoomlevel=8,imagepath=Output_filepath)
    del ie_step_3
    #  bounds[(11.25, 51.618016548773696, 14.0625, 53.330872983017045)]
#----------------------------------------------------------------------------------

