#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging,os
from gdal2mbtiles import MBTilesBuilder

logging.basicConfig(level=logging.DEBUG)
output_directory="output"
# i_parm=0: run programm ; = 1 show info only and return
i_parm=0
i_min_level=0
i_max_level=2
request_url_parms=""
Output_filepath="%s/Request_Url_Sample.mbtiles" % output_directory
#----------------------------------------------------------------------------------
for i_loop in range(i_min_level,i_max_level+1):
 print  "Loop: ",i_loop,"\nOutput: ",Output_filepath
 if i_loop == 0:
  s_wms_server="http://fbinter.stadt-berlin.de/fb/wms/senstadt/berlinzoom"
  s_wms_layers="0"
  request_url_parms="fill"
  parm_wms_options=dict(format="image/jpeg")
  mb_step1 = MBTilesBuilder(mbtiles_output=Output_filepath,request_url=request_url_parms,wms_server=s_wms_server,wms_layers=[s_wms_layers],wms_options=parm_wms_options)
  # 13d21'21.57"E, 52d32'32.85"N
  # 13d24'22.28"E, 52d30'44.83"N
  # mb_step1.add_coverage(bbox=(13.355992, 52.482242, 13.456086, 52.542458), zoomlevels=[16,17])
  mb_step1.add_coverage(bbox=(13.364714,52.502194,13.394683,52.528694), zoomlevels=[11,12,13,14,15,16,17])
  # use u'äöüßÄÖÜ' syntax in insure unicode when inserting into sqlite3
  metadata_list = [
          (u'name',u'Request_Url-Sample'),
          (u'description',u'Zoom 16+17 from different Serers'),
          (u'center',u'13.3777065575123,52.5162690144797,17'),
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
 if i_loop == 1:
  s_wms_server="http://fbinter.stadt-berlin.de/fb/wms/senstadt/k_vms_detailnetz_wms_spatial"
  s_wms_layers="0,1,2,3,4"
  request_url_parms="fill,replace"
  parm_wms_options=dict(format="image/jpeg")
  mb_step2 = MBTilesBuilder(mbtiles_output=Output_filepath,request_url=request_url_parms,wms_server=s_wms_server,wms_layers=[s_wms_layers],wms_options=parm_wms_options)
  # up to zoom 15 (inclusive) : blank images
  mb_step2.add_coverage(bbox=(13.373846,52.513583,13.386250,52.521147), zoomlevels=[17])
  # with 'False' we are telling MBTilesBuilder to exit if the file exists
  # with 'True' we are telling MBTilesBuilder to use an existing mbtiles file
  # i_parm=0: run programm ; = 1 show info only and return
  # - if no parmameters are given: 'False,0' is assumed
  # i_parm=1
  mb_step2.run(True,i_parm)
  del mb_step2
#----------------------------------------------------------------------------------
 if i_loop == 2:
  request_url_parms="load"
  mb_step3 = MBTilesBuilder(mbtiles_output=Output_filepath,request_url=request_url_parms)
  i_parm=0
  mb_step3.run(True,i_parm)
  del mb_step3
 #----------------------------------------------------------------------------------

