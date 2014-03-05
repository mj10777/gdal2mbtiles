#!/usr/bin/env python
# -*- coding: utf-8 -*-
#******************************************************************************
#  $Id: gdal2mbtiles.py 15748 2008-11-17 16:30:54Z klokan $
#
# Project:  Google Summer of Code 2007, 2008 (http://code.google.com/soc/)
# - adapted for mbtiles - 2014 Mark Johnson
# Support:  BRGM (http://www.brgm.fr)
# Purpose:  Convert a raster into TMS (Tile Map Service) tiles in a directory.
#           - generate Google Earth metadata (KML SuperOverlay)
#           - generate simple HTML viewer based on Google Maps and OpenLayers
#           - support of global tiles (Spherical Mercator) for compatibility
#               with interactive web maps a la Google Maps
# Author:   Klokan Petr Pridal, klokan at klokan dot cz
# Web:      http://www.klokan.cz/projects/gdal2mbtiles/
# GUI:      https://github.com/mj10777/mapmbtiles
#
###############################################################################
# Copyright (c) 2008, Klokan Petr Pridal
# - adapted for mbtiles - 2014 Mark Johnson
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#******************************************************************************

from osgeo import gdal,osr

import collections
import json
import logging
import mimetypes
import math
import operator
import os
from pkg_resources import parse_version
import re
import shutil
import sqlite3
import sys
import tempfile
import urllib
import urllib2
from urlparse import urlparse

from gettext import gettext as _

from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED

from io import BytesIO

has_pil = False

try:
 from PIL import Image, ImageEnhance
 has_pil = True
 import numpy
 from xml.etree import ElementTree
 import osgeo.gdal_array as gdalarray
except:
 # 'antialias' resampling is not available
 pass

__version__ = "$Id: gdal2mbtiles.py 15748 2013-12-29 16:30:54Z mj10777 $"

# ---- from globalmercator.py
logger = logging.getLogger(__name__)

MAXZOOMLEVEL = 32
DEG_TO_RAD = math.pi/180
RAD_TO_DEG = 180/math.pi
MAX_LATITUDE = 85.0511287798
EARTH_RADIUS = 6378137

def minmax (a,b,c):
 a = max(a,b)
 a = min(a,c)
 return a

class InvalidCoverageError(Exception):
 """ Raised when coverage bounds are invalid """
 pass

class GlobalMercator(object):
 """
 TMS Global Mercator Profile
 ---------------------------

 Functions necessary for generation of tiles in Spherical Mercator projection,
 EPSG:900913 (EPSG:gOOglE, Google Maps Global Mercator), EPSG:3785, OSGEO:41001.

 Such tiles are compatible with Google Maps, Microsoft Virtual Earth, Yahoo Maps,
 UK Ordnance Survey OpenSpace API, ...
 and you can overlay them on top of base maps of those web mapping applications.

 Pixel and tile coordinates are in TMS notation (origin [0,0] in bottom-left).

 What coordinate conversions do we need for TMS Global Mercator tiles::

      LatLon      <->       Meters      <->     Pixels    <->       Tile

  WGS84 coordinates   Spherical Mercator  Pixels in pyramid  Tiles in pyramid
      lat/lon            XY in metres     XY pixels Z zoom      XYZ from TMS
     EPSG:4326           EPSG:900913
      .----.              ---------               --                TMS
     /      \     <->     |       |     <->     /----/    <->      Google
     \      /             |       |           /--------/          QuadTree
      -----               ---------         /------------/
    KML, public         WebMapService         Web Clients      TileMapService

 What is the coordinate extent of Earth in EPSG:900913?

   [-20037508.342789244, -20037508.342789244, 20037508.342789244, 20037508.342789244]
   Constant 20037508.342789244 comes from the circumference of the Earth in meters,
   which is 40 thousand kilometers, the coordinate origin is in the middle of extent.
      In fact you can calculate the constant as: 2 * math.pi * 6378137 / 2.0
   $ echo 180 85 | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:900913
   Polar areas with abs(latitude) bigger then 85.05112878 are clipped off.

 What are zoom level constants (pixels/meter) for pyramid with EPSG:900913?

   whole region is on top of pyramid (zoom=0) covered by 256x256 pixels tile,
   every lower zoom level resolution is always divided by two
   initialResolution = 20037508.342789244 * 2 / 256 = 156543.03392804062

 What is the difference between TMS and Google Maps/QuadTree tile name convention?

   The tile raster itself is the same (equal extent, projection, pixel size),
   there is just different identification of the same raster tile.
   Tiles in TMS are counted from [0,0] in the bottom-left corner, id is XYZ.
   Google placed the origin [0,0] to the top-left corner, reference is XYZ.
   Microsoft is referencing tiles by a QuadTree name, defined on the website:
   http://msdn2.microsoft.com/en-us/library/bb259689.aspx

 The lat/lon coordinates are using WGS84 datum, yeh?

   Yes, all lat/lon we are mentioning should use WGS84 Geodetic Datum.
   Well, the web clients like Google Maps are projecting those coordinates by
   Spherical Mercator, so in fact lat/lon coordinates on sphere are treated as if
   the were on the WGS84 ellipsoid.

   From MSDN documentation:
   To simplify the calculations, we use the spherical form of projection, not
   the ellipsoidal form. Since the projection is used only for map display,
   and not for displaying numeric coordinates, we don't need the extra precision
   of an ellipsoidal projection. The spherical projection causes approximately
   0.33 percent scale distortion in the Y direction, which is not visually noticable.

 How do I create a raster in EPSG:900913 and convert coordinates with PROJ.4?

   You can use standard GIS tools like gdalwarp, cs2cs or gdaltransform.
   All of the tools supports -t_srs 'epsg:900913'.

   For other GIS programs check the exact definition of the projection:
   More info at http://spatialreference.org/ref/user/google-projection/
   The same projection is degined as EPSG:3785. WKT definition is in the official
   EPSG database.

   Proj4 Text:
     +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
     +k=1.0 +units=m +nadgrids=@null +no_defs

   Human readable WKT format of EPGS:900913:
      PROJCS["Google Maps Global Mercator",
          GEOGCS["WGS 84",
              DATUM["WGS_1984",
                  SPHEROID["WGS 84",6378137,298.2572235630016,
                      AUTHORITY["EPSG","7030"]],
                  AUTHORITY["EPSG","6326"]],
              PRIMEM["Greenwich",0],
              UNIT["degree",0.0174532925199433],
              AUTHORITY["EPSG","4326"]],
          PROJECTION["Mercator_1SP"],
          PARAMETER["central_meridian",0],
          PARAMETER["scale_factor",1],
          PARAMETER["false_easting",0],
          PARAMETER["false_northing",0],
          UNIT["metre",1,
              AUTHORITY["EPSG","9001"]]]
 """

 WSG84 = 'EPSG:4326'
 Spherical_Mercator  = 'EPSG:3857'

 def __init__(self, tms_osm=False,tileSize=256,levels = [0]):
  "Initialize the TMS Global Mercator pyramid"
  self.tileSize = tileSize
  self.tms_osm=tms_osm
  self.initialResolution = 2 * math.pi * 6378137 / self.tileSize
  # 156543.03392804062 for tileSize 256 pixels
  self.originShift = 2 * math.pi * 6378137 / 2.0
  # 20037508.342789244
  if levels:
   self.init_levels(levels)

 def LatLonToMeters(self, lat, lon ):
  "Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
  mx = lon * self.originShift / 180.0
  my = math.log( math.tan((90 + lat) * math.pi / 360.0 )) / (math.pi / 180.0)
  my = my * self.originShift / 180.0
  return mx, my

 def BoundsToMeters(self, tile_bounds ):
  "Converts given lat/lon Bounds in WGS84 Datum to XY in Spherical Mercator EPSG:900913 / 3857 / 3395"
  xmin,ymin,xmax,ymax=tile_bounds
  xmin,ymin=self.LatLonToMeters(ymin,xmin)
  xmax,ymax=self.LatLonToMeters(ymax,xmax)
  return (xmin,ymin,xmax,ymax)

 def MetersToLatLon(self, mx, my ):
  "Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum"
  lon = (mx / self.originShift) * 180.0
  lat = (my / self.originShift) * 180.0
  lat = 180 / math.pi * (2 * math.atan( math.exp( lat * math.pi / 180.0)) - math.pi / 2.0)
  return lat, lon

 def PixelsToMeters(self, px, py, zoom):
  "Converts pixel coordinates in given zoom level of pyramid to EPSG:900913"
  res = self.Resolution( zoom )
  mx = px * res - self.originShift
  my = py * res - self.originShift
  return mx, my

 def MetersToPixels(self, mx, my, zoom):
  "Converts EPSG:900913 to pyramid pixel coordinates in given zoom level"
  res = self.Resolution( zoom )
  px = (mx + self.originShift) / res
  py = (my + self.originShift) / res
  return px, py

 def PixelsToTile(self, px, py):
  "Returns a tile covering region in given pixel coordinates"
  tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
  ty_tms = int( math.ceil( py / float(self.tileSize) ) - 1 )
  return tx, ty_tms

 def PixelsToRaster(self, px, py, zoom):
  "Move the origin of pixel coordinates to top-left corner"
  mapSize = self.tileSize << zoom
  return px, mapSize - py

 def MetersToTile(self, mx, my, zoom):
  "Returns tile for given mercator coordinates"
  px, py = self.MetersToPixels( mx, my, zoom)
  return self.PixelsToTile( px, py)

 def TileBounds(self, tx, ty_tms, zoom):
  "Returns bounds of the given tile in EPSG:900913 coordinates"
  minx, miny = self.PixelsToMeters( tx*self.tileSize, ty_tms*self.tileSize, zoom )
  maxx, maxy = self.PixelsToMeters( (tx+1)*self.tileSize, (ty_tms+1)*self.tileSize, zoom )
  return ( minx, miny, maxx, maxy )

 def TileLatLonBounds(self, tx, ty_tms, zoom ):
  "Returns bounds of the given tile in latutude/longitude using WGS84 datum"
  bounds = self.TileBounds( tx, ty_tms, zoom)
  minLat, minLon = self.MetersToLatLon(bounds[0], bounds[1])
  maxLat, maxLon = self.MetersToLatLon(bounds[2], bounds[3])
  return ( minLat, minLon, maxLat, maxLon )

 def Resolution(self, zoom ):
  "Resolution (meters/pixel) for given zoom level (measured at Equator)"
  # return (2 * math.pi * 6378137) / (self.tileSize * 2**zoom)
  return self.initialResolution / (2**zoom)

 def ZoomForPixelSize(self, pixelSize ):
  "Maximal scaledown zoom of the pyramid closest to the pixelSize."
  for i in range(MAXZOOMLEVEL):
   if pixelSize > self.Resolution(i):
    if i!=0:
     return i-1
    else:
     return 0 # We don't want to scale up

 def GoogleTile(self, tx, ty_tms, zoom):
  "Converts TMS tile coordinates to Google Tile coordinates"
  # coordinate origin is moved from bottom-left to top-left corner of the extent
  return tx, (2**zoom - 1) - ty_tms

 def QuadTree(self, tx, ty_tms, zoom ):
  "Converts TMS tile coordinates to Microsoft QuadTree"
  quadKey = ""
  ty_osm = (2**zoom - 1) - ty_tms
  for i in range(zoom, 0, -1):
   digit = 0
   mask = 1 << (i-1)
   if (tx & mask) != 0:
    digit += 1
   if (ty_osm & mask) != 0:
    digit += 2
   quadKey += str(digit)
  return quadKey

 # needed Landez project functions
 def init_levels(self,levels = [0]):
  self.Bc = []
  self.Cc = []
  self.zc = []
  self.Ac = []
  if not levels:
   raise InvalidCoverageError(_("Wrong zoom levels."))
  self.levels = levels
  self.maxlevel = max(levels) + 1
  c = self.tileSize
  for d in range(self.maxlevel):
   e = c/2;
   self.Bc.append(c/360.0)
   self.Cc.append(c/(2 * math.pi))
   self.zc.append((e,e))
   self.Ac.append(c)
   c *= 2

 # from Landez project
 # https://github.com/makinacorpus/landez
 def project_pixels(self,ll,zoom):
  d = self.zc[zoom]
  e = round(d[0] + ll[0] * self.Bc[zoom])
  f = minmax(math.sin(DEG_TO_RAD * ll[1]),-0.9999,0.9999)
  g = round(d[1] + 0.5*log((1+f)/(1-f))*-self.Cc[zoom])
  return (e,g)

 # from Landez project
 # https://github.com/makinacorpus/landez
 def unproject_pixels(self,px,zoom):
  e = self.zc[zoom]
  f = (px[0] - e[0])/self.Bc[zoom]
  g = (px[1] - e[1])/-self.Cc[zoom]
  h = RAD_TO_DEG * ( 2 * math.atan(math.exp(g)) - 0.5 * math.pi)
  if not self.tms_osm:
   # return tms numbering
   h = - h
  return (f,h)

 def tile_at(self, zoom, position):
  """
  Returns a tuple of (z, x, y)
  """
  x, y = self.project_pixels(position, zoom)
  return (zoom, int(x/self.tileSize), int(y/self.tileSize))

 # from Landez project
 # https://github.com/makinacorpus/landez
 def tile_bbox(self, (z, x, y_tms)):
  """
  Returns the WGS84 bbox of the specified tile [tms notation] west/north/east/south
  """
  topleft = (x * self.tileSize, (y_tms + 1) * self.tileSize)
  bottomright = ((x + 1) * self.tileSize, y_tms * self.tileSize)
  wn = self.unproject_pixels(topleft, z)
  es = self.unproject_pixels(bottomright, z)
  return wn + es

 # from Landez project
 # https://github.com/makinacorpus/landez
 def project(self, (lng, lat)):
  """
  Returns the coordinates in meters from WGS84
  """
  x = lng * DEG_TO_RAD
  lat = max(min(MAX_LATITUDE, lat), -MAX_LATITUDE)
  y = lat * DEG_TO_RAD
  y = math.log(math.tan((math.pi / 4) + (y / 2)))
  return (x*EARTH_RADIUS, y*EARTH_RADIUS)

 # from Landez project
 # https://github.com/makinacorpus/landez
 def unproject(self, (x, y)):
  """
  Returns the coordinates from position in meters
  """
  lng = x/EARTH_RADIUS * RAD_TO_DEG
  lat = 2 * math.atan(math.exp(y/EARTH_RADIUS)) - math.pi/2 * RAD_TO_DEG
  return (lng, lat)

 # from Landez project
 # https://github.com/makinacorpus/landez
 # this return list from north to south and west to east and precise bounds
 def tileslist(self, bbox):
  if not self.levels:
   raise InvalidCoverageError(_("Wrong zoom levels. [init_levels called?]"))
  if len(bbox) != 4:
   raise InvalidCoverageError(_("Wrong format of bounding box."))
  xmin, ymin, xmax, ymax = bbox
  if abs(xmin) > 180 or abs(xmax) > 180 or abs(ymin) > 90 or abs(ymax) > 90:
   s_error="";
   if abs(xmin) > 180:
    s_error+="abs(xmin) (%d) > 180 ; " % abs(xmin)
   if abs(xmax) > 180:
    s_error+="abs(xmax) (%d) > 180 ; " % abs(xmax)
   if abs(ymin) > 90:
    s_error+="abs(ymin) (%d) > 90 ; " % abs(ymin)
   if abs(ymax) > 90:
    s_error+="abs(ymax) (%d) > 90 ; " % abs(ymax)
   raise InvalidCoverageError(_("-E-> Some coordinates exceed [-180,+180], [-90, 90]. [%s]" % (s_error)))
  if xmin >= xmax or ymin >= ymax:
   s_error="";
   if xmin >= xmax:
    s_error+="xmin(%f) >= xmax(%f) ; " % (xmin,xmax)
   if ymin >= ymax:
    s_error+="ymin(%f) >= ymax(%f) ; " % (ymin,ymax)
   raise InvalidCoverageError(_("-E-> Bounding box format is (xmin, ymin, xmax, ymax) [%s]" % (s_error)))
  ll0 = (xmin, ymax)  # left top
  ll1 = (xmax, ymin)  # right bottom
  l = []
  for z in self.levels:
   px0 = self.project_pixels(ll0,z)
   px1 = self.project_pixels(ll1,z)
   for x in range(int(px0[0]/self.tileSize),int(px1[0]/self.tileSize)+1):
    if (x < 0) or (x >= 2**z):
     continue
    for y in range(int(px0[1]/self.tileSize),int(px1[1]/self.tileSize)+1):
     if (y < 0) or (y >= 2**z):
      continue
     if not self.tms_osm:
      # return tms numbering
      y = ((2**z-1) - y)
     lng_min,lat_max,lng_max, lat_min=self.tile_bbox((z,x,y))
     # logger.info(_("Save resulting image xy'%s' - bounds[%s]") % ((x,y),(lng_min,lat_max,lng_max, lat_min)))
     if lng_max > xmax:
      xmax=lng_max
     if lng_min < xmin:
      xmin=lng_min
     if lat_max > ymax:
      ymax=lat_max
     if lat_min < ymin:
      ymin=lat_min
     l.append((z, x, y))
  tile_bounds=(xmin,ymin,xmax,ymax)
  return (l,tile_bounds)

 # from Landez project
 # https://github.com/makinacorpus/landez
 def project_pixels(self,ll,zoom):
  d = self.zc[zoom]
  e = round(d[0] + ll[0] * self.Bc[zoom])
  f = minmax(math.sin(DEG_TO_RAD * ll[1]),-0.9999,0.9999)
  g = round(d[1] + 0.5*math.log((1+f)/(1-f))*-self.Cc[zoom])
  return (e,g)

#---------------------

class GlobalGeodetic(object):
 """
 TMS Global Geodetic Profile
 ---------------------------

 Functions necessary for generation of global tiles in Plate Carre projection,
 EPSG:4326, "unprojected profile".

 Such tiles are compatible with Google Earth (as any other EPSG:4326 rasters)
 and you can overlay the tiles on top of OpenLayers base map.

 Pixel and tile coordinates are in TMS notation (origin [0,0] in bottom-left).

 What coordinate conversions do we need for TMS Global Geodetic tiles?

   Global Geodetic tiles are using geodetic coordinates (latitude,longitude)
   directly as planar coordinates XY (it is also called Unprojected or Plate
   Carre). We need only scaling to pixel pyramid and cutting to tiles.
   Pyramid has on top level two tiles, so it is not square but rectangle.
   Area [-180,-90,180,90] is scaled to 512x256 pixels.
   TMS has coordinate origin (for pixels and tiles) in bottom-left corner.
   Rasters are in EPSG:4326 and therefore are compatible with Google Earth.

      LatLon      <->      Pixels      <->     Tiles

  WGS84 coordinates   Pixels in pyramid  Tiles in pyramid
      lat/lon         XY pixels Z zoom      XYZ from TMS
     EPSG:4326
      .----.                ----
     /      \     <->    /--------/    <->      TMS
     \      /         /--------------/
      -----        /--------------------/
    WMS, KML    Web Clients, Google Earth  TileMapService
 """

 def __init__(self, tileSize = 256):
  self.tileSize = tileSize

 def LatLonToPixels(self, lat, lon, zoom):
  "Converts lat/lon to pixel coordinates in given zoom of the EPSG:4326 pyramid"

  res = 180.0 / self.tileSize / 2**zoom
  px = (180 + lat) / res
  py = (90 + lon) / res
  return px, py

 def PixelsToTile(self, px, py):
  "Returns coordinates of the tile [tms] covering region in pixel coordinates"

  tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
  ty_tms = int( math.ceil( py / float(self.tileSize) ) - 1 )
  return tx, ty_tms

 def LatLonToTile(self, lat, lon, zoom):
  "Returns the tile for zoom which covers given lat/lon coordinates"

  px, py = self.LatLonToPixels( lat, lon, zoom)
  return self.PixelsToTile(px,py)

 def Resolution(self, zoom ):
  "Resolution (arc/pixel) for given zoom level (measured at Equator)"

  return 180.0 / self.tileSize / 2**zoom
  #return 180 / float( 1 << (8+zoom) )

 def ZoomForPixelSize(self, pixelSize ):
  "Maximal scaledown zoom of the pyramid closest to the pixelSize."

  for i in range(MAXZOOMLEVEL):
   if pixelSize > self.Resolution(i):
    if i!=0:
     return i-1
    else:
     return 0 # We don't want to scale up

 def TileBounds(self, tx, ty_tms, zoom):
  "Returns bounds of the given tile [tms]"
  res = 180.0 / self.tileSize / 2**zoom
  return (
   tx*self.tileSize*res - 180,
   ty_tms*self.tileSize*res - 90,
   (tx+1)*self.tileSize*res - 180,
   (ty_tms+1)*self.tileSize*res - 90
  )

 def TileLatLonBounds(self, tx, ty_tms, zoom):
  "Returns bounds of the given tile [tms] in the SWNE form"
  b = self.TileBounds(tx, ty_tms, zoom)
  return (b[1],b[0],b[3],b[2])

#---------------------
# TODO: Finish Zoomify implemtentation!!!
class Zoomify(object):
 """
 Tiles compatible with the Zoomify viewer
 ----------------------------------------
 """

 def __init__(self, width, height, tilesize = 256, tileformat='jpg'):
  """Initialization of the Zoomify tile tree"""

  self.tilesize = tilesize
  self.tileformat = tileformat
  imagesize = (width, height)
  tiles = ( math.ceil( width / tilesize ), math.ceil( height / tilesize ) )

  # Size (in tiles) for each tier of pyramid.
  self.tierSizeInTiles = []
  self.tierSizeInTiles.push( tiles )

  # Image size in pixels for each pyramid tierself
  self.tierImageSize = []
  self.tierImageSize.append( imagesize );

  while (imagesize[0] > tilesize or imageSize[1] > tilesize ):
   imagesize = (math.floor( imagesize[0] / 2 ), math.floor( imagesize[1] / 2) )
   tiles = ( math.ceil( imagesize[0] / tilesize ), math.ceil( imagesize[1] / tilesize ) )
   self.tierSizeInTiles.append( tiles )
   self.tierImageSize.append( imagesize )

  self.tierSizeInTiles.reverse()
  self.tierImageSize.reverse()

  # Depth of the Zoomify pyramid, number of tiers (zoom levels)
  self.numberOfTiers = len(self.tierSizeInTiles)

  # Number of tiles up to the given tier of pyramid.
  self.tileCountUpToTier = []
  self.tileCountUpToTier[0] = 0
  for i in range(1, self.numberOfTiers+1):
   self.tileCountUpToTier.append(
    self.tierSizeInTiles[i-1][0] * self.tierSizeInTiles[i-1][1] + self.tileCountUpToTier[i-1]
   )

 def tilefilename(self, x, y, z):
  """Returns filename for tile with given coordinates"""

  tileIndex = x + y * self.tierSizeInTiles[z][0] + self.tileCountUpToTier[z]
  return os.path.join("TileGroup%.0f" % math.floor( tileIndex / 256 ),
   "%s-%s-%s.%s" % ( z, x, y, self.tileformat))

# ----- end globalmercator.py
# ---- from mbtiles.py
""" Default tiles URL """
DEFAULT_TILES_URL = "http://{s}.tile.openstreetmap.org/{z}/{x}/{y_osm}.png"
""" Default tiles subdomains """
DEFAULT_TILES_SUBDOMAINS = list("abc")
""" Base temporary folder """
DEFAULT_TMP_DIR = os.path.join(tempfile.gettempdir(), 'mapmbtiles')
""" Default output MBTiles file """
DEFAULT_MBTILES_OUTPUT = os.path.join(os.getcwd(), "mbtiles_output.mbtiles")
""" Default tile size in pixels (*useless* in remote rendering) """
DEFAULT_TILE_SIZE = 256
""" Default tile format (mime-type) """
DEFAULT_TILE_FORMAT = 'image/jpeg'
""" Number of retries for remove tiles downloading """
DOWNLOAD_RETRIES = 10
""" Path to fonts for Mapnik rendering """
TRUETYPE_FONTS_PATH = '/usr/share/fonts/truetype/'

# =============================================================================
# UPDATE map SET tile_id = replace(tile_id, '.None', '.tms');
# UPDATE images SET tile_id = replace(tile_id, '.None', '.tms');
# =============================================================================
# base on code faound at. https://github.com/mapbox/mbutil/blob/master/mbutil/util.py
# - creation of mbtiles based on logic used in the geopaparazzi project
# - functions added to be used by gdalt2tiles.py [gdal2mbtiles.py]
# =============================================================================
class MbTiles(object):
 def __init__(self):
  self.verbose=False
  self.tilesize=256
  self.bounds_west=-180.0
  self.bounds_east=180.0
  self.bounds_north=85.05113
  self.bounds_south=-85.05113
  self.min_zoom=0
  self.max_zoom=22
  self.default_zoom=1
  self.mbtiles_name=""
  self.mbtiles_description=""
  self.mbtiles_type="baselayer"
  self.mbtiles_version='1.1'
  self.s_y_type = 'tms'
  self.tms_osm=False
  self.mbtiles_format='jpg'
  self.jpg_quality=75
  self.pil_format='JPEG'
  self.mbtiles_minzoom='0'
  self.mbtiles_maxzoom='22'
  self.mbtiles_bounds="%f,%f,%f,%f"% (self.bounds_west,self.bounds_south,self.bounds_east,self.bounds_north)
  self.center_x=(self.bounds_east+self.bounds_west)/2
  self.center_y=(self.bounds_north+self.bounds_south)/2
  self.mbtiles_center="%f,%f,%s"%(self.center_x,self.center_y,self.default_zoom)
  self.i_request_url_read_value=0
  self.i_request_url_read_db=1
  self.i_request_url_count_create=2
  self.i_request_url_count_drop=3
  self.i_request_url_count=-1

 def open_db(self,s_path_db,mbtiles_dir,mbtiles_format,s_y_type,verbose=False):
  self.s_path_db = s_path_db
  self.mbtiles_dir=mbtiles_dir.strip()
  if self.mbtiles_dir == "":
   self.mbtiles_dir=os.path.dirname(self.mbtiles_file)+ '/'
  if not os.path.exists(self.mbtiles_dir):
   os.makedirs(self.mbtiles_dir)
  if s_y_type == "osm":
   self.s_y_type=s_y_type
   self.tms_osm=True
  if mbtiles_format.find("image/") == -1:
   mbtiles_format.replace("image/ ","")
  if mbtiles_format == "jpeg":
   mbtiles_format="jpg"
  if mbtiles_format == "png":
   self.mbtiles_format=mbtiles_format
   self.pil_format='PNG'
  self.verbose=verbose
  # self.verbose=True
  # setting a default value
  self.mbtiles_name=os.path.splitext(os.path.basename( self.s_path_db ))[0]
  self.mbtiles_description=self.mbtiles_name.replace("."," ")
  self.mbtiles_description=self.mbtiles_description.replace("_"," ")
  db_create=os.path.exists(self.s_path_db)
  self.sqlite3_connection=self.mbtiles_connect(s_path_db)
  self.mbtiles_cursor = self.sqlite3_connection.cursor()
  # self.optimize_connection()
  if not db_create:
   if self.verbose:
    logger.info(_("MbTiles : [open_db] : creating: [%s]") % self.s_path_db)
   self.mbtiles_create()
  else:
   if self.verbose:
    logger.info(_("MbTiles : [open_db] : opening: [%s]") % self.s_path_db)
   self.fetch_metadata()
   self.i_request_url_count=self.get_request_url_count(self.i_request_url_read_db)

 def close_db(self):
  if self.verbose:
   logger.info(_("MbTiles : [close_db] : closing: [%s]") % self.s_path_db)
  # if self.mbtiles_cursor:
  #  self.mbtiles_cursor.close()
  if self.sqlite3_connection:
   self.sqlite3_connection.close()
   if os.path.exists("%s-journal" % self.s_path_db):
    os.remove("%s-journal" % self.s_path_db)

 def flip_y(self,zoom, y):
  # z=19: y[ 352336 ]  y_osm[ 171951 ]  y_tms[ 352336 ]
  # SELECT  CastToInteger(pow(2,19) - 1 - 352336) 'y_osm',CastToInteger(pow(2,19) - 1 - 171951) 'y_tms';
  return (2**zoom-1) - y

 def mbtiles_create(self):
  self.mbtiles_cursor.execute("""CREATE TABLE android_metadata (locale text);""")
  self.mbtiles_cursor.execute("""CREATE TABLE metadata (name text, value text);""")
  self.mbtiles_cursor.execute("""CREATE TABLE grid_key (grid_id TEXT,key_name TEXT);""")
  self.mbtiles_cursor.execute("""CREATE TABLE grid_utfgrid (grid_id TEXT,grid_utfgrid BLOB);""")
  self.mbtiles_cursor.execute("""CREATE TABLE keymap (key_name TEXT,key_json TEXT);""")
  self.mbtiles_cursor.execute("""CREATE TABLE images (tile_data blob,tile_id text);""")
  self.mbtiles_cursor.execute("""CREATE TABLE map (zoom_level INTEGER,tile_column INTEGER,tile_row INTEGER,tile_id TEXT,grid_id TEXT);""")
  self.mbtiles_cursor.execute("""CREATE VIEW tiles AS SELECT map.zoom_level AS zoom_level,map.tile_column AS tile_column,map.tile_row AS tile_row,images.tile_data AS tile_data FROM map JOIN images ON images.tile_id = map.tile_id ORDER BY zoom_level,tile_column,tile_row;""")
  self.mbtiles_cursor.execute("""CREATE VIEW grids AS SELECT map.zoom_level AS zoom_level,map.tile_column AS tile_column,map.tile_row AS tile_row,grid_utfgrid.grid_utfgrid AS grid FROM map JOIN grid_utfgrid ON grid_utfgrid.grid_id = map.grid_id;""")
  self.mbtiles_cursor.execute("""CREATE VIEW grid_data AS SELECT map.zoom_level AS zoom_level,map.tile_column AS tile_column,map.tile_row AS tile_row,keymap.key_name AS key_name,keymap.key_json AS key_json FROM map JOIN grid_key ON map.grid_id = grid_key.grid_id JOIN keymap ON grid_key.key_name = keymap.key_name;""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX name ON metadata (name);""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX grid_key_lookup ON grid_key (grid_id,key_name);""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX grid_utfgrid_lookup ON grid_utfgrid (grid_id);""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX keymap_lookup ON keymap (key_name);""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX images_id ON images (tile_id);""")
  self.mbtiles_cursor.execute("""CREATE UNIQUE INDEX map_index ON map (zoom_level, tile_column, tile_row);""")

 def mbtiles_connect(self,mbtiles_file):
  try:
   con = sqlite3.connect(mbtiles_file)
   return con
  except Exception, e:
   logger.error(_("MbTiles : Could not connect to database[%s]: Error %s:") % (mbtiles_file,e.args[0]))
   logger.exception(e)
   sys.exit(1)

 def fetch_metadata(self):
  self.metadata_return=self.mbtiles_cursor.execute('SELECT LOWER(name), value FROM metadata;').fetchall()
  self.metadata=dict(self.metadata_return)
  if self.metadata:
   self.mbtiles_name=self.metadata.get('name',self.mbtiles_name)
   self.mbtiles_description=self.metadata.get('description',self.mbtiles_description)
   self.mbtiles_type=self.metadata.get('type',self.mbtiles_type)
   self.mbtiles_version=self.metadata.get('version',self.mbtiles_version)
   self.s_y_type = self.metadata.get('tile_row_type',self.s_y_type)
   self.mbtiles_format=self.metadata.get('format',self.mbtiles_format)
   self.mbtiles_bounds=self.metadata.get('bounds',self.mbtiles_bounds)
   sa_bounds=self.mbtiles_bounds.split(",")
   if len(sa_bounds) == 4:
    self.bounds_west=float(sa_bounds[0])
    self.bounds_east=float(sa_bounds[2])
    self.bounds_north=float(sa_bounds[3])
    self.bounds_south=float(sa_bounds[1])
   self.mbtiles_minzoom=self.metadata.get('minzoom',self.mbtiles_minzoom)
   self.mbtiles_maxzoom=self.metadata.get('maxzoom',self.mbtiles_maxzoom)
   self.center_x=(self.bounds_east+self.bounds_west)/2
   self.center_y=(self.bounds_north+self.bounds_south)/2
   self.mbtiles_center="%f,%f,%s"%(self.center_x,self.center_y,self.mbtiles_minzoom)
   self.mbtiles_center=self.metadata.get('center',self.mbtiles_center)
   sa_center=self.mbtiles_center.split(",")
   if len(sa_center) == 3:
    self.center_x=float(sa_center[0])
    self.center_y=float(sa_center[1])
    self.default_zoom=int(sa_center[2])
  return self.metadata_return

 def optimize_connection(self):
  self.mbtiles_cursor.execute("""PRAGMA synchronous=0""")
  self.mbtiles_cursor.execute("""PRAGMA locking_mode=EXCLUSIVE""")
  self.mbtiles_cursor.execute("""PRAGMA journal_mode=DELETE""")

 def optimize_database(self):
  if self.verbose:
   logger.info(_("MbTiles : optimize_database: analyzing db [%s]") % "ANALYZE;")
  self.mbtiles_cursor.execute("""ANALYZE;""")
  if self.verbose:
   logger.info(_("MbTiles : optimize_database: cleaning db [%s]") % "VACUUM;")
  self.mbtiles_cursor.execute("""VACUUM;""")
  if self.verbose:
   logger.info(_("MbTiles : optimize_database: [%s]") % self.s_path_db)

 def insert_metadata(self,metadata_list):
  if metadata_list:
   if self.verbose:
    # logger.info(_("MbTiles : insert_metadata: [%s]") % metadata_list)
    pass
   try:
    # ERROR:mapmbtiles.mbtiles:MbTiles : insert_metadata: Error You must not use 8-bit bytestrings unless you use a text_factory that can interpret 8-bit bytestrings (like text_factory = str). It is highly recommended that you instead just switch your application to Unicode strings.:
    # repr(metadata_list)
    self.mbtiles_cursor.executemany("INSERT OR REPLACE INTO metadata VALUES(?,?)",metadata_list)
    self.sqlite3_connection.commit()
   except sqlite3.Error, e:
    self.sqlite3_connection.rollback()
    logger.error(_("MbTiles : insert_metadata: Error %s:") % e.args[0])
   self.fetch_metadata()

 def save_metadata(self):
  if not self.s_y_type:
   self.s_y_type="tms"
  values_list = [
         ('name', self.mbtiles_name),
         ('description',self.mbtiles_description ),
         ('type', self.mbtiles_type),
         ('version', self.mbtiles_version),
         ('tile_row_type',self.s_y_type),
         ('format', self.mbtiles_format),
         ('bounds', self.mbtiles_bounds),
         ('center', self.mbtiles_center),
         ('minzoom', self.mbtiles_minzoom),
         ('maxzoom', self.mbtiles_maxzoom)
        ]
  self.insert_metadata(values_list)

 def tile_id_to_zxy(self,s_tile_id):
  """
  tile_id format:
  ---------------------------
  there are 2 type of formats: for tiles or rgb-images
  ---
  format-systax: z/r-x/g-y/b.tms/osm/tms
  ---
  'r-' or 'z-' : the red value of the rgb or the zoom-level
  ---
  'g-' or 'x-' : the green value of the rgb or the x-tile
  ---
  'b.' or 'y.' : the blue value of the rgb or the y-tile
  ---
  extention: '.rgb' or '.tms' or '.osm'
  - 'rgb' will be blank images where all pixels have the same rgb-value
  - 'osm' will be tile images where y position is in Open-Street-Map notation (North to South)
  - 'tms' will be tile images where y position is in TMS notation (South to North)
  ---
  """
  sa_split=s_tile_id.split(".")
  value_z=0
  value_x=0
  value_y=0
  value_type=""
  if len(sa_split) == 2:
   value_type=sa_split[1]
   sa_split=sa_split[0].split("-")
  if len(sa_split) == 3:
   value_z=int(sa_split[0])
   value_x=int(sa_split[1])
   value_y=int(sa_split[2])
  return value_z,value_x,value_y,value_type

 def get_request_url_count(self,i_parm):
  """
  i_parm values:
  ---------------------------
  0: return existing value [set when database was opended,not reading the table] [i_request_url_read_value]
  1 : return existing value return existing value [reading the table with count after checking if it exits] [i_request_url_read_db]
  2: create table (if it does not exist) [i_request_url_count_create]
  3: delete table (if it does exist) [i_request_url_count_drop]
  """
  if i_parm == self.i_request_url_count_drop:
   s_sql_request_url = "DROP TABLE IF EXISTS request_url"
   try:
    self.mbtiles_cursor.executemany(s_sql_request_url)
   except sqlite3.Error, e:
    s_sql_request_url=""
    logger.error(_("MbTiles : drop_request_url: Error %s:") % e.args[0])
   if s_sql_request_url != "":
    # no error [DROP was compleated correctly] table has been deleted and is empty
    self.i_request_url_count = -1
  elif i_parm == self.i_request_url_count_create:
   s_sql_request_url = "CREATE TABLE IF NOT EXISTS request_url (tile_id TEXT PRIMARY KEY,tile_url TEXT)"
   try:
    self.mbtiles_cursor.execute(s_sql_request_url)
   except sqlite3.Error, e:
    s_sql_request_url=""
    logger.error(_("MbTiles : create_request_url: Error %s:") % e.args[0])
   if s_sql_request_url != "":
    # no error [CREATE was compleated correctly] table has been created and is empty
    self.i_request_url_count = 0
  elif i_parm == self.i_request_url_read_db:
   s_sql_request_url = "pragma table_info(request_url)"
   try:
    self.mbtiles_cursor.execute(s_sql_request_url)
    pragma_result = self.mbtiles_cursor.fetchone()
    if not pragma_result is None:
     i_field_count=0
     for i_field in range(len(pragma_result)):
      s_field=pragma_result[i_field]
      if s_field == "tile_id" or s_field == "tile_url":
       i_field_count+=1
     if i_field_count > 0:
      s_sql_request_url = "SELECT count(tile_id) AS count_id FROM request_url"
      self.mbtiles_cursor.execute(s_sql_request_url)
      count_result = self.mbtiles_cursor.fetchone()
      if not count_result is None:
       self.i_request_url_count=int(count_result[0])
   except sqlite3.Error, e:
    s_sql_request_url=""
    logger.error(_("MbTiles : count_request_url: Error %s:") % e.args[0])
  elif i_parm == self.i_request_url_count_read_value:
   # return existing value
   pass
  return  self.i_request_url_count

 def fill_request_url(self,bounds,zoomlevels,s_request_url_source,i_parm):
  """
  i_parm values
  ---------------------------
  0: if tile_id is already exists, do not overwrite
  1: replace tile_id
  ---------------------------
  fill the request_url table
  ---------------------------
  checking is done to neither the tile or request exist's
  uses the same logic as in geopaparazzi
  '{z}' - tms - zoom level
  '{x}' - tms - x tile
  '{y_oms}' - tms - y tile - oms notation
  '{y_tms}' - tms - y tile - tms notation
  if '{z}' is not found, wms logic will be used
  """
  b_tms=False
  if s_request_url_source.find('{z}') != -1:
   b_tms=True
  mercator = GlobalMercator(self.tms_osm,self.tilesize,zoomlevels)
  bboxlist,tile_bounds = mercator.tileslist(bounds)
  #if self.verbose:
  logger.debug(_("Found %s tiles for this area.") % len(bboxlist))
  i_count=0
  request_url_list = []
  i_count_tiles=0
  for (z, x, y) in sorted(bboxlist,key=operator.itemgetter(0,1,2)):
   i_count_tiles=self.count_tiles(z,x,y,1)
   if i_count_tiles == 0:
    # tile is not in database
    if i_count_tiles == 0:
     # only if 'replace' has not been used
     i_count_tiles=self.count_tiles(z,x,y,11)
    # tile is not in request_url
    if i_count_tiles == 0:
     if not self.tms_osm:
      y_tms=y
      y_osm=self.flip_y(z,y_tms)
     else:
      y_osm=y
      y_tms=self.flip_y(z,y_osm)
     s_tile_id="{0}-{1}-{2}.{3}".format(str(z),str(x),str(y),self.s_y_type)
     if b_tms:
      s_tile_url = s_request_url_source.format(**locals())
     else:
      mercator = GlobalMercator(True,self.tilesize,[z])
      bbox = mercator.tile_bbox((z, x, y_osm))
      bbox = ','.join(map(str, bbox))
      s_tile_url = urllib.unquote(s_request_url_source)
      s_tile_url += "&bbox=%s" % bbox   # commas are not encoded
     request_url_list.append((s_tile_id,s_tile_url))
     i_count+=1
     if i_count >= 100:
      # save after 100 entries
      self.insert_request_url(request_url_list)
      request_url_list = []
      i_count=0
   # save what has not been saved
   if i_count >= 0:
    self.insert_request_url(request_url_list)

 def insert_request_url(self,request_url_list):
  if request_url_list:
   if self.i_request_url_count < 0:
     # create the request_url table
     self.get_request_url_count(2)
   if self.verbose:
    # logger.info(_("MbTiles : insert_metadata:: [%s]") % metadata_list)
    pass
   try:
    self.mbtiles_cursor.executemany("INSERT OR REPLACE INTO request_url VALUES(?,?)",request_url_list)
    self.sqlite3_connection.commit()
   except sqlite3.Error, e:
    self.sqlite3_connection.rollback()
    logger.error(_("MbTiles : insert_request_url: Error %s:") % e.args[0])

 def delete_request_url(self,s_tile_id):
  if s_tile_id:
   if self.verbose:
    # logger.info(_("MbTiles : delete_request_url: [%s]") % metadata_list)
    pass
   tile_id = (s_tile_id,)
   try:
    self.mbtiles_cursor.execute("DELETE FROM request_url WHERE (tile_id = ?)",tile_id)
   except sqlite3.Error, e:
    logger.error(_("MbTiles : delete_request_url  (%s): Error %s:") % (s_tile_id,e.args[0]))
    s_tile_id=""
   if s_tile_id != "":
    self.i_request_url_count-=1
    if self.i_request_url_count <=0:
     # delete the request_url table
     self.get_request_url_count(3)

 def select_request_url(self,s_tile_id,limit_offset):
  if s_tile_id:
   if self.verbose:
    # logger.info(_("MbTiles : delete_request_url: [%s]") % metadata_list)
    pass
   try:
    result_requst_url=self.mbtiles_cursor.execute("SELECT tile_id,tile_url FROM 'request_url' WHERE ( tile_id = '?')",s_tile_id).fetchall()
   except sqlite3.Error, e:
    logger.error(_("MbTiles : select_request_url: Error %s:") % e.args[0])
  else:
   try:
    result_requst_url=self.mbtiles_cursor.execute("SELECT tile_id,tile_url FROM 'request_url' LIMIT ? OFFSET ?",limit_offset).fetchall()
   except sqlite3.Error, e:
    logger.error(_("MbTiles : select_request_url: Error %s:") % e.args[0])
  return result_requst_url

 def load_request_url(self,i_parm):
  """
  i_parm values
  ---------------------------
  0: not used (yet)
  ---------------------------
  """
  i_limit=100
  i_offset=0
  i_request_url_amount=self.i_request_url_count
  while i_request_url_amount > 0:
   if i_request_url_amount < i_limit:
    i_limit=i_request_url_amount
   result_requst_url=self.select_request_url("",(i_limit,i_offset))
   i_offset+=i_limit
   for i in range(len(result_requst_url)):
    s_tile_id = result_requst_url[i][0]
    s_tile_url = result_requst_url[i][1]
    request = urllib2.Request(s_tile_url )
    try:
     response = urllib2.urlopen(request)
    except urllib2.URLError, e:
     if e.code == 404:
      logger.error(_("MbTiles : dload_request_url: Error %s:") % e)
    else:
      # 200
     s_content_type=response.headers['content-type']
     if s_content_type.find('image/') != -1:
      image_data=response.read()
      if not image_data is None:
       z,x,y,s_type=self.tile_id_to_zxy(s_tile_id)
       self.insert_image(z,x,y,image_data)
       # self.i_request_url_count will be reduced and the table deleted when empty
       self.delete_request_url(s_tile_id)
     else:
      logger.error(_("MbTiles : load_request_url: image not returned: [%s]") % s_content_type)

 def insert_image(self,tz,tx,ty,image_data):
  if not self.s_y_type:
   self.s_y_type="tms"
  if not self.mbtiles_cursor:
   self.mbtiles_cursor = self.sqlite3_connection.cursor()
  s_tile_id="{0}-{1}-{2}.{3}".format(str(tz),str(tx),str(ty),self.s_y_type)
  s_tile_id,output_image=self.check_image(s_tile_id,image_data)
  if output_image:
   image_data=output_image
  sql_insert_map="INSERT OR REPLACE INTO map (tile_id,zoom_level,tile_column,tile_row,grid_id) VALUES(?,?,?,?,?);";
  map_values = [(s_tile_id,tz,tx,ty,'')]
  sql_insert_image="INSERT OR REPLACE INTO images (tile_id,tile_data) VALUES(?,?);"
  image_values = [(s_tile_id,buffer(image_data))]
  # sqlite3.Binary(image_data)
  if self.verbose:
   logger.info(_("MbTiles : insert_image: %d,%d,%d id[%s]") % (tz,tx,ty,s_tile_id))
  try:
   self.mbtiles_cursor.executemany(sql_insert_map,map_values)
   self.mbtiles_cursor.executemany(sql_insert_image,image_values)
   self.sqlite3_connection.commit()
  except sqlite3.Error, e:
   self.sqlite3_connection.rollback()
   logger.error(_("MbTiles : insert_image: Error %s:") % e.args[0])

 def check_image(self,s_tile_id,image_data):
  # a list of (count, color) tuples or None, max amount [we only want information about a blank image]
  output_data=None
  input_image = Image.open(BytesIO(image_data))
  colors = input_image.getcolors(1)
  if self.pil_format != input_image.format:
   if self.pil_format == "JPEG":
    if input_image.mode != "RGB":
     input_image=input_image.convert('RGB')
    # http://effbot.org/imagingbook/pil-index.htm#appendixes
    input_image.save(s_tile_id, format="JPEG", quality=self.jpg_quality, optimize=True, progressive=False)
   else:
    input_image.save(s_tile_id, format="PNG",optimize=True)
   f = open(s_tile_id,'rb')
   output_data = f.read()
   f.close()
   os.remove(s_tile_id)
   input_image = Image.open(BytesIO(output_data))
   colors = input_image.getcolors(1)
   # TypeError: 'buffer' does not have the buffer interface
  if colors:
   # MbTiles : check_image: tile_id[ 18-140789-176144.tms ] colors[ [(65536, (255, 255, 255, 255))] ]
   # color_values[ (65536, (255, 255, 255, 255)) ]
   color_values=colors[0]
   # rgb_values[ (255, 255, 255, 255) ]
   rgb_values=color_values[1]
   # r[ 255 ] g[ 255 ] b[ 255 ]
   r_value=rgb_values[0]
   g_value=rgb_values[1]
   b_value=rgb_values[2]
   s_tile_orig = s_tile_id
   # exception because of hex, avoid this type of formatting - use .format(..)
   s_tile_id = "%2x-%2x-%2x.rgb"%(int(r_value),int(g_value),int(b_value))
   # MbTiles : check_image: [(65536, (255, 255, 255))] 18-140785-176149.tms
   # colors = len(filter(None,image_img.histogram()))
  return s_tile_id,output_data

 def retrieve_blank_image(self,r,g,b):
  s_tile_id="{0}-{1}-{2}.{3}".format(str(r), str(g),str(b),"rgb")
  tile_id = (s_tile_id,)
  self.mbtiles_cursor.execute("SELECT tile_data FROM images WHERE tile_id = ?",tile_id)
  image_data = self.mbtiles_cursor.fetchone()
  if image_data is None:
   image = Image.new("RGB", (self.tilesize, self.tilesize), (r, g, b))
   s_tile_id="{0}-{1}-{2}.{3}".format(str(r), str(g),str(b),self.mbtiles_format)
   image.save(s_tile_id)
   input_file = open(s_tile_id, 'rb')
   if not input_file.closed:
    image_file = input_file.read()
    input_file.close()
    image_data = (image_file,)
    os.remove(s_tile_id)
  return image_data

 def retrieve_bounds(self):
  min_zoom=22
  max_zoom=0
  bounds_west=180.0
  bounds_east=-180.0
  bounds_north=-85.05113
  bounds_south=85.05113
  if not self.mbtiles_cursor:
   self.mbtiles_cursor = self.sqlite3_connection.cursor()
  mercator = GlobalMercator(self.tms_osm)
  zoom_levels = self.mbtiles_cursor.execute('SELECT DISTINCT(zoom_level) FROM map ORDER BY zoom_level;').fetchall()
  for i in range(len(zoom_levels)):
   i_zoom = int(zoom_levels[i][0])
   if i_zoom > max_zoom:
    max_zoom=i_zoom
   if i_zoom < min_zoom:
    min_zoom=i_zoom
   zoom_id = (str(i_zoom),)
   bounds_minmax = self.mbtiles_cursor.execute('SELECT min(tile_column),min(tile_row),max(tile_column),max(tile_row) FROM map WHERE (zoom_level = ?);',zoom_id).fetchone()
   i_x_min = int(bounds_minmax[0])
   i_y_min = int(bounds_minmax[1])
   i_x_max = int(bounds_minmax[2])
   i_y_max = int(bounds_minmax[3])
   # logger.info(_("MbTiles : retrieve_bounds: i_zoom %d min(%d,%d) ; max(%d,%d)")% (i_zoom,i_x_min,i_y_min,i_x_max,i_y_max))
   tile_bounds= mercator.TileLatLonBounds(i_x_min,i_y_min,i_zoom)
   if tile_bounds[0] < bounds_south:
    bounds_south=tile_bounds[0]
   if tile_bounds[1] < bounds_west:
    bounds_west=tile_bounds[1]
   tile_bounds= mercator.TileLatLonBounds(i_x_max,i_y_max,i_zoom)
   if tile_bounds[2] > bounds_north:
    bounds_north=tile_bounds[2]
   if tile_bounds[3] > bounds_east:
    bounds_east=tile_bounds[3]
  self.mbtiles_bounds="%f,%f,%f,%f"% (bounds_west,bounds_south,bounds_east,bounds_north)
  mbtiles_center_x=(bounds_east+bounds_west)/2
  mbtiles_center_y=(bounds_north+bounds_south)/2
  self.mbtiles_center="%f,%f,%s"%(mbtiles_center_x,mbtiles_center_y,min_zoom)
  self.mbtiles_minzoom=min_zoom
  self.mbtiles_maxzoom=max_zoom
  self.save_metadata()
  self.optimize_database()

 def retrieve_image(self,tz,tx,ty):
  if not self.s_y_type:
   self.s_y_type="tms"
  if not self.mbtiles_cursor:
   self.mbtiles_cursor = self.sqlite3_connection.cursor()
  tile_zxy = (str(tz), str(tx),str(ty))
  # SELECT tile_data FROM tiles WHERE ((zoom_level = 1) AND (tile_column = 1) AND (tile_row = 1))
  self.mbtiles_cursor.execute("SELECT tile_data FROM tiles WHERE ((zoom_level = ?) AND (tile_column = ?) AND (tile_row = ?))",tile_zxy)
  image_data = self.mbtiles_cursor.fetchone()
  if image_data is None:
   return  None
  return bytes(image_data[0])

 def retrieve_zoom_images(self,tz,tx,ty):
  if not self.s_y_type:
   self.s_y_type="tms"
  if not self.mbtiles_cursor:
   self.mbtiles_cursor = self.sqlite3_connection.cursor()
  s_tile_source="{0}-{1}-{2}.{3}".format(str(tz), str(tx),str(ty),self.s_y_type)
  tz=tz+1
  image_list = list() # empty list
  for y in range(2*ty,2*ty + 2):
   for x in range(2*tx, 2*tx + 2):
    s_tile_id="{0}-{1}-{2}.{3}".format(str(tz), str(x),str(y),self.s_y_type)
    s_file_id="{0}{1}-{2}-{3}.{4}".format(self.mbtiles_dir,str(tz), str(x),str(y),self.mbtiles_format)
    tile_id = (s_tile_id,)
    self.mbtiles_cursor.execute("SELECT tile_data FROM images WHERE tile_id = ?",tile_id)
    image_data = self.mbtiles_cursor.fetchone()
    if self.verbose:
     logger.info(_("MbTiles : retrieve_zoom_images: source[%s] : fetching[%s]") % (s_tile_source,s_tile_id))
    if image_data is None:
     # 1 / 8051 /media/gb_1500/maps/geo_tiff/rd_Berlin_Schmettau/18-140798-176204.jpg
     # [19-281597-352408.jpg] istilie_id '0-0-0.rgb'
     s_tile_id_orig=s_tile_id
     s_tile_id=self.count_tiles(tz,x,y,10)
     if self.verbose:
      logger.info(_("MbTiles : retrieve_zoom_images: fetching[%s] failed ; may be a rgb-image ; attempting [%s]") % (s_tile_id_orig,s_tile_id))
     tile_id = (s_tile_id,)
     self.mbtiles_cursor.execute("SELECT tile_data FROM images WHERE tile_id = ?",tile_id)
     image_data = self.mbtiles_cursor.fetchone()
     if image_data is None:
      # retireve an blank image fromdatabase, if does not exist, create it
      image_data = self.retrieve_blank_image(0,0,0)
    if image_data:
     output_file = open(s_file_id, 'wb')
     if not output_file.closed:
      output_file.write(image_data[0])
      output_file.close()
      image_list.append(s_file_id)
  return image_list

 def count_tiles(self,tz,tx,ty,i_parm):
  """
   i_parm values:
   ---------------------------
   10: will return string tile_id if thile does not exist
   11: count of tiles based on z,x,y values in request_url
     1: count of tiles based on z,x,y values
  """
  # even when empty, 0 will be returned
  if not self.s_y_type:
   self.s_y_type="tms"
  if not self.mbtiles_cursor:
   self.mbtiles_cursor = self.sqlite3_connection.cursor()
  if i_parm == 10:
   s_sql_command="SELECT tile_id FROM map WHERE zoom_level = ? AND tile_column = ? AND tile_row = ?"
   tile_id = (str(tz),str(tx),str(ty))
  if i_parm == 11:
   if self.i_request_url_count <=0:
     return 0
   s_tile_id="{0}-{1}-{2}.{3}".format(str(tz), str(tx),str(ty),self.s_y_type)
   s_sql_command="SELECT count(tile_id) FROM request_url WHERE tile_id = ?"
   tile_id = (s_tile_id,)
  if i_parm == 0:
   s_tile_id="{0}-{1}-{2}.{3}".format(str(tz), str(tx),str(ty),self.s_y_type)
   s_sql_command="SELECT count(tile_id) FROM map WHERE tile_id = ?"
   tile_id = (s_tile_id,)
  elif i_parm < 4:
   if i_parm == 1:
    # '%-%-%.%' : all tiles
    s_tile_id="%-%-%.%"
    s_sql_command="SELECT count(tile_id) FROM map WHERE tile_id LIKE ?"
   if i_parm == 2:
    # '15-%-%.%' all tiles of z ; Note: there may be .rgb entries
    s_tile_id="{0}".format(str(tz))
    s_sql_command="SELECT count(tile_id)  FROM map WHERE zoom_level = ?"
   if i_parm == 3:
    # '%-17600-%.%' all tiles of x ; Note: there may be .rgb entries
    s_tile_id="{0}".format(str(tx))
    s_sql_command="SELECT count(tile_id)  FROM map WHERE tile_column = ?"
   tile_id = (s_tile_id,)
  self.mbtiles_cursor.execute(s_sql_command,tile_id)
  i_count = self.mbtiles_cursor.fetchone()
  if i_count is None:
   if i_parm == 10:
    s_tile_id="{0}-{1}-{2}.{3}".format(str(0), str(0),str(0),"rgb")
    i_count = (s_tile_id,)
   else:
    i_count = (0,)
  return i_count[0]

 def get_tile_dirs(self,path):
  return [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]

 def mbtiles_from_disk(self,directory_path):
  if not os.path.exists(directory_path):
   logger.error(_("MbTiles : mbtiles_from_disk: directory does not exist : [%s] ") % directory_path)
   return
  if self.verbose:
   logger.info(_("MbTiles : mbtiles_from_disk: fetching[%s] ") % directory_path)
  image_format = ""
  min_zoom=22
  max_zoom=0
  bounds_west=180.0
  bounds_east=-180.0
  bounds_north=-85.05113
  bounds_south=85.05113
  mercator = GlobalMercator(self.tms_osm)
  if self.mbtiles_description == '':
   self.mbtiles_description=os.path.splitext(os.path.basename(directory_path))[0]
   if self.mbtiles_name == '':
    self.mbtiles_name=self.mbtiles_description.replace("."," ")
    self.mbtiles_name=self.mbtiles_name.replace("_"," ")
  for zoomDir in self.get_tile_dirs(directory_path):
   z = int(zoomDir)
   if z > max_zoom:
    max_zoom=z
   if z < min_zoom:
    min_zoom=z
   for rowDir in self.get_tile_dirs(os.path.join(directory_path, zoomDir)):
    x = int(rowDir)
    for current_file in os.listdir(os.path.join(directory_path, zoomDir, rowDir)):
     file_name, ext = current_file.split('.', 1)
     if image_format == "":
      image_format=ext
     f = open(os.path.join(directory_path, zoomDir, rowDir, current_file), 'rb')
     image_data = f.read()
     f.close()
     y = int(file_name)
     self.insert_image(z,x,y,image_data)
     tile_bounds= mercator.TileLatLonBounds(x,y,z)
     if tile_bounds[0] < bounds_south:
      bounds_south=tile_bounds[0]
     if tile_bounds[1] < bounds_west:
      bounds_west=tile_bounds[1]
     if tile_bounds[2] > bounds_north:
      bounds_north=tile_bounds[2]
     if tile_bounds[3] > bounds_east:
      bounds_east=tile_bounds[3]

  self.mbtiles_format=image_format
  self.mbtiles_bounds="%f,%f,%f,%f"% (bounds_west,bounds_south,bounds_east,bounds_north)
  mbtiles_center_x=(bounds_east+bounds_west)/2
  mbtiles_center_y=(bounds_north+bounds_south)/2
  self.mbtiles_center="%f,%f,%s"%(mbtiles_center_x,mbtiles_center_y,min_zoom)
  self.mbtiles_minzoom=min_zoom
  self.mbtiles_maxzoom=max_zoom
  self.save_metadata()
  self.optimize_database()
  if not  os.path.exists(os.path.join(directory_path, "tilemapresource.xml")):
   s_xml=self.mbtiles_create_tilemapresource();
   f = open(os.path.join(directory_path, "tilemapresource.xml"), 'w')
   f.write(s_xml)
   f.close()
  if self.verbose:
   logger.info(_("MbTiles : mbtiles_from_disk: [%s] - Habe fertig") % self.s_path_db)

 def mbtiles_to_disk(self,directory_path):
  if self.verbose:
   logger.info(_("MbTiles : mbtiles_to_disk: reading [%s]]") % self.s_path_db)
  if not os.path.exists(directory_path):
   os.mkdir("%s" % directory_path)
  count = self.mbtiles_cursor.execute('SELECT count(zoom_level) FROM tiles;').fetchone()[0]
  tiles = self.mbtiles_cursor.execute('SELECT zoom_level, tile_column, tile_row, tile_data FROM tiles;')
  t = tiles.fetchone()
  while t:
   z = t[0]
   x = t[1]
   y = t[2]
   tile_dir = os.path.join(directory_path, str(z), str(x))
   if not os.path.isdir(tile_dir):
    os.makedirs(tile_dir)
   tile = os.path.join(tile_dir,'%s.%s' % (y, self.mbtiles_format))
   f = open(tile, 'wb')
   f.write(t[3])
   f.close()
   t = tiles.fetchone()
  s_xml=self.mbtiles_create_tilemapresource();
  f = open(os.path.join(directory_path, "tilemapresource.xml"), 'w')
  f.write(s_xml)
  f.close()
  if self.verbose:
   logger.info(_("MbTiles : mbtiles_to_disk: created directory[%s] with %d tiles - Habe fertig") % directory_path,count)

 # -------------------------------------------------------------------------
 # http://docs.python.org/2/library/xml.etree.elementtree.html
 def mbtiles_read_tilemapresource(self,s_xml_file):
  tile_map = ElementTree.parse(s_xml_file).getroot()
  if tile_map:
   min_zoom=22
   max_zoom=0
   default_zoom=-1
   s_tile_row_type=""
   if tile_map.tag == 'TileMap':
    tile_srs = tile_map.find('SRS')
    if tile_srs:
     if tile_srs.text != 'EPSG:900913':
      # we support only mercator, assume invalid and return
      return
    tile_format = tile_map.find('TileFormat')
    if tile_format:
     tile_format_width = tile_format.get('width')
     tile_format_height = tile_format.get('height')
     if tile_format_width != "256" or tile_format_height != "256":
      # we support tiles of a size of 256, assume invalid and return
      return
     tile_format_mime = tile_format.get('mime-type')
     if tile_format_mime == "image/png" or tile_format_mime == "image/jpeg":
      self.mbtiles_format = tile_format.get('extension')
     else:
      # we support only jpg or png files, assume invalid and return
      return
     tile_format_tile_row_type = tile_format.get('tile_row_type')
     if tile_format_tile_row_type:
      # this is an unofficial parameter
      if tile_format_tile_row_type != "" and tile_format_tile_row_type != "tms":
       if tile_format_tile_row_type == "osm":
        # 'tms' is default - only to reset of osm
        s_tile_row_type=tile_format_tile_row_type
    tile_sets = tile_map.find('TileSets')
    if tile_sets:
     tile_sets_profile = tile_sets.get('profile')
     if tile_sets_profile != 'mercator':
      # we support only mercator, assume invalid and return
      return
     else:
      for zoom_level in tile_sets.getiterator('TileSet'):
       i_zoom = int(zoom_level.get("order"))
       if i_zoom > max_zoom:
        max_zoom=i_zoom
       if i_zoom < mim_zoom:
        min_zoom=i_zoom
      if min_zoom >= 0 and min_zoom <= 22 and max_zoom >= 0 and max_zoom <= 22:
       if min_zoom > max_zoom:
        i_zoom=max_zoom
        max_zoom=min_zoom
        min_zoom=i_zoom
       self.mbtiles_minzoom=min_zoom
       self.mbtiles_maxzoom=max_zoom
      else:
       return
    tile_title = tile_map.find('Title')
    if tile_title:
     self.mbtiles_name=tile_title.text
    tile_description = tile_map.find('Abstract')
    if tile_description:
     self.mbtiles_description=tile_description.text
    tile_origin = tile_map.find('Origin')
    if tile_origin:
     tile_origin_default_z = tile_origin.get('default_z')
     tile_origin_x = tile_origin.get('x')
     tile_origin_y = tile_origin.get('y')
     center_x=float(tile_origin_x)
     center_y=float(tile_origin_y)
     if tile_origin_default_z:
      # this is an unofficial parameter: Original will be interpeded as center of intrest with a default zoom
      default_zoom=int(tile_origin_default_z)
    tile_boundingbox = tile_map.find('BoundingBox')
    if tile_boundingbox:
     tile_minx = tile_boundingbox.get('minx')
     tile_miny = tile_boundingbox.get('miny')
     tile_maxx = tile_boundingbox.get('maxx')
     tile_maxy = tile_boundingbox.get('maxy')
     if default_zoom >= 0:
      self.mbtiles_bounds="%f,%f,%f,%f"% (float(tile_minx),float(tile_miny),float(tile_maxx),float(tile_maxy))
      self.mbtiles_center="%f,%f,%s"%(center_x,center_y,default_zoom)
     else:
      if tile_minx == tile_origin_x and tile_maxy == tile_origin_y:
       self.mbtiles_bounds="%f,%f,%f,%f"% (float(tile_minx),float(tile_miny),float(tile_maxx),float(tile_maxy))
       mbtiles_center_x=(float(tile_maxx)+float(tile_minx))/2
       mbtiles_center_y=(float(tile_maxy)+float(tile_miny))/2
       self.mbtiles_center="%f,%f,%s"%(mbtiles_center_x,mbtiles_center_y,min_zoom)

 # -------------------------------------------------------------------------
 # from Landez project
 # https://github.com/makinacorpus/landez
 def zoomlevels(self):
  rows = self.mbtiles_cursor.execute('SELECT DISTINCT(zoom_level) FROM tiles ORDER BY zoom_level')
  return [int(row[0]) for row in rows]

 # -------------------------------------------------------------------------
 # from Landez project
 # https://github.com/makinacorpus/landez
 def metadata(self,i_parm=0):
  rows = self.mbtiles_cursor.execute('SELECT name, value FROM metadata ORDER BY name  COLLATE NOCASE ASC')
  rows = [(row[0], row[1]) for row in rows]
  if i_parm == 1:
   return rows
  return dict(rows)

 # -------------------------------------------------------------------------
 # from Landez project
 # https://github.com/makinacorpus/landez
 def find_coverage(self, zoom):
  """
  Returns the bounding box (minx, miny, maxx, maxy) of an adjacent
  group of tiles at this zoom level.
  """
  # Find a group of adjacent available tiles at this zoom level
  rows = self.mbtiles_cursor.execute('''SELECT tile_column, tile_row FROM tiles  WHERE zoom_level=? ORDER BY tile_column, tile_row;''', (zoom,))
  tile = rows.fetchone()
  xmin, ymin = tile
  tile_prev = tile
  while tile and tile[0] - tile_prev[0] <= 1:
   # adjacent, go on
   tile_prev = tile
   tile = rows.fetchone()
   xmax, ymax = tile_prev
   # Transform (xmin, ymin) (xmax, ymax) to pixels
   tile_size = self.tilesize
   bottomleft = (xmin * tile_size, (ymax + 1) * tile_size)
   topright = ((xmax + 1) * tile_size, ymin * tile_size)
   # Convert center to (lon, lat)
   mercator = GlobalMercator(self.tms_osm,tile_size,[zoom])
  return mercator.unproject_pixels(bottomleft, zoom) + mercator.unproject_pixels(topright, zoom)

 # -------------------------------------------------------------------------
 # from Landez project
 # https://github.com/makinacorpus/landez
 def grid(self, z, x, y, callback=None):
  y_mercator = (2**int(z) - 1) - int(y)
  rows = self.mbtiles_cursor.execute('''SELECT grid FROM grids WHERE zoom_level=? AND tile_column=? AND tile_row=?;''', (z, x, y_mercator))
  grid = rows.fetchone()
  if not grid:
   raise ExtractionError(_("Could not extract grid %s from %s") % ((z, x, y),self.s_path_db))
  grid_json = json.loads(zlib.decompress(grid[0]))
  rows = self.mbtiles_cursor.execute('''SELECT key_name, key_json FROM grid_data  WHERE zoom_level=? AND tile_column=? AND tile_row=?;''', (z, x, y_mercator))
  # join up with the grid 'data' which is in pieces when stored in mbtiles file
  grid_json['data'] = {}
  grid_data = rows.fetchone()
  while grid_data:
   grid_json['data'][grid_data[0]] = json.loads(grid_data[1])
   grid_data = rows.fetchone()
   serialized = json.dumps(grid_json)
   if callback is not None:
    return '%s(%s);' % (callback, serialized)
  return serialized

 # -------------------------------------------------------------------------
 # from gdal2tiles project
 def mbtiles_create_tilemapresource(self):
  """
     Template for tilemapresource.xml. Returns filled string. Expected variables:
       title, north, south, east, west, isepsg4326, projection, publishurl,
       zoompixels, tilesize, tileformat, profile
       http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
  """
  args = {}
  args['title'] = self.mbtiles_name
  args['description'] = self.mbtiles_description
  args['south'] = self.bounds_south
  args['west'] = self.bounds_west
  args['north'] = self.bounds_north
  args['east'] = self.bounds_east
  args['center_x'] = self.center_x
  args['center_y'] = self.center_y
  args['default_z'] = self.default_zoom
  args['tilesize'] = self.tilesize
  args['tileformat'] = self.mbtiles_format
  s_mime=self.mbtiles_format
  if s_mime == "jpg":
   s_mime="jpeg"
  args['mime'] = "{0}/{1}".format("image",s_mime)
  args['publishurl'] = ""
  args['profile'] = "mercator"
  args['srs'] = "EPSG:900913"
  args['tile_row_type'] =self.s_y_type
  args['default_zoom'] =self.s_y_type

  s_xml = """<?xml version="1.0" encoding="utf-8"?>
<TileMap version="1.0.0" tilemapservice="http://tms.osgeo.org/1.0.0">
 <Title>%(title)s</Title>
 <Abstract>%(description)s</Abstract>
 <SRS>%(srs)s</SRS>
 <BoundingBox minx="%(west).14f" miny="%(south).14f" maxx="%(east).14f" maxy="%(north).14f"/>
 <Origin x="%(center_x).14f" y="%(center_y).14f" default_z="%(default_z)d"/>
 <TileFormat width="%(tilesize)d" height="%(tilesize)d" mime-type="%(mime)s" extension="%(tileformat)s" tile_row_type="%(tile_row_type)s"/>
 <TileSets profile="%(profile)s">
""" % args
  for z in range(int(self.mbtiles_minzoom), int(self.mbtiles_maxzoom)+1):
   s_xml += """  <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, 156543.0339/2**z, z)
  s_xml += """ </TileSets>
</TileMap>
 """
  return s_xml

# from Landez project
# https://github.com/makinacorpus/landez
class EmptyCoverageError(Exception):
 """ Raised when coverage (tiles list) is empty """
 pass

# from Landez project
# https://github.com/makinacorpus/landez
class DownloadError(Exception):
 """ Raised when download at tiles URL fails DOWNLOAD_RETRIES times """
 pass

# from Landez project
# https://github.com/makinacorpus/landez
class ExtractionError(Exception):
 """ Raised when extraction of tiles from specified MBTiles has failed """
 pass

# from Landez project
# https://github.com/makinacorpus/landez
class InvalidFormatError(Exception):
 """ Raised when reading of MBTiles content has failed """
 pass
# from Landez project
# https://github.com/makinacorpus/landez
class TileSource(object):
 def __init__(self, tilesize=None):
  if tilesize is None:
   tilesize = 256
  self.tilesize = tilesize
  self.basename = ''
  self.metadata_input=None
  self.tms_osm=False
  self.s_y_type = 'tms'
  self.mbtiles_format='jpg'
  self.mbtiles_verbose=False
  self.mbtiles_bounds="-180.00000,-85.05113,180.00000,85.05113"
  self.mbtiles_minzoom="0"
  self.mbtiles_maxzoom="22"

  def tile(self, z, x, y):
   raise NotImplementedError

  def metadata(self):
   return dict()

# from Landez project
# https://github.com/makinacorpus/landez
class MBTilesReader(TileSource):
 def __init__(self, mbtiles_input, tilesize=None):
  super(MBTilesReader, self).__init__(tilesize)
  self.mbtiles_input = mbtiles_input.strip()
  self.basename = os.path.basename(self.mbtiles_input)
  self.mbtiles_input_dir=os.path.dirname(self.mbtiles_input)+ '/'
  self.mbtiles_db_input=MbTiles()
  self.mbtiles_db_input.open_db(self.mbtiles_input,self.mbtiles_input_dir,self.mbtiles_format,self.s_y_type,self.mbtiles_verbose)
  self.metadata_input=self.mbtiles_db_input.fetch_metadata()
  self.tms_osm=self.mbtiles_db_input.tms_osm
  self.mbtiles_format=self.mbtiles_db_input.mbtiles_format
  self.s_y_type=self.mbtiles_db_input.s_y_type
  self.mbtiles_bounds=self.mbtiles_db_input.mbtiles_bounds
  self.mbtiles_minzoom=self.mbtiles_db_input.mbtiles_minzoom
  self.mbtiles_maxzoom=self.mbtiles_db_input.mbtiles_maxzoom

 def metadata(self,i_parm=0):
  return self.mbtiles_db_input.metadata(i_parm)

 def zoomlevels(self):
  return self.mbtiles_db_input.zoomlevels()

 def tile(self, z, x, y):
  logger.debug(_("MBTilesReader.Extract tile %s") % ((z, x, y),))
  # y_mercator = (2**int(z) - 1) - int(y)
  return self.mbtiles_db_input.retrieve_image(z,x,y)

 def grid(self, z, x, y, callback=None):
  return self.mbtiles_db_input.grid(z,x,y, callback)

 def find_coverage(self, zoom):
  """
  Returns the bounding box (minx, miny, maxx, maxy) of an adjacent
  group of tiles at this zoom level.
  """
  # Find a group of adjacent available tiles at this zoom level
  return self.mbtiles_db_input.find_coverage(zoom)

# from Landez project
# https://github.com/makinacorpus/landez
class TilesManager(object):
 def __init__(self, **kwargs):
  """
  Manipulates tiles in general. Gives ability to list required tiles on a
  bounding box, download them, render them, extract them from other mbtiles...
  Keyword arguments:
  cache -- use a local cache to share tiles between runs (default True)

  tiles_dir -- Local folder containing existing tiles if cache is
  True, or where temporary tiles will be written otherwise
    (default DEFAULT_TMP_DIR)
  tiles_url -- remote URL to download tiles (*default DEFAULT_TILES_URL*)
  tiles_headers -- HTTP headers to send (*default empty*)
  stylefile -- mapnik stylesheet file (*to render tiles locally*)

  mbtiles_file -- A MBTiles file providing tiles (*to extract its tiles*)

  wms_server -- A WMS server url (*to request tiles*)
  wms_layers -- The list of layers to be requested
  wms_options -- WMS parameters to be requested (see ``landez.reader.WMSReader``)
  tile_size -- default tile size (default DEFAULT_TILE_SIZE)
  tile_format -- default tile format (default DEFAULT_TILE_FORMAT)
  """
  self.tile_size = kwargs.get('tile_size', DEFAULT_TILE_SIZE)
  self.tile_format = kwargs.get('tile_format', DEFAULT_TILE_FORMAT)
  # Tiles Download
  self.tiles_url = kwargs.get('tiles_url', DEFAULT_TILES_URL)
  self.tiles_subdomains = kwargs.get('tiles_subdomains', DEFAULT_TILES_SUBDOMAINS)
  self.tiles_headers = kwargs.get('tiles_headers')
  # Tiles rendering
  self.stylefile = kwargs.get('stylefile')
  # Grids rendering
  self.grid_fields = kwargs.get('grid_fields', [])
  self.grid_layer = kwargs.get('grid_layer', 0)
  # MBTiles reading
  self.mbtiles_input = kwargs.get('mbtiles_input')
   # request_url
  self.request_url = kwargs.get('request_url')
  # WMS requesting
  self.wms_server = kwargs.get('wms_server')
  self.wms_layers = kwargs.get('wms_layers', [])
  self.wms_options = kwargs.get('wms_options', {})
  if self.mbtiles_input:
   self.reader = MBTilesReader(self.mbtiles_input, self.tile_size)
  elif self.wms_server:
   assert self.wms_layers, _("Requires at least one layer (see ``wms_layers`` parameter)")
   self.reader = WMSReader(self.wms_server, self.wms_layers,self.tile_size, **self.wms_options)
   if 'format' in self.wms_options:
    self.tile_format = self.wms_options['format']
    logger.info(_("Tile format set to %s") % self.tile_format)
   self.tiles_url=self.reader.request_url;
  elif self.stylefile:
   self.reader = MapnikRenderer(self.stylefile, self.tile_size)
  else:
   mimetype, encoding = mimetypes.guess_type(self.tiles_url)
   if mimetype and mimetype != self.tile_format:
    self.tile_format = mimetype
    logger.info(_("Tile format set to %s") % self.tile_format)
   self.reader = TileDownloader(self.tiles_url, headers=self.tiles_headers,subdomains=self.tiles_subdomains, tilesize=self.tile_size)
   self.reader.mbtiles_format=self.tile_format
  # Tile files extensions
  self._tile_extension = mimetypes.guess_extension(self.tile_format, strict=False)
  assert self._tile_extension, _("Unknown format %s") % self.tile_format
  if self._tile_extension == '.jpe':
   self._tile_extension = '.jpeg'
  # Cache
  tiles_dir = kwargs.get('tiles_dir', DEFAULT_TMP_DIR)
  if kwargs.get('cache', True):
   self.cache = Disk(self.reader.basename, tiles_dir, extension=self._tile_extension)
  else:
   self.cache = Dummy(extension=self._tile_extension)
  # Overlays
  self._layers = []
  # Filters
  self._filters = []
  # Number of tiles rendered/downloaded here
  self.rendered = 0

 def tileslist(self, bbox, zoomlevels, tms_osm=False):
  """
  Build the tiles list within the bottom-left/top-right bounding
  box (minx, miny, maxx, maxy) at the specified zoom levels.
  Return a list of tuples (z,x,y)
  """
  mercator = GlobalMercator(tms_osm,self.tile_size,zoomlevels)
  return mercator.tileslist(bbox)

 def add_layer(self, tilemanager, opacity=1.0):
  """
  Add a layer to be blended (alpha-composite) on top of the tile.
  tilemanager -- a `TileManager` instance
  opacity -- transparency factor for compositing
  """
  assert has_pil, _("Cannot blend layers without python PIL")
  assert self.tile_size == tilemanager.tile_size, _("Cannot blend layers whose tile size differs")
  assert 0 <= opacity <= 1, _("Opacity should be between 0.0 (transparent) and 1.0 (opaque)")
  self.cache.basename += '%s%.1f' % (tilemanager.cache.basename, opacity)
  self._layers.append((tilemanager, opacity))

 def add_filter(self, filter_):
  """ Add an image filter for post-processing """
  assert has_pil, _("Cannot add filters without python PIL")
  self.cache.basename += filter_.basename
  self._filters.append(filter_)

 def get_metadata(self,i_parm=0):
  metadata_list = self.reader.metadata(i_parm)
  return metadata_list

 def tile(self, (z, x, y)):
  """
  Return the tile (binary) content of the tile and seed the cache.
  """
  output = self.cache.read((z, x, y))
  if output is None:
   # logger.info(_("TilesManager.tile calling sources.tile: ") )
   pass
  output = self.reader.tile(z, x, y)
  if output is None:
   return None
  # Blend layers
  if len(self._layers) > 0:
   logger.debug(_("Will blend %s layer(s)") % len(self._layers))
   output = self._blend_layers(output, (z, x, y))
   # Apply filters
   for f in self._filters:
    image = f.process(self._tile_image(output))
    output = self._image_tile(image)
    # Save result to cache
    self.cache.save(output, (z, x, y))
    self.rendered += 1
  return output

 def grid(self, (z, x, y)):
  """ Return the UTFGrid content """
  # sources.py -> MapnikRenderer -> grid
  content = self.reader.grid(z, x, y, self.grid_fields, self.grid_layer)
  return content

 def _blend_layers(self, imagecontent, (z, x, y)):
  """
  Merge tiles of all layers into the specified tile path
  """
  result = self._tile_image(imagecontent)
  # Paste each layer
  for (layer, opacity) in self._layers:
   try:
    # Prepare tile of overlay, if available
    overlay = self._tile_image(layer.tile((z, x, y)))
   except (DownloadError, ExtractionError), e:
    logger.warn(e)
    continue
   # Extract alpha mask
   overlay = overlay.convert("RGBA")
   r, g, b, a = overlay.split()
   overlay = Image.merge("RGB", (r, g, b))
   a = ImageEnhance.Brightness(a).enhance(opacity)
   overlay.putalpha(a)
   mask = Image.merge("L", (a,))
   result.paste(overlay, (0, 0), mask)
   # Read result
  return self._image_tile(result)

 def _tile_image(self, data):
  """
  Tile binary content as PIL Image.
  """
  image = Image.open(BytesIO(data))
  return image.convert('RGBA')

 def _image_tile(self, image):
  out_image = BytesIO()
  image.save(out_image, self._tile_extension[1:])
  return out_image.getvalue()

# from Landez project
# https://github.com/makinacorpus/landez
class MBTilesBuilder(TilesManager):
 def __init__(self, **kwargs):
  """
  A MBTiles builder for a list of bounding boxes and zoom levels.
  mbtiles_output -- output MBTiles file (default DEFAULT_MBTILES_OUTPUT)
  tmp_dir -- temporary folder for gathering tiles (default DEFAULT_TMP_DIR/mbtiles_output)
  """
  super(MBTilesBuilder, self).__init__(**kwargs)
  self.mbtiles_output = kwargs.get('mbtiles_output', DEFAULT_MBTILES_OUTPUT)
  # Gather tiles for mbutil
  basename, ext = os.path.splitext(os.path.basename(self.mbtiles_output))
  self.tmp_dir = kwargs.get('tmp_dir', DEFAULT_TMP_DIR)
  self.tmp_dir = os.path.join(self.tmp_dir, basename)
  self.tile_format = kwargs.get('tile_format', DEFAULT_TILE_FORMAT)
  # Number of tiles in total
  self.nbtiles = 0
  self._bboxes = []
  self._metadata = []

 def add_coverage(self, bbox, zoomlevels):
  """
  Add a coverage to be included in the resulting mbtiles file.
  """
  self._bboxes.append((bbox, zoomlevels))

 def add_metadata(self, metdatadata_list):
  """
  Add metadata to be included in the resulting mbtiles file.
  """
  self._metadata.append((metdatadata_list, ))

 @property
 def zoomlevels(self):
  """
  Return the list of covered zoom levels
  """
  return self._bboxes[0][1]  #TODO: merge all coverages

 @property
 def bounds(self):
  """
  Return the bounding box of covered areas
  """
  return self._bboxes[0][0]  #TODO: merge all coverages

 def run(self, add=False, i_parm=0):
  """
  Build a MBTile file.
  add -- add if MBTiles file already exists.
  i_parm -- 1 : show info only
  """
  if os.path.exists(self.mbtiles_output):
   if add:
    logger.warn(_("%s already exists. Will be added to.") % self.mbtiles_output)
    # os.remove(self.mbtiles_output)
   else:
    # Already built, do not do anything.
    logger.info(_("%s already exists. Nothing to do.") % self.mbtiles_output)
    return
  if len(self._bboxes) == 0:
   bbox = map(float, self.reader.mbtiles_bounds.split(','))
   zoomlevels = range(int(self.reader.mbtiles_minzoom), int(self.reader.mbtiles_maxzoom)+1)
   self.add_coverage(bbox=bbox, zoomlevels=zoomlevels)
  if self.request_url != "":
    if self.request_url.find('fill') == -1 and self.request_url.find('load') == -1 and self.request_url.find('replace') == -1:
     self.request_url=""
  if self.request_url == "":
   # Compute list of tiles
   tileslist = set()
   for bbox, levels in self._bboxes:
    logger.debug(_("MBTilesBuilder.run: Compute list of tiles for bbox %s on zooms %s.") % (bbox, levels))
    bboxlist,tile_bounds = self.tileslist(bbox, levels,self.reader.tms_osm)
    logger.debug(_("Add %s tiles.") % len(bboxlist))
    tileslist = tileslist.union(bboxlist)
    logger.debug(_("MBTilesBuilder.run: %s tiles in total.") % len(tileslist))
   self.nbtiles = len(tileslist)
   if not self.nbtiles:
    raise EmptyCoverageError(_("No tiles are covered by bounding boxes : %s") % self._bboxes)
   if i_parm == 1:
    logger.info(_("-I-> Computed list of [%d] tiles for bbox %s on zooms %s.") % (self.nbtiles,bbox, levels))
    logger.info(_("-I-> MBTilesBuilder.run: i_parm[%d] ; set i_parm=0 to run ; will now exit.") % i_parm)
    return
  self.mbtiles_output=self.mbtiles_output.strip()
  self.mbtiles_output_dir=os.path.dirname(self.mbtiles_output)+ '/'
  self.mbtiles_db_output=MbTiles()
  self.mbtiles_db_output.open_db(self.mbtiles_output,self.mbtiles_output_dir,self.reader.mbtiles_format,self.reader.s_y_type,self.reader.mbtiles_verbose)
  if self.reader.metadata_input:
   self.mbtiles_db_output.insert_metadata(self.reader.metadata_input)
  if self.request_url == "":
   # Go through whole list of tiles and read from input_db and store in output_db
   self.rendered = 0
   # something is 'unsorting' this
   for (z, x, y) in sorted(tileslist,key=operator.itemgetter(0,1,2)):
    image_data=self.tile((z,x,y))
    if not image_data is None:
     self.mbtiles_db_output.insert_image(z,x,y,image_data)
   # calculate the min/max zoom_levels and bounds
   self.mbtiles_db_output.retrieve_bounds()
   logger.debug(_("MBTilesBuilder.run: %s tiles were missing.") % self.rendered)
   # Package it!
   logger.info(_("Build MBTiles output file '%s'.") % self.mbtiles_output)
   for metadata_list in self._metadata:
    logger.debug(_("MBTilesBuilder.run: adding metadata %s .") % (metadata_list))
    self.mbtiles_db_output.insert_metadata(metadata_list[0])
  else:
   if self.request_url.find('fill') != -1:
    i_parm=0
    if self.request_url.find('replace') != -1:
     i_parm=1
    for bbox, levels in self._bboxes:
     self.mbtiles_db_output.fill_request_url(bbox,levels,self.tiles_url,i_parm)
   if self.request_url.find('load') != -1:
    i_parm=0
    self.mbtiles_db_output.load_request_url(i_parm)
  self.mbtiles_db_output.close_db()

# from Landez project
# https://github.com/makinacorpus/landez
# http://effbot.org/imagingbook/pil-index.htm#appendixes
class ImageExporter(TilesManager):
 def __init__(self, **kwargs):
  """
  Arrange the tiles and join them together to build a single big image.
  Sample Image: 1861_Mercator_Europe zoom_level 8
  - given       bbox: -8.0,36.0,80.0,77.0
  - calulated bbox: -8.4375, 35.46066995149529, 80.15625, 77.157162522661
  Image Size is 16128, 15872 ; Pixel Size = (0.005493164062500,-0.002627047163002)
  Origin = (-8.437500000000000,77.157162522660997)
  PNG:    370   MB
  TIF : 1.000   MB [NONE]
  TIF : 1.000   MB [PACKBITS]  - took 1 minute to create
  TIF :   404.7 MB [DEFLATE,PREDICTOR=2,ZLEVEL=8] - took 13 minutes to create
  TIF :   404.7 MB [DEFLATE,PREDICTOR=2,ZLEVEL=9] - took 13 minutes to create
  TIF :   472.5 MB [LZW,PREDICTOR=2] - took 1 minute to create [default if nothing is supplid]
  TIF :   735.8 MB [LZW,PREDICTOR=1] - took 1 minute to create
  """
  super(ImageExporter, self).__init__(**kwargs)
  # COMPRESS=PACKBITS
  # PIL: TIFF uncompressed, or Packbits, LZW, or JPEG compressed images. In the current version, PIL always writes uncompressed TIFF files
  # http://linfiniti.com/2011/05/gdal-efficiency-of-various-compression-algorithms/
  # predictor for 'DEFLATE' or 'LZW' : 1 or 2
  i_tiff_compression_predictor=2
  # zlevel for 'DEFLATE'  : 1 to 9
  i_tiff_compression_zlevel=8
  self.jpg_quality=75
  self.tiff_compression=[]
  self._metadata = []
  if self.reader.metadata_input:
   self.metadata_input=self.reader.metadata_input
  self.tiff_compress = kwargs.get('tiff_compression', "LZW")
  self.tiff_compress =self.tiff_compress.upper()
  self.jpg_quality = kwargs.get('jpg_quality', self.jpg_quality)
  if self.jpg_quality < 1 or self.jpg_quality > 95:
   self.jpg_quality=75
  i_tiff_compression_predictor = kwargs.get('tiff_predictor', i_tiff_compression_predictor)
  if i_tiff_compression_predictor < 1 or i_tiff_compression_predictor > 2:
   i_tiff_compression_predictor=2
  i_tiff_compression_zlevel = kwargs.get('tiff_zlevel', i_tiff_compression_zlevel)
  if i_tiff_compression_zlevel < 1 or i_tiff_compression_zlevel > 9:
   i_tiff_compression_predictor=8
  if self.tiff_compress == "PACKBITS" :
   self.tiff_compression.append('COMPRESS=PACKBITS')
  elif self.tiff_compress == "DEFLATE":
   self.tiff_compression.append('COMPRESS=%s' % 'DEFLATE')
   self.tiff_compression.append('PREDICTOR=%d' % i_tiff_compression_predictor)
   self.tiff_compression.append('ZLEVEL=%d' % i_tiff_compression_zlevel)
  elif self.tiff_compress == "LZW":
   self.tiff_compression.append('COMPRESS=%s' % 'LZW')
   self.tiff_compression.append('PREDICTOR=%d' % i_tiff_compression_predictor)
  elif self.tiff_compress == "NONE":
   self.tiff_compression.append('COMPRESS=NONE')

 def add_metadata(self, metdatadata_list):
  """
  Add metadata to be included in the resulting mbtiles file.
  """
  self._metadata.append((metdatadata_list, ))

 def grid_tiles(self, bbox, zoomlevel):
  """
  Return a grid of (x, y) tuples representing the juxtaposition
  of tiles on the specified ``bbox`` at the specified ``zoomlevel``.
  """
  tiles,tile_bounds = self.tileslist(bbox, [zoomlevel],self.reader.tms_osm)
  grid = {}
  # for (z, x, y) in sorted(tiles,key=operator.itemgetter(0,1,2),reverse=True):
  for (z, x, y) in tiles:
   if not grid.get(y):
    grid[y] = []
   grid[y].append(x)
  sortedgrid = []
  for y in sorted(grid.keys(),reverse=not self.reader.tms_osm):
   sortedgrid.append([(x, y) for x in sorted(grid[y])])
  return sortedgrid,tile_bounds

 def export_image(self, bbox, zoomlevel, imagepath):
  """
  Writes to ``imagepath`` the tiles for the specified bounding box and zoomlevel.
  """
  assert has_pil, _("Cannot export image without python PIL")
  grid,tile_bounds = self.grid_tiles(bbox, zoomlevel)
  width = len(grid[0])
  height = len(grid)
  widthpix = width * self.tile_size
  heightpix = height * self.tile_size
  result = Image.new("RGBA", (widthpix, heightpix))
  offset = (0, 0)
  for i, row in enumerate(grid):
   for j, (x, y) in enumerate(row):
     offset = (j * self.tile_size, i * self.tile_size)
     img = self._tile_image(self.tile((zoomlevel, x, y)))
     result.paste(img, offset)
  if imagepath.endswith(".tif") or imagepath.endswith(".tiff"):
   # http://effbot.org/imagingbook/pil-index.htm#appendixes
   # Packbits, LZW, or JPEG
   # In the current version, PIL always writes uncompressed TIFF files.
   image_gdal="gdal_input.tif"
   result.save(image_gdal, format="TIFF", compression="JPEG")
   self.geotif_image(tile_bounds,(widthpix,heightpix),imagepath,image_gdal)
  else:
   if imagepath.endswith(".jpg") or imagepath.endswith(".jpeg"):
    # IOError: encoder error -2 when writing image file
    result.save(imagepath, format="JPEG", quality=int(self.jpg_quality), optimize=True, progressive=False)
   elif  imagepath.endswith(".png") :
    result.save(imagepath, format="PNG",optimize=True)
   else:
    result.save(imagepath)
   logger.info(_("-I-> export_image: Save resulting image to '%s' - bounds[%s]") % (imagepath,tile_bounds))

 def geotif_image(self, tile_bounds, image_bounds, imagepath,image_gdal):
  """
  Writes to ``imagepath`` the tiles for the specified bounding box and zoomlevel.
  """
  i_srid=3857
  s_srid="WGS 84 / Pseudo-Mercator"
  # i_srid=3395
  # s_srid="WGS 84 / World Mercator"
  # 4326 Wsg84
  # Upper Left  (  -8.4375000,  77.1571625) (  8d26'15.00"W, 77d 9'25.79"N)
  # Lower Left  (  -8.4375000,  35.4606700) (  8d26'15.00"W, 35d27'38.41"N)
  # Upper Right (  80.1562500,  77.1571625) ( 80d 9'22.50"E, 77d 9'25.79"N)
  # Lower Right (  80.1562500,  35.4606700) ( 80d 9'22.50"E, 35d27'38.41"N)
  # Center      (  35.8593750,  56.3089162) ( 35d51'33.75"E, 56d18'32.10"N)
  # 3857 'WGS 84 / Pseudo-Mercator'
  # Upper Left  ( -939258.204,13932330.020) (  8d26'15.00"W, 77d 9'25.79"N)
  # Lower Left  ( -939258.204, 4226661.916) (  8d26'15.00"W, 35d27'38.41"N)
  # Upper Right ( 8922952.934,13932330.020) ( 80d 9'22.50"E, 77d 9'25.79"N)
  # Lower Right ( 8922952.934, 4226661.916) ( 80d 9'22.50"E, 35d27'38.41"N)
  # Center      ( 3991847.365, 9079495.968) ( 35d51'33.75"E, 62d54'54.84"N)
  # 3395 'WGS 84 / World Mercator'
  # Upper Left  ( -939258.204,13932330.020) (  8d26'15.00"W, 77d14'24.81"N)
  # Lower Left  ( -939258.204, 4226661.916) (  8d26'15.00"W, 35d38'33.56"N)
  # Upper Right ( 8922952.934,13932330.020) ( 80d 9'22.50"E, 77d14'24.81"N)
  # Lower Right ( 8922952.934, 4226661.916) ( 80d 9'22.50"E, 35d38'33.56"N)
  # Center      ( 3991847.365, 9079495.968) ( 35d51'33.75"E, 63d 4'14.87"N)
  bounds_west,bounds_south,bounds_east,bounds_north=tile_bounds
  bounds_wsg84="bounds_wsg84: %f,%f,%f,%f"% (bounds_west,bounds_south,bounds_east,bounds_north)
  mercator = GlobalMercator()
  tile_bounds=mercator.BoundsToMeters(tile_bounds)
  mbtiles_name="";
  mbtiles_description=""
  s_TIFFTAG_DOCUMENTNAME=""
  s_TIFFTAG_IMAGEDESCRIPTION=""
  s_TIFFTAG_SOFTWARE=""
  s_TIFFTAG_DATETIME=""
  s_TIFFTAG_ARTIST=""
  s_TIFFTAG_HOSTCOMPUTER=""
  s_TIFFTAG_COPYRIGHT=""
  if self.metadata_input:
   metadata=dict(self.metadata_input)
   mbtiles_name=metadata.get('name','')
   mbtiles_description=metadata.get('description','')
  if self._metadata:
   for metadata_list in self._metadata:
    metadata=dict(metadata_list[0])
    mbtiles_name=metadata.get('name',mbtiles_name)
    mbtiles_description=metadata.get('description',mbtiles_description)
    s_TIFFTAG_DOCUMENTNAME=metadata.get('TIFFTAG_DOCUMENTNAME',mbtiles_name)
    s_TIFFTAG_IMAGEDESCRIPTION=metadata.get('TIFFTAG_IMAGEDESCRIPTION',mbtiles_description)
    s_TIFFTAG_SOFTWARE=metadata.get('TIFFTAG_SOFTWARE','')
    s_TIFFTAG_DATETIME=metadata.get('TIFFTAG_DATETIME','')
    s_TIFFTAG_ARTIST=metadata.get('TIFFTAG_ARTIST','')
    s_TIFFTAG_HOSTCOMPUTER=metadata.get('TIFFTAG_HOSTCOMPUTER','')
    s_TIFFTAG_COPYRIGHT=metadata.get('TIFFTAG_COPYRIGHT','')
  if s_TIFFTAG_DOCUMENTNAME == "":
   s_TIFFTAG_DOCUMENTNAME=mbtiles_name
  if s_TIFFTAG_IMAGEDESCRIPTION == "":
   s_TIFFTAG_IMAGEDESCRIPTION=mbtiles_description
  tiff_metadata=[]
  if s_TIFFTAG_DOCUMENTNAME != "":
   tiff_metadata.append(('TIFFTAG_DOCUMENTNAME',s_TIFFTAG_DOCUMENTNAME))
  if s_TIFFTAG_IMAGEDESCRIPTION != "":
   tiff_metadata.append(('TIFFTAG_IMAGEDESCRIPTION',s_TIFFTAG_IMAGEDESCRIPTION))
  if s_TIFFTAG_SOFTWARE != "":
   tiff_metadata.append(('TIFFTAG_SOFTWARE',s_TIFFTAG_SOFTWARE))
  else:
   tiff_metadata.append(('TIFFTAG_SOFTWARE',bounds_wsg84))
  if s_TIFFTAG_DATETIME != "":
   tiff_metadata.append(('TIFFTAG_DATETIME',s_TIFFTAG_DATETIME))
  if s_TIFFTAG_ARTIST != "":
   tiff_metadata.append(('TIFFTAG_ARTIST',s_TIFFTAG_ARTIST))
  if s_TIFFTAG_HOSTCOMPUTER != "":
   tiff_metadata.append(('TIFFTAG_HOSTCOMPUTER',s_TIFFTAG_HOSTCOMPUTER))
  if s_TIFFTAG_COPYRIGHT != "":
   tiff_metadata.append(('TIFFTAG_COPYRIGHT',s_TIFFTAG_COPYRIGHT))
  # this assumes the projection is Geographic lat/lon WGS 84
  xmin,ymin,xmax,ymax=tile_bounds
  image_width,image_height=image_bounds
  # Upper Left  (   20800.000,   22000.000)
  # Lower Right (   24000.000,   19600.000)
  # Size is 15118, 11339
  # (24000-20800)/15118 = 3200 = 0,21166821 [xres]
  # (19600-22000)/11339 = 2400 =  −0,211658876 [yres]
  # geo_transform = (20800.0, 0.2116682100807, 0.0, 22000.0, 0.0, -0.21165887644413)
  geo_transform = [xmin, (xmax-xmin)/image_width, 0, ymax, 0, (ymin-ymax)/image_height ]
  spatial_projection = osr.SpatialReference()
  spatial_projection.ImportFromEPSG(i_srid)
  logger.info(_("-I-> geotif_image: Saving as GeoTiff - image[%s] compression[%s]") % (imagepath,self.tiff_compression))
  image_dataset = gdal.Open(image_gdal, gdal.GA_Update )
  image_dataset.SetProjection(spatial_projection.ExportToWkt())
  image_dataset.SetGeoTransform(geo_transform)
  driver = gdal.GetDriverByName("GTiff")
  output_dataset = driver.CreateCopy(imagepath,image_dataset, 0, self.tiff_compression )
  if tiff_metadata:
   logger.info(_("-I-> geotif_image: tiff_metadata[%s]") % tiff_metadata)
   output_dataset.SetMetadata(dict(tiff_metadata))
  # Once we're done, close properly the dataset
  output_dataset = None
  image_dataset = None
  os.remove(image_gdal)
  logger.info(_("-I-> geotif_image: Saved resulting image to '%s' as GeoTiff- bounds[%s]") % (imagepath,tile_bounds))


# from Landez project
# https://github.com/makinacorpus/landez
class TileDownloader(TileSource):
 def __init__(self, url, headers=None, subdomains=None, tilesize=None):
  super(TileDownloader, self).__init__(tilesize)
  self.tiles_url = url
  self.tiles_subdomains = subdomains or ['a', 'b', 'c']
  parsed = urlparse(self.tiles_url)
  self.basename = parsed.netloc
  self.headers = headers or {}

 def tile(self, z, x, y_tms):
  """
  Download the specified tile from `tiles_url`
  """
  logger.debug(_("Download tile %s") % ((z, x, y_tms),))
  # Render each keyword in URL ({s}, {x}, {y}, {z}, {size} ... )
  size = self.tilesize
  s = self.tiles_subdomains[(x + y_tms) % len(self.tiles_subdomains)];
  y_osm = (2**int(z) - 1) - int(y_tms)
  try:
   url = self.tiles_url.format(**locals())
  except KeyError, e:
   raise DownloadError(_("Unknown keyword %s in URL") % e)
  logger.debug(_("Retrieve tile at %s") % url)
  r = DOWNLOAD_RETRIES
  sleeptime = 1
  while r > 0:
   try:
    request = urllib2.Request(url)
    for header, value in self.headers.items():
     request.add_header(header, value)
    stream = urllib2.urlopen(request)
    assert stream.getcode() == 200
    return stream.read()
   except (AssertionError, IOError), e:
    logger.debug(_("Download error, retry (%s left). (%s)") % (r, e))
    r -= 1
    time.sleep(sleeptime)
    # progressivly sleep longer to wait for this tile
    if (sleeptime <= 10) and (r % 2 == 0):
     sleeptime += 1  # increase wait
  raise DownloadError(_("Cannot download URL %s") % url)

# from Landez project
# https://github.com/makinacorpus/landez
class WMSReader(TileSource):
 def __init__(self, url, layers, tilesize=None, **kwargs):
  super(WMSReader, self).__init__(tilesize)
  self.basename = '-'.join(layers)
  self.url = url
  self.wmsParams = dict(
   service='WMS',
   request='GetMap',
   version='1.1.1',
   styles='',
   format=DEFAULT_TILE_FORMAT,
   transparent=False,
   layers=','.join(layers),
   width=self.tilesize,
   height=self.tilesize,
  )
  self.wmsParams.update(**kwargs)
  projectionKey = 'srs'
  if parse_version(self.wmsParams['version']) >= parse_version('1.3'):
   projectionKey = 'crs'
  self.wmsParams[projectionKey] = GlobalMercator.WSG84
  encodedparams = urllib.urlencode(self.wmsParams)
  # this can be sent to the mbtiles request_url logic
  self.request_url = "%s?%s" % (self.url.encode("utf8"), encodedparams)

 def tile(self, z, x, y_tms):
  logger.debug(_("Request WMS tile %s") % ((z, x, y_tms),))
  y_osm = (2**int(z) - 1) - int(y_tms)
  mercator = GlobalMercator(True,self.tilesize,[z])
  bbox = mercator.tile_bbox((z, x, y_osm))
  # bbox = mercator.project(bbox[:2]) + mercator.project(bbox[2:])
  bbox = ','.join(map(str, bbox))
  # Build WMS request URL
  # encodedparams = urllib.urlencode(self.wmsParams)
  # url = "%s?%s" % (self.url, encodedparams)
  url = self.request_url
  url += "&bbox=%s" % bbox   # commas are not encoded
  try:
   logger.debug(_("Download '%s'") % url)
   f = urllib2.urlopen(url)
   header = f.info().typeheader
   assert header == self.wmsParams['format'], "Invalid WMS response type : %s" % header
   return f.read()
  except (AssertionError, IOError),e:
   # BBOX parameter's minimum Y is greater than the maximum Y
   # srs is not permitted: EPSG:3857
   logger.error(_("WMS request URL: Error %s:") % e.args[0])
   raise ExtractionError

# from Landez project
# https://github.com/makinacorpus/landez
class MapnikRenderer(TileSource):
 def __init__(self, stylefile, tilesize=None):
  super(MapnikRenderer, self).__init__(tilesize)
  assert has_mapnik, _("Cannot render tiles without mapnik !")
  self.stylefile = stylefile
  self.basename = os.path.basename(self.stylefile)
  self._mapnik = None
  self._prj = None

 def tile(self, z, x, y):
  """
  Render the specified tile with Mapnik
  """
  logger.debug(_("Render tile %s") % ((z, x, y),))
  mercator = GlobalMercator(False,tilesize,[z])
  return self.render(mercator.tile_bbox((z, x, y)))

 def _prepare_rendering(self, bbox, width=None, height=None):
  if not self._mapnik:
   self._mapnik = mapnik.Map(width, height)
  # Load style XML
  mapnik.load_map(self._mapnik, self.stylefile, True)
  # Obtain <Map> projection
  self._prj = mapnik.Projection(self._mapnik.srs)
  # Convert to map projection
  assert len(bbox) == 4, _("Provide a bounding box tuple (minx, miny, maxx, maxy)")
  c0 = self._prj.forward(mapnik.Coord(bbox[0], bbox[1]))
  c1 = self._prj.forward(mapnik.Coord(bbox[2], bbox[3]))
  # Bounding box for the tile
  bbox = mapnik.Box2d(c0.x, c0.y, c1.x, c1.y)
  self._mapnik.resize(width, height)
  self._mapnik.zoom_to_box(bbox)
  self._mapnik.buffer_size = 128

 def render(self, bbox, width=None, height=None):
  """
  Render the specified tile with Mapnik
  """
  width = width or self.tilesize
  height = height or self.tilesize
  self._prepare_rendering(bbox, width=width, height=height)
  # Render image with default Agg renderer
  tmpfile = NamedTemporaryFile(delete=False)
  im = mapnik.Image(width, height)
  mapnik.render(self._mapnik, im)
  im.save(tmpfile.name, 'png256')  # TODO: mapnik output only to file?
  tmpfile.close()
  content = open(tmpfile.name).read()
  os.unlink(tmpfile.name)
  return content

 def grid(self, z, x, y, fields, layer):
  """
  Render the specified grid with Mapnik
  """
  logger.debug(_("Render grid %s") % ((z, x, y),))
  mercator = GlobalMercator(False,self.tilesize,[z])
  return self.render_grid(mercator.tile_bbox((z, x, y)), fields, layer)

 def render_grid(self, bbox, grid_fields, layer, width=None, height=None):
  """
  Render the specified grid with Mapnik
  """
  width = width or self.tilesize
  height = height or self.tilesize
  self._prepare_rendering(bbox, width=width, height=height)
  grid = mapnik.Grid(width, height)
  mapnik.render_layer(self._mapnik, grid, layer=layer, fields=grid_fields)
  grid = grid.encode()
  return json.dumps(grid)

# from Landez project
# https://github.com/makinacorpus/landez
class Cache(object):
 def __init__(self, **kwargs):
  self.extension = kwargs.get('extension', '.png')

 def tile_file(self, (z, x, y)):
  tile_dir = os.path.join("%s" % z, "%s" % x)
  tile_name = "%s%s" % (y, self.extension)
  return tile_dir, tile_name

 def read(self, (z, x, y)):
  raise NotImplementedError

 def save(self, body, (z, x, y)):
  raise NotImplementedError

 def remove(self, (z, x, y)):
  raise NotImplementedError

 def clean(self):
  raise NotImplementedError

# from Landez project
# https://github.com/makinacorpus/landez
class Dummy(Cache):
 def read(self, (z, x, y)):
  return None

 def save(self, body, (z, x, y)):
  pass

 def remove(self, (z, x, y)):
  pass

 def clean(self):
  pass

# from Landez project
# https://github.com/makinacorpus/landez
class Disk(Cache):
 def __init__(self, basename, folder, **kwargs):
  super(Disk, self).__init__(**kwargs)
  self._basename = None
  self._basefolder = folder
  self.folder = folder
  self.basename = basename

 @property
 def basename(self):
  return self._basename

 @basename.setter
 def basename(self, basename):
  self._basename = basename
  subfolder = re.sub(r'[^a-z^A-Z^0-9]+', '', basename.lower())
  self.folder = os.path.join(self._basefolder, subfolder)

 def tile_fullpath(self, (z, x, y)):
  tile_dir, tile_name = self.tile_file((z, x, y))
  tile_abs_dir = os.path.join(self.folder, tile_dir)
  return os.path.join(tile_abs_dir, tile_name)

 def remove(self, (z, x, y)):
  tile_abs_uri = self.tile_fullpath((z, x, y))
  os.remove(tile_abs_uri)
  parent = os.path.dirname(tile_abs_uri)
  i = 0
  while i <= 3:  # try to remove 3 levels (cache/z/x/)
   try:
    os.rmdir(parent)
    parent = os.path.dirname(parent)
    i += 1
   except OSError:
    break

 def read(self, (z, x, y)):
  tile_abs_uri = self.tile_fullpath((z, x, y))
  if os.path.exists(tile_abs_uri):
   logger.debug(_("Found %s") % tile_abs_uri)
   return open(tile_abs_uri, 'rb').read()
  return None

 def save(self, body, (z, x, y)):
  tile_abs_uri = self.tile_fullpath((z, x, y))
  tile_abs_dir = os.path.dirname(tile_abs_uri)
  if not os.path.isdir(tile_abs_dir):
   os.makedirs(tile_abs_dir)
  logger.debug(_("Save %s bytes to %s") % (len(body), tile_abs_uri))
  open(tile_abs_uri, 'wb').write(body)

 def clean(self):
  logger.debug(_("Clean-up %s") % self.folder)
  try:
   shutil.rmtree(self.folder)
  except OSError:
   logger.warn(_("%s was missing or read-only.") % self.folder)
# ----- end mbtiles.py
# ---- original gdal2mbtiles.py

resampling_list = ('average','near','bilinear','cubic','cubicspline','lanczos','antialias')
tile_formats_list = ('png', 'jpeg', 'hybrid')
profile_list = ('mercator','geodetic','raster','gearth','garmin') #,'zoomify')
webviewer_list = ('all','google','openlayers','none')

format_extension = {
 "PNG" : "png",
 "JPEG" : "jpg"
}

format_mime = {
 "PNG" : "image/png",
 "JPEG" : "image/jpeg"
}

jpeg_quality = 85
jpeg_gdal_options = ["QUALITY=%d" % jpeg_quality]

# =============================================================================
# =============================================================================
# =============================================================================

__doc__globalmaptiles = """
globalmaptiles.py

Global Map Tiles as defined in Tile Map Service (TMS) Profiles
==============================================================

Functions necessary for generation of global tiles used on the web.
It contains classes implementing coordinate conversions for:

  - GlobalMercator (based on EPSG:900913 = EPSG:3785)
       for Google Maps, Yahoo Maps, Microsoft Maps compatible tiles
  - GlobalGeodetic (based on EPSG:4326)
       for OpenLayers Base Map and Google Earth compatible tiles

More info at:

http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
http://wiki.osgeo.org/wiki/WMS_Tiling_Client_Recommendation
http://msdn.microsoft.com/en-us/library/bb259689.aspx
http://code.google.com/apis/maps/documentation/overlays.html#Google_Maps_Coordinates
- adapted for mbtiles - 2014 Mark Johnson
Created by Klokan Petr Pridal on 2008-07-03.
Google Summer of Code 2008, project GDAL2MbTiles for OSGEO.

In case you use this class in your product, translate it to another language
or find it usefull for your project please let me know.
My email: klokan at klokan dot cz.
I would like to know where it was used.

Class is available under the open-source GDAL license (www.gdal.org).
"""

import math

class GDAL2MbTiles(object):

 def flip_y(self,zoom, y):
  # z=19: y[ 352336 ]  y_osm[ 171951 ]  y_tms[ 352336 ]
  return (2**zoom-1) - y
 # -------------------------------------------------------------------------
 def process(self):
  """The main processing function, runs all the main steps of processing"""
  # Opening and preprocessing of the input file
  if self.options.mbtiles_fromdisk or self.options.mbtiles_todisk:
   if self.options.mbtiles_fromdisk:
    i_parm=10
   if self.options.mbtiles_todisk:
    i_parm=11
   if self.options.verbose:
    print "GDAL2MbTiles :mbtiles from/to disk [",i_parm,"] mbtiles_fromdisk[",self.options.mbtiles_fromdisk,"] mbtiles_todisk[",self.options.mbtiles_todisk,"]"
   self.mbtiles_setup(i_parm)
   return
  else:
   if self.options.verbose:
    print "GDAL2MbTiles :tile creation mbtiles[",self.options.mbtiles,"]"
   self.open_input()
   # Generation of main metadata files and HTML viewers
   self.generate_metadata()
   # Generation of the lowest tiles
   self.generate_base_tiles()
   # Generation of the overview tiles (higher in the pyramid)
   self.generate_overview_tiles()
   # Generating of KML
   self.generate_kml()

 # -------------------------------------------------------------------------
 def error(self, msg, details = "" ):
  """Print an error message and stop the processing"""

  if details is not None:
   msg += "\n\n" + details

  if not self.is_subprocess:
   self.parser.error(msg)
  else:
   raise Exception(msg)

 # -------------------------------------------------------------------------
 def progressbar(self, complete = 0.0):
  """Print progressbar for float value 0..1"""

  if self.is_subprocess:
   sys.stderr.write("%f\n" % complete)
   sys.stderr.flush()
  else:
   gdal.TermProgress_nocb(complete)

 # -------------------------------------------------------------------------
 def stop(self):
  """Stop the rendering immediately"""
  self.stopped = True

 # -------------------------------------------------------------------------
 def __init__(self, arguments, is_subprocess=False, gdalcache=None):
  """Constructor function - initialization"""
  self.mbtiles_db=None
  self.stopped = False
  self.input = None
  self.output = None
  self.is_subprocess = is_subprocess

  # Should we read bigger window of the input raster and scale it down?
  # Note: Modified leter by open_input()
  # Not for 'near' resampling
  # Not for Wavelet based drivers (JPEG2000, ECW, MrSID)
  # Not for 'raster' profile
  self.scaledquery = True

  # Should we use Read on the input file for generating overview tiles?
  # Note: Modified later by open_input()
  # Otherwise the overview tiles are generated from existing underlying tiles
  self.overviewquery = False

  # RUN THE ARGUMENT PARSER:

  self.optparse_init()
  self.options, self.args = self.parser.parse_args(args=arguments)
  if not self.args:
   self.error("No input file specified")

  # POSTPROCESSING OF PARSED ARGUMENTS:

  # Tile size
  try:
   self.tilesize = int(self.options.tilesize)
  except:
   self.tilesize = 256
   if self.options.profile == 'garmin':
    self.tilesize = 512

  # How big should be query window be for scaling down
  # Later on reset according the chosen resampling algorightm
  self.querysize = 4 * self.tilesize

  # Workaround for old versions of GDAL
  try:
   if (self.options.verbose and self.options.resampling == 'near') or gdal.TermProgress_nocb:
    pass
  except:
   self.error("This version of GDAL is not supported. Please upgrade to 1.6+.")
   #,"You can try run crippled version of gdal2mbtiles with parameters: -v -r 'near'")

  # Is output directory the last argument?
  # Test output directory, if it doesn't exist
  if os.path.isdir(self.args[-1]) or ( len(self.args) > 1 and not os.path.exists(self.args[-1])):
   self.output = self.args[-1]
   self.args = self.args[:-1]

  self.mbtiles_file=""
  # More files on the input not directly supported yet

  if (len(self.args) > 2):
   self.error("Processing of several input files is not supported.",
   """Please first use a tool like gdal_vrtmerge.py or gdal_merge.py on the files:
gdal_vrtmerge.py -o merged.vrt %s""" % " ".join(self.args))
   # TODO: Call functions from gdal_vrtmerge.py directly

  self.input = self.args[0]

  # Default values for not given options

  if not self.output:
   # Directory with input filename without extension in actual directory
   self.output = os.path.splitext(os.path.basename( self.input ))[0]

  if not self.options.title:
   self.options.title = os.path.basename( self.input )

  if self.options.url and not self.options.url.endswith('/'):
   self.options.url += '/'
  if self.options.url:
   self.options.url += os.path.basename( self.output ) + '/'

  # Supported options

  if self.options.resampling == 'average':
   try:
    if gdal.RegenerateOverview:
     pass
   except:
    self.error("'average' resampling algorithm is not available.", "Please use -r 'near' argument or upgrade to newer version of GDAL.")

  elif self.options.resampling == 'antialias':
   try:
    if numpy:
     pass
   except:
    self.error("'antialias' resampling algorithm is not available.", "Install PIL (Python Imaging Library) and numpy.")

  elif self.options.resampling == 'near':
   self.querysize = self.tilesize
  elif self.options.resampling == 'bilinear':
   self.querysize = self.tilesize * 2

  # Tile format.
  if self.options.tile_format is None:
   if self.options.profile == 'gearth':
    self.options.tile_format = 'hybrid'
   elif self.options.profile == 'garmin':
    self.options.tile_format = 'jpeg'
   else:
    self.options.tile_format = 'png'

  # Tile size and no webviewer for garmin profile
  if self.options.profile == 'garmin':
   self.options.webviewer = 'none'
   self.options.url = ''

  # Webviewer default depends on tile format.
  if self.options.webviewer is None:
   if self.options.tile_format == 'hybrid':
    self.options.webviewer = 'none'
   else:
    self.options.webviewer = 'all'
  # We don't support webviewers with hybrid trees yet.
  elif self.options.tile_format == 'hybrid' and self.options.webviewer != 'none':
   print ("WARNING: hybrid tile format is incompatible with webviewers you selected (%s), " +
       "so they will not be created.") % self.options.webviewer
   self.options.webviewer = "none"

  # No need support webviewers with mbtiles.
  if self.options.mbtiles:
   self.options.webviewer = "none"
  # User specified zoom levels
  self.tminz = None
  self.tmaxz = None
  if self.options.zoom:
   minmax = self.options.zoom.split('-',1)
   minmax.extend([''])
   min, max = minmax[:2]
   self.tminz = int(min)
   if max:
    self.tmaxz = int(max)
   else:
    self.tmaxz = int(min)

   self.te_minx=0
   self.te_miny=0
   self.te_maxx=0
   self.te_maxy=0
   if self.options.te_bounds != '':
    if self.options.te_bounds.find(',') != -1:
     sa_bounds=self.options.te_bounds.split(",")
    else:
     sa_bounds=self.options.te_bounds.split(" ")
    if len(sa_bounds) == 4:
     self.te_minx=float(sa_bounds[0])
     self.te_maxx=float(sa_bounds[2])
     self.te_maxy=float(sa_bounds[3])
     self.te_miny=float(sa_bounds[1])

  # KML generation
  self.kml = self.options.kml
  if self.options.kml_depth is not None:
   self.kml_depth = int(self.options.kml_depth)
   assert self.kml_depth > 0
  else:
   if self.options.profile == 'gearth':
    self.kml_depth = 3
   else:
    self.kml_depth = 1

  if self.options.kmz is None:
   self.options.kmz = self.options.profile == 'gearth'

  # GDAL Cache
  if gdalcache is not None: # default gdal.GetCacheMax() == 40*1024*1024:
   gdal.SetCacheMax(gdalcache)

  # Output the results
  if self.options.verbose:
   print "Options:", self.options
   print "Input:", self.input
   print "Output:", self.output
   print "Cache: %s MB" % (gdal.GetCacheMax() / 1024 / 1024)
   print

 # -------------------------------------------------------------------------
 def optparse_init(self):
  """Prepare the option parser for input (argv)"""

  from optparse import OptionParser, OptionGroup
  usage = "Usage: %prog [options] input_file(s) [output]"
  p = OptionParser(usage, version="%prog "+ __version__)
  p.add_option("-p", "--profile", dest='profile', type='choice', choices=profile_list,
        help="Tile cutting profile (%s) - default 'mercator' (Google Maps compatible)" % ",".join(profile_list))
  p.add_option("-r", "--resampling", dest="resampling", type='choice', choices=resampling_list,
      help="Resampling method (%s) - default 'average'" % ",".join(resampling_list))
  p.add_option("-f", "--tile-format", dest="tile_format", type='choice', choices=tile_formats_list,
      help="Image format of generated tiles (%s) - default 'png'" % ",".join(tile_formats_list))
  p.add_option('-s', '--s_srs', dest="s_srs", metavar="SRS",
        help="The spatial reference system used for the source input data")
  p.add_option('-z', '--zoom', dest="zoom",
        help="Zoom levels to render (format:'2-5' or '10').")
  p.add_option('-e', '--resume', dest="resume", action="store_true",
        help="Resume mode. Generate only missing files.")
  p.add_option('-a', '--srcnodata', dest="srcnodata", metavar="NODATA",
        help="NODATA transparency value to assign to the input data")
  p.add_option('-i', '--init-dest', dest="init_dest",
        help="Colour used to initialize output, only for 'jpeg' tile format")
  p.add_option('', '--tilesize', dest="tilesize",
        help="Size of the tiles - default 256")
  p.add_option('', '--osm', dest="tms_osm", action="store_true",
        help="tms or osm numbering - default tms")
  p.add_option('', '--mbtiles', dest="mbtiles", action="store_true",
        help="mbtiles - tiles creation to mbtiles file")
  p.add_option('', '--mbtiles_to_disk', dest="mbtiles_todisk", action="store_true",
        help="mbtiles tiles- write mbtiles tiles to a directory")
  p.add_option('', '--mbtiles_from_disk', dest="mbtiles_fromdisk", action="store_true",
        help="mbtiles tiles- create mbtiles file from tile directory")
  p.add_option('', "--te", dest="te_bounds",
        help="bounds to extract (coordinates in the output SRS): xmin ymin xmax ymax OR xmin,ymin,xmax,ymax")
  p.add_option("-v", "--verbose", dest="verbose",action="store_true",
        help="Print status messages to stdout")
  # KML options
  g = OptionGroup(p, "KML (Google Earth) options", "Options for generated Google Earth SuperOverlay metadata")
  g.add_option("-k", "--force-kml", dest='kml', action="store_true",
        help="Generate KML for Google Earth - default for 'geodetic' profile and 'raster' in EPSG:4326. For a dataset with different projection use with caution!")
  g.add_option("-n", "--no-kml", dest='kml', action="store_false",
        help="Avoid automatic generation of KML files for EPSG:4326")
  g.add_option("-u", "--url", dest='url',
        help="URL address where the generated tiles are going to be published")
  g.add_option('-d', '--kml-depth', dest="kml_depth",
        help="How many levels to store before linking, default 1")
  g.add_option('--kmz', dest="kmz", action="store_true",
        help="Compress KML files into KMZ format, default for 'gearth' profile")
  g.add_option('--no-kmz', dest="kmz", action="store_false",
        help="Do not compress KML files into KMZ format, default for 'mercator', 'geodetic' and 'raster' profiles")
  p.add_option_group(g)

  # HTML options
  g = OptionGroup(p, "Web viewer options", "Options for generated HTML viewers a la Google Maps")
  g.add_option("-w", "--webviewer", dest='webviewer', type='choice', choices=webviewer_list,
        help="Web viewer to generate (%s) - default 'all'" % ",".join(webviewer_list))
  g.add_option("-t", "--title", dest='title',
        help="Title of the map")
  g.add_option("-c", "--copyright", dest='copyright',
        help="Copyright for the map")
  g.add_option("-g", "--googlekey", dest='googlekey',
        help="Google Maps API key from http://code.google.com/apis/maps/signup.html")
  g.add_option("-y", "--yahookey", dest='yahookey',
        help="Yahoo Application ID from http://developer.yahoo.com/wsregapp/")
  p.add_option_group(g)

  # TODO: MapFile + TileIndexes per zoom level for efficient MapServer WMS
  #g = OptionGroup(p, "WMS MapServer metadata", "Options for generated mapfile and tileindexes for MapServer")
  #g.add_option("-i", "--tileindex", dest='wms', action="store_true"
  #      help="Generate tileindex and mapfile for MapServer (WMS)")
  # p.add_option_group(g)

  p.set_defaults(verbose=False, profile="mercator", kml=False, url=None,
  copyright='', resampling='average', resume=False, tilesize=None,mbtiles=False,tms_osm=False,
  mbtiles_todisk=False,mbtiles_fromdisk=False,te_bounds='',
  googlekey='INSERT_YOUR_KEY_HERE', yahookey='INSERT_YOUR_YAHOO_APP_ID_HERE')

  self.parser = p

 # -------------------------------------------------------------------------
 def open_input(self):
  """Initialization of the input raster, reprojection if necessary"""
  gdal.SetConfigOption("GDAL_PAM_ENABLED", "YES")
  gdal.AllRegister()
  # self.options.verbose=True
  if self.options.tms_osm:
   self.s_y_type="osm"
  else:
   self.s_y_type="tms"
  if self.options.verbose:
   print "open_input :", self.input," osm[",self.options.tms_osm,",",self.s_y_type,"] mbtiles[",self.options.mbtiles,"] mbtiles_todisk[",self.options.mbtiles_todisk,"] mbtiles_fromdisk[",self.options.mbtiles_fromdisk,"]";
  # Open the input file
  if self.input:
   self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
  else:
   raise Exception("No input file was specified")

  if self.options.verbose:
   print "Input file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount)

  if not self.in_ds:
   # Note: GDAL prints the ERROR message too
   self.error("It is not possible to open the input file '%s'." % self.input )

  # Read metadata from the input file
  if self.in_ds.RasterCount == 0:
   self.error( "Input file '%s' has no raster band" % self.input )

  if self.in_ds.GetRasterBand(1).GetRasterColorTable():
   # TODO: Process directly paletted dataset by generating VRT in memory
   self.error( "Please convert this file to RGB/RGBA and run gdal2mbtiles on the result.",
   """From paletted file you can create RGBA file (temp.vrt) by:
gdal_translate -of vrt -expand rgba %s temp.vrt
then run:
gdal2mbtiles temp.vrt""" % self.input )

  # Get NODATA value
  # User supplied values overwrite everything else.
  if self.options.srcnodata is not None:
   nds = map(float, self.options.srcnodata.split(','))
   if len(nds) < self.in_ds.RasterCount:
    self.in_nodata = (nds * self.in_ds.RasterCount)[:self.in_ds.RasterCount]
   else:
    self.in_nodata = nds
  else:
   # If the source dataset has NODATA, use it.
   self.in_nodata = []
   for i in range(1, self.in_ds.RasterCount+1):
    if self.in_ds.GetRasterBand(i).GetNoDataValue() != None:
     self.in_nodata.append( self.in_ds.GetRasterBand(i).GetNoDataValue() )

  if self.options.verbose:
   print "NODATA: %s" % self.in_nodata

  # INIT DEST
  if self.options.init_dest is not None:
   if self.options.tile_format == "jpeg":
    if self.in_ds.RasterCount == 4:
     nbands = 3
    else:
     nbands = self.in_ds.RasterCount

    nds = map(float, self.options.init_dest.split(','))

    if len(nds) == 1:
     init_dest = nds * nbands
    elif len(nds) == nbands:
     init_dest = nds
    else:
     print "WARNING: you suplied %d '--init-dest' values but the dataset has %d data bands" % (len(nds), nbands)
     init_dest = None
   else:
    init_dest = None
    print "WARNING: --init-dest can be used only with 'jpeg' tile format"
  else:
   if self.options.tile_format == "jpeg":
    init_dest = [255,255,255]
   else:
    init_dest = None

  #
  # Here we should have RGBA input dataset opened in self.in_ds
  #

  if self.options.verbose:
   print "Preprocessed file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount)

  # Spatial Reference System of the input raster


  self.in_srs = None

  if self.options.s_srs:
   self.in_srs = osr.SpatialReference()
   self.in_srs.SetFromUserInput(self.options.s_srs)
   self.in_srs_wkt = self.in_srs.ExportToWkt()
  else:
   self.in_srs_wkt = self.in_ds.GetProjection()
   if not self.in_srs_wkt and self.in_ds.GetGCPCount() != 0:
    self.in_srs_wkt = self.in_ds.GetGCPProjection()
   if self.in_srs_wkt:
    self.in_srs = osr.SpatialReference()
    self.in_srs.ImportFromWkt(self.in_srs_wkt)
   #elif self.options.profile != 'raster':
   # self.error("There is no spatial reference system info included in the input file.","You should run gdal2mbtiles with --s_srs EPSG:XXXX or similar.")

  # Spatial Reference System of tiles

  self.out_srs = osr.SpatialReference()

  if self.options.profile == 'mercator':
   self.out_srs.ImportFromEPSG(900913)
  elif self.options.profile in ('geodetic', 'gearth', 'garmin'):
   self.out_srs.ImportFromEPSG(4326)
  else:
   self.out_srs = self.in_srs

  # Are the reference systems the same? Reproject if necessary.

  self.out_ds = None

  if self.options.profile in ('mercator', 'geodetic', 'gearth', 'garmin'):

   if (self.in_ds.GetGeoTransform() == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)) and (self.in_ds.GetGCPCount() == 0):
    self.error("There is no georeference - neither affine transformation (worldfile) nor GCPs. You can generate only 'raster' profile tiles.",
    "Either gdal2mbtiles with parameter -p 'raster' or use another GIS software for georeference e.g. gdal_transform -gcp / -a_ullr / -a_srs")

   if self.in_srs:

    if (self.in_srs.ExportToProj4() != self.out_srs.ExportToProj4()) or (self.in_ds.GetGCPCount() != 0):

     # Generation of VRT dataset in tile projection, default 'nearest neighbour' warping
     self.out_ds = gdal.AutoCreateWarpedVRT( self.in_ds, self.in_srs_wkt, self.out_srs.ExportToWkt() )

     # TODO: HIGH PRIORITY: Correction of AutoCreateWarpedVRT according the max zoomlevel for correct direct warping!!!

     if self.options.verbose:
      print "Warping of the raster by AutoCreateWarpedVRT (result saved into 'tiles.vrt')"
      self.out_ds.GetDriver().CreateCopy("tiles.vrt", self.out_ds)

     # Note: self.in_srs and self.in_srs_wkt contain still the non-warped reference system!!!

     # Correction of AutoCreateWarpedVRT for NODATA values
     if self.in_nodata != []:
      import tempfile
      tempfilename = tempfile.mktemp('-gdal2mbtiles.vrt')
      self.out_ds.GetDriver().CreateCopy(tempfilename, self.out_ds)
      # open as a text file
      s = open(tempfilename).read()
      # Add the warping options
      s = s.replace("""<GDALWarpOptions>""","""<GDALWarpOptions>
   <Option name="UNIFIED_SRC_NODATA">YES</Option>
   <Option name="INIT_DEST">NO_DATA</Option>""")
      # replace BandMapping tag for NODATA bands....
      if init_dest is None:
       dstnodata = self.in_nodata
      else:
       dstnodata = init_dest
      for i in range(len(self.in_nodata)):
       s = s.replace("""<BandMapping src="%i" dst="%i"/>""" % ((i+1),(i+1)),"""<BandMapping src="%i" dst="%i">
       <SrcNoDataReal>%i</SrcNoDataReal>
       <SrcNoDataImag>0</SrcNoDataImag>
       <DstNoDataReal>%i</DstNoDataReal>
       <DstNoDataImag>0</DstNoDataImag>
     </BandMapping>""" % ((i+1), (i+1), self.in_nodata[i], dstnodata[i]))
      # save the corrected VRT
      open(tempfilename,"w").write(s)
      # open by GDAL as self.out_ds
      self.out_ds = gdal.Open(tempfilename) #, gdal.GA_ReadOnly)
      # delete the temporary file
      os.unlink(tempfilename)

      # set NODATA_VALUE metadata
      self.out_ds.SetMetadataItem('NODATA_VALUES','%s' % " ".join(str(int(f)) for f in self.in_nodata))

      if self.options.verbose:
       print "Modified warping result saved into 'tiles1.vrt'"
       open("tiles1.vrt","w").write(s)

     # -----------------------------------
     # Correction of AutoCreateWarpedVRT for Mono (1 band) and RGB (3 bands) files without NODATA:
     # equivalent of gdalwarp -dstalpha
     elif self.in_nodata == [] and self.out_ds.RasterCount in (1,3):
      import tempfile
      tempfilename = tempfile.mktemp('-gdal2mbtiles.vrt')
      self.out_ds.GetDriver().CreateCopy(tempfilename, self.out_ds)
      # open as a text file
      s = open(tempfilename).read()
      # Add the warping options
      s = s.replace("""<BlockXSize>""","""<VRTRasterBand dataType="Byte" band="%i" subClass="VRTWarpedRasterBand">
    <ColorInterp>Alpha</ColorInterp>
  </VRTRasterBand>
  <BlockXSize>""" % (self.out_ds.RasterCount + 1))
      s = s.replace("""</GDALWarpOptions>""", """<DstAlphaBand>%i</DstAlphaBand>
  </GDALWarpOptions>""" % (self.out_ds.RasterCount + 1))
      if init_dest is None:
       init_dest_str = "0"
      else:
       init_dest_str = ",".join(str(f) for f in init_dest)
      s = s.replace("""</WorkingDataType>""", """</WorkingDataType>
    <Option name="INIT_DEST">%s</Option>""" % init_dest_str)
      # save the corrected VRT
      open(tempfilename,"w").write(s)
      # open by GDAL as self.out_ds
      self.out_ds = gdal.Open(tempfilename) #, gdal.GA_ReadOnly)
      # delete the temporary file
      os.unlink(tempfilename)

      if self.options.verbose:
       print "Modified -dstalpha warping result saved into 'tiles1.vrt'"
       open("tiles1.vrt","w").write(s)

     elif init_dest is not None:
      import tempfile
      tempfilename = tempfile.mktemp('-gdal2mbtiles.vrt')
      self.out_ds.GetDriver().CreateCopy(tempfilename, self.out_ds)
      # open as a text file
      s = open(tempfilename).read()
      # Add the warping options
      s = s.replace("""</WorkingDataType>""", """</WorkingDataType>
    <Option name="INIT_DEST">%s</Option>""" % ",".join(str(f) for f in init_dest))
      # save the corrected VRT
      open(tempfilename,"w").write(s)
      # open by GDAL as self.out_ds
      self.out_ds = gdal.Open(tempfilename) #, gdal.GA_ReadOnly)
      # delete the temporary file
      os.unlink(tempfilename)

      if self.options.verbose:
       print "Modified warping result saved into 'tiles1.vrt'"
       open("tiles1.vrt","w").write(s)

     # For raster with 4-bands: 4th unknown band set to alpha
     if (self.out_ds.RasterCount == 4
      and self.out_ds.GetRasterBand(4).GetRasterColorInterpretation() == gdal.GCI_Undefined):
      self.out_ds.GetRasterBand(4).SetRasterColorInterpretation(gdal.GCI_AlphaBand)

     s = '''
     '''

   else:
    self.error("Input file has unknown SRS.", "Use --s_srs ESPG:xyz (or similar) to provide source reference system." )

   if self.out_ds and self.options.verbose:
    print "Projected file:", "tiles.vrt", "( %sP x %sL - %s bands)" % (self.out_ds.RasterXSize, self.out_ds.RasterYSize, self.out_ds.RasterCount)

  if not self.out_ds:
   self.out_ds = self.in_ds

  #
  # Here we should have a raster (out_ds) in the correct Spatial Reference system
  #

  # KML test
  self.isepsg4326 = False
  srs4326 = osr.SpatialReference()
  srs4326.ImportFromEPSG(4326)
  if self.out_srs and srs4326.ExportToProj4() == self.out_srs.ExportToProj4():
   self.kml = True
   self.isepsg4326 = True
   if self.options.verbose:
    print "KML autotest OK!"

  # Read the georeference

  self.out_gt = self.out_ds.GetGeoTransform()

  #originX, originY = self.out_gt[0], self.out_gt[3]
  #pixelSize = self.out_gt[1] # = self.out_gt[5]

  # Test the size of the pixel

  # MAPTILER - COMMENTED
  #if self.out_gt[1] != (-1 * self.out_gt[5]) and self.options.profile != 'raster':
   # TODO: Process corectly coordinates with are have swichted Y axis (display in OpenLayers too)
   #self.error("Size of the pixel in the output differ for X and Y axes.")

  # Report error in case rotation/skew is in geotransform (possible only in 'raster' profile)
  if (self.out_gt[2], self.out_gt[4]) != (0,0):
   self.error("Georeference of the raster contains rotation or skew. Such raster is not supported. Please use gdalwarp first.")
   # TODO: Do the warping in this case automaticaly

  #
  # Here we expect: pixel is square, no rotation on the raster
  #

  # Output Bounds - coordinates in the output SRS
  self.ominx = self.out_gt[0]
  self.omaxx = self.out_gt[0]+self.out_ds.RasterXSize*self.out_gt[1]
  self.omaxy = self.out_gt[3]
  self.ominy = self.out_gt[3]-self.out_ds.RasterYSize*self.out_gt[1]
  # Note: maybe round(x, 14) to avoid the gdal_translate behaviour, when 0 becomes -1e-15
  # user defined bounds to extract  - coordinates in the output SRS
  if self.options.te_bounds != '':
   if self.te_minx >= self.ominx and self.te_minx <= self.omaxx:
    if self.te_maxx >= self.ominx and self.te_maxx <= self.omaxx:
     if self.te_miny >= self.ominy and self.te_miny <= self.omaxy:
      if self.te_maxy >= self.ominy and self.te_maxy <= self.omaxy:
       # replace only if inside the read bounds
       self.ominx = self.te_minx
       self.omaxx = self.te_maxx
       self.ominy = self.te_miny
       self.omaxy = self.te_maxy
       if self.options.verbose:
        print "User defined Bounds (output srs) have been set:", round(self.ominx, 13), self.ominy, self.omaxx, self.omaxy

  if self.options.verbose:
   print "Bounds (output srs):", round(self.ominx, 13), self.ominy, self.omaxx, self.omaxy

  if self.options.mbtiles:
   self.options.profile = 'mercator'
  if self.options.profile == 'mercator':
   self.mercator = GlobalMercator(self.options.tms_osm) # from globalmaptiles.py

  #
  # Calculating ranges for tiles in different zoom levels
  #

   # Function which generates SWNE in LatLong for given tile
   self.tileswne = self.mercator.TileLatLonBounds

   # Generate table with min max tile coordinates for all zoomlevels
   self.tminmax = range(0,32)
   for tz in range(0, 32):
    tminx, tminy = self.mercator.MetersToTile( self.ominx, self.ominy, tz )
    tmaxx, tmaxy = self.mercator.MetersToTile( self.omaxx, self.omaxy, tz )
    # crop tiles extending world limits (+-180,+-90)
    tminx, tminy = max(0, tminx), max(0, tminy)
    tmaxx, tmaxy = min(2**tz-1, tmaxx), min(2**tz-1, tmaxy)
    self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

   # TODO: Maps crossing 180E (Alaska?)

   # Get the minimal zoom level (map covers area equivalent to one tile)
   if self.tminz == None:
    self.tminz = self.mercator.ZoomForPixelSize( self.out_gt[1] * max( self.out_ds.RasterXSize, self.out_ds.RasterYSize) / float(self.tilesize) )

   # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
   if self.tmaxz == None:
    self.tmaxz = self.mercator.ZoomForPixelSize( self.out_gt[1] )

   if self.options.verbose:
    print "Bounds (latlong):", self.mercator.MetersToLatLon( self.ominx, self.ominy), self.mercator.MetersToLatLon( self.omaxx, self.omaxy)
    print 'MinZoomLevel:', self.tminz
    print "MaxZoomLevel:", self.tmaxz, "(", self.mercator.Resolution( self.tmaxz ),")"

  # this must be call befor ImageOutput is called (self.output may be changed)
  if self.options.mbtiles:
   if not self.mbtiles_db:
    self.mbtiles_setup(1);

  # Instantiate image output.
  self.image_output = ImageOutput(self.options.tile_format, self.out_ds, self.tilesize,
          self.options.resampling, init_dest, self.output,
          self.options.verbose,self.options.mbtiles)
  if self.options.profile == 'geodetic':

   self.geodetic = GlobalGeodetic() # from globalmaptiles.py

   # Function which generates SWNE in LatLong for given tile
   self.tileswne = self.geodetic.TileLatLonBounds

   # Generate table with min max tile coordinates for all zoomlevels
   self.tminmax = range(0,32)
   for tz in range(0, 32):
    tminx, tminy = self.geodetic.LatLonToTile( self.ominx, self.ominy, tz )
    tmaxx, tmaxy = self.geodetic.LatLonToTile( self.omaxx, self.omaxy, tz )
    # crop tiles extending world limits (+-180,+-90)
    tminx, tminy = max(0, tminx), max(0, tminy)
    tmaxx, tmaxy = min(2**(tz+1)-1, tmaxx), min(2**tz-1, tmaxy)
    self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

   # TODO: Maps crossing 180E (Alaska?)

   # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
   if self.tminz == None:
    self.tminz = self.geodetic.ZoomForPixelSize( self.out_gt[1] * max( self.out_ds.RasterXSize, self.out_ds.RasterYSize) / float(self.tilesize) )

   # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
   if self.tmaxz == None:
    self.tmaxz = self.geodetic.ZoomForPixelSize( self.out_gt[1] )

   if self.options.verbose:
    print "Bounds (latlong):", self.ominx, self.ominy, self.omaxx, self.omaxy

  if self.options.profile in ('raster', 'gearth', 'garmin'):

   log2 = lambda x: math.log10(x) / math.log10(2) # log2 (base 2 logarithm)

   self.nativezoom = int(max( math.ceil(log2(self.out_ds.RasterXSize/float(self.tilesize))),
                              math.ceil(log2(self.out_ds.RasterYSize/float(self.tilesize)))))

   if self.options.verbose:
    print "Native zoom of the raster:", self.nativezoom

   # Get the minimal zoom level (whole raster in one tile)
   if self.tminz == None:
    self.tminz = 0

   # Get the maximal zoom level (native resolution of the raster)
   if self.tmaxz == None:
    self.tmaxz = self.nativezoom

   # Garmin has maximally 100 tiles - lower the tmaxz if necessary
   if self.options.profile == 'garmin':
    tno = math.ceil(self.out_ds.RasterXSize / self.tilesize) * math.ceil(self.out_ds.RasterYSize / self.tilesize)
    for tz in range(self.tmaxz, 1, -1):
     if tno > 100:
      tno /= 4
      self.tmaxz -= 1
      print "Warning: GARMIN has a limit 100 tiles per device: lowering the max zoom level to:", self.tmaxz
     else:
      continue

   # Force only one zoom level for the 'garmin' tile profile
   if self.options.profile == 'garmin':
    self.tminz = self.tmaxz

   # Generate table with min max tile coordinates for all zoomlevels
   self.tminmax = range(0, self.tmaxz+1)
   self.tsize = range(0, self.tmaxz+1)
   for tz in range(0, self.tmaxz+1):
    tsize = 2.0**(self.nativezoom-tz)*self.tilesize
    tminx, tminy = 0, 0
    tmaxx = int(math.ceil( self.out_ds.RasterXSize / tsize )) - 1
    tmaxy = int(math.ceil( self.out_ds.RasterYSize / tsize )) - 1
    self.tsize[tz] = math.ceil(tsize)
    self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

   # Function which generates SWNE in LatLong for given tile
   if self.kml and self.in_srs_wkt:
    self.ct = osr.CoordinateTransformation(self.in_srs, srs4326)
    def rastertileswne(x,y,z):
     pixelsizex = (2**(self.nativezoom-z) * self.out_gt[1]) # X-pixel size in level
     pixelsizey = (2**(self.nativezoom-z) * self.out_gt[5]) # Y-pixel size in level (usually -1*pixelsizex)
     west = self.out_gt[0] + x*self.tilesize*pixelsizex
     east = west + self.tilesize*pixelsizex
     south = self.ominy + y*self.tilesize*pixelsizex
     north = south + self.tilesize*pixelsizex
     if not self.isepsg4326:
      # Transformation to EPSG:4326 (WGS84 datum)
      west, south = self.ct.TransformPoint(west, south)[:2]
      east, north = self.ct.TransformPoint(east, north)[:2]
     return south, west, north, east

    self.tileswne = rastertileswne
   else:
    self.tileswne = lambda x, y, z: (0,0,0,0)

 # -------------------------------------------------------------------------
 def generate_metadata(self):
  """Generation of main metadata files and HTML viewers (metadata related to particular tiles are generated during the tile processing)."""
  if self.options.mbtiles:
   return
  if not os.path.exists(self.output):
   os.makedirs(self.output)

  if self.options.profile == 'mercator':

   south, west = self.mercator.MetersToLatLon( self.ominx, self.ominy)
   north, east = self.mercator.MetersToLatLon( self.omaxx, self.omaxy)
   south, west = max(-85.05112878, south), max(-180.0, west)
   north, east = min(85.05112878, north), min(180.0, east)
   self.swne = (south, west, north, east)

   # Generate googlemaps.html
   if self.options.webviewer in ('all','google') and self.options.profile == 'mercator':
    if not self.options.resume or not os.path.exists(os.path.join(self.output, 'googlemaps.html')):
     f = open(os.path.join(self.output, 'googlemaps.html'), 'w')
     f.write( self.generate_googlemaps() )
     f.close()

   # Generate openlayers.html
   if self.options.webviewer in ('all','openlayers'):
    if not self.options.resume or not os.path.exists(os.path.join(self.output, 'openlayers.html')):
     f = open(os.path.join(self.output, 'openlayers.html'), 'w')
     f.write( self.generate_openlayers() )
     f.close()

  elif self.options.profile == 'geodetic':

   west, south = self.ominx, self.ominy
   east, north = self.omaxx, self.omaxy
   south, west = max(-90.0, south), max(-180.0, west)
   north, east = min(90.0, north), min(180.0, east)
   self.swne = (south, west, north, east)

   # Generate openlayers.html
   if self.options.webviewer in ('all','openlayers'):
    if not self.options.resume or not os.path.exists(os.path.join(self.output, 'openlayers.html')):
     f = open(os.path.join(self.output, 'openlayers.html'), 'w')
     f.write( self.generate_openlayers() )
     f.close()

  elif self.options.profile in ['raster','gearth','garmin']:

   west, south = self.ominx, self.ominy
   east, north = self.omaxx, self.omaxy

   self.swne = (south, west, north, east)

   # Generate openlayers.html
   if self.options.webviewer in ('all','openlayers'):
    if not self.options.resume or not os.path.exists(os.path.join(self.output, 'openlayers.html')):
     f = open(os.path.join(self.output, 'openlayers.html'), 'w')
     f.write( self.generate_openlayers() )
     f.close()


  # Generate tilemapresource.xml.
  if (self.options.tile_format != 'hybrid' and self.options.profile != 'garmin'
   and (not self.options.resume or not os.path.exists(os.path.join(self.output, 'tilemapresource.xml')))):
   f = open(os.path.join(self.output, 'tilemapresource.xml'), 'w')
   f.write( self.generate_tilemapresource())
   f.close()

 # -------------------------------------------------------------------------
 def generate_base_tiles(self):
  """Generation of the base tiles (the lowest in the pyramid) directly from the input raster"""

  gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")

  print "Generating Base Tiles:"
  if self.options.verbose:
   #mx, my = self.out_gt[0], self.out_gt[3] # OriginX, OriginY
   #px, py = self.mercator.MetersToPixels( mx, my, self.tmaxz)
   #print "Pixel coordinates:", px, py, (mx, my)
   print
   print "Tiles generated from the max zoom level:"
   print "----------------------------------------"
   print


  # Set the bounds
  tminx, tminy, tmaxx, tmaxy = self.tminmax[self.tmaxz]
  querysize = self.querysize

  # Just the center tile
  #tminx = tminx+ (tmaxx - tminx)/2
  #tminy = tminy+ (tmaxy - tminy)/2
  #tmaxx = tminx
  #tmaxy = tminy

  #print tminx, tminy, tmaxx, tmaxy
  tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
  #print tcount
  ti = 0
  i_y_column_count=((tmaxy-tminy)+1)
  ds = self.out_ds
  tz = self.tmaxz
  if self.options.verbose:
   # tx in range(tminx, tmaxx+1) tminx[ 281596 ] tmaxx[ 281744 ] ; ((tmaxx-tmaxy)+1) x_tiles[ 23393 ]
   print "\ttz=[",tz,"] : tx in range(tminx, tmaxx+1) tminx[",tminx,"] tmaxx[",tmaxx,"] ; ((tmaxx-tmaxy)+1) x_tiles[",tcount,"]"
   #  ty_tms in range(tmaxy, tminy-1, -1) tmaxy[ 352409 ] tminy[ 352253 ] ; ((tmaxy-tminy)) y_tiles[ 157 ] 352409-(352253-1)
   print "\ttz=[",tz,"] : ty_tms in range(tmaxy, tminy-1, -1) tmaxy[",tmaxy,"] tminy[",tminy,"] ; ((tmaxy-tminy+1)) y_tiles[",i_y_column_count,"]"
  if self.options.resume:
   i_count = self.tile_exists(0, 0, tz,2)
   if i_count == tcount:
    if self.options.verbose:
     print "\tTile generation skipped because of --resume ;  x/y-tiles of z[",tz,"]  y_tiles[",tcount,"]"
    return
  for tx in range(tminx, tmaxx+1):
   tmaxy_work=tmaxy
   if self.options.resume:
    i_count = self.tile_exists(tx, 0, tz,3)
    if i_count == i_y_column_count:
     if self.options.verbose:
      print "\tTile generation skipped because of --resume ;  z =",tz," ; y-tiles of x[",tx,"]  y_tiles[",i_y_column_count,"]"
     break
    else:
     if i_count > 0:
      # this assums the rows are compleate, which may NOT be true
      tmaxy_work-=i_count
      if self.options.verbose:
       print "\tTile generation skipped to tmaxy[",tmaxy_work,"] because of --resume ;  z =",tz," ; y-tiles of x[",tx,"]  y_tiles[",i_y_column_count,"]"
   for ty_tms in range(tmaxy_work, tminy-1, -1): #range(tminy, tmaxy+1):
    ty_osm=self.flip_y(tz,ty_tms)
    ty=ty_tms
    if self.options.tms_osm:
     ty=ty_osm
    if self.stopped:
     if self.options.mbtiles:
      if self.mbtiles_db:
       self.mbtiles_db.close_db()
     break
    ti += 1

    if self.options.resume:
     exists = self.tile_exists(tx, ty, tz,0)
     if exists and self.options.verbose:
      print "\tTile generation skipped because of --resume ;  z =",tz," ; x =",tx," ; y_tms =",ty_tms, "; y_osm =",ty_osm
    else:
     exists = False

    if not exists:
     if self.options.verbose:
      print ti, '/', tcount, self.get_verbose_tile_name(tx, ty, tz)
     # Don't scale up by nearest neighbour, better change the querysize
     # to the native resolution (and return smaller query tile) for scaling
     if self.options.profile in ('mercator','geodetic'):
      if self.options.profile == 'mercator':
       # Tile bounds in EPSG:900913
       b = self.mercator.TileBounds(tx, ty_tms, tz)
      elif self.options.profile == 'geodetic':
       b = self.geodetic.TileBounds(tx, ty_tms, tz)

      rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1])
      nativesize = wb[0]+wb[2] # Pixel size in the raster covering query geo extent
      if self.options.verbose:
       print "\tNative Extent (querysize",nativesize,"): ", rb, wb

      querysize = self.querysize
      # Tile bounds in raster coordinates for ReadRaster query
      rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1], querysize=querysize)

      rx, ry, rxsize, rysize = rb
      wx, wy, wxsize, wysize = wb
     else: # 'raster' or 'gearth' or 'garmin' profile:
      tsize = int(self.tsize[tz]) # tilesize in raster coordinates for actual zoom
      xsize = self.out_ds.RasterXSize # size of the raster in pixels
      ysize = self.out_ds.RasterYSize
      if tz >= self.nativezoom:
       querysize = self.tilesize # int(2**(self.nativezoom-tz) * self.tilesize)

      rx = (tx) * tsize
      rxsize = 0
      if tx == tmaxx:
       rxsize = xsize % tsize
      if rxsize == 0:
       rxsize = tsize

      rysize = 0
      if ty_tms == tmaxy:
       rysize = ysize % tsize
      if rysize == 0:
       rysize = tsize
      ry = ysize - (ty_tms * tsize) - rysize

      wx, wy = 0, 0

      wxsize, wysize = int(rxsize/float(tsize) * querysize), int(rysize/float(tsize) * querysize)
      if wysize != querysize:
       wy = querysize - wysize
     xyzzy = Xyzzy(querysize, rx, ry, rxsize, rysize, wx, wy, wxsize, wysize)
     try:
      if self.options.verbose:
       print ti,'/',tcount,' total ; z =',tz,' ; x =',tx,' ; y_tms =',ty_tms,' ; y_osm =',ty_osm
       print "\tReadRaster Extent: ", (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)
      self.write_base_tile(tx, ty, tz, xyzzy)
     except ImageOutputException, e:
      self.error("'%d/%d/%d': %s" % (tz, tx, ty, e.message))

    if not self.options.verbose or self.is_subprocess:
     self.progressbar( ti / float(tcount) )
  if self.options.mbtiles:
   if self.mbtiles_db:
    self.mbtiles_db.close_db()
    self.mbtiles_db=None

 # -------------------------------------------------------------------------
 # MBTiles support -begin -
 def mbtiles_setup(self,i_parm):
  self.mbtiles_format=self.options.tile_format
  if self.mbtiles_format == "jpeg":
   self.mbtiles_format="jpg"
  if self.options.tms_osm:
   self.s_y_type="osm"
  else:
   self.s_y_type="tms"
  if not self.mbtiles_db:
   if i_parm == 1:
    if self.output == "":
     self.mbtiles_dir=os.path.dirname(self.input)+ '/'
     self.mbtiles_file=self.mbtiles_dir+os.path.splitext(os.path.basename( self.input ))[0] + '.mbtiles'
    else:
     if self.mbtiles_file == '':
      self.mbtiles_file=self.output
      self.mbtiles_dir=os.path.dirname(self.mbtiles_file)+ '/'
      self.output=self.mbtiles_dir
   if i_parm == 10:
    self.mbtiles_file=self.output
    self.mbtiles_dir=os.path.dirname(self.mbtiles_file)+ '/'
   if i_parm == 11:
    self.mbtiles_file=self.input
    self.mbtiles_dir=os.path.dirname(self.mbtiles_file)+ '/'
   self.mbtiles_db=MbTiles()
   # if self.options.verbose:
   self.mbtiles_db.open_db(self.mbtiles_file.strip(),self.mbtiles_dir,self.mbtiles_format,self.s_y_type,self.options.verbose)
   if i_parm == 1:
    minLat, minLon = self.mercator.MetersToLatLon(self.ominx,self.ominy)
    maxLat, maxLon = self.mercator.MetersToLatLon(self.omaxx, self.omaxy)
    mbtiles_center_x=(maxLon+minLon)/2
    mbtiles_center_y=(maxLat+minLat)/2
    mbtiles_name=self.options.title
    mbtiles_description=self.options.copyright
    mbtiles_minzoom="%d"%(self.tminz)
    mbtiles_maxzoom="%d"%(self.tmaxz)
    mbtiles_bounds="%f,%f,%f,%f"% (minLon,minLat,maxLon,maxLat)
    mbtiles_center="%f,%f,%s"%(mbtiles_center_x,mbtiles_center_y,mbtiles_minzoom)
    values_list = [
         ('name', mbtiles_name),
         ('description',mbtiles_description ),
         ('type', 'baselayer'),
         ('version', '1.1'),
         ('tile_row_type',self.s_y_type),
         ('format', self.mbtiles_format),
         ('bounds', mbtiles_bounds),
         ('center', mbtiles_center),
         ('minzoom', mbtiles_minzoom),
         ('maxzoom', mbtiles_maxzoom)
        ]
    # print "MbTiles :", self.options.mbtiles,"[",self.mbtiles_file,"]",values_list
    self.mbtiles_db.insert_metadata(values_list)
   if i_parm == 10:
    print "input_dir[",self.input,"]"
    self.mbtiles_db.mbtiles_from_disk(self.input)
   if i_parm == 11:
    self.mbtiles_db.mbtiles_to_disk(self.output)

 # -------------------------------------------------------------------------
 def tile_exists(self,tx, ty, tz, i_parm):
  if self.options.mbtiles:
   if not self.mbtiles_db:
    self.mbtiles_setup(1);
   i_count = self.mbtiles_db.count_tiles(tz,tx,ty,i_parm)
  else:
   if i_parm == 0:
    i_count = self.image_output.tile_exists(tx, ty, tz)
   else:
    i_count=0
  return i_count
 # -------------------------------------------------------------------------
 def write_base_tile(self, tx, ty, tz, xyzzy):
  if self.options.mbtiles:
   if not self.mbtiles_db:
    self.mbtiles_setup(1);
   if self.mbtiles_db:
    self.image_output.write_base_tile(tx, ty, tz, xyzzy)
    if self.image_output.tile_exists(tx, ty, tz):
     filename = self.image_output.get_full_path(tx, ty, tz, format_extension[self.image_output.format])
     input_file = open(filename, 'rb')
     if not input_file.closed:
      self.image_data = input_file.read()
      input_file.close()
      self.mbtiles_db.insert_image(tz,tx,ty,self.image_data)
      os.remove(filename)
   else:
    pass
    # print "write_base_tile: self.mbtiles_db None"
  else:
   self.image_output.write_base_tile(tx, ty, tz, xyzzy)
 # -------------------------------------------------------------------------
 def write_overview_tile(self, tx, ty, tz,tms_osm):
  if self.options.mbtiles:
   if not self.mbtiles_db:
    self.mbtiles_setup(1);
   if self.mbtiles_db:
    # retrieve the 4 images and write to disk
    image_list=self.mbtiles_db.retrieve_zoom_images(tz,tx,ty)
    self.image_output.write_overview_tile(tx, ty, tz,tms_osm)
    if self.image_output.tile_exists(tx, ty, tz):
     filename = self.image_output.get_full_path(tx, ty, tz, format_extension[self.image_output.format])
     input_file = open(filename, 'rb')
     if not input_file.closed:
      image_data = input_file.read()
      input_file.close()
      self.mbtiles_db.insert_image(tz,tx,ty,image_data)
      # exit();
      os.remove(filename)
      # delete the 4 images
      for image_file in image_list:
        os.remove(image_file)
   else:
    pass
    # print "write_overview_tile: self.mbtiles_db None"
  else:
   self.image_output.write_overview_tile(tx, ty, tz,tms_osm)
  # MBTiles support -end -
 # -------------------------------------------------------------------------
 def generate_overview_tiles(self):
  """Generation of the overview tiles (higher in the pyramid) based on existing tiles"""

  gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")

  print "Generating Overview Tiles:"

  if self.options.profile == 'garmin': # no overview tiles for 'garmin'
   return
  # Usage of existing tiles: from 4 underlying tiles generate one as overview.

  tcount = 0
  zcount = 0
  for tz in range(self.tmaxz-1, self.tminz-1, -1):
   tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
   tcount += (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
   zcount+=1
  if self.options.resume:
   count_tiles=tcount
   zcount+=1
   tminx, tminy, tmaxx, tmaxy = self.tminmax[self.tmaxz]
   count_tiles += (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
   i_count = self.tile_exists(0, 0, 0,1)
   if i_count == count_tiles:
    if self.options.verbose:
     print "\tTile generation skipped because of --resume ;  all-tiles [",zcount,"] zoom-levels with  tiles[",count_tiles,"]"
    return
  ti = 0

  # querysize = tilesize * 2

  for tz in range(self.tmaxz-1, self.tminz-1, -1):
   tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
   i_x_column_count=((tmaxx-tminx)+1)
   i_y_column_count=((tmaxy-tminy)+1)
   if self.options.verbose:
    #  tx in range(tminx, tmaxx+1) tminx[ 140798 ] tmaxx[ 140872 ] ; ((tmaxx-tmaxy)+1) x_tiles[ -35331 ]
    print "\ttz=[",tz,"] : tx in range(tminx, tmaxx+1) tminx[",tminx,"] tmaxx[",tmaxx,"] ; ((tmaxx-tminx)+1) x_tiles[",i_x_column_count,"]"
    # ty_tms in range(tmaxy, tminy-1, -1) tmaxy[ 176204 ] tminy[ 176126 ] ; ((tmaxy-tminy)) y_tiles[ 78 ]
    print "\ttz=[",tz,"] :ty_tms in range(tmaxy, tminy-1, -1) tmaxy[",tmaxy,"] tminy[",tminy,"] ; ((tmaxy-tminy)) y_tiles[",i_y_column_count,"]"
   if self.options.resume:
    i_count = self.tile_exists(0, 0, tz,2)
    print "\tTile generation skipped because of --??? ;  x/y-tiles of z[",tz,"]  x/y_tiles[",tcount,"] i_count[",i_count,"]"
    if i_count == tcount:
     if self.options.verbose:
      print "\tTile generation skipped because of --resume ;  x/y-tiles of z[",tz,"]  x/y_tiles[",tcount,"]"
     break
   for tx in range(tminx, tmaxx+1):
    tmaxy_work=tmaxy
    if self.options.resume:
     i_count = self.tile_exists(tx, 0, tz,3)
     print "\tTile generation skipped because of --??? ;  z =",tz," ; y-tiles of x[",tx,"]  y_tiles[",i_y_column_count,"] i_count[",i_count,"]"
     if i_count == i_y_column_count:
      if self.options.verbose:
       print "\tTile generation skipped because of --resume ;  z =",tz," ; y-tiles of x[",tx,"]  y_tiles[",i_y_column_count,"]"
      break
     else:
      if i_count > 0:
       # this assums the rows are compleate, which may NOT be true 18-140798-176204.jpg
       tmaxy_work-=i_count
       if self.options.verbose:
        print "\tTile generation skipped to tmaxy[",tmaxy_work,"] because of --resume ;  z =",tz," ; y-tiles of x[",tx,"]  y_tiles[",i_y_column_count,"]"
    for ty_tms in range(tmaxy_work, tminy-1, -1): #range(tminy, tmaxy+1):
     ty_osm=self.flip_y(tz,ty_tms)
     ty=ty_tms
     if self.options.tms_osm:
      ty=ty_osm
     if self.stopped:
      if self.options.mbtiles:
       if self.mbtiles_db:
        self.mbtiles_db.close_db()
        self.mbtiles_db=None
      break

     ti += 1

     if self.options.resume:
      exists = self.tile_exists(tx, ty, tz,0)
      if exists and self.options.verbose:
       print "\tTile generation skipped because of --resume"
     else:
      exists = False

     if not exists:
      if self.options.verbose:
       print ti, '/', tcount, self.get_verbose_tile_name(tx, ty, tz)
      try:
       self.write_overview_tile(tx, ty, tz,self.options.tms_osm)
      except ImageOutputException, e:
       self.error("'%d/%d/%d': %s" % (tz, tx, ty, e.message))

     if not self.options.verbose or self.is_subprocess:
      self.progressbar( ti / float(tcount) )
  if self.options.mbtiles:
   if self.mbtiles_db:
    self.mbtiles_db.close_db()
    self.mbtiles_db=None

 # -------------------------------------------------------------------------
 def generate_kml(self):
  if not self.kml:
   return

  # The KMZ specific to 'garmin' profile
  if self.options.profile == 'garmin':
   if os.path.basename( self.output ):
    zipname = os.path.basename( self.output ) + '.kmz'
   else:
    zipname = os.path.basename( self.output[:-1] + '.kmz' )
   f = ZipFile(os.path.join( self.output, zipname), "w", ZIP_DEFLATED)
   kml_tiles = {}
   children_kml = []
   tz = self.tmaxz
   tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
   for ty_tms in range(tminy, tmaxy+1):
    for tx in range(tminx, tmaxx+1):
     image_format = self.image_output.try_to_use_existing_tile(tx, ty_tms, tz)
     if image_format is None:
      continue
     filename = self.image_output.get_full_path(tx, ty_tms, tz, format_extension[image_format])
     f.write(filename, get_tile_filename(tx, ty_tms, tz, format_extension[image_format],False), ZIP_STORED)
     os.unlink(filename)

     d = self.get_kml_dict(tx, ty_tms, tz, image_format, draworder = 50)
     children_kml.append( self.generate_garmin_kml(d) )

   f.writestr("doc.kml", self.generate_document_kml(self.options.title, "".join(children_kml)))
   f.close()
   import shutil
   shutil.rmtree( os.path.join(self.output, str(tz)) )
   return

  # Base level KML.
  kml_tiles = {}
  tz = self.tmaxz
  tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
  for ty_tms in range(tminy, tmaxy+1):
   for tx in range(tminx, tmaxx+1):
    image_format = self.image_output.try_to_use_existing_tile(tx, ty_tms, tz)
    if image_format is None:
     continue

    d = self.get_kml_dict(tx, ty_tms, tz, image_format)

    if self.kml_depth == 1 or self.tmaxz == self.tminz:
     self.write_kml_tile(tx, ty_tms, tz, self.generate_node_kml(d, []))
     kml_tiles[tx,ty_tms,tz] = self.generate_link_kml(d)
    else:
     kml_tiles[tx,ty_tms,tz] = self.generate_leaf_kml(d)

  # Overviews KML.
  for tz in range(self.tmaxz-1, self.tminz-1, -1):
   tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
   for ty_tms in range(tminy, tmaxy+1):
    for tx in range(tminx, tmaxx+1):
     image_format = self.image_output.try_to_use_existing_tile(tx, ty_tms, tz)
     if image_format is None:
      continue

     d = self.get_kml_dict(tx, ty_tms, tz, image_format)

     children = [kml_tiles[x,y,tz+1]
        for y in range(2*ty_tms, 2*ty_tms + 2)
        for x in range(2*tx, 2*tx + 2)
        if (x,y,tz+1) in kml_tiles]

     node_kml = self.generate_node_kml(d, children)
     if (self.tmaxz-tz + 1) % self.kml_depth == 0 or tz == self.tminz:
      self.write_kml_tile(tx, ty_tms, tz, node_kml)
      kml_tiles[tx,ty_tms,tz] = self.generate_link_kml(d)
     else:
      kml_tiles[tx,ty_tms,tz] = node_kml

  # Root KML.
  tminx, tminy, tmaxx, tmaxy = self.tminmax[self.tminz]
  children_kml = [kml_tiles[x,y,self.tminz]
      for y in range(tminy, tmaxy+1)
      for x in range(tminx, tmaxx+1)
      if (x,y,self.tminz) in kml_tiles]

  lookat = computeFitLookAt(self.swne[1], self.swne[0], self.swne[3], self.swne[2])
  lookat_kml = self.generate_lookat_kml_block(lookat[0], lookat[1], lookat[2])
  self.write_kml_file("doc.kml", self.options.title, lookat_kml+"\n".join(children_kml))

 # -------------------------------------------------------------------------
 def geo_query(self, ds, ulx, uly, lrx, lry, querysize = 0):
  """For given dataset and query in cartographic coordinates
  returns parameters for ReadRaster() in raster coordinates and
  x/y shifts (for border tiles). If the querysize is not given, the
  extent is returned in the native resolution of dataset ds."""

  geotran = ds.GetGeoTransform()
  rx= int((ulx - geotran[0]) / geotran[1] + 0.001)
  ry= int((uly - geotran[3]) / geotran[5] + 0.001)
  rxsize= int((lrx - ulx) / geotran[1] + 0.5)
  rysize= int((lry - uly) / geotran[5] + 0.5)

  if not querysize:
   wxsize, wysize = rxsize, rysize
  else:
   wxsize, wysize = querysize, querysize

  # Coordinates should not go out of the bounds of the raster
  wx = 0
  if rx < 0:
   rxshift = abs(rx)
   wx = int( wxsize * (float(rxshift) / rxsize) )
   wxsize = wxsize - wx
   rxsize = rxsize - int( rxsize * (float(rxshift) / rxsize) )
   rx = 0
  if rx+rxsize > ds.RasterXSize:
   wxsize = int( wxsize * (float(ds.RasterXSize - rx) / rxsize) )
   rxsize = ds.RasterXSize - rx

  wy = 0
  if ry < 0:
   ryshift = abs(ry)
   wy = int( wysize * (float(ryshift) / rysize) )
   wysize = wysize - wy
   rysize = rysize - int( rysize * (float(ryshift) / rysize) )
   ry = 0
  if ry+rysize > ds.RasterYSize:
   wysize = int( wysize * (float(ds.RasterYSize - ry) / rysize) )
   rysize = ds.RasterYSize - ry

  return (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)

 # -------------------------------------------------------------------------
 def write_kml_tile(self, tx, ty_tms, tz, kml):
  if self.options.kmz:
   filename = get_tile_filename(tx, ty_tms, tz, "kmz",False)
   self.write_kmz_file(filename, filename, kml)
  else:
   filename = get_tile_filename(tx, ty_tms, tz, "kml",False)
   self.write_kml_file(filename, filename, kml)

 def write_kml_file(self, filename, title, content):
  f = open(os.path.join(self.output, filename), 'w')
  f.write(self.generate_document_kml(title, content))
  f.close()

 def write_kmz_file(self, filename, title, content):
  f = ZipFile(os.path.join(self.output, filename), "w", ZIP_DEFLATED)
  f.writestr("doc.kml", self.generate_document_kml(title, content))
  f.close()

 def generate_node_kml(self, d, children):
  """Return KML describing tile image and its children."""
  return self.generate_leaf_kml(d, "\n".join(children))

 def generate_garmin_kml(self, d ):
  """Return Garmin KML block describing an tile image."""
  return ("""
      <GroundOverlay>
        <Icon>
          <href>%(image_url)s</href>
          <DrawOrder>%(draw_order)d</DrawOrder>
        </Icon>
        <LatLonBox>
          <north>%(north).14f</north>
          <south>%(south).14f</south>
          <east>%(east).14f</east>
          <west>%(west).14f</west>
        </LatLonBox>
      </GroundOverlay>""" % d )

 def generate_leaf_kml(self, d, content=""):
  """Return KML describing tile image and insert content."""
  return ("""\
    <Folder>
      <Region>
        <Lod>
          <minLodPixels>%(minlodpixels)d</minLodPixels>
          <maxLodPixels>%(maxlodpixels)d</maxLodPixels>
        </Lod>
        <LatLonAltBox>
          <north>%(north).14f</north>
          <south>%(south).14f</south>
          <east>%(east).14f</east>
          <west>%(west).14f</west>
        </LatLonAltBox>
      </Region>
      <GroundOverlay>
        <drawOrder>%(draw_order)d</drawOrder>
        <Icon>
          <href>%(image_url)s</href>
        </Icon>
        <LatLonBox>
          <north>%(north).14f</north>
          <south>%(south).14f</south>
          <east>%(east).14f</east>
          <west>%(west).14f</west>
        </LatLonBox>
      </GroundOverlay>""" % d
   + """\
%s
    </Folder>""" % content)

 def generate_link_kml(self, d):
  """Return KML linking to the tile."""
  return """\
    <NetworkLink>
      <name>%(image_filename)s</name>
      <Region>
        <Lod>
          <minLodPixels>%(minlodpixels)d</minLodPixels>
          <maxLodPixels>-1</maxLodPixels>
        </Lod>
        <LatLonAltBox>
          <north>%(north).14f</north>
          <south>%(south).14f</south>
          <east>%(east).14f</east>
          <west>%(west).14f</west>
        </LatLonAltBox>
      </Region>
      <Link>
        <href>%(link_url)s</href>
        <viewRefreshMode>onRegion</viewRefreshMode>
      </Link>
    </NetworkLink>""" % d

 def generate_document_kml(self, title, content):
  """Return full KML document with given title and content."""
  return """\
<?xml version="1.0" encoding="utf-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <name>%s</name>
    <description></description>
    <Style>
      <ListStyle id="hideChildren">
        <listItemType>checkHideChildren</listItemType>
      </ListStyle>
    </Style>
%s
  </Document>
</kml>""" % (title.replace('\\','/'), content)

 def generate_lookat_kml_block(self, lng, lat, viewrange):
  """Return the KML string containing correct <LookAt> block"""
  return """
    <LookAt>
      <longitude>%.14f</longitude>
      <latitude>%.14f</latitude>
      <altitude>0</altitude>
      <range>%.f</range>
      <tilt>0</tilt>
      <heading>0</heading>
    </LookAt>
""" % (lng, lat, viewrange)

 def get_kml_dict(self, tx, ty_tms, tz, image_format, draworder = 0):
  """Return dictionary describing KML info about tile to be used in templates."""
  d = {}

  d["south"], d["west"], d["north"], d["east"] = self.tileswne(tx, ty_tms, tz)

  image_filename = get_tile_filename(tx, ty_tms, tz, format_extension[image_format],False)
  d["image_filename"] = image_filename
  d["image_filename"] = d["image_filename"].replace("\\","/")

  if self.options.url is None:
   d["image_url"] = "../../%s" % image_filename
  else:
   d["image_url"] = "%s%s" % (self.options.url, image_filename)
  d["image_url"] = d["image_url"].replace("\\","/")

  url = self.options.url
  if url is None:
   # Top level KML is linked from `doc.kml' and it needs different path.
   if tz == self.tminz:
    url = ""
   else:
    url = "../../"

  if self.options.kmz:
   extension = "kmz"
  else:
   extension = "kml"

  d["link_url"] = "%s%s" % (url, get_tile_filename(tx, ty_tms, tz, extension,False))
  d["link_url"] = d["link_url"].replace("\\","/")

  d["minlodpixels"] = int(self.tilesize / 2)
  d["maxlodpixels"] = -1 # int(self.tilesize * 8)

  if tx == 0:
   d["draw_order"] = draworder + 2 * tz + 1
  else:
   d["draw_order"] = draworder + 2 * tz

  return d

 def get_verbose_tile_name(self, tx, ty, tz):
  if self.options.tile_format == "hybrid":
   extension = "?"
  else:
   extension = format_extension[self.image_output.format]

  return self.image_output.get_full_path(tx, ty, tz, extension)

 # -------------------------------------------------------------------------
 def generate_tilemapresource(self):
  """
     Template for tilemapresource.xml. Returns filled string. Expected variables:
       title, north, south, east, west, isepsg4326, projection, publishurl,
       zoompixels, tilesize, tileformat, profile
       http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
  """

  args = {}
  args['title'] = self.options.title
  args['south'], args['west'], args['north'], args['east'] = self.swne
  args['tilesize'] = self.tilesize
  args['tileformat'] = format_extension[self.image_output.format]
  args['mime'] = format_mime[self.image_output.format]
  args['publishurl'] = "" if self.options.url is None else self.options.url
  args['profile'] = self.options.profile

  if self.options.profile == 'mercator':
   args['srs'] = "EPSG:900913"
  elif self.options.profile in ('geodetic', 'gearth'):
   args['srs'] = "EPSG:4326"
  elif self.options.s_srs:
   args['srs'] = self.options.s_srs
  elif self.out_srs:
   args['srs'] = self.out_srs.ExportToWkt()
  else:
   args['srs'] = ""

  s = """<?xml version="1.0" encoding="utf-8"?>
 <TileMap version="1.0.0" tilemapservice="http://tms.osgeo.org/1.0.0">
   <Title>%(title)s</Title>
   <Abstract></Abstract>
   <SRS>%(srs)s</SRS>
   <BoundingBox minx="%(south).14f" miny="%(west).14f" maxx="%(north).14f" maxy="%(east).14f"/>
   <Origin x="%(west).14f" y="%(south).14f"/>
   <TileFormat width="%(tilesize)d" height="%(tilesize)d" mime-type="%(mime)s" extension="%(tileformat)s"/>
   <TileSets profile="%(profile)s">
""" % args
  for z in range(self.tminz, self.tmaxz+1):
   if self.options.profile == 'raster':
    s += """     <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, (2**(self.nativezoom-z) * self.out_gt[1]), z)
   elif self.options.profile == 'mercator':
    s += """     <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, 156543.0339/2**z, z)
   elif self.options.profile == 'geodetic':
    s += """     <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, 0.703125/2**z, z)
  s += """   </TileSets>
 </TileMap>
 """
  return s

 # -------------------------------------------------------------------------
 def generate_googlemaps(self):
  """
  Template for googlemaps.html implementing Overlay of tiles for 'mercator' profile.
  It returns filled string. Expected variables:
  title, googlemapskey, north, south, east, west, minzoom, maxzoom, tilesize, tileformat, publishurl
  """
  args = {}
  args['title'] = self.options.title
  args['googlemapskey'] = self.options.googlekey
  args['south'], args['west'], args['north'], args['east'] = self.swne
  args['minzoom'] = self.tminz
  args['maxzoom'] = self.tmaxz
  args['tilesize'] = self.tilesize
  args['tileformat'] = format_extension[self.image_output.format]
  args['publishurl'] = "" if self.options.url is None else self.options.url
  args['copyright'] = self.options.copyright

  s = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
   <html xmlns="http://www.w3.org/1999/xhtml" xmlns:v="urn:schemas-microsoft-com:vml">
     <head>
       <title>%(title)s</title>
       <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
       <meta http-equiv='imagetoolbar' content='no'/>
       <style type="text/css"> v\:* {behavior:url(#default#VML);}
           html, body { overflow: hidden; padding: 0; height: 100%%; width: 100%%; font-family: 'Lucida Grande',Geneva,Arial,Verdana,sans-serif; }
           body { margin: 10px; background: #fff; }
           h1 { margin: 0; padding: 6px; border:0; font-size: 20pt; }
           #header { height: 43px; padding: 0; background-color: #eee; border: 1px solid #888; }
           #subheader { height: 12px; text-align: right; font-size: 10px; color: #555;}
           #map { height: 95%%; border: 1px solid #888; }
       </style>
       <script src='http://maps.google.com/maps?file=api&amp;v=2&amp;key=%(googlemapskey)s' type='text/javascript'></script>
       <script type="text/javascript">
       //<![CDATA[

       /*
        * Constants for given map
        * TODO: read it from tilemapresource.xml
        */

       var mapBounds = new GLatLngBounds(new GLatLng(%(south)s, %(west)s), new GLatLng(%(north)s, %(east)s));
       var mapMinZoom = %(minzoom)s;
       var mapMaxZoom = %(maxzoom)s;

       var opacity = 0.75;
       var map;
       var ge;
       var hybridOverlay;

       /*
        * Create a Custom Opacity GControl
        * https://github.com/mj10777/mapmbtilesgoogle-maps-overlay-opacity-control/
        */

       var CTransparencyLENGTH = 58;
       // maximum width that the knob can move (slide width minus knob width)

       function CTransparencyControl( overlay ) {
           this.overlay = overlay;
           this.opacity = overlay.getTileLayer().getOpacity();
       }
       CTransparencyControl.prototype = new GControl();

       // This function positions the slider to match the specified opacity
       CTransparencyControl.prototype.setSlider = function(pos) {
           var left = Math.round((CTransparencyLENGTH*pos));
           this.slide.left = left;
           this.knob.style.left = left+"px";
           this.knob.style.top = "0px";
       }

       // This function reads the slider and sets the overlay opacity level
       CTransparencyControl.prototype.setOpacity = function() {
        // set the global variable
           opacity = this.slide.left/CTransparencyLENGTH;
           this.map.clearOverlays();
           this.map.addOverlay(this.overlay, { zPriority: 0 });
           if (this.map.getCurrentMapType() == G_HYBRID_MAP) {
               this.map.addOverlay(hybridOverlay);
           }
       }

       // This gets called by the API when addControl(new CTransparencyControl())
       CTransparencyControl.prototype.initialize = function(map) {
           var that=this;
           this.map = map;

           // Is this MSIE, if so we need to use AlphaImageLoader
           var agent = navigator.userAgent.toLowerCase();
           if ((agent.indexOf("msie") > -1) && (agent.indexOf("opera") < 1)){this.ie = true} else {this.ie = false}

           // create the background graphic as a <div> containing an image
           var container = document.createElement("div");
           container.style.width="70px";
           container.style.height="21px";

           // Handle transparent PNG files in MSIE
           if (this.ie) {
             var loader = "filter:progid:DXImageTransform.Microsoft.AlphaImageLoader(src='https://github.com/mj10777/mapmbtiles/img/opacity-slider.png', sizingMethod='crop');";
             container.innerHTML = '<div style="height:21px; width:70px; ' +loader+ '" ></div>';
           } else {
             container.innerHTML = '<div style="height:21px; width:70px; background-image: url(https://github.com/mj10777/mapmbtiles/img/opacity-slider.png)" ></div>';
           }

           // create the knob as a GDraggableObject
           // Handle transparent PNG files in MSIE
           if (this.ie) {
             var loader = "progid:DXImageTransform.Microsoft.AlphaImageLoader(src='https://github.com/mj10777/mapmbtiles/img/opacity-slider.png', sizingMethod='crop');";
             this.knob = document.createElement("div");
             this.knob.style.height="21px";
             this.knob.style.width="13px";
      this.knob.style.overflow="hidden";
             this.knob_img = document.createElement("div");
             this.knob_img.style.height="21px";
             this.knob_img.style.width="83px";
             this.knob_img.style.filter=loader;
      this.knob_img.style.position="relative";
      this.knob_img.style.left="-70px";
             this.knob.appendChild(this.knob_img);
           } else {
             this.knob = document.createElement("div");
             this.knob.style.height="21px";
             this.knob.style.width="13px";
             this.knob.style.backgroundImage="url(https://github.com/mj10777/mapmbtiles/img/opacity-slider.png)";
             this.knob.style.backgroundPosition="-70px 0px";
           }
           container.appendChild(this.knob);
           this.slide=new GDraggableObject(this.knob, {container:container});
           this.slide.setDraggableCursor('pointer');
           this.slide.setDraggingCursor('pointer');
           this.container = container;

           // attach the control to the map
           map.getContainer().appendChild(container);

           // init slider
           this.setSlider(this.opacity);

           // Listen for the slider being moved and set the opacity
           GEvent.addListener(this.slide, "dragend", function() {that.setOpacity()});
           //GEvent.addListener(this.container, "click", function( x, y ) { alert(x, y) });

           return container;
         }

         // Set the default position for the control
         CTransparencyControl.prototype.getDefaultPosition = function() {
           return new GControlPosition(G_ANCHOR_TOP_RIGHT, new GSize(7, 47));
         }

       /*
        * Full-screen Window Resize
        */

       function getWindowHeight() {
           if (self.innerHeight) return self.innerHeight;
           if (document.documentElement && document.documentElement.clientHeight)
               return document.documentElement.clientHeight;
           if (document.body) return document.body.clientHeight;
           return 0;
       }

       function getWindowWidth() {
           if (self.innerWidth) return self.innerWidth;
           if (document.documentElement && document.documentElement.clientWidth)
               return document.documentElement.clientWidth;
           if (document.body) return document.body.clientWidth;
           return 0;
       }

       function resize() {
           var map = document.getElementById("map");
           var header = document.getElementById("header");
           var subheader = document.getElementById("subheader");
           map.style.height = (getWindowHeight()-80) + "px";
           map.style.width = (getWindowWidth()-20) + "px";
           header.style.width = (getWindowWidth()-20) + "px";
           subheader.style.width = (getWindowWidth()-20) + "px";
           // map.checkResize();
       }


       /*
        * Main load function:
        */

       function load() {

          if (GBrowserIsCompatible()) {

             // Bug in the Google Maps: Copyright for Overlay is not correctly displayed
             var gcr = GMapType.prototype.getCopyrights;
             GMapType.prototype.getCopyrights = function(bounds,zoom) {
                 return ["%(copyright)s"].concat(gcr.call(this,bounds,zoom));
             }

             map = new GMap2( document.getElementById("map"), { backgroundColor: '#fff' } );

             map.addMapType(G_PHYSICAL_MAP);
             map.setMapType(G_PHYSICAL_MAP);

             map.setCenter( mapBounds.getCenter(), map.getBoundsZoomLevel( mapBounds ));

             hybridOverlay = new GTileLayerOverlay( G_HYBRID_MAP.getTileLayers()[1] );
             GEvent.addListener(map, "maptypechanged", function() {
               if (map.getCurrentMapType() == G_HYBRID_MAP) {
                   map.addOverlay(hybridOverlay);""" % args
  if self.kml:
   s += """
               } else if (map.getCurrentMapType() == G_SATELLITE_3D_MAP) {
                   var url = document.location.toString();
                   if (url.substr(0,4) != 'http') alert('You have to upload the tiles to a webserver to see the overlay in Google Earth Plugin');
                   if (!ge) map.getEarthInstance(getEarthInstanceCB);"""
  s += """
               } else {
                  map.removeOverlay(hybridOverlay);
               }
             } );

             var tilelayer = new GTileLayer(GCopyrightCollection(''), mapMinZoom, mapMaxZoom);
             var mercator = new GMercatorProjection(mapMaxZoom+1);
             tilelayer.getTileUrl = function(tile,zoom) {
                 if ((zoom < mapMinZoom) || (zoom > mapMaxZoom)) {
                     return "https://github.com/mj10777/mapmbtiles/img/none.png";
                 }
                 var ymax = 1 << zoom;
                 var y = ymax - tile.y -1;
                 var tileBounds = new GLatLngBounds(
                     mercator.fromPixelToLatLng( new GPoint( (tile.x)*256, (tile.y+1)*256 ) , zoom ),
                     mercator.fromPixelToLatLng( new GPoint( (tile.x+1)*256, (tile.y)*256 ) , zoom )
                 );
                 if (mapBounds.intersects(tileBounds)) {
                     return zoom+"/"+tile.x+"/"+y+".%(tileformat)s";
                 } else {
                     return "https://github.com/mj10777/mapmbtiles/img/none.png";
                 }
             }
             // IE 7-: support for PNG alpha channel
             // Unfortunately, the opacity for whole overlay is then not changeable, either or...
             tilelayer.isPng = function() { return true;};
             tilelayer.getOpacity = function() { return opacity; }

             overlay = new GTileLayerOverlay( tilelayer );
             map.addOverlay(overlay);

             map.addControl(new GLargeMapControl3D());
             map.addControl(new GHierarchicalMapTypeControl());
             map.addControl(new CTransparencyControl( overlay ));
  """ % args
  if self.kml:
   s += """
             map.addMapType(G_SATELLITE_3D_MAP);
  """
  s += """

             map.enableContinuousZoom();
             map.enableScrollWheelZoom();

             map.setMapType(G_HYBRID_MAP);
          }
          resize();
       }
  """
  if self.kml:
   s += """
       function getEarthInstanceCB(object) {
          ge = object;
          var url = document.location.toString();
          var newurl = url.substr(0,url.lastIndexOf('/'))+'/doc.kml';
          if (ge) {
              var link = ge.createLink("");
              if ("%(publishurl)s") { link.setHref("%(publishurl)s/doc.kml") }
              else { link.setHref(newurl) };
              var networkLink = ge.createNetworkLink("");
              networkLink.set(link, false, false);
              ge.getFeatures().appendChild(networkLink);
          } else {
              // alert("Initialization of the Google Earth Plugin failed. You can still open the KML file in normal Google Earth.");
              // window.location = newurl; // JavaScript redirect to the URL of KML
          }
       }
  """ % args
  s += """
       onresize=function(){ resize(); };

       //]]>
       </script>
     </head>
     <body onload="load()">
         <div id="header"><h1>%(title)s</h1></div>
         <div id="subheader">Generated by <a href="https://github.com/mj10777/mapmbtiles">MapMbTiles</a>/<a href="http://www.klokan.cz/projects/gdal2mbtiles/">GDAL2MbTiles</a>, Copyright &copy; 2008 <a href="http://www.klokan.cz/">Klokan Petr Pridal</a>,  <a href="http://www.gdal.org/">GDAL</a> &amp; <a href="http://www.osgeo.org/">OSGeo</a> <a href="http://code.google.com/soc/">GSoC</a>
   <!-- PLEASE, LET THIS NOTE ABOUT AUTHOR AND PROJECT SOMEWHERE ON YOUR WEBSITE, OR AT LEAST IN THE COMMENT IN HTML. THANK YOU -->
         </div>
          <div id="map"></div>
     </body>
   </html>
  """ % args

  return s


 # -------------------------------------------------------------------------
 def generate_openlayers( self ):
  """
  Template for openlayers.html implementing overlay of available Spherical Mercator layers.

  It returns filled string. Expected variables:
  title, googlemapskey, yahooappid, north, south, east, west, minzoom, maxzoom, tilesize, tileformat, publishurl
  """

  args = {}
  args['title'] = self.options.title
  args['googlemapskey'] = self.options.googlekey
  args['yahooappid'] = self.options.yahookey
  args['south'], args['west'], args['north'], args['east'] = self.swne
  args['minzoom'] = self.tminz
  args['maxzoom'] = self.tmaxz
  args['tilesize'] = self.tilesize
  args['tileformat'] = format_extension[self.image_output.format]
  if self.image_output.format == "PNG":
   args['has_alpha'] = 'true'
  else:
   args['has_alpha'] = 'false'
  args['publishurl'] = "" if self.options.url is None else self.options.url
  args['copyright'] = self.options.copyright
  if self.options.profile in ('raster', 'gearth'):
   args['rasterzoomlevels'] = self.tmaxz+1
   args['rastermaxresolution'] = 2**(self.nativezoom) * self.out_gt[1]

  s = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
  <html xmlns="http://www.w3.org/1999/xhtml>"
    <head>
      <title>%(title)s</title>
      <meta http-equiv='imagetoolbar' content='no'/>
      <style type="text/css"> v\:* {behavior:url(#default#VML);}
          html, body { overflow: hidden; padding: 0; height: 100%%; width: 100%%; font-family: 'Lucida Grande',Geneva,Arial,Verdana,sans-serif; }
          body { margin: 10px; background: #fff; }
          h1 { margin: 0; padding: 6px; border:0; font-size: 20pt; }
          #header { height: 43px; padding: 0; background-color: #eee; border: 1px solid #888; }
          #subheader { height: 12px; text-align: right; font-size: 10px; color: #555;}
          #map { height: 95%%; border: 1px solid #888; }
      </style>""" % args

  if self.options.profile == 'mercator':
   s += """
      <script src='http://dev.virtualearth.net/mapcontrol/mapcontrol.ashx?v=6.1'></script>
      <script src='http://maps.google.com/maps?file=api&amp;v=2&amp;key=%(googlemapskey)s' type='text/javascript'></script>
      <script src="http://api.maps.yahoo.com/ajaxymap?v=3.0&amp;appid=%(yahooappid)s"></script>""" % args

  s += """
      <script src="http://www.openlayers.org/api/2.7/OpenLayers.js" type="text/javascript"></script>
      <script type="text/javascript">
          var map;
       var mapBounds = new OpenLayers.Bounds( %(west)s, %(south)s, %(east)s, %(north)s);
       var mapMinZoom = %(minzoom)s;
       var mapMaxZoom = %(maxzoom)s;

          // avoid pink tiles
          OpenLayers.IMAGE_RELOAD_ATTEMPTS = 3;
          OpenLayers.Util.onImageLoadErrorColor = "transparent";

          function init(){""" % args

  if self.options.profile == 'mercator':
   s += """
             var options = {
                 controls: [],
                 projection: new OpenLayers.Projection("EPSG:900913"),
                 displayProjection: new OpenLayers.Projection("EPSG:4326"),
                 units: "m",
                 maxResolution: 156543.0339,
                 maxExtent: new OpenLayers.Bounds(-20037508, -20037508, 20037508, 20037508.34)
              };
             map = new OpenLayers.Map('map', options);

             // create Google Mercator layers
             var gmap = new OpenLayers.Layer.Google("Google Streets",
     { sphericalMercator: true, numZoomLevels: 20} );
             var gsat = new OpenLayers.Layer.Google("Google Satellite",
     {type: G_SATELLITE_MAP, sphericalMercator: true, numZoomLevels: 20} );
             var ghyb = new OpenLayers.Layer.Google("Google Hybrid",
     {type: G_HYBRID_MAP, sphericalMercator: true, numZoomLevels: 20});
             var gter = new OpenLayers.Layer.Google("Google Terrain",
     {type: G_PHYSICAL_MAP, sphericalMercator: true, numZoomLevels: 20 });

             // create Virtual Earth layers
    OpenLayers.Layer.VirtualEarth.prototype.MAX_ZOOM_LEVEL=19;
    OpenLayers.Layer.VirtualEarth.prototype.RESOLUTIONS=OpenLayers.Layer.Google.prototype.RESOLUTIONS
             var veroad = new OpenLayers.Layer.VirtualEarth("Virtual Earth Roads",
     {'type': VEMapStyle.Road, 'sphericalMercator': true, numZoomLevels: 20});
             var veaer = new OpenLayers.Layer.VirtualEarth("Virtual Earth Aerial",
     {'type': VEMapStyle.Aerial, 'sphericalMercator': true, numZoomLevels: 20 });
             var vehyb = new OpenLayers.Layer.VirtualEarth("Virtual Earth Hybrid",
                 {'type': VEMapStyle.Hybrid, 'sphericalMercator': true});

             // create Yahoo layer
             var yahoo = new OpenLayers.Layer.Yahoo("Yahoo Street",
                 {'sphericalMercator': true});
             var yahoosat = new OpenLayers.Layer.Yahoo("Yahoo Satellite",
                 {'type': YAHOO_MAP_SAT, 'sphericalMercator': true});
             var yahoohyb = new OpenLayers.Layer.Yahoo("Yahoo Hybrid",
                 {'type': YAHOO_MAP_HYB, 'sphericalMercator': true});

             // create OSM/OAM layer
             var osm = new OpenLayers.Layer.TMS( "OpenStreetMap",
                 "http://tile.openstreetmap.org/",
                 { type: 'png', getURL: osm_getTileURL, displayOutsideMaxExtent: true,
       attribution: '<a href="http://www.openstreetmap.org/">OpenStreetMap</a>'} );
             var oam = new OpenLayers.Layer.TMS( "OpenAerialMap",
                 "http://tile.openaerialmap.org/tiles/1.0.0/openaerialmap-900913/",
                 { type: 'png', getURL: osm_getTileURL } );

             // create TMS Overlay layer
             var tmsoverlay = new OpenLayers.Layer.TMS( "TMS Overlay", "",
                 {   // url: '', serviceVersion: '.', layername: '.',
      type: '%(tileformat)s', getURL: overlay_getTileURL, alpha: %(has_alpha)s,
      isBaseLayer: false
                 });
    if (OpenLayers.Util.alphaHack() == false) { tmsoverlay.setOpacity(0.7); }

             map.addLayers([gmap, gsat, ghyb, gter, veroad, veaer, vehyb,
                            yahoo, yahoosat, yahoohyb, osm, oam,
                            tmsoverlay]);

             var switcherControl = new OpenLayers.Control.LayerSwitcher();
             map.addControl(switcherControl);
             switcherControl.maximizeControl();

             map.zoomToExtent( mapBounds.transform(map.displayProjection, map.projection ) );
   """ % args

  elif self.options.profile == 'geodetic':
   s += """
             var options = {
                 controls: [],
              projection: new OpenLayers.Projection("EPSG:4326"),
              maxResolution: 0.703125,
              maxExtent: new OpenLayers.Bounds(-180, -90, 180, 90)
              };
             map = new OpenLayers.Map('map', options);

             layer = new OpenLayers.Layer.WMS( "Blue Marble",
                     "http://labs.metacarta.com/wms-c/Basic.py?", {layers: 'satellite' } );
             map.addLayer(layer);
             wms = new OpenLayers.Layer.WMS( "VMap0",
                     "http://labs.metacarta.com/wms-c/Basic.py?", {layers: 'basic', format: 'image/png' } );
             map.addLayer(wms);

             var tmsoverlay = new OpenLayers.Layer.TMS( "TMS Overlay", "",
                 {
                     serviceVersion: '.', layername: '.', alpha: %(has_alpha)s,
      type: '%(tileformat)s', getURL: overlay_getTileURL,
      isBaseLayer: false
                 });
             map.addLayer(tmsoverlay);
    if (OpenLayers.Util.alphaHack() == false) { tmsoverlay.setOpacity(0.7); }

             var switcherControl = new OpenLayers.Control.LayerSwitcher();
             map.addControl(switcherControl);
             switcherControl.maximizeControl();

             map.zoomToExtent( mapBounds );
   """ % args

  elif self.options.profile in ('raster', 'gearth'):
   s += """
             var options = {
                 controls: [],
              maxExtent: new OpenLayers.Bounds(  %(west)s, %(south)s, %(east)s, %(north)s ),
              maxResolution: %(rastermaxresolution)f,
              numZoomLevels: %(rasterzoomlevels)d
              };
             map = new OpenLayers.Map('map', options);

          var layer = new OpenLayers.Layer.TMS( "TMS Layer","",
              {  url: '', serviceVersion: '.', layername: '.', alpha: %(has_alpha)s,
      type: '%(tileformat)s', getURL: overlay_getTileURL
     });
          map.addLayer(layer);
    map.zoomToExtent( mapBounds );
  """ % args


  s += """
             map.addControl(new OpenLayers.Control.PanZoomBar());
             map.addControl(new OpenLayers.Control.MousePosition());
             map.addControl(new OpenLayers.Control.MouseDefaults());
             map.addControl(new OpenLayers.Control.KeyboardDefaults());
         }
   """ % args

  if self.options.profile == 'mercator':
   s += """
         function osm_getTileURL(bounds) {
             var res = this.map.getResolution();
             var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
             var y = Math.round((this.maxExtent.top - bounds.top) / (res * this.tileSize.h));
             var z = this.map.getZoom();
             var limit = Math.pow(2, z);

             if (y < 0 || y >= limit) {
                 return "https://github.com/mj10777/mapmbtiles/img/none.png";
             } else {
                 x = ((x %% limit) + limit) %% limit;
                 return this.url + z + "/" + x + "/" + y + "." + this.type;
             }
         }

         function overlay_getTileURL(bounds) {
             var res = this.map.getResolution();
             var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
             var y = Math.round((bounds.bottom - this.tileOrigin.lat) / (res * this.tileSize.h));
             var z = this.map.getZoom();
             if (this.map.baseLayer.name == 'Virtual Earth Roads' || this.map.baseLayer.name == 'Virtual Earth Aerial' || this.map.baseLayer.name == 'Virtual Earth Hybrid') {
                z = z + 1;
             }
          if (mapBounds.intersectsBounds( bounds ) && z >= mapMinZoom && z <= mapMaxZoom ) {
                //console.log( this.url + z + "/" + x + "/" + y + "." + this.type);
                return this.url + z + "/" + x + "/" + y + "." + this.type;
                } else {
                   return "https://github.com/mj10777/mapmbtiles/img/none.png";
                }
         }
   """ % args

  elif self.options.profile == 'geodetic':
   s += """
         function overlay_getTileURL(bounds) {
    bounds = this.adjustBounds(bounds);
             var res = this.map.getResolution();
             var x = Math.round((bounds.left - this.tileOrigin.lon) / (res * this.tileSize.w));
             var y = Math.round((bounds.bottom - this.tileOrigin.lat) / (res * this.tileSize.h));
             var z = this.map.getZoom();
    var path = this.serviceVersion + "/" + this.layername + "/" + z + "/" + x + "/" + y + "." + this.type;
    var url = this.url;
          if (mapBounds.intersectsBounds( bounds ) && z >= mapMinZoom && z <= mapMaxZoom) {
                // console.log( this.url + z + "/" + x + "/" + y + "." + this.type);
                return this.url + z + "/" + x + "/" + y + "." + this.type;
                } else {
                   return "https://github.com/mj10777/mapmbtiles/img/none.png";
                }
         }
   """ % args

  elif self.options.profile in ('raster','gearth'):
   s += """
         function overlay_getTileURL(bounds) {
             var res = this.map.getResolution();
             var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
             var y = Math.round((bounds.bottom - this.maxExtent.bottom) / (res * this.tileSize.h));
             var z = this.map.getZoom();
    if (x >= 0 && y >= 0) {
              return this.url + z + "/" + x + "/" + y + "." + this.type;
    } else {
                 return "https://github.com/mj10777/mapmbtiles/img/none.png";
    }
   }
   """ % args

  s += """
     function getWindowHeight() {
          if (self.innerHeight) return self.innerHeight;
          if (document.documentElement && document.documentElement.clientHeight)
              return document.documentElement.clientHeight;
          if (document.body) return document.body.clientHeight;
           return 0;
      }

      function getWindowWidth() {
       if (self.innerWidth) return self.innerWidth;
       if (document.documentElement && document.documentElement.clientWidth)
           return document.documentElement.clientWidth;
       if (document.body) return document.body.clientWidth;
           return 0;
      }

      function resize() {
       var map = document.getElementById("map");
       var header = document.getElementById("header");
       var subheader = document.getElementById("subheader");
       map.style.height = (getWindowHeight()-80) + "px";
       map.style.width = (getWindowWidth()-20) + "px";
       header.style.width = (getWindowWidth()-20) + "px";
       subheader.style.width = (getWindowWidth()-20) + "px";
    if (map.updateSize) { map.updateSize(); };
      }

      onresize=function(){ resize(); };

      </script>
    </head>
    <body onload="init()">
   <div id="header"><h1>%(title)s</h1></div>
   <div id="subheader">Generated by <a href="https://github.com/mj10777/mapmbtiles">MapMbTiles</a>/<a href="http://www.klokan.cz/projects/gdal2mbtiles/">GDAL2MbTiles</a>, Copyright &copy; 2008 <a href="http://www.klokan.cz/">Klokan Petr Pridal</a>,  <a href="http://www.gdal.org/">GDAL</a> &amp; <a href="http://www.osgeo.org/">OSGeo</a> <a href="http://code.google.com/soc/">GSoC</a>
   <!-- PLEASE, LET THIS NOTE ABOUT AUTHOR AND PROJECT SOMEWHERE ON YOUR WEBSITE, OR AT LEAST IN THE COMMENT IN HTML. THANK YOU -->
   </div>
      <div id="map"></div>
      <script type="text/javascript" >resize()</script>
    </body>
  </html>""" % args

  return s

# =============================================================================
# =============================================================================
# =============================================================================


def ImageOutput(name, out_ds, tile_size, resampling, init_dest, output_dir, verbose,mbtiles):

 """Return object representing tile image output implementing given parameters."""

 resampler = Resampler(resampling)

 if name == "hybrid":
  return HybridImageOutput(out_ds, tile_size, resampler, init_dest, output_dir, verbose)

 if name == "png":
  image_format = "PNG"
 elif name == "jpeg":
  image_format = "JPEG"

 return SimpleImageOutput(out_ds, tile_size, resampler, init_dest, output_dir, verbose, [image_format],mbtiles)


class ImageOutputException(Exception):

 """Raised when the tile image can't be saved to disk."""


class BaseImageOutput(object):

 """Base class for image output.

 Child classes are supposed to provide two methods `write_base_tile' and
 `write_overview_tile'. These will call `create_base_tile' and `create_overview_tile'
 with arguments appropriate to their output strategy.

 When this class is instantiated with only one image format, it is stored in
 a member field `format'.
 """

 def __init__(self, out_ds, tile_size, resampler, init_dest, output_dir, verbose, image_formats,mbtiles):
  self.out_ds = out_ds
  self.tile_size = tile_size
  self.resampler = resampler
  self.init_dest = init_dest
  self.output_dir = output_dir
  self.verbose = verbose
  self.image_formats = image_formats
  self.mbtiles=mbtiles;
  if len(self.image_formats) == 1:
   self.format = self.image_formats[0]

  self.mem_drv = get_gdal_driver("MEM")
  self.alpha = None
  self.alpha_filler = None

  # Get alpha band (either directly or from NODATA value)
  self.alpha_band = self.out_ds.GetRasterBand(1).GetMaskBand()
  if (self.alpha_band.GetMaskFlags() & gdal.GMF_ALPHA) or self.out_ds.RasterCount in (2, 4):
   # TODO: Better test for alpha band in the dataset
   self.data_bands_count = self.out_ds.RasterCount - 1
  else:
   self.data_bands_count = self.out_ds.RasterCount

 def write_base_tile(self, tx, ty, tz, xyzzy):

  """Create image of a base level tile and write it to disk."""

  data_bands = range(1, self.data_bands_count+1)
  data = self.out_ds.ReadRaster(xyzzy.rx, xyzzy.ry, xyzzy.rxsize, xyzzy.rysize,
           xyzzy.wxsize, xyzzy.wysize, band_list=data_bands)

  image_format = self.get_base_tile_format(tx, ty, tz, xyzzy)

  if image_format is None:
   return
  else:
   num_bands = self.get_num_bands(image_format)

  if self.verbose:
   print "\tReadRaster Extent: ", (xyzzy.rx, xyzzy.ry, xyzzy.rxsize, xyzzy.rysize),
   print  'z =',tz,' ; x =',tx,' ; y =',ty, (xyzzy.wx, xyzzy.wy, xyzzy.wxsize, xyzzy.wysize)

  dstile = self.mem_drv.Create('', self.tile_size, self.tile_size, num_bands)

  path = self.get_full_path(tx, ty, tz, format_extension[image_format])

  # Query is in 'nearest neighbour' but can be bigger in then the tilesize
  # We scale down the query to the tilesize by supplied algorithm.
  if self.tile_size == xyzzy.querysize:
   self.fill_init_dest(dstile)

   # Use the ReadRaster result directly in tiles ('nearest neighbour' query)
   dstile.WriteRaster(xyzzy.wx, xyzzy.wy, xyzzy.wxsize, xyzzy.wysize, data, band_list=data_bands)
   if image_format == "PNG":
    dstile.WriteRaster(xyzzy.wx, xyzzy.wy, xyzzy.wxsize, xyzzy.wysize, self.alpha, band_list=[num_bands])

   gdal_write(path, dstile, image_format)

   # Note: For source drivers based on WaveLet compression (JPEG2000, ECW, MrSID)
   # the ReadRaster function returns high-quality raster (not ugly nearest neighbour)
   # TODO: Use directly 'near' for WaveLet files
  else:
   # Big ReadRaster query in memory scaled to the tilesize - all but 'near' algo
   dsquery = self.mem_drv.Create('', xyzzy.querysize, xyzzy.querysize, num_bands)
   self.fill_init_dest(dsquery)

   dsquery.WriteRaster(xyzzy.wx, xyzzy.wy, xyzzy.wxsize, xyzzy.wysize, data, band_list=data_bands)
   if image_format == "PNG":
    dsquery.WriteRaster(xyzzy.wx, xyzzy.wy, xyzzy.wxsize, xyzzy.wysize, self.alpha,band_list=[num_bands])

   self.resampler(path, dsquery, dstile, image_format)

  self.alpha = None

 def write_overview_tile(self, tx, ty, tz,tms_osm):

  """Create image of a overview level tile and write it to disk."""

  image_format = self.get_overview_tile_format(tx, ty, tz)

  if image_format is None:
   return
  else:
   num_bands = self.get_num_bands(image_format)

  dsquery = self.mem_drv.Create('', 2*self.tile_size, 2*self.tile_size, num_bands)
  self.fill_init_dest(dsquery)
  # tms: z=19: 281626
  # -z=18-140813 176168*2=352336; 176168*2+1=352337
  # -- 352336,352337
  y_from=2*ty
  y_to=2*ty + 1
  ty_tms=ty;
  s_y_type="tms"
  if tms_osm:
   # osm: z=19: 281626
   # -z=18-140813 85975*2+1=171951; 85975*2=171950
   # -- 171951,171950 [in range: last/end not used]
   y_from=2*ty + 1
   y_to=2*ty
   ty_tms=(2**tz-1) - ty
   s_y_type="osm"
  s_tile_id="{0}-{1}-{2}.{3}".format(str(tz), str(tx),str(ty),s_y_type)
  if self.verbose:
   #  Build from zoom 19  tiles: (281626, 171951) (281627, 171951) (281626, 171950) (281627, 171950)
   print "\tBuild  [",s_tile_id,"] from [",self.output_dir,"] zoom", tz+1," tiles [",s_y_type,"]: ", (2*tx, y_from), (2*tx+1, y_from),(2*tx, y_to), (2*tx+1, y_to)

  for cx, cy, child_image_format in self.iter_children(tx, ty, tz):
   if (ty_tms==0 and cy==1) or (ty_tms!=0 and (cy % (y_from)) != 0):
    tileposy = 0
   else:
    tileposy = self.tile_size
   if tx:
    tileposx = cx % (2*tx) * self.tile_size
   elif tx==0 and cx==1:
    tileposx = self.tile_size
   else:
    tileposx = 0

   path = self.get_full_path(cx, cy, tz+1, format_extension[child_image_format])

   dsquerytile = gdal.Open(path, gdal.GA_ReadOnly)

   dsquery.WriteRaster(tileposx, tileposy, self.tile_size, self.tile_size,
    dsquerytile.ReadRaster(0, 0, self.tile_size, self.tile_size),
    band_list=range(1, dsquerytile.RasterCount+1))

   if image_format == "PNG" and dsquerytile.RasterCount != num_bands:
    dsquery.WriteRaster(tileposx, tileposy, self.tile_size, self.tile_size,
     self.get_alpha_filler(), band_list=[num_bands])

  dstile = self.mem_drv.Create('', self.tile_size, self.tile_size, num_bands)
  path = self.get_full_path(tx, ty, tz, format_extension[image_format])
  self.resampler(path, dsquery, dstile, image_format)

 def iter_children(self, tx, ty, tz,):
  """Generate all children of the given tile produced on the lower level."""
  # no next to bake changes due to different tms/osm logic
  for y in range(2*ty,2*ty + 2):
   for x in range(2*tx, 2*tx + 2):
    image_format = self.try_to_use_existing_tile(x, y, tz+1)
    if image_format is not None:
     yield x, y, image_format

 def read_alpha(self, xyzzy):
  self.alpha = self.alpha_band.ReadRaster(xyzzy.rx, xyzzy.ry, xyzzy.rxsize, xyzzy.rysize, xyzzy.wxsize, xyzzy.wysize)

 def fill_init_dest(self, image):
  if self.init_dest is not None:
   for i,v in enumerate(self.init_dest[:image.RasterCount]):
    image.GetRasterBand(i+1).Fill(v)

 def get_num_bands(self, image_format):
  if image_format == "JPEG":
   return self.data_bands_count
  else:
   return self.data_bands_count + 1

 def get_alpha_filler(self):
  if self.alpha_filler is None:
   self.alpha_filler = "\xff" * (self.tile_size * self.tile_size)
  return self.alpha_filler

 def try_to_use_existing_tile(self, tx, ty, tz):
  """Return image format of the tile if it exists already on disk."""
  for image_format in self.image_formats:
   if os.path.exists(self.get_full_path(tx, ty, tz, format_extension[image_format])):
    return image_format
  return None

 def tile_exists(self, tx, ty, tz):
  return self.try_to_use_existing_tile(tx, ty, tz) != None

 def get_full_path(self, tx, ty, tz, extension):
  return os.path.join(self.output_dir, get_tile_filename(tx, ty, tz, extension,self.mbtiles))


class SimpleImageOutput(BaseImageOutput):

 """Image output using only one image format."""

 def get_base_tile_format(self, tx, ty, tz, xyzzy):
  if self.format == "PNG":
   self.read_alpha(xyzzy)

  return self.format

 def get_overview_tile_format(self, tx, ty, tz,):
  return self.format


class HybridImageOutput(BaseImageOutput):

 """Image output which skips fully transparent tiles, saves the fully opaque
 as JPEG and the rest as PNG.
 """

 def __init__(self, out_ds, tile_size, resampler, init_dest, output_dir, verbose):
  BaseImageOutput.__init__(self, out_ds, tile_size, resampler, init_dest, output_dir, verbose, ["JPEG", "PNG"],False)

  img = self.mem_drv.Create("", self.tile_size, self.tile_size, 1)
  rb = img.GetRasterBand(1)
  rb.Fill(0)
  self.transparent_checksum = rb.Checksum(0, 0, self.tile_size, self.tile_size)
  rb.Fill(255)
  self.opaque_checksum = rb.Checksum(0, 0, self.tile_size, self.tile_size)

 def get_base_tile_format(self, tx, ty, tz, xyzzy):
  if xyzzy.rxsize == self.tile_size and xyzzy.rysize == self.tile_size:
   c = self.alpha_band.Checksum(xyzzy.rx, xyzzy.ry, self.tile_size, self.tile_size)

   if c == self.transparent_checksum:
    if self.verbose:
     print "\tTile generation skipped because it is fully transparent"
    return None
   elif c == self.opaque_checksum:
    image_format = "JPEG"
   else:
    image_format = "PNG"
    self.read_alpha(xyzzy)
  else:
   self.read_alpha(xyzzy)
   transparent, opaque = self.transparent_or_opaque(self.alpha)

   if transparent:
    if self.verbose:
     print "\tTile generation skipped because it is fully transparent"
    return None
   elif opaque:
    image_format = "JPEG"
   else:
    image_format = "PNG"

  if self.verbose:
   print "\tSaving tile in %s format" % image_format

  return image_format

 def get_overview_tile_format(self, tx, ty, tz):
  children = list(self.iter_children(tx, ty, tz))

  if len(children) == 0:
   if self.verbose:
    print "\tTile generation skipped because it is fully transparent"
   return None

  if any(image_format == "PNG" for x, y, image_format in children) or len(children) < 4:
   image_format = "PNG"
  else:
   image_format = "JPEG"

  if self.verbose:
   print "\tSaving tile in %s format" % image_format

  return image_format

 def transparent_or_opaque(self, alpha):
  transparent = opaque = True
  for c in alpha:
   transparent = transparent and c == '\x00'
   opaque = opaque and c == '\xff'
  assert not (transparent and opaque)
  return transparent, opaque


def Resampler(name):

 """Return a function performing given resampling algorithm."""

 def resample_average(path, dsquery, dstile, image_format):
  for i in range(1, dstile.RasterCount+1):
   res = gdal.RegenerateOverview(dsquery.GetRasterBand(i), dstile.GetRasterBand(i), "average")
   if res != 0:
       raise ImageOutputException("RegenerateOverview() failed with error %d" % res)

  gdal_write(path, dstile, image_format)

 def resample_antialias(path, dsquery, dstile, image_format):
  querysize = dsquery.RasterXSize
  tilesize = dstile.RasterXSize

  array = numpy.zeros((querysize, querysize, 4), numpy.uint8)
  for i in range(dstile.RasterCount):
   array[:,:,i] = gdalarray.BandReadAsArray(dsquery.GetRasterBand(i+1), 0, 0, querysize, querysize)
  im = Image.fromarray(array, 'RGBA') # Always four bands
  im1 = im.resize((tilesize,tilesize), Image.ANTIALIAS)

  if os.path.exists(path):
   im0 = Image.open(path)
   im1 = Image.composite(im1, im0, im1)

  ensure_dir_exists(path)

  if image_format == "JPEG":
   im1.save(path, image_format, quality=jpeg_quality)
  else:
   im1.save(path, image_format)


 if name == "average":
  return resample_average
 elif name == "antialias":
  return resample_antialias

 resampling_methods = {
  "near"        : gdal.GRA_NearestNeighbour,
  "bilinear"    : gdal.GRA_Bilinear,
  "cubic"       : gdal.GRA_Cubic,
  "cubicspline" : gdal.GRA_CubicSpline,
  "lanczos"     : gdal.GRA_Lanczos
 }

 resampling_method = resampling_methods[name]

 def resample_gdal(path, dsquery, dstile, image_format):
  querysize = dsquery.RasterXSize
  tilesize = dstile.RasterXSize

  dsquery.SetGeoTransform( (0.0, tilesize / float(querysize), 0.0, 0.0, 0.0, tilesize / float(querysize)) )
  dstile.SetGeoTransform( (0.0, 1.0, 0.0, 0.0, 0.0, 1.0) )

  res = gdal.ReprojectImage(dsquery, dstile, None, None, resampling_method)
  if res != 0:
      raise ImageOutputException("ReprojectImage() failed with error %d" % res)

  gdal_write(path, dstile, image_format)

 return resample_gdal


def gdal_write(path, dstile, image_format):
 ensure_dir_exists(path)
 driver = get_gdal_driver(image_format)

 if image_format == "JPEG":
  driver.CreateCopy(path, dstile, strict=0, options=jpeg_gdal_options)
 else:
  driver.CreateCopy(path, dstile, strict=0)


def get_gdal_driver(name):
 driver = gdal.GetDriverByName(name)
 if driver is None:
  raise Exception("The '%s' driver was not found, is it available in this GDAL build?" % name)
 else:
  return driver


def get_tile_filename(tx, ty, tz, extension,mbtiles):
 if mbtiles:
  return "%s-%s-%s.%s" % (str(tz), str(tx),ty, extension)
 else:
  return os.path.join(str(tz), str(tx), "%s.%s" % (ty, extension))


def ensure_dir_exists(path):
 dirname = os.path.dirname(path)
 if not os.path.exists(dirname):
  os.makedirs(dirname)


class Xyzzy(object):

 """Collection of coordinates describing what to read where for the given tile at the base level."""

 def __init__(self, querysize, rx, ry, rxsize, rysize, wx, wy, wxsize, wysize):
  self.querysize = querysize
  self.rx = rx
  self.ry = ry
  self.rxsize = rxsize
  self.rysize = rysize
  self.wx = wx
  self.wy = wy
  self.wxsize = wxsize
  self.wysize = wysize

# ------------------------------------------------
# Functions for calculating KML <LookAt> tag:

def angular_distance(lng1, lat1, lng2, lat2):
 """
 Haversine formula on sphere - for higher precission on ellipsoid use the Vincenty formula
 http://www.movable-type.co.uk/scripts/latlong.html
 """

 phi1 = math.radians(lat1)
 phi2 = math.radians(lat2)

 d_phi = math.radians(lat2 - lat1)
 d_lmd = math.radians(lng2 - lng1)

 A = math.pow(math.sin(d_phi / 2), 2) + \
  math.cos(phi1) * math.cos(phi2) * \
  math.pow(math.sin(d_lmd / 2), 2)

 return 2 * math.atan2(math.sqrt(A), math.sqrt(1 - A))

def computeFitLookAt(minlng, minlat, maxlng, maxlat):
 "Ported from http://code.google.com/p/earth-api-samples/source/browse/trunk/lib/ge-poly-fit-hack.js"

 DEGREES = math.pi / 180.0
 EARTH_RADIUS = 6378137

 viewrange = 0.0
 center = [ minlng+(maxlng-minlng)/2.0, minlat+(maxlat-minlat)/2.0 ]

 lngSpan = EARTH_RADIUS * angular_distance(center[0], minlat, center[0], maxlat)
 latSpan = EARTH_RADIUS * angular_distance(minlng, center[1], maxlng, center[1])

 aspectRatio = 1.0
 PAD_FACTOR = 1.5 # add 50% to the computed range for padding

 aspectUse = max(aspectRatio, min(1.0, lngSpan / latSpan))
 alpha = (45.0 / (aspectUse + 0.4) - 2.0) * DEGREES # computed experimentally

 # create LookAt using distance formula
 if (lngSpan > latSpan): # polygon is wide
  beta = min(90 * DEGREES, alpha + lngSpan / 2 / EARTH_RADIUS)
 else: # polygon is taller
  beta = min(90 * DEGREES, alpha + latSpan / 2 / EARTH_RADIUS)

 viewrange = PAD_FACTOR * EARTH_RADIUS * (math.sin(beta) * \
  math.sqrt(1 / math.pow(math.tan(alpha),2) + 1) - 1)

 return center[0], center[1], viewrange

# =============================================================================


if __name__=='__main__':
 argv = gdal.GeneralCmdLineProcessor( sys.argv )
 if argv:
  gdal2mbtiles = GDAL2MbTiles( argv[1:], gdalcache=128*1024*1024 )
  gdal2mbtiles.process()
