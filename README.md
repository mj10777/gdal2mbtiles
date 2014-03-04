gdal2mbtiles
============

gdal2mbtiles - gdal2tiles with mbtiles support

***

* this is an adapted version from [https://github.com/mj10777/mapmbtiles] where
   * `gdal2mbtiles.py`,`globalmercator.py` and `mbtiles.py`
      * have been combinded into one file `gdal2mbtiles.py`

* goal is that `gdal2mbtiles.py` will run 
   * in the same way as `gdal2tiles.py`, but with mbtiles support
   * that mbtyles task can be run from small pythons scripts
      * using the functionality taken from the `Mapbox` and `Landez` projects

* final goal:
   * to submit `gdal2mbtiles.py` to the GDAL-Project

***

* the `Map Tile Cutter` is based on the original work of:
    *  Klokan Petr Pridal `klokan.petr.pridal@gmail.com`
    *  the original project source can be found at [http://code.google.com/p/maptiler]
       * the source code used here was taken from the `maptiler_1.0.beta2_all.deb`
          * this source is newer than the source in the git-repository
       * to my knowlage the support for this project has been disscontinued in its Open-Source form.

* some of the Mbtiles creation logic came from the `Mapbox` project:
    *  the original project source can be found at [https://github.com/mapbox/mbutil]
       * much of this code however differs

* some of the Mbtiles funtionality came from the `Landez` project:
    *  the original project source can be found at [https://github.com/makinacorpus/landez]
       * much of this code as been adaped to read/write directly to the mbtiles files
          * beforhand it used a cache that was filled by another class


***

* The main goal was to adapt the `gdal2tiles.py` to also support the creation of a mbtiles Databases
    * `gdal2tiles.py` has been renamed to `gdal2mbtiles.py` to avoid conficts with the original gdal version
    * the basic functionality has otherwise not been changed
       * `tilemapresource.xml` : is now written correctly (x and y values were switched)
       * all of the y tiles are created first (each x directory is filled compleatly, before the next is created)
          * this makes it quicker to resume after an interuption

* The created mbtiles Databases are based on the same logic used in the geopaparrazi project:
    * [https://github.com/geopaparazzi/geopaparazzi/wiki/mbtiles-Implementation]
    * it uses `tiles` as a view and not a table
    * it will check for `blank` images (all pixels have the same RGB value) and store this only once

* When install on linux, a `soft-link`  called `gdal2mbtiles` will be created
    * this can be called with the same paramaters as `gdal2tiles.py`
    * when called with `--mbtiles` the [output] parameter will be ignored
       * in the directory as the image: a `image_name.mbtiles` will be created
       * the created tiles will then be stored in a mbtiles Database and NOT as a tile-directory
    * when called with `--mbtiles_from_disk` 
       * the [input_file] parameter will be used as the tile-directory
       * the [output] parameter will be used as the file-name for mbtiles
          * the tiles of the tile-directory will be imported into the mbtiles file
    * when called with `--mbtiles_to_disk` 
       * the [input_file] parameter will be used as the file-name for mbtiles
       * the [output] parameter will be used as the tile-directory
          * the tiles of the mbtiles will be exported into the tile-directory

***

In the `samples` directory, there are some python scripts that use this project when installed

* they are based on the functionality taken from the `Landez` project
   * importing all or a portion of one mbtiles to another
      * this can be used to `convert` a table based mbtiles to a `view` base mbtiles
         * this will also check for `blank` images
   * filling a mbtiles from a WMS-Server
   * exporting all or a portion of one mbtiles to a image
      * when `tif` is used, the result will be a geotif

A sample mbtiles Database has been included, which some of the samples use

* due to size limitations, other mbtiles files cannot be included in this project
   * where used in the sample sripts, a link has been supplided where this file can be downloaded:
      * [http://www.mj10777.de/public/download/mbtiles/]

The `Landez` project also supports other functions, not yet tested:

* `Blend tiles together`
* `Merge multiple sources of tiles (URL, WMS, MBTiles, Mapnik stylesheet) together.`
* `Composite a WMS layer with OpenStreetMap using transparency`
* `Add post-processing filters`
* `Replace a specific color by transparent pixels`


---

2014-02-28

---
