from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
import qgis.utils
from Ui_GeopapaTile import Ui_GeopapaTile

import os, sys, time,  math, subprocess, platform, datetime
from osgeo import gdal, ogr
from osgeo.gdalconst import *

from xml.dom.minidom import parseString

import numpy as np
import fnmatch



class GeopapaTileDialog(QDialog, Ui_GeopapaTile):

    def which(self,name, flags=os.X_OK):
        """Search PATH for executable files with the given name.

        On newer versions of MS-Windows, the PATHEXT environment variable will be
        set to the list of file extensions for files considered executable. This
        will normally include things like ".EXE". This fuction will also find files
        with the given name ending with any of these extensions.

        On MS-Windows the only flag that has any meaning is os.F_OK. Any other
        flags will be ignored.

        @type name: C{str}
        @param name: The name for which to search.

        @type flags: C{int}
        @param flags: Arguments to L{os.access}.

        @rtype: C{list}
        @param: A list of the full paths to files found, in the
        order in which they were found.
        mj10777: 2014-03-12 added from: http://twistedmatrix.com/trac/browser/tags/releases/twisted-8.2.0/twisted/python/procutils.py
        mj10777: 2014-03-12 adapted: return result as in linux : first found result as string otherwsise empty
    """
        result = []
        exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
        path = os.environ.get('PATH', None)
        if path is None:
            return []
        s_path=""
        for p in os.environ.get('PATH', '').split(os.pathsep):
            p = os.path.join(p, name)
            if os.access(p, flags):
                result.append(p)
            for e in exts:
                pext = p + e
                if os.access(pext, flags):
                    result.append(pext)
        if len(result) > 0:
         s_path=result[0]
        return s_path

    def __init__(self, iface):
        QDialog.__init__(self)
        self.iface = iface
        self.setupUi(self)
        self.NOVALUE=-9999
        self.connect(self.btnOutputDir, SIGNAL("clicked()"), self.outDirFun)
        self.connect(self.buttonBox, SIGNAL("accepted()"),self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"),self, SLOT("reject()"))
        QObject.connect(self.buttonBox, SIGNAL("helpRequested()"),self.call_help)

        mapCanvas = self.iface.mapCanvas()
        # get projecr epsg
        self.mapRenderer = iface.mapCanvas().mapRenderer()
        if QGis.QGIS_VERSION_INT < 10900:
                projectCRS=self.mapRenderer.destinationSrs().epsg()
        else:
                projectCRS=int(self.mapRenderer.destinationCrs().postgisSrid())
        self.lineEditSourceSrs.clear()
        self.lineEditSourceSrs.insert(str(projectCRS))

        # get map rectangle
        mapRect = self.mapRenderer.extent()
        self.xStart = mapRect.xMinimum()
        self.xEnd = mapRect.xMaximum()
        self.yStart = mapRect.yMinimum()
        self.yEnd = mapRect.yMaximum()

#       width = mapRenderer.width()
#       height = mapRenderer.height()

        #printed raster will have (width * metersPerPixel) meters in X and
        #(height * metersPerPixel) meters in Y
#       xDiff = width * metersPerPixel
#       yDiff = height * metersPerPixel

        #Generate first rectangle
#       newRect = QgsRectangle(xStart, yStart, xStart + xDiff, yStart + yDiff)
#       mapRenderer.setExtent(newRect)


    def call_help(self):
        qgis.utils.showPluginHelp()

    def createWorldFile(self,fileName, mainScale):
        #>>Create World File
        mapRect = self.mapRenderer.extent()
        f = open(fileName + ".tfw", 'w')
        f.write(str(mainScale) + '\n')
        f.write(str(0) + '\n')
        f.write(str(0) + '\n')
        f.write('-' + str(mainScale) + '\n')
        f.write(str(mapRect.xMinimum() + mainScale/2) + '\n')
        f.write(str(mapRect.yMaximum() - mainScale/2))
        f.close()


    def outDirFun(self):
        "Display file dialog for output directory"
        self.lineOutput.clear()
        dirName = QFileDialog.getExistingDirectory( self,
                                                      self.tr("Select save directory"),
                                                      QDir.currentPath(),
                                                      QFileDialog.ShowDirsOnly|
                                                      QFileDialog.ReadOnly )
        if QGis.QGIS_VERSION_INT < 10900:
                if not dirName.isEmpty():
                        self.lineOutput.clear()
                        self.lineOutput.insert(dirName)
        else:
                if dirName:
                        self.lineOutput.clear()
                        self.lineOutput.insert(dirName)
        return dirName


    def doMapurlfile(self,outDir,sourceName):

        sourceDir=outDir+"/"+sourceName+"/"
        path_xml_file=sourceDir+"tilemapresource.xml"
        if not os.path.isfile(path_xml_file):
                self.textEdit.append( "Error  : [" + path_xml_file+ "] Output tilemapresource.xml was not created  : will exit")
                return
        # read data from tilemapresource.xml"
        xmlfile = open(path_xml_file);
        data = xmlfile.read();
        xmlfile.close()
        dom = parseString(data)
        tagName = "BoundingBox"
        BoundingBox = dom.getElementsByTagName(tagName)[0] #.replace("<"+tagName+">","").replace("</"+tagName+">","")
        # get center coords
        ## attenzione! gdal2tiles.py usa x ed y invertite!
        miny = float(BoundingBox.attributes["minx"].value)
        minx = float(BoundingBox.attributes["miny"].value)
        maxy = float(BoundingBox.attributes["maxx"].value)
        maxx = float(BoundingBox.attributes["maxy"].value)

        xcenter = minx + ((maxx - minx)/2)
        ycenter = miny + ((maxy - miny)/2)

        orders = []
        # get minzoom and max zoom fro TileSet
        TileSets = dom.getElementsByTagName("TileSet")
        for TileSet in TileSets:
                #print TileSet.attributes["order"].value
                orders.append(int(TileSet.attributes["order"].value))

        minzoom = min(orders)
        maxzoom = max(orders)

        ### Create mapurl file
        mapurlfile = str(sourceDir).rstrip("/")+".mapurl"
        f = open(mapurlfile, 'w')
        f.write("url=%s/ZZZ/XXX/YYY.png\n" % (sourceName))
        f.write("minzoom=%s\n" % (minzoom))
        f.write("maxzoom=%s\n" % (maxzoom))
        f.write("bounds=%s %s %s %s\n" % (minx,miny,maxx,maxy))
        f.write("center=%s %s\n" % (xcenter,ycenter))
        f.write("type=tms\n")
        f.close()

        return mapurlfile


    def doImage(self,myImagePath,myImagepathNoExt):

        dpi = self.spinBoxDpi.value()
        myWidth = self.spinBoxWidth.value() #width in pixels

        metersPerPixel = (self.xEnd-self.xStart)/(myWidth)
        myHeight = int((self.yEnd - self.yStart)/metersPerPixel)
        #print myWidth
        #print myHeight
        xPaperSize = myWidth/(dpi/25.4)
        yPaperSize = myHeight/(dpi/25.4)
        #self.textEdit.append("xPaperSize %f" % (xPaperSize))
        #self.textEdit.append("yPaperSize %f" % (yPaperSize))

################################################
        # create image
#       img = QImage(QSize(myWidth,myHeight), QImage.Format_ARGB32_Premultiplied)

        # set image's background color
#       color = QColor(255,255,255)
#       img.fill(color.rgb())

        # create painter
#       p = QPainter()
#       p.begin(img)
#       p.setRenderHint(QPainter.Antialiasing)
#       mapRenderer2 = mapRenderer
#       mapRenderer2.setOutputSize(img.size(), dpi)
        # do the rendering
#       mapRenderer2.render(p)
#       p.end()

        # save image
#       img.save("/tmp/amap.tif","tif")
#       createWorldFile("/tmp/amap", metersPerPixel)
############################################

        composition = QgsComposition(self.mapRenderer)
        composition.setPrintResolution(dpi)
        composition.setPaperSize(xPaperSize, yPaperSize)
        composition.setPlotStyle(QgsComposition.Print)
        #dpi = composition.printResolution()
        dpmm = dpi / 25.4 #get dots per mm
        # add a map to the composition
        composerMap = QgsComposerMap(composition,0,0,composition.paperWidth(),composition.paperHeight())
        composition.addItem(composerMap)
        # create output image and initialize it
        image = QImage(QSize(myWidth, myHeight), QImage.Format_ARGB32) #output image size
        image.setDotsPerMeterX(dpmm * 1000) #mm to meters
        image.setDotsPerMeterY(dpmm * 1000) #mm to meters
        image.fill(0)
        # render the composition
        imagePainter = QPainter(image)
        sourceArea = QRectF(0, 0, composition.paperWidth(), composition.paperHeight() )
        targetArea = QRectF(0, 0, myWidth, myHeight)
        composition.render(imagePainter, targetArea, sourceArea)
        imagePainter.end()
        self.textEdit.append("width:"+str(myWidth)+" height: "+str(myHeight))
        image.save(myImagePath , "tif")
        self.createWorldFile(myImagepathNoExt, metersPerPixel)



    def accept(self):
        # Called when "OK" button pressed

        self.textEdit.clear()
        inizio = datetime.datetime.now()
        self.textEdit.append("Starting...")
        self.textEdit.append("QGis.QGIS_VERSION_INT: - "+str(QGis.QGIS_VERSION_INT) )
        #  20100

        outDir=self.lineOutput.text()
        if (outDir == ''):
                QMessageBox.critical(None,"Exiting gracefully","Output directory not defined %s!" % (outDir))
                return
        if (len(outDir.split(' ')) > 1):
                QMessageBox.critical(None,"Exiting gracefully","Output directory contains spaces %s!" % (outDir))
                return
        gdalwarp="gdalwarp"
        gdal2mbtiles="gdal2mbtiles.py"
        gdal2mbtiles_parm="--mbtiles"
        if platform.system() == "Windows":
                gdal2tiles="gdal2tiles.bat"
        else:
                gdal2tiles="gdal2tiles.py"
        path_gdalwarp=self.which(gdalwarp)
        path_gdal2tiles=self.which(gdal2tiles)
        path_gdal2mbtiles=self.which(gdal2mbtiles)
        if path_gdalwarp == "":
         self.textEdit.append( "Error  : [" + gdalwarp+ "] not found : will exit")
         return
        if path_gdal2tiles == "" and path_gdal2mbtiles == "":
         self.textEdit.append( "Error  : [" + gdal2tiles+ "] not found : will exit")
         return
        if path_gdal2mbtiles != "":
         self.textEdit.append( "Info  : [" + gdal2mbtiles+ "] was found and will be used")
         gdal2tiles=gdal2mbtiles
        else:
         gdal2mbtiles_parm=""
        outImage=self.lineEditImageName.text()
        ssrs=self.lineEditSourceSrs.text()
        zoom=self.lineEditZoom.text()
        zoom=zoom.replace(" ","")
        sourceName=self.lineEditSourceName.text()
        if gdal2mbtiles_parm == "":
         sourceDir=str(outDir+"/"+sourceName+"/")
         if not os.path.isdir(sourceDir):
                 os.makedirs(sourceDir)
        else:
         sourceDir=str(outDir+"/"+outImage[:-4]+".mbtiles")
        imagepath=outDir+"/"+outImage
        imagepathNoExt=imagepath[:-4]

        self.doImage(imagepath,imagepathNoExt)
        if not os.path.isfile(imagepath):
                self.textEdit.append( "Error  : [" + imagepath+ "] Output image was not created  : will exit")
                return

        # do warp...
        gdalwarpCMD=str("-s_srs EPSG:"+ssrs+" -t_srs EPSG:3785 -r bilinear "+imagepath+" "+imagepathNoExt+"_3785.tif")
        self.textEdit.append("gdalwarp "+gdalwarpCMD)
#       os.system(gdalwarpCMD)
#       p=subprocess.Popen(gdalwarpCMD,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#       output, errors = p.communicate()
#       self.textEdit.append("gdalwarp output:\n - "+output)
#       self.textEdit.append("gdalwarp errors:\n - "+errors)

        processWarp = QProcess( parent=None )
        if QGis.QGIS_VERSION_INT < 10900:
                processWarp.start( gdalwarp,QStringList() << gdalwarpCMD.split(" "), QIODevice.ReadOnly )
        else:
                processWarp.start( gdalwarp,gdalwarpCMD.split(" "), QIODevice.ReadOnly )
        arr = QByteArray()
        if processWarp.waitForFinished(-1):
                arr = processWarp.readAllStandardOutput()
                processWarp.close()
        self.textEdit.append("gdalwarp:\n - "+ str(arr) )

        # do tiles...
        # mj10777; -r bilinear does not always work with some images: the default -r avarage will be used
        # gdal2tilesCMD=str("-z "+zoom+" -r bilinear "+imagepathNoExt+"_3785.tif "+sourceDir)
        gdal2tilesCMD=str(gdal2mbtiles_parm+" -z "+zoom+" "+imagepathNoExt+"_3785.tif "+sourceDir)
        self.textEdit.append("")

        processTiles = QProcess( parent=None )
        self.textEdit.append(gdal2tiles + " " + gdal2tilesCMD)
        if QGis.QGIS_VERSION_INT < 10900:
                processTiles.start( gdal2tiles,QStringList() << gdal2tilesCMD.split(" "), QIODevice.ReadOnly )
        else:
                processTiles.start( gdal2tiles,gdal2tilesCMD.split(" "), QIODevice.ReadOnly )

        arrT = QByteArray()
        if processTiles.waitForFinished(-1):
                arrT = processTiles.readAllStandardOutput()
                processTiles.close()

        self.textEdit.append(gdal2tiles+":\n - "+str(arrT) )
        if gdal2mbtiles_parm == "":
         # should not run when mbtiles is used
         self.doMapurlfile(outDir,sourceName)

        self.textEdit.append("")
        fine = datetime.datetime.now()
        self.textEdit.append( 'Starting time: '+ str(inizio) )
        self.textEdit.append( 'Ending time: '+ str(fine) )
        diff = fine-inizio
        self.textEdit.append( "Elapsed time: " + str(diff) )

        self.textEdit.append( "\nDone!\n")
        self.textEdit.append( "Go to " + outDir)
        self.textEdit.append( "then copy " + sourceName + " directory and " +sourceName + ".mapurl file")
        self.textEdit.append( "in the  /mnt/sdcard/maps directory on your Android device.")
        self.textEdit.append("")

