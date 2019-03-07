#!/usr/bin/python
"""Deal with COSMOS catalogs"""

import re
import sys
import urllib2
import numpy as np
import xml.etree.ElementTree as ET

import lsst.afw.image as afwImage
from   lsst.afw.fits.fitsLib import MemFileManager, memmove
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from . import utils

import lsst.afw.display as afwDisplay

def readAlexieFITS(fileName="/Users/rhl/Dropbox/Robert/cosmos_forhsc_feb20_2012.fits"):
    """Read Alexie's COSMOS table, converting it to a minimal form that can be used with afwTable.matchRaDec
    """
    import pyfits

    tbl = pyfits.open(fileName)
    data = tbl[1].data

    schema = afwTable.SimpleTable.makeMinimalSchema()
    schema.addField(afwTable.Field["I"]("mu_class", "S/G classification. 1: galaxy; 2: star; 3: other"))
    schema.addField(afwTable.Field["F"]("mag_auto", "SExtractor's mag_auto"))

    cat = afwTable.SimpleCatalog(schema)

    n = len(data)
    # There has to be a better way! https://jira.lsstcorp.org/browse/DM-347
    cat.reserve(n)
    for i in range(n):
        cat.addNew()

    cat["id"][:]        = data["IDENT"]
    cat["coord.ra"][:]  = data["ALPHA_J2000"]
    cat["coord.dec"][:] = data["DELTA_J2000"]
    cat["mu_class"][:]  = data["MU_CLASS"]
    cat["mag_auto"][:]  = data["MAG_AUTO"]

    cat.writeFits("cosmos_sg.fits")

def readAlexieMASKED_reg(fileName="/Users/rhl/Dropbox/Robert/MASKED.reg"):
    """Read Alexie's COSMOS table, converting it to a minimal form that can be used with afwTable.matchRaDec
    """
    with open(fileName, "r") as fd:
        while True:
            line = fd.readline()
            if not line:
                break
            mat = re.search(r"fk5;circle\(\s*(\d+\.\d+),\s*(\d+\.\d+)\s*,\s*(\d+\.\d+)\"\)# color=\w+\s*$",
                            line)
            if mat:
                ra, dec, rad = [float(_) for _ in mat.groups()]
                print ra, dec, rad
            return

    schema = afwTable.SimpleTable.makeMinimalSchema()
    schema.addField(afwTable.Field["I"]("mu_class", "S/G classification. 1: galaxy; 2: star; 3: other"))
    schema.addField(afwTable.Field["F"]("mag_auto", "SExtractor's mag_auto"))

    cat = afwTable.SimpleCatalog(schema)

    n = len(data)
    # There has to be a better way! https://jira.lsstcorp.org/browse/DM-347
    cat.reserve(n)
    for i in range(n):
        cat.addNew()

    cat["id"][:]        = data["IDENT"]
    cat["coord.ra"][:]  = data["ALPHA_J2000"]
    cat["coord.dec"][:] = data["DELTA_J2000"]
    cat["mu_class"][:]  = data["MU_CLASS"]
    cat["mag_auto"][:]  = data["MAG_AUTO"]

    cat.writeFits("cosmos_sg.fits")


def getCosmosCutout(ra=150.23983, dec=+2.56283, sizeX=3):
    """Return an ACS COSMOS cutout from IRSA

    \param ra   Right ascension in decimal degrees (J2000)
    \param dec  Declination in decimal degrees (J2000)
    \param size Size of cutout (arcsec)

    See http://irsa.ipac.caltech.edu/applications/Cutouts/docs/CutoutsProgramInterface.html
    but beware that (as of 2014-03-31) the locstr in the example isn't quite right (needs
    a %20 between the ra and dec)
    """
    url = "http://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts?" + \
        "mission=COSMOS&max_size=180&locstr=%(ra).5f%%20%(dec).5f&sizeX=%(sizeX)d&ntable_cutouts=1&cutouttbl1=acs_mosaic_2.0&mode=PI" % \
        dict(ra=ra, dec=dec, sizeX=sizeX)

    response = urllib2.urlopen(url)
    html = response.read()

    root = ET.fromstring(html)

    status = root.get("status")
    if status != "ok":
        raise RuntimeError("Failed to open cutout URL: %s %s" % (status, root.get("message")))
    #
    # Find the cutouts
    #
    cutouts = root.find("images").find("cutouts")
    if not len(cutouts):
        html2 = urllib2.urlopen(root.find("summary").find("resultHtml").text.strip()).read()
        # fix things for ET
        html2 = re.sub(r"&(nbsp;|micro)", "", html2)

        root2 = ET.fromstring(html2)

        body = root2[1]
        try:
            msg = body[-1][0][3][0][0][-1][0][0].text # aaarghhhh
            msg = re.sub(r":\s*$", "", msg)
        except:
            msg = "Unable to retrieve cutout" + ":"
        print >> sys.stderr, "Failed to retrieve cutout %s" % re.sub(r"\s+", " ", msg)
        return afwImage.ExposureF(30,30)
    #
    # Read the fits data
    #
    try:
        fitsFD = urllib2.urlopen(cutouts[0].text.strip()) # 0: fits, 1: jpg
        fitsData = fitsFD.read()
    except Exception as e:
        fitsFD = None
    finally:
        del fitsFD
    #
    # and create an afwExposure
    #
    manager = MemFileManager(len(fitsData))
    memmove(manager.getData(), fitsData)

    return afwImage.ExposureF(manager)

def acsEventCallback(key, source, im, frame):
    """Callback for event handlers to find COSMOS ACS cutout.

    \param key     Key struck
    \param source  The Source under the cursor
    \param im      The (HSC) image cutout displayed in frame
    \param frame   The frame that the HSC data's displayed in (we'll use the next one)

    We also use the following static members of utils.EventHandler (if set):
    sizeCutout   The size of the HSC cutout (arcsec; default: 4.0)
    scale   Make the COSMOS image with pixel size scale*HSC's pixel size (default 0.25 => 0.42mas)

    Use as e.g. utils.eventCallbacks['c'] = cosmos.acsEventCallback
    """
    sizeCutout = utils.EventHandler.sizeCutout if hasattr(utils.EventHandler, "sizeCutout") else 4.0 # arcsec
    scale = utils.EventHandler.scale if hasattr(utils.EventHandler, "scale") else 0.25   # Pixel size scaling

    pos = source.get("coord")
    exp = getCosmosCutout(*pos.getPosition(), sizeX=sizeCutout)

    if im and exp and exp.getWcs():
        #
        # Resample and rotate to the HSC orientation
        #
        warpingControl = afwMath.WarpingControl("lanczos3")
        rat = im.getWcs().pixelScale().asArcseconds()/exp.getWcs().pixelScale().asArcseconds()
        hsize = int(0.5*exp.getWidth()/(scale*rat))
        rexp = afwImage.ExposureF(2*hsize + 1, 2*hsize + 1)
        rexp.setWcs(afwImage.Wcs(pos.getPosition(), afwGeom.Point2D(hsize, hsize),
                                im.getWcs().getCDMatrix()*scale))
        afwMath.warpExposure(rexp, exp, warpingControl)
    else:
        print "\nI'm unable to remap the cosmos image to your coordinates, sorry"
        rexp = exp.getMaskedImage().getImage()

    frame += 1
    rim = rexp
    if hasattr(rim, "getMaskedImage"):
        rim = rim.getMaskedImage().getImage()
    if hasattr(rim, "getImage"):
        rim = rim.getImage()
    disp = afwDisplay.Display(frame=frame)
    disp.mtv(rim)

    if hasattr(rexp, "getWcs"):
        cen = rexp.getWcs().skyToPixel(pos) - afwGeom.PointD(rexp.getXY0())
        disp.pan(*cen)
        disp.dot('+', *cen)

if __name__ == "__main__":
    readAlexieMASKED_reg(fileName="/Users/rhl/Dropbox/Robert/MASKED.reg")
