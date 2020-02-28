"""Utilities to look at deblender outputs

Based on Dustin Lang's code in meas_deblender/examples, but heavily rewritten since then, so don't blame him

E.g.
import lsst.daf.persistence as dafPersist
import lsst.afw.display as afwDisplay
import lsst.analysis.utils as utils
import lsst.analysis.deblender as deblender

disp = afwDisplay.Display()

butler = dafPersist.Butler(os.path.join(os.environ["SUPRIME_DATA_DIR"], "rerun", "rhl", "realFlats"))
dataId = dict(visit=905518, ccd=31)
calexp = butler.get("calexp", **dataId)
ss = butler.get('src', **dataId)
families = deblender.Families(ss, butler, nChildMin=0)

frame=1            # We'll use 0 for the children
disp.mtv(calexp, title="%(visit)s %(ccd)s" % dataId, frame=frame)
utils.showSourceSet(ss, frame=frame)
utils.showSourceSet(ss, mapperInfo=families.mapperInfo, frame=frame, symb='id', size=1.0, fontFamily="times", ctype=afwDisplay.GREEN)
#
# This puts us into a loop waiting on afwDisplay.  Sigh.
# Use d on an object in frame 1 to show the children in frame 0; q to quit the interactive loop
#
deblender.showBlend(calexp, families, display=disp)

"""

import math, re, sys
import numpy as np

from . import utils
import lsst.geom as geom
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.display as afwDisplay
import lsst.afw.display.rgb as afwRgb
import lsst.afw.display.utils as displayUtils

def getEllipses(src, nsigs=[1.], **kwargs):
    from matplotlib.patches import Ellipse

    xc = src.getX()
    yc = src.getY()
    x2 = src.getIxx()
    y2 = src.getIyy()
    xy = src.getIxy()
    # SExtractor manual v2.5, pg 29.
    a2 = (x2 + y2)/2. + math.sqrt(((x2 - y2)/2.)**2 + xy**2)
    b2 = (x2 + y2)/2. - math.sqrt(((x2 - y2)/2.)**2 + xy**2)
    theta = math.rad2deg(math.arctan2(2.*xy, (x2 - y2)) / 2.)
    a = math.sqrt(a2)
    b = math.sqrt(b2)
    ells = []
    for nsig in nsigs:
        ells.append(Ellipse([xc,yc], 2.*a*nsig, 2.*b*nsig, angle=theta, **kwargs))
    return ells

def drawEllipses(plt, src, **kwargs):
    return

    els = getEllipses(src, **kwargs)
    for el in els:
        plt.gca().add_artist(el)
    return els

def makeDeblendFamilyMosaic(mi, parent, kids, mapperInfo=None,
                            background=-10, maskbit=False, imBbox=None):
    """Create a mosaic of an object's children
    """

    aa = {}
    if maskbit:
        aa.update(mask=True)
    parent_im = footprintToImage(parent.getFootprint(), mi, **aa)
    bbox = geom.BoxD(parent.getFootprint().getBBox())
    pext = (bbox.getMinX(), bbox.getMaxX(), bbox.getMinY(), bbox.getMaxY())

    pks = parent.getFootprint().getPeaks()
    pix = [pk.getIx() for pk in pks]
    piy = [pk.getIy() for pk in pks]
    pfx = [pk.getFx() for pk in pks]
    pfy = [pk.getFy() for pk in pks]

    N = 1 + len(kids)
    S = np.ceil(np.sqrt(N))
    C = S
    R = np.ceil(float(N)/C)

    Rx,Ry = [],[]
    tts = []
    stys = []
    xys = []
    #
    # Find how large an image we need to display the parent and all the children
    #
    kidImages, kim = {}, None
    for kid in kids:
        kim = footprintToImage(kid.getFootprint(), mi, **aa)
        kidImages[kid] = kim

    if not kim:
        kim = parent_im.clone()

    if not imBbox:
        imBbox = parent_im.getBBox(afwImage.PARENT)
        for kid in kids:
            imBbox.include(kidImages[kid].getBBox(afwImage.PARENT))

    mos = displayUtils.Mosaic(background=background)
        
    bbox = geom.Box2I(geom.Point2I(kim.getX0() - imBbox.getMinX(),
                                         kim.getY0() - imBbox.getMinY()), kim.getDimensions())

    kidImages[parent] = parent_im       # not strictly a kid

    for kid in [parent] + kids:
        kim = kidImages[kid]
        #
        # Put the child into the correct place in the parent image.  We have to do this for
        # the parent too if some of the children extended outside its BBox
        #
        bbox = geom.Box2I(geom.Point2I(kim.getX0() - imBbox.getMinX(),
                                             kim.getY0() - imBbox.getMinY()), kim.getDimensions())

        _kim = parent_im.Factory(imBbox)
        if kid.getFootprint().getBBox().isEmpty():
            _kim[:] = 0
            pass
        else:
            _kim[bbox, afwImage.LOCAL] = kim

        mos.append(_kim, '%d%s' % (mapperInfo.getId(kid) if mapperInfo else (kid.getId() & 0xfff),
                                   "P" if kid == parent else "C"))
        del _kim

    return mos

def plotDeblendFamily(mi, parent, kids, mapperInfo=None, dkids=[],
                      background=-10, symbolSize=2,
                      plotb=False,
                      arcsinh=True, maskbit=False, display=afwDisplay.getDisplay(0)):
    """Display a deblend using afwDisplay

Each child is marked with a + at its centre (green if deblended-as-psf else red)
all the other peaks in its footprint are marked with x (cyan if deblended-as-psf else magenta)
    """

    if mi:
        try:
            mi = mi.getMaskedImage()        # maybe it's an Exposure?
        except AttributeError:
            pass
    
    mos = makeDeblendFamilyMosaic(mi, parent, kids, mapperInfo, background, maskbit)

    if mapperInfo:
        # some displays, e.g. ds9, don't handle those chars well
        title = re.sub(r"[{}']", "", str(mapperInfo.getId(parent, None)))
    else:
        title = "0x%x == %d" % (parent.getId(), (parent.getId() & 0xffff))

    mosaicImage = mos.makeMosaic(display=display, title=title)

    if display is not None:
        display.dot("%s  (%.1f, %1.f)" % (title, parent.getX(), parent.getY()),
                0.5*mosaicImage.getWidth(), 1.03*mosaicImage.getHeight(),
                ctype=afwDisplay.BLACK, fontFamily="times", size=3)

        px0, py0 = footprintToImage(parent.getFootprint(), mi).getXY0()

        with display.Buffering():
            for i, src in enumerate([parent] + kids):    
                x0, y0 = mos.getBBox(i).getMin()
                x0 -= px0; y0 -= py0

                if src.get("deblend_deblendedAsPsf"):
                    centroid_ctype = afwDisplay.GREEN
                    peak_ctype = afwDisplay.CYAN
                else:
                    centroid_ctype = afwDisplay.RED
                    peak_ctype = afwDisplay.MAGENTA

                display.dot("+", src.getX() + x0, src.getY() + y0,
                        size=symbolSize, ctype=centroid_ctype)
                for p in src.getFootprint().getPeaks():
                    display.dot("x", p.getFx() + x0, p.getFy() + y0,
                            size=0.5*symbolSize if i == 0 else symbolSize,
                            ctype=afwDisplay.YELLOW if i == 0 else peak_ctype)


    return mosaicImage

def footprintToImage(fp, mi=None, mask=False):
    if fp.isHeavy():
        try:
            fp = afwDet.cast_HeavyFootprintF(fp) # not needed with pybind11
        except AttributeError:
            pass
        pass
    elif mi is None:
        print("Unable to make a HeavyFootprint as image is None", file=sys.stderr)
    else:
        fp = afwDet.makeHeavyFootprint(fp, mi)
    bb = fp.getBBox()
    if mask:
        im = afwImage.MaskedImageF(bb.getWidth(), bb.getHeight())
    else:
        im = afwImage.ImageF(bb.getWidth(), bb.getHeight())
    im.setXY0(bb.getMinX(), bb.getMinY())

    try:
        fp.insert(im)
    except AttributeError:              # we failed to make it heavy
        assert not mi
        pass
    
    if mask:
        im = im.getMask()
    return im

class Families(list):
    def __init__(self, cat, butler=None, nChildMin=None):
        '''
        Returns [ (parent0, [child0, child1]), (parent1, [child0, ...]), ...]
        where parents are sorted by ID.  Only objects that are deblended are included (unless nChildMin == 0)

        if nChildMin is not None, include only deblends with at least than many children.  As a
        special case, nChildMin == 0 includes all objects, even those that aren't blended
        '''
        self.cat = cat
        self.mapperInfo = utils.makeMapperInfo(butler)

        # parent -> [children] map.
        children = {}
        for src in cat:
            pid = src.getParent()
            if not pid:
                continue

            if pid in children:
                children[pid].append(src)
            else:
                children[pid] = [src]

        parentIds = children.keys()
        #
        # Impose nChildMin
        #
        if nChildMin == 0:                  # include all objects:
            for src in cat:
                if src.getParent():         # already accounted for
                    continue

                sid = src.getId()
                if sid not in parentIds:
                    children[sid] = []
        elif nChildMin is not None:
            for pid in parentIds:
                if len(children[pid]) < nChildMin:
                    del children[pid]


        parentIds = sorted(children.keys())

        for pid in parentIds:
            self.append((cat.find(pid), children[pid]))

    def find(self, objId, matchRadius=10):
        """Return the object's family (you may specify either the ID for the parent or a child)"""
        x, y = None, None
        try:
            x, y = objId
        except TypeError:
            pass

        if x is not None:
            oneObjCatalog = afwTable.SourceCatalog(self.cat.getSchema())
            centroidName = self.cat.table.getSchema().getAliasMap().get("slot_Centroid")
            oneObjCatalog.table.defineCentroid(centroidName)

            s = oneObjCatalog.addNew()
            s.set("%s_x" % centroidName, x)
            s.set("%s_y" % centroidName, y)

            matched = afwTable.matchXy(self.cat, oneObjCatalog, matchRadius)

            if len(matched) == 0:
                print("Unable to find object at (%.2f, %.2f)" % (x, y), file=sys.stderr)
                return None

            if False:
                objId = self.mapperInfo.splitId(matched[0][0].getId(), asDict=True)["objId"]
            else:
                objId = matched[0][0].getId()

        if False:
            family = [f for f in self if self.mapperInfo.getId(f[0]) == objId]
        else:
            family = [f for f in self if f[0].getId() == objId]
        if family:
            return family[0]

        for family in self:
            for child in family[1]:
                if False:
                    if self.mapperInfo.getId(child) == objId:
                        return family
                else:
                    if child.getId() == objId:
                        return family

        return None

# backwards compatibility
    
def getFamilies(cat, butler, nChildMin=None):
    return Families(cat, butler, nChildMin)

def findFamily(families, objId):
    """Return the object's family (you may specify either the ID for the parent or a child)"""

    return families.find(objId)

def makeDisplayFamily(calexp, families, matchRadius=20, background=-0.1, display=None,
                      rgb=False, obeyXY0=True):
    """Factory function for callback function implementing showBlend"""
    def display_family(k, x, y, obeyXY0=obeyXY0):
        if calexp is not None:
            if obeyXY0:
                xy0 = calexp.getXY0()
                x += xy0[0]
                y += xy0[1]
        fam = families.find((x, y), matchRadius=matchRadius)
        if fam:
            if rgb:
                plotDeblendFamilyRGB(fam[0], rgbFileFmt="AI-%.0f,%.0f-%%s.png" % (x, y))
            else:
                plotDeblendFamily(calexp, *fam, mapperInfo=None, background=background, 
                                  display=display)

    return display_family

#
# plotDeblendFamilyRGB uses a couple of dicts to avoid repeated I/O
#
try:
    coaddDict
except NameError:
    coaddDict = {}
try:
    familiesDict
except NameError:
    familiesDict = {}

def plotDeblendFamilyRGB(parent, bands=['g', 'r', 'i'],
                         min=0.01, max=0.5, Q=8, rgbFileFmt=None):
    x, y = parent.getX(), parent.getY()

    fams = {}
    imBbox = geom.BoxI()
    for bandName in "GRI".upper():
        filterName = "HSC-%s" % bandName

        x0, y0 = coaddDict[filterName].getXY0()
        x0, y0 = 0, 0
        fams[filterName] = familiesDict[filterName].find((x + x0, y + y0), matchRadius=20)
        if not fams[filterName]:
            return
        parent, kids = fams[filterName]

        bbox = parent.getFootprint().getBBox()
        # Can children extend outside parent BBox?
        for kid in kids:
            kim = footprintToImage(kid.getFootprint(), coaddDict[filterName].getMaskedImage())
            bbox.include(kim.getBBox(afwImage.PARENT))

        imBbox.include(bbox)

    images = {} 
    for bandName in bands:
        filterName = "HSC-%s" % bandName.upper()

        images[bandName] = makeDeblendFamilyMosaic(coaddDict[filterName].getMaskedImage(),
                                                   *fams[filterName],
                                                    background=-0.1, imBbox=imBbox).makeMosaic(display=None)

    for bands in [bands]:
        B, G, R = bands
        rgb = afwRgb.makeRGB(images[R], images[G], images[B], min, max - min, Q)

        afwRgb.displayRGB(rgb, show=True)

        if rgbFileFmt:
            afwRgb.writeRGB(rgbFileFmt % "".join(bands), rgb)

def showBlend(calexp, families, key='d', background=0.0, display=afwDisplay.getDisplay(0),
              imageDisplay=None, mtv=False, cleanupCallbacks=True):
    """Show blends interactively on an afwDisplay

    \param calexp   Exposure containing objects of interest
    \param families A Families object
    \param key      Key to display the family under the cursor in Display display
    \param display The afwDisplay.Display to display the families
    \param imageDisplay  The afwDisplay.Display displaying calexp (see mtv)
    \param mtv      If true, display calexp on display
    \param cleanupCallbacks If True reset afwDisplay to its initial state upon exit
                            (False is only useful for debugging errors in callbacks)

E.g.
import lsst.daf.persistence as dafPersist
import lsst.analysis.deblender as deblender

butler = dafPersist.Butler("/home/astro/hsc/hsc/HSC/rerun/rhl/tmp")
did = dict(visit=905518, ccd=31)
calexp = butler.get("calexp", **did)
ss = butler.get("src", **did)
families = deblender.Families(ss, butler, nChildMin=0)
deblender.showBlend(calexp, families, display=afwDisplay.Display(1))

Then hit 'key' (default: d) on objects of interest; 'r' for an rgb image
"""
    if display is not None:
        if mtv:
            imageDisplay.mtv(calexp)

    old = {}
    try:
        for k in set(list("ha") + [key]):
            old[k] = display.setCallback(k)

        display.setCallback(key, makeDisplayFamily(calexp, families, display=display))
        def new_h(*args):
            old['h'](*args)
            print("   1,2,4,8: Zoom to specified scale")
            print("   a:       show All the pixels")
            print("   %s:      show family under the cursor and return to python prompt" % key)
            print("   l:       cycle through stretch types")
        display.setCallback('h', new_h)
        display.setCallback('a', lambda k, x, y: display.zoom("to fit"))
        for z in [1, 2, 4, 8]:
            def _zoom(k, x, y, z=z):
                """Zoom by %d""" % z
                display.zoom(z)
            display.setCallback('%d' % z, _zoom)
        
        def callbackLog(k, x, y, i=[0]):
            """Cycle through stretches"""
            i[0] = (i[0] + 1)%3
            if i[0] == 0:
                display.scale("log", "minmax")
            elif i[0] == 1:
                display.scale("linear", "minmax")
            elif i[0] == 2:
                display.scale("linear", "zscale")
        display.setCallback('l', callbackLog)

        display.setCallback('r', makeDisplayFamily(calexp, families, rgb=True))

        display.interact()
    except Exception as e:
        print("Error in callback: %s" % e)
    finally:
        # Reset callbacks
        if cleanupCallbacks:
            for k, func in old.items():
                display.setCallback(k, func)

