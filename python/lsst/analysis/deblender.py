"""Utilities to look at deblender outputs

Based on Dustin Lang's code in meas_deblender/examples, but heavily rewritten since then, so don't blame him

E.g.
import lsst.daf.persistence as dafPersist
import lsst.afw.display.ds9 as ds9
import lsst.analysis.utils as utils
import lsst.analysis.deblender as deblender

butler = dafPersist.Butler(os.path.join(os.environ["SUPRIME_DATA_DIR"], "rerun", "rhl", "realFlats"))
dataId = dict(visit=905518, ccd=31)
calexp = butler.get("calexp", **dataId)
ss = butler.get('src', **dataId)
families = deblender.Families(ss, butler, nChildMin=0)

frame=1            # We'll use 0 for the children
ds9.mtv(calexp, title="%(visit)s %(ccd)s" % dataId, frame=frame)
utils.showSourceSet(ss, frame=frame)
utils.showSourceSet(ss, mapperInfo=families.mapperInfo, frame=frame, symb='id', size=1.0, fontFamily="times", ctype=ds9.GREEN)
#
# This puts us into a loop waiting on ds9.  Sigh.
# Use d on an object in frame 1 to show the children in frame 0; q to quit the ds9 interactive loop
#
deblender.showBlend(calexp, families, frame=0)

"""

import math, re, sys

import utils
import lsst.pex.logging as pexLog
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.display.ds9 as ds9
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

def plotDeblendFamily(mi, parent, kids, mapperInfo=None, dkids=[],
                      background=-10, symbolSize=2,
                      plotb=False, ellipses=True,
                      arcsinh=True, maskbit=False, frame=0):
    """Display a deblend on ds9

Each child is marked with a + at its centre (green if deblended-as-psf else red)
all the other peaks in its footprint are marked with x (cyan if deblended-as-psf else magenta)
    """

    if mi:
        try:
            mi = mi.getMaskedImage()        # maybe it's an Exposure?
        except AttributeError:
            pass

    aa = {}
    if maskbit:
        aa.update(mask=True)
    parent_im = footprintToImage(parent.getFootprint(), mi, **aa)
    bbox = afwGeom.BoxD(parent.getFootprint().getBBox())
    pext = (bbox.getMinX(), bbox.getMaxX(), bbox.getMinY(), bbox.getMaxY())

    pks = parent.getFootprint().getPeaks()
    pix = [pk.getIx() for pk in pks]
    piy = [pk.getIy() for pk in pks]
    pfx = [pk.getFx() for pk in pks]
    pfy = [pk.getFy() for pk in pks]
    if ellipses:
        ell = (parent.getX(), parent.getY(), parent.getIxx(), parent.getIyy(), parent.getIxy())

    N = 1 + len(kids)
    S = math.ceil(math.sqrt(N))
    C = S
    R = math.ceil(float(N)/C)

    Rx,Ry = [],[]
    tts = []
    stys = []
    xys = []
    #
    # Find how large an image we need to display the parent and all the children
    #
    imBbox = parent_im.getBBox(afwImage.PARENT)
    kidImages, kim = {}, None
    for kid in kids:
        kim = footprintToImage(kid.getFootprint(), mi, **aa)
        kidImages[kid] = kim
        
        imBbox.include(kim.getBBox(afwImage.PARENT))

    mos = displayUtils.Mosaic(background=background)

    if not kim:
        kim = parent_im.clone()
        
    bbox = afwGeom.Box2I(afwGeom.Point2I(kim.getX0() - imBbox.getMinX(),
                                         kim.getY0() - imBbox.getMinY()), kim.getDimensions())

    kidImages[parent] = parent_im       # not strictly a kid

    for kid in [parent] + kids:
        kim = kidImages[kid]
        #
        # Put the child into the correct place in the parent image.  We have to do this for
        # the parent too if some of the children extended outside its BBox
        #
        bbox = afwGeom.Box2I(afwGeom.Point2I(kim.getX0() - imBbox.getMinX(),
                                             kim.getY0() - imBbox.getMinY()), kim.getDimensions())

        _kim = parent_im.Factory(imBbox)
        _kim[bbox] <<= kim
        mos.append(_kim, '%d%s' % (mapperInfo.getId(kid) if mapperInfo else (kid.getId() & 0xfff),
                                   "P" if kid == parent else "C"))
        del _kim

    if mapperInfo:
        title = re.sub(r"[{}']", "", str(mapperInfo.getId(parent, None))) # ds9 doesn't handle those chars well
    else:
        title = "0x%x == %d" % (parent.getId(), (parent.getId() & 0xffff))
    mosaicImage = mos.makeMosaic(frame=frame, title=title)
    ds9.dot("%s  (%.1f, %1.f)" % (title, parent.getX(), parent.getY()),
            0.5*mosaicImage.getWidth(), 1.03*mosaicImage.getHeight(), frame=frame,
            ctype=ds9.BLACK, fontFamily="times", size=3)

    with ds9.Buffering():
        for i, src in enumerate([parent] + kids):    
            x0, y0 = mos.getBBox(i).getMin()
            x0 -= parent_im.getX0(); y0 -= parent_im.getY0()          

            if src.get("deblend.deblended-as-psf"):
                centroid_ctype = ds9.GREEN
                peak_ctype = ds9.CYAN
            else:
                centroid_ctype = ds9.RED
                peak_ctype = ds9.MAGENTA
            
            ds9.dot("+", src.getX() + x0, src.getY() + y0, frame=frame,
                    size=symbolSize, ctype=centroid_ctype)
            for p in src.getFootprint().getPeaks():
                ds9.dot("x", p.getFx() + x0, p.getFy() + y0, frame=frame,
                        size=0.5*symbolSize if i == 0 else symbolSize,
                        ctype=ds9.YELLOW if i == 0 else peak_ctype)

        if False:
            if len(kid.flags):
                tt += ', ' + ', '.join(kid.flags)

        if False:
            # peak(s)
            plt.plot(kid.pfx, kid.pfy, 'x', **sty2)
            xys.append((kid.pfx, kid.pfy, sty2))
            # centroid
            plt.plot([kid.cx], [kid.cy], 'x', **sty1)
            xys.append(([kid.cx], [kid.cy], sty1))
            # ellipse
            if ellipses and not kid.ispsf:
                drawEllipses(plt, kid, ec=sty1['color'], fc='none', alpha=0.7)
            if plotb:
                plt.axis(ext)
            else:
                plt.axis(pax)


    # add child centers and ellipses...
    if False:
        for x,y,sty in xys:
            plt.plot(x, y, 'x', **sty)
    if ellipses:
        for kid,sty in zip(kids,stys):
            if kid.ispsf:
                continue
            drawEllipses(plt, kid, ec=sty['color'], fc='none', alpha=0.7)
    if False:
        plt.plot([parent.getX()], [parent.getY()], 'x', color='b')
        if ellipses:
            drawEllipses(plt, parent, ec='b', fc='none', alpha=0.7)

    # Plot dropped kids
    for kid in dkids:
        ext = kid.ext
        # bounding box
        xx = [ext[0],ext[1],ext[1],ext[0],ext[0]]
        yy = [ext[2],ext[2],ext[3],ext[3],ext[2]]
        plt.plot(xx, yy, 'y-')
        # peak(s)
        plt.plot(kid.pfx, kid.pfy, 'yx')
    
    #plt.axis(pax)
    

def footprintToImage(fp, mi=None, mask=False):
    if fp.isHeavy():
        fp = afwDet.cast_HeavyFootprintF(fp)
    elif mi is None:
        print >> sys.stderr, "Unable to make a HeavyFootprint as image is None"
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
            centroidName = self.cat.table.getCentroidDefinition()
            oneObjCatalog.table.defineCentroid(centroidName)

            s = oneObjCatalog.addNew()
            s.set("%s.x" % centroidName, x)
            s.set("%s.y" % centroidName, y)

            matched = afwTable.matchXy(self.cat, oneObjCatalog, matchRadius)

            if len(matched) == 0:
                print >> sys.stderr, "Unable to find object at (%.2f, %.2f)" % (x, y)
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

def makeDisplayFamily(calexp, families, matchRadius=20, background=-10, frame=None):
    """Factory function for callback function implementing showBlend"""
    def display_family(k, x, y):
        x0, y0 = calexp.getXY0()
        fam = families.find((x + x0, y + y0), matchRadius=matchRadius)
        if fam:
            plotDeblendFamily(calexp, *fam, mapperInfo=families.mapperInfo, background=background, frame=frame)

    return display_family

def showBlend(calexp, families, key='d', frame0=0, frame=None, mtv=False):
    """Show blends interactively on ds9

    \param calexp   Exposure containing objects of interest
    \param families A Families object
    \param key      Key to display the family under the cursor in frame frame
    \param frame0   The ds9 frame to display the families
    \param frame    The ds9 frame displaying calexp (see mtv)
    \param mtv      If true, display calexp in ds9 frame frame

E.g.
import lsst.daf.persistence as dafPersist
import lsst.analysis.deblender as deblender

butler = dafPersist.Butler("/home/astro/hsc/hsc/HSC/rerun/rhl/tmp")
did = dict(visit=905518, ccd=31)
calexp = butler.get("calexp", **did)
ss = butler.get("src", **did)
families = deblender.Families(ss, butler, nChildMin=0)
deblender.showBlend(calexp, families, frame=1)

Then hit 'key' (default: d) on objects of interest
"""
    if frame is not None:
        if mtv:
            ds9.mtv(calexp, frame=frame)
        ds9.ds9Cmd(ds9.selectFrame(frame))

    old = {}
    try:
        for k in set(list("ha") + [key]):
            old[k] = ds9.setCallback(k)

        ds9.setCallback(key, makeDisplayFamily(calexp, families, frame=frame0))
        def new_h(*args):
            old['h'](*args)
            print "   1,2,4,8: Zoom to specified scale"
            print "   a:       show All the pixels"
            print "   %s:      show family under the cursor and return to python prompt" % key
            print "   l:       cycle through sretch types"
        ds9.setCallback('h', new_h)
        ds9.setCallback('a', lambda k, x, y: ds9.ds9Cmd("zoom to fit", frame=frame0))
        for z in [1, 2, 4, 8]:
            def _zoom(k, x, y, z=z):
                ds9.zoom(z, frame=frame0)
            ds9.setCallback('%d' % z, _zoom)
        
        def callbackLog(k, x, y, i=[0]):
            """Cycle through stretches"""
            i[0] = (i[0] + 1)%3
            if i[0] == 0:
                ds9.ds9Cmd("scale log; scale minmax")
            elif i[0] == 1:
                ds9.ds9Cmd("scale linear; scale minmax")
            elif i[0] == 2:
                ds9.ds9Cmd("scale linear; scale zscale")

        ds9.setCallback('l', callbackLog)

        ds9.interact()
    except Exception, e:
        print "RHL", e
    finally:
        print "Cleaning up"
        for k, func in old.items():
            ds9.setCallback(k, func)
