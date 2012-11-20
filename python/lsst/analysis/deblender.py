"""Utilities to look at deblender outputs

Based on Dustin Lang's code in meas_deblender/examples, but heavily rewritten since then, so don't blame him

E.g.
import lsst.afw.display.ds9 as ds9
import lsst.analysis.utils as utils

da = utils.Data(dataRoot=os.path.join(os.environ["SUPRIME_DATA_DIR"], "SUPA"),
                Mapper=utils.SuprimecamMapperMit, rerun="rhl/mit");
dataId = utils.DataId(visit=24473, ccd=3);
calexp = da.getDataset("calexp", dataId)[0]
ss = da.getDataset('src', dataId)[0]
families = deblender.Families(ss, nChildMin=0)

frame=1
ds9.mtv(calexp, title=dataId, frame=frame)
utils.showSourceSet(ss, frame=frame)
utils.showSourceSet(ss, frame=frame, symb='id', size=2.5, fontFamily="times", ctype=ds9.GREEN)

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

def plotDeblendFamily(mi, parent, kids, dkids=[],
                      background=-10, symbolSize=2,
                      plotb=False, ellipses=True,
                      arcsinh=True, maskbit=False, frame=0):
    """Display a deblend on ds9"""
    
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

    mos = displayUtils.Mosaic(background=background)
    mos.append(parent_im, '%dP' % utils.butler.mapperInfo.getId(parent))

    for i, kid in enumerate(kids):
        kim = footprintToImage(kid.getFootprint(), mi, **aa)
        peak = kid.getFootprint().getPeaks()[0]
        #
        # Put the child into the correct place in the parent image
        #
        _kim = parent_im.Factory(parent_im.getDimensions())
        _kim.setXY0(parent_im.getXY0())

        bbox = afwGeom.Box2I(afwGeom.Point2I(kim.getX0() - parent_im.getX0(),
                                             kim.getY0() - parent_im.getY0()), kim.getDimensions())
        
        sim = kim.Factory(_kim, bbox)
        sim <<= kim
        mos.append(_kim, '%dC' % utils.butler.mapperInfo.getId(kid)); del _kim

    title = re.sub(r"[{}']", "",
                   str(utils.butler.mapperInfo.getId(parent, None))) # ds9 doesn't handle those chars well
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
                ds9.dot("+", p.getFx() + x0, p.getFy() + y0, frame=frame,
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
    else:
        fp = afwDet.makeHeavyFootprint(fp, mi)
    bb = fp.getBBox()
    if mask:
        im = afwImage.MaskedImageF(bb.getWidth(), bb.getHeight())
    else:
        im = afwImage.ImageF(bb.getWidth(), bb.getHeight())
    im.setXY0(bb.getMinX(), bb.getMinY())
    fp.insert(im)
    if mask:
        im = im.getMask()
    return im

class Families(list):
    def __init__(self, cat, nChildMin=None):
        '''
        Returns [ (parent0, [child0, child1]), (parent1, [child0, ...]), ...]
        where parents are sorted by ID.  Only objects that are deblended are included (unless nChildMin == 0)

        if nChildMin is not None, include only deblends with at least than many children.  As a
        special case, nChildMin == 0 includes all objects, even those that aren't blended
        '''
        self.cat = cat

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

            objId = utils.butler.mapperInfo.splitId(matched[0][0].getId(), asDict=True)["objId"]

        family = [f for f in self if utils.butler.mapperInfo.getId(f[0]) == objId]
        if family:
            return family[0]

        for family in self:
            for child in family[1]:
                if utils.butler.mapperInfo.getId(child) == objId:
                    return family

        return None

# backwards compatibility
    
def getFamilies(cat, nChildMin=None):
    return Families(cat, nChildMin)

def findFamily(families, objId):
    """Return the object's family (you may specify either the ID for the parent or a child)"""

    return families.find(objId)

def makeDisplayFamily(calexp, families, matchRadius=20):
    """Factory function for callback function implementing showBlend"""
    def display_family(k, x, y):
        fam = families.find((x, y), matchRadius=matchRadius)
        if fam:
            plotDeblendFamily(calexp, *fam, background=1000)
            #return True

    return display_family

def showBlend(calexp, families, frame=None, key='d', mtv=False):
    if frame is not None:
        if mtv:
            ds9.mtv(calexp, frame=frame)
        ds9.ds9Cmd(ds9.selectFrame(frame))

    old = {}
    try:
        for k in "hdp":
            old[k] = ds9.setCallback(k)

        ds9.setCallback(key, makeDisplayFamily(calexp, families))
        def new_h(*args):
            old['h'](*args)
            print "   d: show family under the cursor and return to python prompt"
            print "   p: pan to this point"
        ds9.setCallback('h', new_h)
        ds9.setCallback('p', lambda k, x, y: ds9.pan(x, y, frame=frame))
        
        ds9.interact()
    except Exception, e:
        print "RHL", e
    finally:
        print "Cleaning up"
        for k, func in old.items():
            ds9.setCallback(k, func)
