"""Utilities to look at deblender outputs"""

import math, re
import numpy as np

import utils
import lsst.pex.logging as pexLog
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

def getEllipses(src, nsigs=[1.], **kwargs):
    xc = src.getX()
    yc = src.getY()
    x2 = src.getIxx()
    y2 = src.getIyy()
    xy = src.getIxy()
    # SExtractor manual v2.5, pg 29.
    a2 = (x2 + y2)/2. + np.sqrt(((x2 - y2)/2.)**2 + xy**2)
    b2 = (x2 + y2)/2. - np.sqrt(((x2 - y2)/2.)**2 + xy**2)
    theta = np.rad2deg(np.arctan2(2.*xy, (x2 - y2)) / 2.)
    a = np.sqrt(a2)
    b = np.sqrt(b2)
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


def getId(src, field="objId"):
    idDict = utils.butler.mapperInfo.splitId(src.getId(), asDict=True)

    return idDict[field] if field else idDict

def findFamily(families, objId):
    """Return the object's family (either parent or child)"""

    family = [f for f in families if (f[0].getId() & 0xffff) == objId]
    if family:
        return family[0]

    for family in families:
        for child in family[1]:
            if (child.getId() & 0xffff) == objId:
                return family

    return None
    

# Real thing: make plots given the Sources
def plotDeblendFamily(mi, parent, kids, dkids=[],
                      background=-10, symbolSize=2,
                      plotb=False, ellipses=True,
                      arcsinh=True, maskbit=False, frame=0):

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
    mos.append(parent_im, '%dP' % getId(parent))

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
        mos.append(_kim, '%dC' % getId(kid)); del _kim

    title = re.sub(r"[{}']", "", str(getId(parent, None))) # ds9 doesn't handle those chars well
    mosaicImage = mos.makeMosaic(frame=frame, title=title)
    ds9.dot("%s  (%.1f, %1.f)" % (title, parent.getX(), parent.getY()),
            0.5*mosaicImage.getWidth(), 5 + 1.05*mosaicImage.getHeight(), frame=frame,
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

def getFamilies(cat, nChildMin=None):
    '''
    Returns [ (parent0, [child0, child1]), (parent1, [child0, ...]), ...]
    where parents are sorted by ID.  Only objects that are deblended are included (unless nChildMin == 0)

    if nChildMin is not None, include only deblends with at least than many children.  As a
    special case, nChildMin == 0 includes all objects, even those that aren't blended
    '''
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
    return [(cat.find(pid), children[pid]) for pid in parentIds]
