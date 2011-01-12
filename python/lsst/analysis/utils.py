"""
Some utilities for looking at the outputs of processing runs

Currently specialised to handle ImSim outputs --- but this restriction should be lifted!
"""
import array, math, os, re, sys
import numpy as np
try:
    import matplotlib.pyplot as pyplot
except ImportError:
    pyplot = None
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
try:
    log
except:
    log = pexLogging.Log.getDefaultLog()
    log.setThreshold(pexLogging.Log.DEBUG);

import lsst.afw.image as afwImage
from lsst.afw.coord import DEGREES
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms.utils as maUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as ds9Utils

import lsst.daf.persistence as dafPersist
from lsst.obs.lsstSim import LsstSimMapper

class Data(object):
    class Id(object):
        """The identifier for a CCD-level dataset, as stored in PerCcdData"""
        def __init__(self, visit=None, raft=None, sensor=None):
            self.visit = visit
            self.raft = raft
            self.sensor = sensor

        def __repr__(self):
            return "(%d, '%s', '%s')" % (self.visit, self.raft, self.sensor)

    def __init__(self, *args, **kwargs):
        self.setButler(*args, **kwargs)
        self.ZP0 = 25
        self.fig = None
        self.visit = None

    def setButler(self, run=None,
                  dataRoot=None, dataRootFormat="/lsst2/datarel-runs/pt1prod_im%04d/update",
                  registryRoot=None, registryRun=None):
        if run is None:
            self.run = None
        else:
            self.run = run
            self.butler = None

        try:
            if butler is not None:
                return butler
        except NameError:
            butler = None

        if run is None and dataRoot is None and not butler:
            raise RuntimeError("Please specify a run or root")

        if run is None:
            if dataRoot is None:
                return butler
        else:
            if not dataRoot:
                dataRoot = dataRootFormat % run

        if registryRoot:
            assert(registryRun is None)
        else:
            if registryRun is None and run:
                registryRun = run

            if registryRun:
                registryRoot = dataRootFormat % registryRun
            else:
                registryRoot = dataRoot

        registry = None
        if registryRoot:
            registryFileName = "registry.sqlite3"
            for d in (".", "update"):
                f = os.path.join(registryRoot, d, registryFileName)
                if os.path.exists(f):
                    registry = f
                    break
            if not registry:
                raise RuntimeError("I'm unable to find your registry in %s", registryRoot)

        bf = dafPersist.ButlerFactory(mapper=LsstSimMapper(root=dataRoot, registry=registry))
        self.butler = bf.create()

    def lookupDataBySkytile(self, dataType):
        dataSets = {}
        for st, v, f, r, s in self.butler.queryMetadata("raw", "skyTile",
                                                        ["skyTile", "visit", "filter", "raft", "sensor"]):
            if self.butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s):
                if not dataSets.has_key(st):
                    dataSets[st] = {}
                if not dataSets[st].has_key((v, f)):
                    dataSets[st][(v, f)] = []
                dataSets[st][(v, f)].append((r, s))

        self.dataSets = dataSets
        
        return self.dataSets

    def lookupDataByVisit(self, dataType, *args, **kwargs):
        """Lookup all the available data of a given type (e.g. "psf") and maybe visit/raft/sensor.

        See also getDataset()
        """

        kwargs = dict([(k, v) for k, v in kwargs.items() if v is not None])

        dataSets = {}
        for v, f, r, s in self.butler.queryMetadata("raw", "visit",
                                                    ["visit", "filter", "raft", "sensor",],
                                                    *args, **kwargs):
            if self.butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s):
                if not dataSets.has_key(v):
                    dataSets[v] = []
                dataSets[v].append((r, s))

        self.dataSets = dataSets
        return self.dataSets

    def getDataset(self, dataType, visit=None, raft=None, sensor=None, ids=True):
        """Get all the data of the given type (e.g. "psf"); visit may be None (meaning use default);
raft or sensor may be None (meaning get all)

N.b. This routine resets the self.ids list unless ids is False; it is assumed that you retrieve all data items from the same set of sensors, but this is not checked.
"""
        if visit:
            self.visit = visit

        dataSets = self.lookupDataByVisit(dataType, visit=visit, raft=raft, sensor=sensor)
        if not self.visit and len(dataSets) == 1:
            self.visit = dataSets.keys()[0]

        dataSets = dataSets.get(self.visit)
        if not dataSets:
            raise RuntimeError("I cannot find your data" +
                               (" for visit %d" % self.visit if self.visit else ""))

        data = []
        if ids:
            self.ids = []
        for raft, sensor in dataSets:
            data.append(self.butler.get(dataType, visit=self.visit, raft=raft, sensor=sensor))
            if ids:
                self.ids.append(Data.Id(visit=self.visit, raft=raft, sensor=sensor))

        return data

    def getEimages(self, raft, sensor):
        pass

    def getPsfs(self, *args, **kwargs):
        return self.getDataset("psf", *args, **kwargs)
    
    def getSources(self, *args, **kwargs):
        return [pss.getSources() for pss in self.getDataset("src", *args, **kwargs)]

    def getMags(self, ss, calculateApCorr=False):
        """Return numpy arrays constructed from SourceSet ss"""
        apMags = np.empty(len(ss))
        psfMags = np.empty(len(ss))
        flags = np.empty(len(ss))

        for i in range(len(ss)):
            apMags[i]  = ss[i].getApFlux()
            psfMags[i] = ss[i].getPsfFlux()
            flags[i] = ss[i].getFlagForDetection()

        apMags  = self.ZP0 - 2.5*np.log10(apMags)
        psfMags = self.ZP0 - 2.5*np.log10(psfMags)

        good = np.logical_and(np.isfinite(apMags), np.isfinite(psfMags))

        apMags = apMags[good]
        psfMags = psfMags[good]
        flags = flags[good]

        if calculateApCorr:
            delta = psfMags - apMags
            apCorr = np.median(delta[psfMags < 12])
            print "RHL", apCorr
            apMags += apCorr

        return apMags, psfMags, flags

    def _getMags(self, ss, apMags=None, psfMags=None, flags=None):
        """Return python arrays constructed from SourceSet ss, and possibly extending apMags/psfMags"""
        if not apMags:
            apMags = array.array('f')
        if not psfMags:
            psfMags = array.array('f')
        if not flags:
            flags = array.array('L')

        _apMags, _psfMags, _flags = self.getMags(ss, True)

        apMags.extend(_apMags)
        psfMags.extend(_psfMags)
        flags.extend(_flags)

        return apMags, psfMags, flags

    def getMagsByVisit(self, visit=None, nSensor=0, sensor=None, raft=None):
        """Set the "self.magnitudes" for a visit"""
        d = self.lookupDataByVisit("src", visit=visit, sensor=sensor, raft=raft)

        if visit:
            self.visit = visit

        if not self.visit and len(d) > 0:
            self.visit = d.keys()[0]

        if not self.visit:
            raise RuntimeError("Please specify a visit")

        apMags, psfMags, flags = None, None, None

        if nSensor == 0:
            nSensor = len(d[self.visit])     # i.e. all

        for r, s in d[self.visit][0:nSensor]:
            if raft and r != raft or sensor and sensor != s:
                continue
            ss = self.getSources(raft=r, sensor=s)[0]
            apMags, psfMags, flags = self._getMags(ss, apMags=apMags, psfMags=psfMags, flags=flags)

        self.apMags  = np.ndarray(len(apMags),  dtype='f', buffer=apMags)
        self.psfMags = np.ndarray(len(psfMags), dtype='f', buffer=psfMags)
        self.flags   = np.ndarray(len(flags),   dtype=long, buffer=flags)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, fixup=False, **kwargs):
        if kwargs.has_key("visit"):
            self.visit = kwargs["visit"]
            del kwargs["visit"]

        try:
            self.visit = int(re.sub(r"^v", "", self.visit))
        except:
            pass

        self.name = 'imsim-v%d-r%s-s%s' % (self.visit,
                                          kwargs["raft"].replace(',' ,''), kwargs["sensor"].replace(',', ''))

        try:
            filters = self.butler.queryMetadata('raw', 'filter', visit=self.visit, **kwargs)
            filterName = filters[0]
        except IndexError:
            raise RuntimeError("No filters are available for visit %d" % self.visit)

        psources = self.butler.get('icSrc', visit=self.visit, **kwargs)
        pmatches = self.butler.get('icMatch', visit=self.visit, **kwargs)
        if True:
            print 'Got sources', psources
            print 'Got matches', pmatches

        matchmeta = pmatches.getSourceMatchMetadata()
        self.matches = pmatches.getSourceMatches()
        if True:
            print 'Match metadata:', matchmeta
        self.sources = psources.getSources()

        if False:                       # only read the metadata; not quite enough unfortunately
            calexp_md = self.butler.get('calexp_md', visit=self.visit, **kwargs)
            wcs = afwImage.makeWcs(calexp_md)
        else:
            calexp = self.butler.get('calexp', visit=self.visit, **kwargs)
            wcs = calexp.getWcs()

        if False:                       # why do this?
            wcs = afwImage.cast_TanWcs(wcs)
        if False:
            print 'Got wcs', wcs.getFitsMetadata().toString()

        self.zp = calexp.getCalib().getMagnitude(1.0)
        print 'Zeropoint is', self.zp

        # ref sources
        W, H = calexp.getWidth(), calexp.getHeight()
        xc, yc = 0.5*W, 0.5*H
        radec = wcs.pixelToSky(xc, yc)
        ra = radec.getLongitude(DEGREES)
        dec = radec.getLatitude(DEGREES)
        radius = wcs.pixelScale()*math.hypot(xc, yc)*1.1
        print 'Image W,H', W,H
        print 'Image center RA,Dec', ra, dec
        print 'Searching radius', radius, 'arcsec'
        pol = pexPolicy.Policy()
        pol.set('matchThreshold', 30)
        solver = measAstrom.createSolver(pol, log)
        idName = 'id'
        # could get this from matchlist meta...
        anid = matchmeta.getInt('ANINDID')
        if False:
            print 'Searching index with ID', anid
            print 'Using ID column name', idName
            print 'Using filter column name', filterName
        X = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
        self.ref = X.first
        inds = X.second

        referrs, stargal = None, None
        print 'Tag-along columns:'
        cols = solver.getTagAlongColumns(anid)
        print cols
        for c in cols:
            print 'column: ', c.name, c.fitstype, c.ctype, c.units, c.arraysize
        colnames = [c.name for c in cols]

        if False:
            col = filterName + '_err'
            if col in colnames:
                referrs = solver.getTagAlongDouble(anid, col, inds)

        col = 'starnotgal'
        if col in colnames:
            stargal1 = solver.getTagAlongBool(anid, col, inds)
            stargal = []
            for i in range(len(stargal1)):
                stargal.append(stargal1[i])

        if False:
            print 'Got', len(self.ref), 'reference catalog sources'

        fdict = maUtils.getDetectionFlags()

        keepref = []
        keepi = []
        for i in xrange(len(self.ref)):
            x, y = wcs.skyToPixel(self.ref[i].getRa(), self.ref[i].getDec())
            if x < 0 or y < 0 or x > W or y > H:
                continue
            self.ref[i].setXAstrom(x)
            self.ref[i].setYAstrom(y)
            if stargal[i]:
                self.ref[i].setFlagForDetection(self.ref[i].getFlagForDetection() | fdict["STAR"])
            keepref.append(self.ref[i])
            keepi.append(i)
            
        print 'Kept', len(keepref), 'reference sources'
        self.ref = keepref

        if referrs is not None:
            referrs = [referrs[i] for i in keepi]
        if stargal is not None:
            stargal = [stargal[i] for i in keepi]

        self.stargal = stargal
        self.referrs = referrs

        if False:
            m0 = self.matches[0]
            f,s = m0.first, m0.second
            print 'match 0: ref %i, source %i' % (f.getSourceId(), s.getSourceId())
            print '  ref x,y,flux = (%.1f, %.1f, %.1f)' % (f.getXAstrom(), f.getYAstrom(), f.getPsfFlux())
            print '  src x,y,flux = (%.1f, %.1f, %.1f)' % (s.getXAstrom(), s.getYAstrom(), s.getPsfFlux())

        measAstrom.joinMatchList(self.matches, self.ref, first=True, log=log)
        args = {}
        if fixup:
            # ugh, mask and offset req'd because source ids are assigned at write-time
            # and match list code made a deep copy before that.
            # (see svn+ssh://svn.lsstcorp.org/DMS/meas/astrom/tickets/1491-b r18027)
            args['mask'] = 0xffff
            args['offset'] = -1
        measAstrom.joinMatchList(self.matches, self.sources, first=False, log=log, **args)

        if False:
            for m in self.matches:
                x0,x1 = m.first.getXAstrom(), m.second.getXAstrom()
                y0,y1 = m.first.getYAstrom(), m.second.getYAstrom()
                print 'x,y, dx,dy', x0, y0, x1-x0, y1-y0

        if False:
            m0 = self.matches[0]
            f,s = m0.first, m0.second
            print 'match 0: ref %i, source %i' % (f.getSourceId(), s.getSourceId())
            print '  ref x,y,flux = (%.1f, %.1f, %.1f)' % (f.getXAstrom(), f.getYAstrom(), f.getPsfFlux())
            print '  src x,y,flux = (%.1f, %.1f, %.1f)' % (s.getXAstrom(), s.getYAstrom(), s.getPsfFlux())
            r,d = 2.31262000000000, 3.16386000000000
            x,y = wcs.skyToPixel(r,d)
            print 'x,y', x,y
            r2d2 = wcs.pixelToSky(x,y)
            r2 = r2d2.getLongitude(DEGREES)
            d2 = r2d2.getLatitude(DEGREES)
            print r,d
            print r2,d2

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmag(data, fig=None, markersize=0.1, color="red"):
    """Plot (aper - psf) v. psf mags"""
    if fig:
        data.fig = fig

    if data.fig:
        data.fig.clf()
    else:
        if not pyplot:
            raise RuntimeError("I am unable to plot as I failed to import matplotlib")
        data.fig = pyplot.figure()

    axes = data.fig.add_axes((0.1, 0.1, 0.85, 0.80));

    fdict = maUtils.getDetectionFlags()

    bad = np.bitwise_and(data.flags, fdict["INTERP_CENTER"])
    good = np.logical_not(bad)
    a = data.apMags[good]
    p = data.psfMags[good]

    axes.plot(p, a - p, "o", markersize=markersize, color=color)
    axes.plot((0, 20), (0, 0), "b-")
    axes.set_ylim(-0.4, 0.6)
    axes.set_xlim(12, 18)
    axes.set_xlabel("psf")
    axes.set_ylabel("aper - psf")
    title = "Visit %d, all CCDs, " % (data.visit)
    if data.run:
        title = "pt1prod_im%04d, " % (data.run)
    title += "per-CCD aperture correction"
    axes.set_title(title)

    data.fig.show()

def plotCalibration(data, plotBand=0.05, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data
    
If plotBand is provided, draw lines at +- plotBand
    """
    if fig:
        data.fig = fig

    if data.fig:
        data.fig.clf()
    else:
        if not pyplot:
            raise RuntimeError("I am unable to plot as I failed to import matplotlib")
        data.fig = pyplot.figure()

    fdict = maUtils.getDetectionFlags()

    mstars = [m for m in data.matches if (m.second.getFlagForDetection() & fdict["STAR"])] # data
    realStars = [(m.first.getFlagForDetection() & fdict["STAR"]) != 0                      # catalogue
                 for m in data.matches if (m.second.getFlagForDetection() & fdict["STAR"])]

    axes = data.fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars])
    instmag = np.array([data.zp - 2.5*math.log10(s.second.getPsfFlux()) for s in mstars])
    realStars = np.array([m for m in realStars])

    delta = refmag - instmag
    #markersize 
    if False:                                                            # plot top/bottom of chip differently
        top = np.array([(m.second.getYAstrom() > 2048) for m in mstars]) # actually, top of chip

        axes.plot(refmag[top], delta[top], "r+")
        axes.plot(refmag[np.logical_not(top)], delta[np.logical_not(top)], "b+")
    else:
        axes.plot(refmag, delta, "r+")

    axes.plot(refmag[realStars], delta[realStars], "mo")
        
    axes.plot((-100, 100), (0, 0), "g-")
    if plotBand:
        for x in (-plotBand, plotBand):
            axes.plot((-100, 100), x*np.ones(2), "g--")

    axes.set_ylim(-1.1, 1.1)
    axes.set_xlim(24, 13)
    axes.set_xlabel("Reference")
    axes.set_ylabel("Reference - Instrumental")
    axes.set_title(data.name)

    data.fig.show()

    if False:
        _mstars = [m for m in mstars if math.fabs(-2.5*math.log10(m.first.getPsfFlux()) - 17.5) < 0.5]
        for m in _mstars:
            print "%.1f %.1f  %.2f %.2f" % (m.second.getXAstrom(), m.second.getYAstrom(),
                                            -2.5*math.log10(m.first.getPsfFlux()),
                                            -2.5*math.log10(m.first.getPsfFlux()) - 
                                            (data.zp - 2.5*math.log10(m.second.getPsfFlux())))
            ds9.dot("*", m.second.getXAstrom(), m.second.getYAstrom(), ctype=ds9.MAGENTA, size=10)

def displayCalibration(data, frame=0):
    """display the calibration objects in Data object data on ds9
Stars are green; galaxies are red based on our processing
    """
    fdict = maUtils.getDetectionFlags()

    for m in data.matches:
        ref, src = m.first, m.second
        x1, y1 = src.getXAstrom(), src.getYAstrom()
        if ref.getFlagForDetection() & fdict["STAR"]:
            ptype = "+"
            if src.getFlagForDetection() & fdict["STAR"]:
                ctype = ds9.GREEN
            else:
                ctype = ds9.YELLOW
        else:
            ptype = "o"
            if not (src.getFlagForDetection() & fdict["STAR"]):
                ctype = ds9.RED
            else:
                ctype = ds9.MAGENTA

        ds9.dot(ptype, x1, y1, ctype=ctype)

def showPsfs(data, psfs, frame=None):
    mos = ds9Utils.Mosaic()
    for i in range(len(psfs)):
        psf = psfs[i]
        mos.append(psf.computeImage(), "%s:%s" % (data.ids[i].raft, data.ids[i].sensor))

    mos.makeMosaic(frame=frame)
