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
import lsst.pex.exceptions as pexExcept
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
try:
    log
except:
    log = pexLogging.Log.getDefaultLog()
    #log.setThreshold(pexLogging.Log.DEBUG);
import MySQLdb
try:
    default_db_host
except NameError:
    default_db_host = "lsst10.ncsa.uiuc.edu"
    default_db_table = "rplante_DC3b_u_pt11final"
    default_db_table = "rplante_DC3b_u_weeklytest2_2011_0214_science"

try:
    butlerDataRoot
except NameError:
    butlerDataRoot = None
    butler = None

import lsst.afw.image as afwImage
from lsst.afw.coord import DEGREES
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as ds9Utils

import lsst.daf.persistence as dafPersist
from lsst.obs.lsstSim import LsstSimMapper

try:
    mpFigures
except NameError:
    mpFigures = {0 : None}              # matplotlib (actually pyplot) figures

def getMpFigure(fig=None, clear=True):
    """Return a pyplot figure(); if fig is supplied save it and make it the default
    fig may also be a bool (make a new figure) or an int (return or make a figure (1-indexed;
    python-list style -n supported)
    """

    if not pyplot:
        raise RuntimeError("I am unable to plot as I failed to import matplotlib")

    if isinstance(fig, bool):       # we want a new one
        fig = len(mpFigures) + 1    # matplotlib is 1-indexed

    if isinstance(fig, int):
        i = fig
        if i == 0:
            raise RuntimeError("I'm sorry, but matplotlib uses 1-indexed figures")
        if i < 0:
            try:
                i = sorted(mpFigures.keys())[i] # simulate list's [-n] syntax
            except IndexError:
                if mpFigures:
                    print >> sys.stderr, "Illegal index: %d" % i
                i = 1

        if not mpFigures.has_key(i):
            for j in range(1, i):
                getMpFigure(j)
                
            mpFigures[i] = pyplot.figure()
        fig = mpFigures[i]

    if not fig:
        i = sorted(mpFigures.keys())[0]
        if i > 0:
            fig = mpFigures[i[-1]]
        else:
            fig = getMpFigure(1)

    if clear:
        fig.clf()

    return fig

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

        self.visit = None
        self.raft = None
        self.sensor = None

        self.matches = None
        self.ZP0 = 31

    def setVRS(self, visit=None, raft=None, sensor=None, reset=False):
        """Save the values of visit, raft, sensor; if None use the pre-set values (which are cleared if reset is True)"""

        if reset:
            self.visit = None
            self.raft = None
            self.sensor = None

        if visit is None:
            visit = self.visit
        if raft is None:
            raft = self.raft
        if sensor is None:
            sensor = self.sensor

        try:
            visit = int(re.sub(r"^v", "", visit))
        except:
            pass

        self.visit = visit
        self.raft = raft
        self.sensor = sensor

    def setButler(self, run=None,
                  dataRoot=None, dataRootFormat="/lsst2/datarel-runs/pt1prod_im%04d/update",
                  registryRoot=None, registryRun=None):
        if run is None:
            self.run = None
        else:
            self.run = run
            self.butler = None

        global butlerDataRoot, butler
        if butlerDataRoot and dataRoot != butlerDataRoot:
            butler = None

        try:
            if butler is not None:
                self.butler = butler
                self.dataRoot = butlerDataRoot
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
        self.dataRoot = dataRoot

        butler, butlerDataRoot = self.butler, self.dataRoot # globals

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

        kwargs = dict([(k, v) for k, v in kwargs.items() if 1 or v is not None])
        extra = {}
        if dataType == "raw":
            for k in ("snap", "channel",):
                extra[k] = kwargs.get(k, 1)
                kwargs[k] = extra[k]

        if False:
            vfsr = self.butler.queryMetadata("raw", "visit",
                                             ["visit", "filter", "raft", "sensor",],
                                             *args, **kwargs)
        else:
            if kwargs.get("raft"):
                if kwargs.get("sensor"):
                    vfsr = self.butler.queryMetadata("raw", "visit",
                                                     ["visit", "filter", "raft", "sensor",],
                                                     visit=kwargs["visit"], raft=kwargs["raft"], sensor=kwargs["sensor"])
                else:
                    vfsr = self.butler.queryMetadata("raw", "visit",
                                                     ["visit", "filter", "raft", "sensor",],
                                                     visit=kwargs["visit"], raft=kwargs["raft"])
            else:
                    vfsr = self.butler.queryMetadata("raw", "visit",
                                                     ["visit", "filter", "raft", "sensor",],
                                                     visit=kwargs["visit"])

        dataSets = {}
        for v, f, r, s in vfsr:
            if self.butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s, **extra):
                if not dataSets.has_key(v):
                    dataSets[v] = []
                dataSets[v].append((r, s))

        self.dataSets = dataSets
        return self.dataSets

    def getDataset(self, dataType, visit=None, raft=None, sensor=None, ids=True, **kwargs):
        """Get all the data of the given type (e.g. "psf"); visit may be None (meaning use default);
raft or sensor may be None (meaning get all)

N.b. This routine resets the self.ids list unless ids is False; it is assumed that you retrieve all data items from the same set of sensors, but this is not checked.
"""
        self.setVRS(visit, raft, sensor, reset=True)

        dataSets = self.lookupDataByVisit(dataType, visit=self.visit, raft=self.raft, sensor=self.sensor)
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
            dataElem = self.butler.get(dataType, visit=self.visit, raft=self.raft, sensor=self.sensor)

            if dataType == "calexp" and kwargs.get("calibrate"):
                mi = dataElem.getMaskedImage();
                if kwargs.get("db"):
                    fluxMag0 = getFluxMag0DB(self.visit, raft, sensor)
                else:
                    fluxMag0 = dataElem.getCalib().getFluxMag0()[0]
                mi *= fluxMag0*1e-12
                del mi

            data.append(dataElem)
            if ids:
                self.ids.append(Data.Id(visit=self.visit, raft=raft, sensor=sensor))

        return data

    def getEimages(self, visit=None, raft=None, sensor=None):
        self.setVRS(visit, raft, sensor)

    def getPsfs(self, *args, **kwargs):
        return self.getDataset("psf", *args, **kwargs)
    
    def getSources(self, *args, **kwargs):
        return [pss.getSources() for pss in self.getDataset("src", *args, **kwargs)]

    def getMags(self, ss, calculateApCorr=False):
        """Return numpy arrays constructed from SourceSet ss"""
        apMags = np.empty(len(ss))
        modelMags = np.empty(len(ss))
        psfMags = np.empty(len(ss))
        flags = np.empty(len(ss), dtype='L')

        for i in range(len(ss)):
            apMags[i]  = ss[i].getApFlux()
            modelMags[i]  = ss[i].getModelFlux()
            psfMags[i] = ss[i].getPsfFlux()
            flags[i] = ss[i].getFlagForDetection()

        if True:
            calexp_md = self.butler.get('calexp_md', visit=self.visit, raft=self.raft, sensor=self.sensor)
            zp = afwImage.Calib(calexp_md).getMagnitude(1.0)
        else:
            zp = self.zp if hasattr(self, "zp") else self.ZP0

        apMags  = zp - 2.5*np.log10(apMags)
        modelMags  = zp - 2.5*np.log10(modelMags)
        psfMags = zp - 2.5*np.log10(psfMags)

        good = np.logical_and(np.isfinite(apMags), np.isfinite(psfMags))

        apMags = apMags[good]
        modelMags = modelMags[good]
        psfMags = psfMags[good]
        flags = flags[good]

        if calculateApCorr:
            delta = psfMags - apMags
            apCorr = np.median(delta[psfMags < 12])
            apMags += apCorr

        return apMags, modelMags, psfMags, flags

    def _getMags(self, ss, apMags=None, modelMags=None, psfMags=None, flags=None):
        """Return python arrays constructed from SourceSet ss, and possibly extending {ap,model,psf}Mags"""
        if not apMags:
            apMags = array.array('f')
        if not modelMags:
            modelMags = array.array('f')
        if not psfMags:
            psfMags = array.array('f')
        if not flags:
            flags = array.array('L')

        _apMags, _modelMags, _psfMags, _flags = self.getMags(ss)

        apMags.extend(_apMags)
        modelMags.extend(_modelMags)
        psfMags.extend(_psfMags)
        flags.extend(_flags)

        return apMags, modelMags, psfMags, flags

    def getMagsByVisit(self, visit=None, nSensor=0, sensor=None, raft=None):
        """Set the "self.magnitudes" for a visit"""

        d = self.lookupDataByVisit("src", visit=visit, sensor=sensor, raft=raft)

        self.setVRS(visit, raft, sensor)

        if visit:
            self.visit = visit

        if not self.visit and len(d) > 0:
            self.visit = d.keys()[0]

        if not self.visit:
            raise RuntimeError("Please specify a visit")

        apMags, modelMags, psfMags, flags = None, None, None, None

        if nSensor == 0:
            nSensor = len(d[self.visit])     # i.e. all

        for r, s in d[self.visit][0:nSensor]:
            if raft and r != raft or sensor and sensor != s:
                continue
            ss = self.getSources(visit=self.visit, raft=r, sensor=s)[0]
            apMags, modelMags, psfMags, flags = \
                    self._getMags(ss, apMags=apMags, modelMags=modelMags, psfMags=psfMags, flags=flags)

        visitStr = ("%d" % visit) if visit else "*"
        raftStr =   raft.replace(',' ,'')   if raft else "*"
        sensorStr = sensor.replace(',' ,'') if sensor else "*"
        self.name = 'imsim-%s-r%s-s%s [%s]' % (visitStr, raftStr, sensorStr,
                                                os.path.basename(os.readlink(os.path.split(self.dataRoot)[0])))

        self.apMags  = np.ndarray(len(apMags),  dtype='f', buffer=apMags)
        self.modelMags  = np.ndarray(len(modelMags),  dtype='f', buffer=modelMags)
        self.psfMags = np.ndarray(len(psfMags), dtype='f', buffer=psfMags)
        self.flags   = np.ndarray(len(flags),   dtype=long, buffer=flags)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, fixup=False, useDB=False, **kwargs):
        self.setVRS(kwargs.get("visit"), kwargs.get("raft"), kwargs.get("sensor"))

        kwargs["visit"] = self.visit
        kwargs["raft"] = self.raft
        kwargs["sensor"] = self.sensor

        try:
            kwargs = dict([(k, v) for k, v in kwargs.items() if v is not None])
            filters = set(); rafts = set(); sensors = set()
            for f, r, s in self.butler.queryMetadata("raw", "filter", ["filter", "raft", "sensor"], **kwargs):
                filters.add(f)
                rafts.add(r)
                sensors.add(s)
            filterName = list(filters)[0]
            if not kwargs.get("raft"):
                kwargs["raft"] = list(rafts)[0]
            if not kwargs.get("sensor"):
                kwargs["sensor"] = list(sensors)[0]
        except IndexError:
            msg = "No filters are available for visit %s" % self.visit
            if self.raft:
                msg += " raft %s" % self.raft
            if self.sensor:
                msg += " sensor %s" % self.sensor
            
            raise RuntimeError(msg)

        visitStr = ("%d" % kwargs["visit"]) if kwargs["visit"] else "*"
        raftStr =   kwargs["raft"].replace(',' ,'')   if kwargs["raft"] else "*"
        sensorStr = kwargs["sensor"].replace(',' ,'') if kwargs["sensor"] else "*"
        self.name = 'imsim-v%s-r%s-s%s [%s]' % (self.visit, raftStr, sensorStr,
                                                os.path.basename(os.readlink(os.path.split(self.dataRoot)[0])))
        psources = self.butler.get('icSrc', **kwargs)
        pmatches = self.butler.get('icMatch', **kwargs)
        if False:
            print 'Got sources', psources
            print 'Got matches', pmatches

        matchmeta = pmatches.getSourceMatchMetadata()
        self.matches = pmatches.getSourceMatches()
        if False:
            print 'Match metadata:', matchmeta
        self.sources = psources.getSources()

        useOutputSrc = True             # use fluxes from the "src", not "icSrc"
        if useOutputSrc:
            srcs = self.butler.get('src', **kwargs).getSources()
            import lsst.afw.detection as afwDetect
            pmMatch = afwDetect.matchXy(self.sources, srcs, 1.0, True)
            for icSrc, src, d in pmMatch:
                icSrc.setPsfFlux(src.getPsfFlux())

        calexp_md = self.butler.get('calexp_md', **kwargs)
        wcs = afwImage.makeWcs(calexp_md)
        W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")

        if useDB:
            self.zp = 2.5*math.log10(getFluxMag0DB(self.visit, self.raft, self.sensor))
        else:
            self.zp = afwImage.Calib(calexp_md).getMagnitude(1.0)

        # ref sources
        xc, yc = 0.5*W, 0.5*H
        radec = wcs.pixelToSky(xc, yc)
        ra = radec.getLongitude(DEGREES)
        dec = radec.getLatitude(DEGREES)
        radius = wcs.pixelScale()*math.hypot(xc, yc)*1.1
        if False:
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
        self.ref = X.refsources
        inds = X.inds

        referrs, stargal = None, None
        cols = solver.getTagAlongColumns(anid)
        colnames = [c.name for c in cols]

        if False:
            print 'Tag-along columns:'
            print cols
            for c in cols:
                print 'column: ', c.name, c.fitstype, c.ctype, c.units, c.arraysize

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

def queryDB(query, host=None, db=None, read_default_file="~/.my.cnf",):
    """Submit a query.  If host/db is None, use default_db_{host,db}
"""
    if host is None:
        host = default_db_host
    if db is None:
        db = default_db_table

    conn = MySQLdb.connect(host=host, db=db, read_default_file=read_default_file)

    conn.query(query)

    r = conn.use_result()
    return r.fetch_row(0)

def getFluxMag0DB(visit, raft, sensor, db=None):
    """Get fluxMag0 from the DB for the specified visit/raft/sensor"""
    
    return queryDB("""
       select
          fluxMag0
       from
          Science_Ccd_Exposure_Mapped_View
       where
          visit = %ld and raftName = '%s' and ccdName = '%s'""" % (visit, raft, sensor), db=db
                 )[0][0]

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmag(data, fig=None, magType="model", markersize=1, color="red"):
    """Plot (aper - psf) v. psf mags"""
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    fdict = maUtils.getDetectionFlags()

    bad = np.bitwise_and(data.flags, fdict["INTERP_CENTER"])
    good = np.logical_not(bad)
    if magType == "ap":
        a = data.apMags[good]
    elif magType == "model":
        a = data.modelMags[good]
    else:
        raise RuntimeError("Unknown magnitude type %s" % magType)
    p = data.psfMags[good]

    axes.plot(p, a - p, "o", markersize=markersize, color=color)
    axes.plot((0, 30), (0, 0), "b-")
    axes.set_ylim(-0.4, 0.6)
    axes.set_xlim(24, 13)
    axes.set_xlabel("psf")
    axes.set_ylabel("%s - psf" % magType)
    title = "Visit %d, all CCDs, " % (data.visit)
    title = data.name
    if data.run:
        title = "pt1prod_im%04d, " % (data.run)
    #title += "per-CCD aperture correction"
    axes.set_title(title)

    fig.show()

    return fig

def plotCalibration(data, plotBand=0.05, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand
    """
    fig = getMpFigure(fig)

    fdict = maUtils.getDetectionFlags()

    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    mstars = [m for m in data.matches if (m.second.getFlagForDetection() & fdict["STAR"])] # data
    realStars = [(m.first.getFlagForDetection() & fdict["STAR"]) != 0                      # catalogue
                 for m in data.matches if (m.second.getFlagForDetection() & fdict["STAR"])]

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars])
    measuredMagType = "psf"
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

    #axes.set_ylim(-1.1, 1.1)
    axes.set_ylim(-2.1, 2.1)
    axes.set_xlim(24, 13)
    axes.set_xlabel("Reference")
    axes.set_ylabel("Reference - %s" % measuredMagType)
    axes.set_title(data.name)

    fig.show()

    if False:
        _mstars = [m for m in mstars if math.fabs(-2.5*math.log10(m.first.getPsfFlux()) - 17.5) < 0.5]
        for m in _mstars:
            print "%.1f %.1f  %.2f %.2f" % (m.second.getXAstrom(), m.second.getYAstrom(),
                                            -2.5*math.log10(m.first.getPsfFlux()),
                                            -2.5*math.log10(m.first.getPsfFlux()) - 
                                            (data.zp - 2.5*math.log10(m.second.getPsfFlux())))
            ds9.dot("*", m.second.getXAstrom(), m.second.getYAstrom(), ctype=ds9.MAGENTA, size=10)

    return fig

def displayCalibration(data, frame=0):
    """display the calibration objects in Data object data on ds9
Stars are green; galaxies are red based on our processing
    """
    fdict = maUtils.getDetectionFlags()

    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

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

def showObject(data, objectId, filterId, frame0=0, raw=False, visits=None, maxDs9=0, showPsfs=False):
    """Show and return all the images for some objectId and filterId (a maximum of maxDs9 are displayed, starting at ds9 frame0).
If you specify visits, only that set of visits are considered for display

If showPsfs is true, include a reconstructed psf mosaic in the lower left corner 
    """

    ims = {}
    kwargs = {}
    if raw:
        imType = 'raw'
        kwargs['snap'] = 0
    else:
        imType = 'calexp'

    frame = frame0                      # frame for ds9
    nDs9 = 0                            # number of frames shown
    for visit, raft, ccd, XAstrom, YAstrom, psfMag in queryDB("""
       select
          visit, raftName, ccdName, XAstrom, YAstrom, dnToABMag(s.psfFlux, exp.fluxMag0) as psfMag
       from
          Source as s
          join Object o on s.objectId = o.objectId
          join Science_Ccd_Exposure_Mapped_Materialized as exp on s.scienceCcdExposureId = exp.scienceCcdExposureId
       where
          o.objectId = %ld and exp.filterId = %d and
          (s.flagForDetection & 0xa01) = 0
          """ % (objectId, filterId)):

        if visits and visit not in visits: # (85563714, 85597925, 85501938)
            continue

        print "visit=%d, raft='%s', ccd='%s'" % (visit, raft, ccd)

        ims[visit] = data.getDataset(imType, visit=visit, raft=raft, sensor=ccd, calibrate=True, db=True, **kwargs)[0]

        if maxDs9 and nDs9 < maxDs9:
            im = ims[visit].getMaskedImage().getImage()

            if showPsfs:
                im = im.Factory(im, True)

                psf = data.getDataset("psf", visit=visit, raft=raft, sensor=ccd, **kwargs)[0]
                nx = 15
                psfMosaic = maUtils.showPsfMosaic(ims[visit], psf, nx=nx, frame=None).makeMosaic(mode=nx)
                sim = im.Factory(im, afwImage.BBox(afwImage.PointI(0, 0), psfMosaic.getWidth(), psfMosaic.getHeight()))
                sim <<= psfMosaic
                sim *= 1000
                del sim
            
            ds9.mtv(im, wcs=ims[visit].getWcs(), title="%ld %s %s" % (visit, raft, ccd), frame=frame)
            ds9.setMaskTransparency(75)
            ds9.dot("+", XAstrom, YAstrom, frame=frame)
            ds9.zoom(4, XAstrom, YAstrom, frame=frame)
            nDs9 += 1

        frame += 1

    return ims

def subtractModels(data, objectId, filterId, frame0=0, raw=False, visits=None, maxDs9=0, showExposure=True):
    """Show and return all the images for some objectId and filterId (a maximum of maxDs9 are displayed, starting at ds9 frame0).
If you specify visits, only that set of visits are considered for display
    """

    ims = {}
    kwargs = {}
    if raw:
        imType = 'raw'
        kwargs['snap'] = 0
    else:
        imType = 'calexp'

    frame = frame0                      # frame for ds9
    nDs9 = 0                            # number of frames shown
    for visit, raft, ccd, XAstrom, YAstrom, psfMag in queryDB("""
       select
          visit, raftName, ccdName, XAstrom, YAstrom, dnToABMag(s.psfFlux, exp.fluxMag0) as psfMag
       from
          Source as s
          join Object o on s.objectId = o.objectId
          join Science_Ccd_Exposure_Mapped_Materialized as exp on s.scienceCcdExposureId = exp.scienceCcdExposureId
       where
          o.objectId = %ld and exp.filterId = %d and
          (s.flagForDetection & 0xa01) = 0
          """ % (objectId, filterId)):

        if visits and visit not in visits: # (85563714, 85597925, 85501938)
            continue

        print "visit=%d, raft='%s', sensor='%s'" % (visit, raft, ccd)

        exp = data.getDataset(imType, visit=visit, raft=raft, sensor=ccd, calibrate=True, db=True, **kwargs)[0]
        psf = data.getDataset("psf", visit=visit, raft=raft, sensor=ccd, **kwargs)[0]

        # We're missing this copy constructor: exp.Factory(exp, True)
        subtracted =  exp.Factory(exp.getMaskedImage().Factory(exp.getMaskedImage(), True), exp.getWcs().clone())

        for s in data.getSources(visit=visit, raft=raft, sensor=ccd)[0]:
            try:
                measAlg.subtractPsf(psf, subtracted.getMaskedImage(), s.getXAstrom(), s.getYAstrom())
            except pexExcept.LsstCppException, e:
                pass

        ims[visit] = (exp, subtracted)

        if maxDs9 and nDs9 < maxDs9*2:
            ds9.setMaskTransparency(75)
            if showExposure:
                ds9.mtv(exp, title="%ld %s %s" % (visit, raft, ccd), frame=frame)
                ds9.mtv(subtracted, title="%ld %s %s" % (visit, raft, ccd), frame=frame + 1)
            else:
                ds9.mtv(exp.getMaskedImage().getImage(), wcs=ims[visit][0].getWcs(), title="%ld %s %s" % (visit, raft, ccd), frame=frame)
                ds9.mtv(subtracted.getMaskedImage().getImage(), title="%ld %s %s" % (visit, raft, ccd), frame=frame + 1)

            for f in (frame, frame + 1):
                ds9.dot("+", XAstrom, YAstrom, frame=f)
                ds9.zoom(4, XAstrom, YAstrom, frame=f)

            nDs9 += 2

        frame += 2

    return ims
