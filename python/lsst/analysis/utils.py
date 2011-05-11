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
import lsst.daf.persistence as dafPersist

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.afw.coord import DEGREES
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as ds9Utils
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
try:
    import lsst.afw.extensions.rgb as afwRgb
except ImportError:
    afwRgb = None

import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils

from lsst.obs.lsstSim import LsstSimMapper

try:
    import lsst.ip.isr as ipIsr
except ImportError:
    ipIsr = None

try:
    skyToPixel_takes_degrees
except:
    skyToPixel_takes_degrees = True     # Aaargh.
try:
    log
except:
    log = pexLogging.Log.getDefaultLog()
    #log.setThreshold(pexLogging.Log.DEBUG);

try:
    import MySQLdb
except ImportError:
    MySQLdb = None

try:
    default_db_host
except NameError:
    default_db_host = "lsst10.ncsa.uiuc.edu"
    default_db_table = "buildbot_PT1_2_u_wp_trunk_2011_0507_160801"

try:
    butlerDataRoot
except NameError:
    butlerDataRoot = None
    butler = None

try:
    mpFigures
except NameError:
    mpFigures = {0 : None}              # matplotlib (actually pyplot) figures

flagsDict = maUtils.getDetectionFlags()

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
            #
            # Modify pyplot.figure().show() to make it raise the plot too
            #
            def show(self, _show=mpFigures[i].show):
                _show(self)
                try:
                    self.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word
                except Exception, e:
                    pass
            # create a bound method
            import types
            mpFigures[i].show = types.MethodType(show, mpFigures[i], mpFigures[i].__class__)

        fig = mpFigures[i]

    if not fig:
        i = sorted(mpFigures.keys())[0]
        if i > 0:
            fig = mpFigures[i[-1]]
        else:
            fig = getMpFigure(1)

    if clear:
        fig.clf()

    pyplot.figure(fig.number)           # make it active

    return fig

def makeSubplots(figure, nx=2, ny=2):
    """Return a generator of a set of subplots"""
    for window in range(nx*ny):  
        yield figure.add_subplot(nx, ny, window + 1) # 1-indexed

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
        self.zp = None

        self.name = "??"

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
        global butlerDataRoot, butler, inButler, rerun

        if run is None:
            self.run = None
        else:
            self.run = run
            butlerDataRoot = None

        if butlerDataRoot and dataRoot != butlerDataRoot:
            butler = None
            rerun = ""

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

        butler = dafPersist.ButlerFactory(mapper=LsstSimMapper(root=dataRoot, registry=registry)).create()
        butlerDataRoot = dataRoot

        inputRoot = os.path.join(os.path.split(dataRoot)[0], "input")
        inButler = dafPersist.ButlerFactory(mapper=LsstSimMapper(root=inputRoot, registry=registry)).create()

        bdr0 = os.path.split(butlerDataRoot)[0] # chop off last element (e.g. "update")
        if os.path.islink(bdr0):
            bdr0 = os.readlink(bdr0)
        rerun = os.path.basename(bdr0)

    def lookupDataBySkytile(self, dataType):
        dataSets = {}
        for st, v, f, r, s in butler.queryMetadata("raw", "skyTile",
                                                        ["skyTile", "visit", "filter", "raft", "sensor"]):
            if butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s):
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
        extra = {}
        if dataType == "raw":
            for k in ("snap", "channel",):
                extra[k] = kwargs.get(k, 1)
                kwargs[k] = extra[k]


        dataSets = {}
        for v, f, r, s in butler.queryMetadata("raw", "visit", ["visit", "filter", "raft", "sensor",],
                                               *args, **kwargs):
            if butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s, **extra):
                if not dataSets.has_key(v):
                    dataSets[v] = []
                dataSets[v].append((r, s))

        self.dataSets = dataSets
        return self.dataSets

    def getDataset(self, dataType, ids=True, calibrate=False, useDb=False, setMask=None, fixOrientation=None,
                   **dataId):
        """Get all the data of the given type (e.g. "psf"); visit may be None (meaning use default);
raft or sensor may be None (meaning get all)
"""
        if dataType == "eimage":
            if fixOrientation is None:
                fixOrientation = True
            if setMask is None:
                setMask = True
        else:
            if setMask:
                raise RuntimeError("setMask only makes sense when reading an eimage")
            if fixOrientation:
                raise RuntimeError("fixOrientation only makes sense when reading an eimage")

        if dataType in ("raw", "flat", "bias", "dark",):
            return [assembleCcd(dataType, inButler, reNorm=False, **dataId)]
        elif dataType in ("eimage",):
            if fixOrientation is None:
                fixOrientation = True
            raw_filename = inButler.get('raw_filename', channel='0,0', snap=0, **dataId)[0]
            dirName, fileName = os.path.split(re.sub(r"/raw/", "/raw/../eimage/", raw_filename))
            fileName = re.sub(r"^imsim", "eimage", fileName)
            fileName = re.sub(r"_C00_", "_", fileName)
            fileName = os.path.join(os.path.split(dirName)[0], fileName)
            eimageExposure = afwImage.ExposureF(fileName) # as read
            #
            # Fix orientation
            #
            eimage = eimageExposure.getMaskedImage().getImage()
            if fixOrientation:
                eimage = afwMath.rotateImageBy90(afwMath.flipImage(eimage, False, True), 1)

            if setMask:
                calexp = self.getDataset("calexp", **dataId)[0]
                mask = calexp.getMaskedImage().getMask()
                if not fixOrientation:
                    mask = afwMath.rotateImageBy90(afwMath.flipImage(mask, False, True), 1)
            else:
                mask = None

            return afwImage.makeExposure(afwImage.makeMaskedImage(eimage, mask), eimageExposure.getWcs())

        self.setVRS(dataId["visit"], dataId["raft"], dataId["sensor"], reset=True)

        dataSets = self.lookupDataByVisit(dataType, visit=self.visit, raft=self.raft, sensor=self.sensor)
        if not self.visit and len(dataSets) == 1:
            self.visit = dataSets.keys()[0]

        dataSets = dataSets.get(self.visit)
        if not dataSets:
            raise RuntimeError("I cannot find your data" +
                               (" for visit %d" % self.visit if self.visit else ""))

        data = []
        for raft, sensor in dataSets:
            dataElem = butler.get(dataType, visit=self.visit, raft=self.raft, sensor=self.sensor)

            if dataType == "calexp":
                psf = butler.get("psf", visit=self.visit, raft=self.raft, sensor=self.sensor)
                dataElem.setPsf(psf)
                del psf
            if dataType == "calexp" and calibrate:
                mi = dataElem.getMaskedImage();
                if useDb:
                    fluxMag0 = getFluxMag0DB(self.visit, raft, sensor)
                else:
                    fluxMag0 = dataElem.getCalib().getFluxMag0()[0]
                mi *= fluxMag0*1e-12
                del mi

            data.append(dataElem)

        return data

    def getCalexp(self, **dataId):
        """Like data.butler.get, but get the PSF too"""
        exp = butler.get("calexp", **dataId)
        psf = butler.get("psf", **dataId)
        exp.setPsf(psf)

        return exp

    def getEimages(self, visit=None, raft=None, sensor=None):
        self.setVRS(visit, raft, sensor)

    def getPsfs(self, *args, **dataId):
        return self.getDataset("psf", *args, **dataId)
    
    def getSources(self, *args, **dataId):
        return [pss.getSources() for pss in self.getDataset("src", *args, **dataId)]

    def getMags(self, ss, resetCalib=True):
        """Return numpy arrays constructed from SourceSet ss"""
        ids = np.empty(len(ss), dtype='L')
        flags = np.empty(len(ss), dtype='L')
        apMags = np.empty(len(ss))
        modelMags = np.empty(len(ss))
        psfMags = np.empty(len(ss))
        xAstrom = np.empty(len(ss))
        yAstrom = np.empty(len(ss))

        for i in range(len(ss)):
            ids[i] = ss[i].getId()
            flags[i] = ss[i].getFlagForDetection()
            apMags[i]  = ss[i].getApFlux()
            modelMags[i]  = ss[i].getModelFlux()
            psfMags[i] = ss[i].getPsfFlux()

            xAstrom[i] = ss[i].getXAstrom()
            yAstrom[i] = ss[i].getYAstrom()

        if resetCalib:
            calexp_md = butler.get('calexp_md', visit=self.visit, raft=self.raft, sensor=self.sensor)
            self.wcs = afwImage.makeWcs(calexp_md)
            self.calib = afwImage.Calib(calexp_md)
            self.zp = self.calib.getMagnitude(1.0)

        for i in range(len(ss)):
            try:
                apMags[i]  = self.calib.getMagnitude(apMags[i])
            except:
                apMags[i]  = np.nan

            try:
                psfMags[i]  = self.calib.getMagnitude(psfMags[i])
            except:
                psfMags[i]  = np.nan

            try:
                modelMags[i]  = self.calib.getMagnitude(modelMags[i])
            except:
                modelMags[i]  = np.nan

        good = np.logical_and(np.isfinite(apMags), np.isfinite(psfMags))

        ids = ids[good]
        flags = flags[good]
        apMags = apMags[good]
        modelMags = modelMags[good]
        psfMags = psfMags[good]
        xAstrom = xAstrom[good]
        yAstrom = yAstrom[good]

        return ids, flags, xAstrom, yAstrom, apMags, modelMags, psfMags

    def _getMags(self, ss, ids=None, flags=None, xAstrom=None, yAstrom=None,
                 apMags=None, modelMags=None, psfMags=None):
        """Return python arrays constructed from SourceSet ss, and possibly extending {ap,model,psf}Mags"""
        if not ids:
            ids = array.array('L')
        if not flags:
            flags = array.array('L')
        if not xAstrom:
            xAstrom = array.array('d')
        if not yAstrom:
            yAstrom = array.array('d')
        if not apMags:
            apMags = array.array('d')
        if not modelMags:
            modelMags = array.array('d')
        if not psfMags:
            psfMags = array.array('d')

        _ids, _flags, _xAstrom, _yAstrom, _apMags, _modelMags, _psfMags = self.getMags(ss)

        ids.extend(_ids)
        flags.extend(_flags)
        xAstrom.extend(_xAstrom)
        yAstrom.extend(_yAstrom)
        apMags.extend(_apMags)
        modelMags.extend(_modelMags)
        psfMags.extend(_psfMags)

        return ids, flags, xAstrom, yAstrom, apMags, modelMags, psfMags

    def getMagsByVisit(self, visit=None, nSensor=0, sensor=None, raft=None):
        """Set the "self.magnitudes" for a visit"""

        d = self.lookupDataByVisit("src", visit=visit, sensor=sensor, raft=raft)
        if not d:
            raise RuntimeError("Unable to find any data that matches your requirements")

        self.setVRS(visit, raft, sensor)

        if visit:
            self.visit = visit

        if not self.visit and len(d) > 0:
            self.visit = d.keys()[0]

        if not self.visit:
            raise RuntimeError("Please specify a visit")

        ids, flags, xAstrom, yAstrom, apMags, modelMags, psfMags = None, None, None, None, None, None, None

        if nSensor == 0:
            nSensor = len(d[self.visit])     # i.e. all

        for r, s in d[self.visit][0:nSensor]:
            if raft and r != raft or sensor and sensor != s:
                continue
            ss = self.getSources(visit=self.visit, raft=r, sensor=s)[0]
            ids, flags, xAstrom, yAstrom, apMags, modelMags, psfMags = \
                    self._getMags(ss, ids=ids, flags=flags, xAstrom=xAstrom, yAstrom=yAstrom,
                                  apMags=apMags, modelMags=modelMags, psfMags=psfMags)

        visitStr = ("%d" % visit) if visit else "*"
        raftStr =   raft.replace(',' ,'')   if raft else "*"
        sensorStr = sensor.replace(',' ,'') if sensor else "*"
        self.name = 'imsim-%s-r%s-s%s [%s]' % (visitStr, raftStr, sensorStr, rerun)

        self.ids       = np.ndarray(len(ids),       dtype='L', buffer=ids)
        self.flags     = np.ndarray(len(flags),     dtype='L', buffer=flags)
        self.xAstrom   = np.ndarray(len(xAstrom),   dtype='d', buffer=xAstrom)
        self.yAstrom   = np.ndarray(len(yAstrom),   dtype='d', buffer=yAstrom)
        self.apMags    = np.ndarray(len(apMags),    dtype='d', buffer=apMags)
        self.modelMags = np.ndarray(len(modelMags), dtype='d', buffer=modelMags)
        self.psfMags   = np.ndarray(len(psfMags),   dtype='d', buffer=psfMags)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, fixup=False, useDB=False, **kwargs):
        self.setVRS(kwargs.get("visit"), kwargs.get("raft"), kwargs.get("sensor"))

        kwargs["visit"] = self.visit
        kwargs["raft"] = self.raft
        kwargs["sensor"] = self.sensor

        try:
            kwargs = dict([(k, v) for k, v in kwargs.items() if v is not None])
            filters = set(); rafts = set(); sensors = set()
            for f, r, s in butler.queryMetadata("raw", "filter", ["filter", "raft", "sensor"], **kwargs):
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
        self.name = 'imsim-v%s-r%s-s%s [%s]' % (self.visit, raftStr, sensorStr, rerun)
        if True:
            self.getCalibObjectsImpl(filterName, fixup=fixup, useDB=False, **kwargs)
        else:
            self.matches, calib, self.ref = getCalibObjectsImpl2(filterName, **kwargs)
            self.zp = calib.getMagnitude(1.0)

    def getCalibObjectsImpl(self, filterName, fixup=False, useDB=False, **kwargs):
        psources = butler.get('icSrc', **kwargs)
        pmatches = butler.get('icMatch', **kwargs)

        matchmeta = pmatches.getSourceMatchMetadata()
        self.matches = pmatches.getSourceMatches()
        self.sources = psources.getSources()

        useOutputSrc = False             # use fluxes from the "src", not "icSrc"
        if useOutputSrc:
            srcs = butler.get('src', **kwargs).getSources()
            import lsst.afw.detection as afwDetect
            pmMatch = afwDetect.matchXy(self.sources, srcs, 1.0, True)
            for icSrc, src, d in pmMatch:
                icSrc.setPsfFlux(src.getPsfFlux())

        calexp_md = butler.get('calexp_md', **kwargs)
        self.wcs = afwImage.makeWcs(calexp_md)
        W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")

        if useDB:
            self.zp = 2.5*math.log10(getFluxMag0DB(self.visit, self.raft, self.sensor))
        else:
            self.calib = afwImage.Calib(calexp_md)
            self.zp = self.calib.getMagnitude(1.0)

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

        keepref = []
        keepi = []
        for i in xrange(len(self.ref)):
            if skyToPixel_takes_degrees:
                x, y = wcs.skyToPixel(math.degrees(self.ref[i].getRa()), math.degrees(self.ref[i].getDec()))
            else:
                x, y = wcs.skyToPixel(self.ref[i].getRa(), self.ref[i].getDec())
            if x < 0 or y < 0 or x > W or y > H:
                continue
            self.ref[i].setXAstrom(x)
            self.ref[i].setYAstrom(y)
            if stargal[i]:
                self.ref[i].setFlagForDetection(self.ref[i].getFlagForDetection() | flagsDict["STAR"])
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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getCalibObjectsImpl2(filterName, useOutputSrc=False, **kwargs):
    """
    A version of getCalibObjectsImpl that isn't a class method, for use by other code
    
    @param useOutputSrc             # use fluxes from the "src", not "icSrc"
    """
    psources = butler.get('icSrc', **kwargs)
    pmatches = butler.get('icMatch', **kwargs)
    calexp_md = butler.get('calexp_md', **kwargs)

    wcs = afwImage.makeWcs(calexp_md)
    W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
    calib = afwImage.Calib(calexp_md)

    matches = pmatches.getSourceMatches()
    sources = psources.getSources()

    anid = pmatches.getSourceMatchMetadata().getInt('ANINDID')

    del psources; del pmatches; del calexp_md # cleanup

    if useOutputSrc:
        srcs = butler.get('src', **kwargs).getSources()
        import lsst.afw.detection as afwDetect
        pmMatch = afwDetect.matchXy(sources, srcs, 1.0, True)
        for icSrc, src, d in pmMatch:
            icSrc.setPsfFlux(src.getPsfFlux())

    # ref sources
    xc, yc = 0.5*W, 0.5*H
    radec = wcs.pixelToSky(xc, yc)
    ra = radec.getLongitude(DEGREES)
    dec = radec.getLatitude(DEGREES)
    radius = wcs.pixelScale()*math.hypot(xc, yc)*1.1

    pol = pexPolicy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAstrom.createSolver(pol, log)
    idName = 'id'

    X = solver.getCatalogue(ra, dec, radius, filterName, idName, anid)
    refsources = X.refsources
    inds = X.inds

    referrs, stargal = None, None
    cols = solver.getTagAlongColumns(anid)
    colnames = [c.name for c in cols]

    col = 'starnotgal'
    if col in colnames:
        stargal1 = solver.getTagAlongBool(anid, col, inds)
        stargal = []
        for i in range(len(stargal1)):
            stargal.append(stargal1[i])

    keepref = []
    keepi = []
    for i in xrange(len(refsources)):
        x, y = wcs.skyToPixel(refsources[i].getRa(), refsources[i].getDec())
        if x < 0 or y < 0 or x > W or y > H:
            continue
        refsources[i].setXAstrom(x)
        refsources[i].setYAstrom(y)
        if stargal[i]:
            refsources[i].setFlagForDetection(refsources[i].getFlagForDetection() | flagsDict["STAR"])
        keepref.append(refsources[i])
        keepi.append(i)

    refsources = keepref

    if referrs is not None:
        referrs = [referrs[i] for i in keepi]
    if stargal is not None:
        stargal = [stargal[i] for i in keepi]

    measAstrom.joinMatchList(matches, refsources, first=True, log=log)
    args = {}
    measAstrom.joinMatchList(matches, sources, first=False, log=log, **args)

    return matches, calib, refsources

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

def plotDmag(data, fig=None, magType="model", maglim=20, markersize=1, color="red"):
    """Plot (aper - psf) v. psf mags"""
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = np.bitwise_and(data.flags, (flagsDict["INTERP_CENTER"] | flagsDict["EDGE"]))
    good = np.logical_not(bad)

    if False:
        amp30 = np.logical_and(data.xAstrom > 1526, data.xAstrom < 2036)
        amp30 = np.logical_and(amp30, data.yAstrom < 2000)
        good = np.logical_and(good, amp30)

    if magType == "ap":
        a = data.apMags[good]
    elif magType == "model":
        a = data.modelMags[good]
    else:
        raise RuntimeError("Unknown magnitude type %s" % magType)
    p = data.psfMags[good]
    delta = a - p

    sgVal = 0.05
    stellar = np.abs(delta) < sgVal
    locus = np.logical_and(a < maglim, stellar)

    try:
        stats = afwMath.makeStatistics(delta[locus], afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
    except:
        mean, stdev = np.nan, np.nan

    if False:
        axes.plot(p[locus], delta[locus], "o", markersize=2*markersize, color="blue")
    else:
        axes.plot((0, maglim, maglim, 0, 0), sgVal*np.array([-1, -1, 1, 1, -1]), "g:")

    axes.plot(a, delta, "o", markersize=markersize, color=color)
    axes.plot((0, 30), (0, 0), "b-")

    if False:
        if False:                           # clump
            magmin, magmax = 20.03, 20.13
            if magType == "model":
                dmin, dmax = -0.72, -0.67
            else:
                dmin, dmax = -0.576, -0.512
        else:                               # "saturated"
            magmin, magmax = 14, 14.9
            dmin, dmax = -0.2, 0.2
            dmin, dmax = -0.2, -0.1

        axes.plot([magmax, magmin, magmin, magmax, magmax], [dmin, dmin, dmax, dmax, dmin], "r-")
        objs = np.logical_and(np.logical_and(a > magmin, a < magmax),
                              np.logical_and(delta > dmin, delta < dmax))
        ids = data.ids[good]
        for i in sorted(ids[objs]):
            print i, splitId(i)

    axes.set_ylim(-0.8, 0.2)
    if False:
        axes.set_ylim(-5.4, 0.6)
    axes.set_xlim(13, 24)
    axes.set_xlabel(magType)
    axes.set_ylabel("%s - psf" % magType)
    title = "Visit %d, all CCDs, " % (data.visit)
    title = data.name
    if data.run:
        title = "pt1prod_im%04d, " % (data.run)
    #title += "per-CCD aperture correction"
    axes.set_title(title)

    fig.text(0.20, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID
    #
    global eventHandler
    eventHandler = EventHandler(axes, a, delta, data.ids[good], data)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmagHistograms(data, fig=None, magType="model", color="red"):
    """Plot (aper - psf) v. psf mags"""
    fig = getMpFigure(fig)

    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = np.bitwise_and(data.flags, (flagsDict["INTERP_CENTER"] | flagsDict["EDGE"]))
    good = np.logical_not(bad)

    if False:
        amp30 = np.logical_and(data.xAstrom > 1526, data.xAstrom < 2036)
        amp30 = np.logical_and(amp30, data.yAstrom < 2000)
        good = np.logical_and(good, amp30)

    if magType == "ap":
        a = data.apMags[good]
    elif magType == "model":
        a = data.modelMags[good]
    else:
        raise RuntimeError("Unknown magnitude type %s" % magType)

    p = data.psfMags[good]
    delta = a - p

    sgVal = 0.10
    stellar = np.abs(delta) < sgVal

    nbin = 15
    minmag, maxmag, dmag = 15, 20, 0.5
    binEdge = np.arange(minmag, maxmag - dmag/2, dmag)

    subplots = makeSubplots(fig, nx=1, ny=len(binEdge))

    for edge in binEdge:
        inBin = np.logical_and(np.logical_and(a >= edge, a < edge + dmag), stellar)

        pp = p[inBin]        
        dd = delta[inBin]

        if len(dd) == 0:
            continue

        stats = afwMath.makeStatistics(dd, afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
        
        axes = subplots.next()
        n, bins, patches = axes.hist(dd, nbin)
        axes.plot((mean, mean), (-100, max(n) + 100), "g:")
        axes.set_xlim(-1.05*sgVal, 1.05*sgVal)
        axes.set_ylim(0, 1.05*max(n) + 1)

        axes.xaxis.set_major_formatter(pyplot.NullFormatter())
        axes.yaxis.set_major_formatter(pyplot.NullFormatter())

        axes.set_xlabel("%.2f" % (edge + 0.5*dmag))
        axes.set_title("%.4f" % (stdev), fontsize="smaller") 

    fig.suptitle("%s - psf: %s" % (magType, data.name))

    fig.show()

    return fig

def plotCalibration(data, plotBand=0.05, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand
    """
    fig = getMpFigure(fig)

    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    mstars = [m for m in data.matches if (m.second.getFlagForDetection() & flagsDict["STAR"])] # data
    realStars = [(m.first.getFlagForDetection() & flagsDict["STAR"]) != 0                      # catalogue
                 for m in data.matches if (m.second.getFlagForDetection() & flagsDict["STAR"])]

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars])
    measuredMagType = "psf"
    instmag = np.array([data.zp - 2.5*math.log10(s.second.getPsfFlux()) for s in mstars])
    realStars = np.array([m for m in realStars])

    delta = refmag - instmag

    stats = afwMath.makeStatistics(delta[refmag < 20], afwMath.STDEVCLIP | afwMath.MEANCLIP)
    mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)

    axes.plot(refmag, delta, "mo")

    #axes.plot(refmag[realStars], delta[realStars], "gx")
        
    axes.plot((-100, 100), (0, 0), "g-")
    if plotBand:
        for x in (-plotBand, plotBand):
            axes.plot((-100, 100), x*np.ones(2), "g--")

    axes.set_ylim(-0.6, 0.6)
    #axes.set_ylim(-2.1, 2.1)
    axes.set_xlim(13, 24)
    axes.set_xlabel("Reference")
    axes.set_ylabel("Reference - %s" % measuredMagType)
    axes.set_title(data.name)

    fig.text(0.75, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")

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
    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    for m in data.matches:
        ref, src = m.first, m.second
        x1, y1 = src.getXAstrom(), src.getYAstrom()
        if ref.getFlagForDetection() & flagsDict["STAR"]:
            ptype = "+"
            if src.getFlagForDetection() & flagsDict["STAR"]:
                ctype = ds9.GREEN
            else:
                ctype = ds9.YELLOW
        else:
            ptype = "o"
            if not (src.getFlagForDetection() & flagsDict["STAR"]):
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

def findSource(src, x, y, radius=2):
    ss = afwDetect.SourceSet()

    s = afwDetect.Source();
    s.setXAstrom(x); s.setYAstrom(y)
    ss.push_back(s)

    matches = [m.first for m in afwDetect.matchXy(src, ss, radius)]
    if len(matches) == 0:
        raise RuntimeError("Unable to find a source within %g pixels of (%g,%g)" % (radius, x, y))

    if len(matches) == 1:
        return matches[0]
    else:
        return matches

def getRefCatalog(data, **dataId):
    """Get the reference catalogue for dataId (a single CCD) as a SourceSet"""
    ss = afwDetect.SourceSet()

    calexp_md = butler.get('calexp_md', **dataId)
    data.wcs = afwImage.makeWcs(calexp_md)
    data.calib = afwImage.Calib(calexp_md)
    data.zp = data.calib.getMagnitude(1.0)

    for refObjectId, ra, dec, mag, isStar, isVar, visit, raftName, ccdName \
            in queryDB("""
    select
       refObjectId, sro.ra, sro.decl,
       CASE WHEN exp.filterId = 0 THEN uMag
            WHEN exp.filterId = 1 THEN gMag
            WHEN exp.filterId = 2 THEN rMag
            WHEN exp.filterId = 3 THEN iMag
            WHEN exp.filterId = 4 THEN zMag
            WHEN exp.filterId = 5 THEN yMag
       END as Mag,
       sro.isStar, sro.isVar,
       exp.visit, exp.raftName, exp.ccdName
    from
       SimRefObject as sro,
       Science_Ccd_Exposure as exp
    where
       visit = %(visit)d and raftName = '%(raft)s' and ccdName = '%(sensor)s' and
       qserv_ptInSphPoly(sro.ra, sro.decl, concat_ws(' ',
                         llcRa, llcDecl, lrcRa, lrcDecl, urcRa, urcDecl, ulcRa, ulcDecl))
          """ % dataId):
        s = afwDetect.Source();

        s.setId(refObjectId)
        ra, dec = math.radians(ra), math.radians(dec)
        s.setRa(ra); s.setDec(dec)

        if skyToPixel_takes_degrees:
            ra, dec = math.degrees(ra), math.degrees(dec)

        x, y = data.wcs.skyToPixel(ra, dec)
        s.setXAstrom(x), s.setYAstrom(y)

        flux = 10**(-0.4*(mag - data.zp))
        s.setApFlux(flux); s.setPsfFlux(flux); s.setModelFlux(flux)

        flag = flagsDict["BINNED1"]
        if isStar:
            flag |= flagsDict["STAR"]
        if isVar:
            flag |= flagsDict["PEAKCENTER"] # XXX
        s.setFlagForDetection(flag)

        ss.push_back(s)

    return ss

def showRefCatalog(data, wcs=None, frame=0, **dataId):
    ss = data.getSources(**dataId)[0]
    refCat = getRefCatalog(data, **dataId)

    if frame is not None:
        if not wcs:
            wcs = data.wcs

        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="+")
        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="x", mask=["PEAKCENTER",])
        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="o", mask=["STAR",])

        showSourceSet(ss, wcs=wcs, frame=frame, ctype=ds9.RED, symb="x")

    return ss, refCat

def showSourceSet(sourceSet, exp=None, wcs=None, raDec=None, magmin=-100, magmax=None, magType="psf",
                  frame=0, ctype=ds9.GREEN, symb="+", size=2, mask=None):
    """Show a SourceSet on ds9.

    If mask, it's a set of bitplane names (e.g. INTERP) which must be set to display the source"""
    
    ds9.cmdBuffer.pushSize()

    if mask:
        bmask = 0L
        for name in mask:
            bmask |= flagsDict[name]
        mask = bmask

    if raDec is None and wcs is not None:
        raDec = True

    if raDec and not wcs:
        if exp:
            wcs = exp.getWcs()
        if not wcs:
            raise RuntimeError("You'll have to provide a wcs if you want me to plot the (ra, dec) positions")

    for s in sourceSet:
        if mask:
            if not (s.getFlagForDetection() & mask):
                continue

        if magmax is not None:
            if magType == "ap":
                flux = s.getApFlux()
            elif magType == "model":
                flux = s.getModelFlux()
            elif magType == "psf":
                flux = s.getPsfFlux()
            else:
                raise RuntimeError("Uknown magnitude type %s" % magType)
                
            mag = exp.getCalib().getMagnitude(flux)

            if mag < magmin or mag > magmax:
                continue

        if raDec:
            ra, dec = s.getRa(), s.getDec()

            if skyToPixel_takes_degrees:
                ra, dec = math.degrees(ra), math.degrees(dec)

            x, y = wcs.skyToPixel(ra, dec)
        else:
            x, y = s.getXAstrom(), s.getYAstrom()

        ds9.dot(symb, x, y, frame=frame, ctype=ctype, size=size)

    ds9.cmdBuffer.popSize()

def plotCompleteness(data, dataId, refCat=None, matchRadius=2, **kwargs):
    ss = data.getSources(**dataId)[0]
    if not refCat:
        refCat = getRefCatalog(data, **dataId)

    matched, detectedOnly, refOnly = getMatches(ss, refCat, matchRadius)

    if not kwargs.has_key("title"):
        kwargs["title"] = str(dataId)

    return (plotCounts(data, matched["src"], refOnly, detectedOnly, **kwargs), ss, refCat)

def plotCounts(data, matched, refOnly, detectedOnly, includeDetOnly=True, magType="model",
               magmin=14, magmax=26, dmag=0.5, log=True, title=None, fig=None):
    
    if not title:
        title = data.name

    bins = np.arange(magmin, magmax + 0.5*dmag, dmag)

    arrays = []                         # our data arrays
    for s in [matched, refOnly, detectedOnly]:
        ids, flags, xAstrom, yAstrom, apMags, modelMags, psfMags = data.getMags(s, resetCalib=False)

        if magType == "ap":
            mags = apMags
        elif magType == "model":
            mags = modelMags
        elif magType == "psf":
            mags = psfMags
        else:
            raise RuntimeError("Uknown magnitude type %s" % magType)

        bad = np.bitwise_and(flags, (flagsDict["INTERP_CENTER"] | flagsDict["EDGE"]))
        mags = mags[np.logical_not(bad)]

        arrays.append(np.histogram(mags, bins)[0])

    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    left = bins[0:-1] - 0.5*dmag        # left side of bins
    axes.bar(left, arrays[0], width=dmag, log=log, label="matched", color="green")
    axes.bar(left, arrays[1], bottom=arrays[0], width=dmag, log=log, label="refOnly", color="blue")
    if includeDetOnly:
        axes.bar(left, arrays[2], width=dmag, log=log, label="detOnly", alpha=0.7, color="red")

    axes.legend(loc=2)
    axes.set_xlim(magmin, magmax)
    axes.set_ylim(0.5 if log else -1, 1.05*max(arrays[0] + arrays[1]))
    axes.set_xlabel(magType)

    axes.set_title(title)

    fig.show()

    return fig

def getMatches(ss, refCat, radius=2):
    """Get the matches between a SourceSet of measured Sources and another, refCat, generated from the DB

    Return three lists:
       matches          As dictionary of two SourceSets (keys ref and src) of Sources from refCat and ss
       detectedOnly     A SourceSet of Sources from ss that don't match
       refOnly          A SourceSet of Sources from refCat that don't match
    """
    matches = afwDetect.matchRaDec(refCat, ss, 2000) # all objects should match

    matchedRefs = {}                    # matched pairs;  reference object is unique
    matchedSrcs = {}                    # list of matched objects from ss
    
    detectedOnly = afwDetect.SourceSet()
    refOnly = afwDetect.SourceSet()

    for m in matches:
        d = m.distance
        ref, src = m.first, m.second

        sid = src.getId()
        if matchedSrcs.has_key(sid):
            oref, osrc, od = matchedSrcs[sid] # other match
            if od < d:                  # the other match was better
                refOnly.push_back(ref)
                continue

            del matchedRefs[oref.getId()]
            refOnly.push_back(oref)

        matchedRefs[ref.getId()] = (ref, m)
        matchedSrcs[sid] = (ref, src, d)

    matched = dict(ref = afwDetect.SourceSet(), src = afwDetect.SourceSet())

    for ref, m in matchedRefs.values():
        d = m.distance
        ref, src = m.first, m.second

        if m.distance < radius:
            # Set the quality flags from the src
            ref.setFlagForDetection( ref.getFlagForDetection() |
                                    (src.getFlagForDetection() & ~flagsDict["STAR"]))

            matched["ref"].push_back(ref)
            matched["src"].push_back(src)
        else:
            refOnly.push_back(ref)

    sids = set()
    for ref, src, distance in matchedSrcs.values():
        sids.add(src.getId())
        if distance >= radius:
            detectedOnly.push_back(src)
    #
    # It's possible for a src to match to no reference object (think of one in the interior of a triangle)
    #
    for src in ss:
        if src.getId() not in sids:
            detectedOnly.push_back(src)
            
    return matched, detectedOnly, refOnly

def dataIdToTitle(**dataId):
    "%(visit)ld %(raft)s %(sensor)s" % dataId

def subtractModels(data, frame=0, showExposure=True, **dataId):
    """Show and return all the images for some objectId and filterId (a maximum of maxDs9 are displayed, starting at ds9 frame0).
If you specify visits, only that set of visits are considered for display
    """

    exp = data.getDataset("calexp", **dataId)[0]
    psf = exp.getPsf()

    # We're missing this copy constructor: exp.Factory(exp, True)
    subtracted =  exp.Factory(exp.getMaskedImage().Factory(exp.getMaskedImage(), True), exp.getWcs().clone())

    for s in data.getSources(**dataId)[0]:
        try:
            measAlg.subtractPsf(psf, subtracted.getMaskedImage(), s.getXAstrom(), s.getYAstrom())
        except pexExcept.LsstCppException, e:
            pass

    ds9.setMaskTransparency(75)

    title = dataIdToTitle(**dataId)
    if showExposure:
        ds9.mtv(exp, title=title, frame=frame)
        ds9.mtv(subtracted, title=title, frame=frame + 1)
    else:
        ds9.mtv(exp.getMaskedImage().getImage(), wcs=ims[visit][0].getWcs(), title=title, frame=frame)
        ds9.mtv(subtracted.getMaskedImage().getImage(), title=title, frame=frame + 1)

def psf(exposure, ngrid=40, fig=None):
    psf = exposure.getPsf()
    k = psf.getKernel()
    funcs = k.getSpatialFunctionList()
    nfunc = k.getNKernelParameters()
    assert nfunc == len(funcs)

    width, height = exposure.getDimensions()
    deltax = width/float(ngrid - 1)
    deltay = height/float(ngrid - 1)
    X, Y = np.meshgrid(deltax*np.arange(ngrid + 1), deltay*np.arange(ngrid + 1))

    Z = []
    for j in range(nfunc):
        Z.append(100*np.ones(X.size))
        
    j = 0
    for _x, _y in zip(X, Y):
        for x, y in zip(_x, _y):
            for k in range(nfunc):
                Z[k][j] = funcs[k](x, y)
            j += 1

    for j in range(nfunc):
        Z[j] = Z[j].reshape(X.shape)
        if j > 0:
            Z[j] /= Z[0]

    Z[0] /= np.max(Z[0])
    #
    # Time for matplotlib
    #
    fig = getMpFigure(fig)

    n = nfunc - 1
    nx = int(math.sqrt(n))
    ny = n/nx
    if ny == 0:
        ny = 1
    while nx*ny < n:
        nx += 1
        
    subplots = makeSubplots(fig, nx=nx, ny=ny)

    ccd = cameraGeom.cast_Ccd(exposure.getDetector())
    visit = exposure.getMetadata().get("OBSID")

    fig.suptitle("PSF eigen components, %s %s [%s]" % (visit, ccd.getId().getName(), rerun), fontsize=14) 

    for i in range(1, nfunc):
        axes = subplots.next()
        CS = axes.contour(X, Y, Z[i], 20, colors="red")
        axes.clabel(CS, inline=1, fontsize=10)
        axes.set_xlim(-5, width - 1 + 5)
        axes.set_ylim(-5, height - 1 + 5)

        axes.text(0.1, 0.1, "Component %d" % i, transform=axes.transAxes)

    fig.show()

    return fig

def grayScale(image, rgbFile, min=0, max=50, Q=8, bin=1):
    """Write a grayscale file (currently tiff) to rgbFile

The default stretch is an asinh stretch; Q == 0 corresponds to a linear.  If Q == 0 then min and max are
the range of the stretch; as Q increases they still define the slope of the linear portion of the stretch

If bin is specified, it should be an integer > 1;  the output file will be binned by that factor.
    """
    
    if not afwRgb:
        raise RuntimeError("I was unable to import lsst.afw.extensions.rgb")

    #
    # Handle Exposures and MaskedImages
    #
    try:
        image = image.getMaskedImage()    # maybe it's an Exposure
    except AttributeError:
        pass

    try:
        image = image.getImage()          # maybe it's (now) a MaskedImage
    except AttributeError:
        pass

    if bin > 1:
        image = afwMath.binImage(image, bin)

    afwRgb.RgbImageF(image, image, image, afwRgb.asinhMappingF(min, max - min, Q)).writeTiff(rgbFile)

def assembleCcd(dataType, butler, snap=0, reNorm=False, **dataId):
    """Return an Exposure of a raw data frame (or related image stored per-amp such as a flat or bias)

    @param dataType	Desired type (e.g. "raw")
    @param butler       An input butler
    @param snap	        The desired snap, if appropriate
    @param reNorm       Renormalise each amp image to unit gain
    @param dataId       Dictionary of visit, ccd, etc.
    """
    if not ipIsr:
        raise RuntimeError("Failed to import lsst.ip.isr -- is it setup?")

    cg = butler.get("camera", **dataId)
    md = butler.get("raw_md",  snap=snap, channel="0,0", **dataId)

    mat = re.search(r"R(\d)(\d)_S(\d)(\d)", md.get("CCDID"))
    r0, r1, s0, s1 = [int(d) for d in mat.groups()]
    id = cameraGeom.Id("R:%d,%d S:%d,%d" % (r0, r1, s0, s1))

    ccd = cameraGeomUtils.findCcd(cg, id)
    
    ampList = []
    for n,a in enumerate(ccd):
        serial = a.getId().getSerial()
        channel = "%d,%d" % (serial//8, serial%8)

        exp = butler.get(dataType, snap=snap, channel=channel, **dataId)
        try:
            exp = exp.convertF()
        except AttributeError:          # already F
            pass

        exp.setDetector(a)
        exp.setFilter(afwImage.Filter(md))
        try:
            exp.setCalib(afwImage.Calib(md))
        except AttributeError:          # pre afw 4.3.2.1
            calib = afwImage.Calib(md)
            exp.getCalib().setExptime(calib.getExptime())
            
        ampList.append(exp)

    return ipIsr.assembleCcd(ampList, ccd, reNorm=reNorm)

def splitId(oid):
    """Split an ObjectId into visit, raft, sensor, and objId"""
    objId = int((oid & 0xffff) - 1)     # Should be the same value as was set by apps code
    oid >>= 16
    raftSensorId = oid & 0x1ff
    oid >>= 9
    visit = int(oid)

    raftId, sensorId = int(raftSensorId//10), int(raftSensorId%10)
    raft = "%d,%d" % (raftId//5, raftId%5)
    sensor = "%d,%d" % (sensorId//3, sensorId%3)

    return visit, raft, sensor, objId

class EventHandler(object):
    """A class to handle key strokes with matplotlib displays"""
    def __init__(self, axes, xs, ys, ids, data, frame=0):
        self.axes = axes
        self.xs = xs
        self.ys = ys
        self.ids = ids
        self.data = data
        self.frame = frame

        self.cid = self.axes.figure.canvas.mpl_connect('key_press_event', self)

    def __call__(self, ev):
        if ev.inaxes != self.axes:
            return
        
        if ev.key in ("i", "I", "p", "P"):
            dist = np.hypot(self.xs - ev.xdata, self.ys - ev.ydata)
            dmin = min(dist)

            dist = np.hypot(self.xs - ev.xdata, self.ys - ev.ydata)
            objId = self.ids[dist == min(dist)][0]

            print "\r>>>",
            print ("%.3f %.3f" % (ev.xdata, ev.ydata)), splitId(objId), \
                  maUtils.explainDetectionFlags(self.data.flags[objId == self.data.ids][0]),
            sys.stdout.flush()

            if ev.key == "I":
                print ""
            elif ev.key in ("p", "P"):
                x = self.data.xAstrom[objId == self.data.ids][0]
                y = self.data.yAstrom[objId == self.data.ids][0]
                ds9.pan(x, y, frame=self.frame)

        else:
            pass
