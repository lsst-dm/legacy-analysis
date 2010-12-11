"""
Some utilities for looking at the outputs of processing runs

Currently specialised to handle ImSim outputs --- but this restriction should be lifted!
"""
import array, math, os, re
import numpy
try:
    import matplotlib.pyplot as pyplot
except ImportError:
    pyplot = None
import lsst.afw.image as afwImage
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
        if run is not None:
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

        return dataSets

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

        return dataSets

    def getDataset(self, dataType, visit=None, raft=None, sensor=None, ids=True):
        """Get all the data of the given type (e.g. "psf"); visit may be None (meaning use default);
raft or sensor may be None (meaning get all)

N.b. This routine resets the self.ids listless ids is False; it is assumed that you retrieve all data items from the same set
of sensors, but this is not checked.
"""
        if visit:
            self.visit = visit

        dataSets = self.lookupDataByVisit(dataType, raft=raft, sensor=sensor)
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
        apMags = numpy.empty(len(ss))
        psfMags = numpy.empty(len(ss))
        flags = numpy.empty(len(ss))

        for i in range(len(ss)):
            apMags[i]  = ss[i].getApFlux()
            psfMags[i] = ss[i].getPsfFlux()
            flags[i] = ss[i].getFlagForDetection()

        apMags  = self.ZP0 - 2.5*numpy.log10(apMags)
        psfMags = self.ZP0 - 2.5*numpy.log10(psfMags)

        good = numpy.logical_and(numpy.isfinite(apMags), numpy.isfinite(psfMags))

        apMags = apMags[good]
        psfMags = psfMags[good]
        flags = flags[good]

        if calculateApCorr:
            delta = psfMags - apMags
            apCorr = numpy.median(delta[psfMags < 12])
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

    def getMagsByVisit(self, visit=None, nSensor=0):
        """Set the "self.magnitudes" for a visit"""
        d = self.lookupDataByVisit("src")

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
            ss = self.getSources(raft=r, sensor=s)[0]
            apMags, psfMags, flags = self._getMags(ss, apMags=apMags, psfMags=psfMags, flags=flags)

        self.apMags  = numpy.ndarray(len(apMags),  dtype='f', buffer=apMags)
        self.psfMags = numpy.ndarray(len(psfMags), dtype='f', buffer=psfMags)
        self.flags   = numpy.ndarray(len(flags),   dtype=long, buffer=flags)

    def plotDmag(self, fig=None, markersize=0.1, color="red"):
        """Plot (aper - psf) v. psf mags"""
        if fig:
            self.fig = fig

        if self.fig:
            self.fig.clf()
        else:
            if not pyplot:
                raise RuntimeError("I am unable to plot as I failed to import matplotlib")
            self.fig = pyplot.figure()

        axes = self.fig.add_axes((0.1, 0.1, 0.85, 0.80));

        fdict = maUtils.getDetectionFlags()

        bad = numpy.bitwise_and(self.flags, fdict["INTERP_CENTER"])
        good = numpy.logical_not(bad)
        a = self.apMags[good]
        p = self.psfMags[good]

        axes.plot(p, a - p, "o", markersize=markersize, color=color)
        axes.plot((0, 20), (0, 0), "b-")
        axes.set_ylim(-0.4, 0.6)
        axes.set_xlim(8, 17)
        axes.set_xlabel("psf")
        axes.set_ylabel("aper - psf")
        axes.set_title("Visit %d, all CCDs, pt1prod_im%04d, per-CCD aperture correction" %
                       (self.visit, self.run))

        self.fig.show()

    def showPsfs(self, psfs, frame=None):
        mos = ds9Utils.Mosaic()
        for i in range(len(psfs)):
            psf = psfs[i]
            mos.append(psf.computeImage(), "%s:%s" % (self.ids[i].raft, self.ids[i].sensor))

        mos.makeMosaic(frame=frame)
