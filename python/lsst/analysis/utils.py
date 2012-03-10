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
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as ds9Utils
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.table as afwTable
try:
    import lsst.afw.extensions.rgb as afwRgb
except ImportError:
    afwRgb = None

import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils

try:
    from lsst.obs.lsstSim import LsstSimMapper
except ImportError:
    class LsstSimMapper(object): pass

try:
    from lsst.obs.suprimecam.suprimecamMapper import SuprimecamMapper
except ImportError:
    class SuprimecamMapper(object): pass

try:
    flagsDict
except:
    flagsDict, i = {}, 0
    for n in ["STAR", "flags.pixel.interpolated.center", "flags.pixel.interpolated.any",
              "flags.pixel.saturated.center", "flags.pixel.saturated.any",
              "flags.pixel.edge", "BINNED1", "flags.badcentroid"]:
        flagsDict[n] = 1 << i; i += 1

try:
    import lsst.ip.isr as ipIsr
except ImportError:
    ipIsr = None

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
    default_db_table = None # "buildbot_PT1_2_u_wp_trunk_2011_0507_160801"

try:
    butlerDataRoot
except NameError:
    butlerDataRoot = None
    butlerRerun = None
    butler = None

try:
    mpFigures
except NameError:
    mpFigures = {0 : None}              # matplotlib (actually pyplot) figures
    eventHandlers = {}                  # event handlers for matplotlib figures


def getMapper(registryRoot, defaultMapper=None, mapperFile="mapper.py"):
    """
    Get proper mapper given a directory registryRoot

    registryRoot should contain a file, mapperFile, which
    defines a variable called Mapper

    E.g. mapper.py might contain:
    from lsst.obs.lsstSim import LsstSimMapper as Mapper
    """
    Mapper = None
    _locals = {}
    mapper_py = os.path.join(registryRoot, mapperFile)
    try:
        execfile(mapper_py, {}, _locals)
        Mapper = _locals.get("Mapper")
    except IOError:
        pass                        # the file doesn't exist; this isn't an error
    except ImportError, e:
        raise

    if not Mapper:
        if defaultMapper:
            return defaultMapper

        raise RuntimeError("I'm unable to find a viable mapper")
    else:
        return Mapper

def getNameOfSet(vals):
    """Convert a list of numbers into a string, merging consecutive values"""
    if not vals:
        return ""

    valName = []
    val0 = vals[0]; val1 = val0
    for val in vals[1:]:
        if val == val1 + 1:
            val1 = val
        else:
            valName.append("%d-%d" % (val0, val1) if val1 != val0 else str(val0))
            val0 = val; val1 = val0

    valName.append("%d-%d" % (val0, val1) if val1 != val0 else str(val0))

    return ", ".join(valName)

def makeMapperInfo(mapper):
    """Return an object with extra per-mapper information (e.g. which fields fully specify an exposure)"""

    class LsstSimMapperInfo(object):
        def __init__(self, Mapper):
            LsstSimMapperInfo.Mapper = Mapper

        @staticmethod
        def getFields(dataType):
            fields = ["visit", "filter", "raft", "sensor",]
            if dataType == "raw":
                fields += ["snap", "channel",]

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw", "flat", "bias", "dark",)

        @staticmethod
        def dataIdToTitle(dataIds):
            try:
                dataId = dataIds[0]

                if len(dataIds) > 1:
                    print >> sys.stderr, "Fit dataIdToTitle for more than one dataId"
            except KeyError:
                dataId = dataIds
                
            filterName = afwImage.Filter(butler.get("calexp_md", **dataId)).getName()
            return "%ld %s %s [%s]" % (dataId["visit"], dataId["raft"], dataId["sensor"], filterName)

        @staticmethod
        def exposureToStr(exposure):
            ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getName()
            visit = exposure.getMetadata().get("OBSID")

            return "%s %s" % (visit, ccdId)

        assembleCcd = staticmethod(assembleCcdLsst)

        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            inputRoot = os.path.join(os.path.split(dataRoot)[0], "input")
            return dafPersist.ButlerFactory(mapper=LsstSimMapperInfo.Mapper(root=inputRoot,
                                                                            registry=registry)).create()

        @staticmethod
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

        @staticmethod
        def getDataIdMask(dataId):
            """Return a mask to | with a Source.id to uniquely identify the object"""
            return 0x0                  # no bits to add

    class SubaruMapperInfo(object):
        def __init__(self, Mapper):
            SubaruMapperInfo.Mapper = Mapper

        @staticmethod
        def getFields(dataType):
            if dataType in ("flat",):
                fields = ["visit", "ccd"]
            else:
                fields = ["visit", "filter", "ccd"]

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw",)

        @staticmethod
        def dataIdToTitle(dataIds):
            filters = set()
            ccds = set()
            visits = set()
            for dataId in dataIds:
                if dataId["ccd"] == None:
                    did = dataId.copy(); did["ccd"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get("calexp_md", **did)).getName())
                    except:
                        filters.add("?")
                    ccds.add("all")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get("calexp_md", **dataId)).getName())
                    except:
                        filters.add("?")
                    ccds.add(dataId["ccd"])
                try:
                    visits.add(dataId["visit"])
                except TypeError:
                    for v in dataId["visit"]:
                        visits.add(v)

            ccds = sorted(list(ccds))
            filters = sorted(list(filters))
            visits = sorted(list(visits))

            if len(visits) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple visits and filters: %s %s" % \
                      (visits, filters)
                visits = visits[0:1]

            title = "%s CCD%s [%s]" % \
                    (", ".join([str(x) for x in visits]), getNameOfSet(ccds), ", ".join(filters))
            title = "%s CCD%s [%s]" % (getNameOfSet(visits), getNameOfSet(ccds), ", ".join(filters))
            if rerunName:
                title += " %s" % rerunName

            return title

        @staticmethod
        def exposureToStr(exposure):
            try:
                ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getSerial()
                visit = re.sub(r"^SUPA", "", exposure.getMetadata().get("FRAMEID"))
            except AttributeError:
                return "??"

            return "%s %s" % (visit, ccdId)

        assembleCcd = staticmethod(assembleCcdSubaru)
    
        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            return butler

        @staticmethod
        def splitId(oid):
            """Split an ObjectId into visit, ccd, and objId"""
            oid = long(oid)
            objId = int(oid & 0xffff)     # Should be the same value as was set by apps code
            oid >>= 16
            ccd = int(oid & 0xff)
            oid >>= 8
            visit = int(oid)

            return visit, ccd, objId

        @staticmethod
        def getDataIdMask(dataId):
            """Return a mask to | with a Source.id to uniquely identify the object"""
            return ((dataId["visit"] << 8) | dataId["ccd"]) << 16        

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if isinstance(mapper, LsstSimMapper):
        return LsstSimMapperInfo(LsstSimMapper)
    elif isinstance(butler.mapper, SuprimecamMapper):
        return SubaruMapperInfo(SuprimecamMapper)
    else:
        raise RuntimeError("Impossible mapper")

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

        def lift(fig):
            fig.canvas._tkcanvas._root().lift() # == Tk's raise, but raise is a python reserved word

        if mpFigures.has_key(i):
            try:
                lift(mpFigures[i])
            except Exception, e:
                del mpFigures[i]

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
                    lift(self)
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

class DataId(dict):
    """A class to represent an dataset, as specified by e.g. visit/ccd

    Note that DataId inherits from dict
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def __str__(self):
        return butler.mapperInfo.dataIdToTitle([self])

    def copy(self, **kwargs):
        """Return a copy, setting any fields specified by kwargs"""
        new = DataId(self)
        for k, v in kwargs.items():
            new[k] = v

        return new

class Data(object):
    def __init__(self, *args, **kwargs):
        self.setButler(*args, **kwargs)

        self.dataId = {}

        self.matches = None
        self.ZP0 = 31
        self.zp = None

        self.name = "??"

    def setVRS(self, reset=False, **dataId):
        """Save the values of visit, raft, sensor; if None use the pre-set values (which are cleared if reset is True)"""

        if reset:
            self.dataId = {}

        self.dataId.update(dataId)
            
    def setButler(self, rerun=None, dataRoot=None, registryRoot=None):
        global butlerDataRoot, butlerRerun, butler, inButler, rerunName

        if \
           (butlerDataRoot and dataRoot != butlerDataRoot) or \
           (rerun and rerun != butlerRerun):
            butler = None
            rerunName = ""

        try:
            if butler is not None:
                return butler
        except NameError:
            butler = None

        if rerun is None and dataRoot is None and not butler:
            raise RuntimeError("Please specify a rerun or root")

        if rerun is None:
            if dataRoot is None:
                if butler:
                    return butler
                raise RuntimeError("Please specify a rerun or root")

        if not registryRoot:
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
                raise RuntimeError("I'm unable to find your registry in %s" % registryRoot)

        Mapper = getMapper(registryRoot, defaultMapper=LsstSimMapper)
        try:
            butler = dafPersist.ButlerFactory(mapper=Mapper(outRoot=dataRoot, rerun=rerun,
                                                            registry=registry)).create()
        except TypeError:               # outRoot/rerun weren't accepted
            butler = dafPersist.ButlerFactory(mapper=Mapper(root=dataRoot, registry=registry)).create()
        butler.mapperInfo = makeMapperInfo(butler.mapper)

        butlerDataRoot = dataRoot
        butlerRerun = rerun

        inButler = butler.mapperInfo.getInButler(dataRoot, registry, butler)

        if rerun:
            rerunName = rerun
        else:
            bdr0 = butlerDataRoot # os.path.split(butlerDataRoot)[0] # chop off last element (e.g. "update")
            if os.path.islink(bdr0):
                bdr0 = os.readlink(bdr0)

            rerunName = os.path.basename(bdr0)

    def lookupDataBySkytile(self, dataType):
        """N.b. not converted to use dataId --- Lsst specific"""
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

    def lookupDataByVisit(self, dataType, dataId):
        """Lookup all the available data of a given type (e.g. "psf") and maybe dataId.

        See also getDataset()
        """

        dataId = dict([(k, v) for k, v in dataId.items() if v is not None])
        fields = butler.mapperInfo.getFields(dataType)

        try:
            dataId["visit"][0]
        except TypeError:
            dataId["visit"] = [dataId["visit"]]
        visits = dataId["visit"]; dataId["visit"] = None

        dataSets = {}
        for v in visits:
            dataId["visit"] = v
            for vals in butler.queryMetadata("raw", "visit", fields, **dataId):
                _dataId = dict(zip(fields, vals))

                if butler.datasetExists(dataType, **_dataId):
                    if not dataSets.has_key(v):
                        dataSets[v] = []
                    dataSets[v].append(_dataId)

            self.dataSets = dataSets

        return self.dataSets

    def getDataset(self, dataType, dataId, ids=True, calibrate=False, setMask=None, fixOrientation=None,
                   fixAmpLevels=False):
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

        if dataType in butler.mapperInfo.getTrimmableData():
            return [butler.mapperInfo.assembleCcd(dataType, inButler, dataId, fixAmpLevels=fixAmpLevels)]
        elif dataType in ("eimage",):
            sdataId = dataId.copy(); sdataId["snap"] = dataId.get("snap", 0)
            raw_filename = inButler.get('raw_filename', channel='0,0', **sdataId)[0]
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
                calexp = self.getDataset("calexp", dataId)[0]
                mask = calexp.getMaskedImage().getMask()
                if not fixOrientation:
                    mask = afwMath.rotateImageBy90(afwMath.flipImage(mask, False, True), 1)
            else:
                mask = None

            return [afwImage.makeExposure(afwImage.makeMaskedImage(eimage, mask), eimageExposure.getWcs())]

        self.setVRS(reset=True, **dataId)

        dataSets = self.lookupDataByVisit(dataType, dataId)
        visit = dataId.get("visit")
        if not visit and len(dataSets) == 1:
            visit = dataSets.keys()[0]

        dataSets = dataSets.get(visit)
        if not dataSets:
            raise RuntimeError("I cannot find your data for visit %d" % visit if visit else "")

        data = []
        for dataId in dataSets:
            dataElem = butler.get(dataType, **dataId)

            if dataType == "calexp":
                psf = butler.get("psf", **dataId)
                dataElem.setPsf(psf)
                del psf
            if dataType == "calexp" and calibrate:
                mi = dataElem.getMaskedImage();
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

    def getPsfs(self, dataId):
        return self.getDataset("psf", dataId)
    
    def getSources(self, dataId):
        dataIdMask = butler.mapperInfo.getDataIdMask(dataId)
        sources = []
        dataSets = self.dataSets
        for ss_file in self.getDataset("src_filename", dataId):
            ss = afwTable.SourceCatalog.readFits(ss_file[0])
            for s in ss:
                s.setId(s.getId() | dataIdMask)

            sources.append(ss)

        if dataSets:
            self.dataSets = dataSets
            
        return sources

    def getMagsByType(self, magType, good=None, magDict=None):
        if magDict:
            mags = magDict[magType]
        else:
            if magType == "ap":
                mags = self.apMags
            elif magType == "inst":
                mags = self.instMags
            elif magType == "model":
                mags = self.modelMags
            elif magType == "psf":
                mags = self.psfMags
            else:
                raise RuntimeError("Unknown magnitude type %s" % magType)

        if good is not None:
            mags = mags[good]

        return mags

    def expandDataId(self, dataId, nSensor=0):
        """Return a list of fully qualified dataIds, given a possibly-ambiguous one

ccd may be a list"""
        
        try:
            ccds = dataId["ccd"]
            dataId["ccd"] = None        # all
        except:
            ccds = None

        dataId = dict([(k, v) for k, v in dataId.items() if v is not None])

        d = self.lookupDataByVisit("src", dataId)
        if not d:
            return []

        self.setVRS(False, **dataId)

        if not dataId["visit"] and len(d) > 0:
            dataId["visit"] = [d.keys()[0]]

        try:
            dataId.get("visit")[0]
        except TypeError:
            dataId["visit"] = [dataId.get("visit")]

        dids = []
        for visit in dataId["visit"]:
            if not visit and len(d) > 0:
                visit = d.keys()[0]

            if not visit:
                raise RuntimeError("Please specify a visit")

            if nSensor > 0:
                d[visit] = d[visit][0:nSensor]

            try:
                dids += [DataId(did) for did in d[visit]]
            except KeyError:
                pass

        if ccds:
            dids = [d for d in dids if d["ccd"] in ccds]
            for v, ds in self.dataSets.items():
                self.dataSets[v] = [d for d in ds if d["ccd"] in ccds]

        return dids

    def getMagsByVisit(self, dataId, nSensor=0, extraApFlux=0.0):
        """Set the "self.magnitudes" for a visit"""

        try:
            dataId.count
        except:
            dataId = [dataId]

        dataIds = sum([self.expandDataId(did, nSensor) for did in dataId], [])
        dataIds = [d for d in dataIds if d]

        ids, flags, stellar, xAstrom, yAstrom, shape, apMags, instMags, modelMags, psfMags = \
             None, None, None, None, None, None, None, None, None, None

        for did in dataIds:
            ss = self.getSources(did)[0]
            ids, flags, stellar, xAstrom, yAstrom, shape, apMags, instMags, modelMags, psfMags = \
                _getMagsFromSS(ss, did, ids=ids, flags=flags, stellar=stellar, extraApFlux=extraApFlux,
                               xAstrom=xAstrom, yAstrom=yAstrom, shape=shape,
                               apMags=apMags, instMags=instMags, modelMags=modelMags, psfMags=psfMags)

        self.name = butler.mapperInfo.dataIdToTitle(dataIds)

        if not ids:
            raise RuntimeError("Failed to read any data for %s" % " ".join([str(d) for d in dataId]))

        self.ids       = np.ndarray(len(ids),       dtype='L', buffer=ids)
        self.flags     = np.ndarray(len(flags),     dtype='L', buffer=flags)
        self.stellar   = np.ndarray(len(stellar),   dtype='d', buffer=stellar)
        self.xAstrom   = np.ndarray(len(xAstrom),   dtype='d', buffer=xAstrom)
        self.yAstrom   = np.ndarray(len(yAstrom),   dtype='d', buffer=yAstrom)
        self.shape = []
        for i in range(len(shape)):
            self.shape.append(np.ndarray(len(shape[i]),   dtype='d', buffer=shape[i]))
        self.apMags    = np.ndarray(len(apMags),    dtype='d', buffer=apMags)
        self.instMags  = np.ndarray(len(instMags),  dtype='d', buffer=instMags)
        self.modelMags = np.ndarray(len(modelMags), dtype='d', buffer=modelMags)
        self.psfMags   = np.ndarray(len(psfMags),   dtype='d', buffer=psfMags)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, dataId, nSensor=0, fixup=False, useDB=False):
        dataIds = self.expandDataId(dataId, nSensor)
        if not dataIds:
            raise RuntimeError("%s doesn't expand to any valid dataIds" % dataId)

        self.setVRS([dataIds])

        self.name = butler.mapperInfo.dataIdToTitle(dataIds)
        filterName = afwImage.Filter(butler.get("calexp_md", **dataIds[0])).getName()

        self.sources = []
        self.matches = []
        self.ref = []
        self.referrs = []
        zpArr = []

        for dataId in dataIds:
            sources, matches, ref, referrs, zp = self._getCalibObjectsImpl(dataId, fixup=fixup)
            
            self.sources += sources
            self.matches += matches
            self.ref += ref
            if referrs:
                self.referrs += referrs

            zpArr.append(zp)

        self.zp = np.mean(zpArr)

    @staticmethod
    def _getCalibObjectsImpl(dataId, fixup=False, verbose=0):
        sources = butler.get('icSrc', **dataId)
        matches = butler.get('icMatch', **dataId)
        
        matchmeta = pmatches.getSourceMatchMetadata()
        
        useOutputSrc = False             # use fluxes from the "src", not "icSrc"
        if useOutputSrc:
            srcs = butler.get('src', **dataId)
            pmMatch = afwDetect.matchXy(sources, srcs, 1.0, True)
            for icSrc, src, d in pmMatch:
                icSrc.setPsfFlux(src.getPsfFlux())

        calexp_md = butler.get('calexp_md', **dataId)
        wcs = afwImage.makeWcs(calexp_md)
        W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")

        filterName = afwImage.Filter(butler.get("calexp_md", **dataId)).getName()

        calib = afwImage.Calib(calexp_md)
        zp = calib.getMagnitude(1.0)

        # ref sources
        xc, yc = 0.5*W, 0.5*H
        radec = wcs.pixelToSky(xc, yc)
        ra = radec.getLongitude()
        dec = radec.getLatitude()
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
        ref = X.refsources
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
            print 'Got', len(ref), 'reference catalog sources'

        keepref = []
        keepi = []
        for i in xrange(len(ref)):
            x, y = wcs.skyToPixel(ref[i].getRa(), ref[i].getDec())
            if x < 0 or y < 0 or x > W or y > H:
                continue
            ref[i].setXAstrom(x)
            ref[i].setYAstrom(y)
            if stargal[i]:
                ref[i].setFlagForDetection(ref[i].getFlagForDetection() | flagsDict["STAR"])
            keepref.append(ref[i])
            keepi.append(i)
            
        if verbose:
            print 'Kept', len(keepref), 'reference sources'
        ref = keepref

        if referrs is not None:
            referrs = [referrs[i] for i in keepi]
        if stargal is not None:
            stargal = [stargal[i] for i in keepi]

        stargal = stargal
        referrs = referrs

        if False:
            m0 = matches[0]
            f,s = m0.first, m0.second
            print 'match 0: ref %i, source %i' % (f.getSourceId(), s.getSourceId())
            print '  ref x,y,flux = (%.1f, %.1f, %.1f)' % (f.getXAstrom(), f.getYAstrom(), f.getPsfFlux())
            print '  src x,y,flux = (%.1f, %.1f, %.1f)' % (s.getXAstrom(), s.getYAstrom(), s.getPsfFlux())

        measAstrom.joinMatchList(matches, ref, first=True, log=log)
        args = {}
        if fixup:
            # ugh, mask and offset req'd because source ids are assigned at write-time
            # and match list code made a deep copy before that.
            # (see svn+ssh://svn.lsstcorp.org/DMS/meas/astrom/tickets/1491-b r18027)
            args['mask'] = 0xffff
            args['offset'] = -1
        measAstrom.joinMatchList(matches, sources, first=False, log=log, **args)

        return sources, matches, ref, referrs, zp
    

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Convert SourceSets to numpy arrays
#
def getMagsFromSS(table, dataId, extraApFlux=0.0):
    """Return numpy arrays constructed from SourceTable table"""

    nrow = len(table)
    ss = table.getColumnView()

    ids = ss["id"]
    flags = np.zeros(nrow, dtype=np.uint64)
    for name, val in flagsDict.items():
        if name == "BINNED1":
            flags |= val
        elif name == "STAR":
            pass
        else:
            flags |= np.where(ss[name], np.uint64(val), np.uint64(0x0))
        
    stellar = ss["classification.extendedness"] < 0.5
    psfMags = np.copy(ss[table.table.getPsfFluxKey()]) # we're going to convert flux to mag
    apMags = np.copy(ss[table.table.getApFluxKey()]) + extraApFlux
    instMags = np.copy(ss[table.table.getInstFluxKey()])
    modelMags = np.copy(ss[table.table.getModelFluxKey()])
    xAstrom = ss[table.table.getCentroidKey().getX()]
    yAstrom = ss[table.table.getCentroidKey().getY()]
    shape = [ss[table.table.getShapeKey().getIxx()],
             ss[table.table.getShapeKey().getIxy()],
             ss[table.table.getShapeKey().getIyy()]]

    calexp_md = butler.get('calexp_md', **dataId)
    calib = afwImage.Calib(calexp_md)
    fluxMag0, fluxMag0Err = calib.getFluxMag0()
    if fluxMag0 <= 0.0:
        fluxMag0 = 1e13
        print >> sys.stderr, "Setting fluxMag0 to %g" % fluxMag0
        calib.setFluxMag0(fluxMag0, fluxMag0Err)

    for i in range(nrow):
        try:
            apMags[i]  = calib.getMagnitude(apMags[i])
        except:
            apMags[i]  = np.nan

        try:
            psfMags[i]  = calib.getMagnitude(psfMags[i])
        except:
            psfMags[i]  = np.nan

        try:
            instMags[i]  = calib.getMagnitude(instMags[i])
        except:
            instMags[i]  = np.nan

        try:
            modelMags[i]  = calib.getMagnitude(modelMags[i])
        except:
            modelMags[i]  = np.nan

    good = np.logical_and(np.isfinite(apMags), np.isfinite(psfMags))

    ids = ids[good]
    flags = flags[good]
    stellar = stellar[good]
    apMags = apMags[good]
    instMags = instMags[good]
    modelMags = modelMags[good]
    psfMags = psfMags[good]
    xAstrom = xAstrom[good]
    yAstrom = yAstrom[good]
    for i in range(len(shape)):
        shape[i] = shape[i][good]

    return ids, flags, stellar, xAstrom, yAstrom, shape, \
        dict(ap=apMags, inst=instMags, model=modelMags, psf=psfMags)

def _getMagsFromSS(ss, dataId, ids=None, flags=None, stellar=None, extraApFlux=0, xAstrom=None, yAstrom=None,
                   shape=None,
                   apMags=None, instMags=None, modelMags=None, psfMags=None):
    """Return python arrays constructed from SourceSet ss, and possibly extending {ap,model,psf}Mags"""
    if not ids:
        ids = array.array('L')
    if not flags:
        flags = array.array('L')
    if not stellar:
        stellar = array.array('d')
    if not xAstrom:
        xAstrom = array.array('d')
    if not yAstrom:
        yAstrom = array.array('d')
    if not shape:
        shape = [array.array('d'), array.array('d'), array.array('d')]
    if not apMags:
        apMags = array.array('d')
    if not modelMags:
        modelMags = array.array('d')
    if not instMags:
        instMags = array.array('d')
    if not psfMags:
        psfMags = array.array('d')

    _ids, _flags, _stellar, _xAstrom, _yAstrom, _shape, _mags = getMagsFromSS(ss, dataId, extraApFlux)
    _apMags, _instMags, _modelMags, _psfMags = _mags["ap"], _mags["inst"], _mags["model"], _mags["psf"]

    ids.extend(_ids)
    flags.extend(_flags)
    stellar.extend(_stellar)
    xAstrom.extend(_xAstrom)
    yAstrom.extend(_yAstrom)
    for i in range(len(shape)):
        shape[i].extend(_shape[i])
    apMags.extend(_apMags)
    instMags.extend(_instMags)
    modelMags.extend(_modelMags)
    psfMags.extend(_psfMags)

    return ids, flags, stellar, xAstrom, yAstrom, shape, apMags, instMags, modelMags, psfMags

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def queryDB(query, host=None, db=None, read_default_file="~/.my.cnf",):
    """Submit a query.  If host/db is None, use default_db_{host,db}
"""
    if host is None:
        host = default_db_host
    if db is None:
        db = default_db_table

    if not db:
        raise RuntimeError("Please set utils.default_db_table")

    conn = MySQLdb.connect(host=host, db=db, read_default_file=read_default_file)

    for q in query.split(";"):
        if q.strip():
            conn.query(q)

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

def plotDmag(data, magType1="model", magType2="psf", maglim=20, magmin=14,
             showMedians=False, sgVal=0.05, meanDelta=0.0, adjustMean=False,
             xmin=None, xmax=None, ymin=None, ymax=None,
             title="+", markersize=1, color="red", frames=[0], fig=None):
    """Plot (magType1 - magType2) v. magType1 mags (e.g. "model" and "psf")

The magnitude limit for the "locus" box used to calculate statistics is maglim, and only objects within
+- sgLim of meanDelta are included in the locus (the axes are also adjusted to include meanDelta in the
visible area).  If adjustMean is True, adjust all the CCDs being plotted to have the same mean within
the "locus" area.

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits (y{min,max} are interpreted relative to meanDelta
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = np.bitwise_and(data.flags, (flagsDict["flags.pixel.interpolated.center"] |
                                      flagsDict["flags.pixel.edge"]))
    suspect = np.bitwise_and(data.flags, (flagsDict["flags.pixel.interpolated.any"] |
                                          flagsDict["flags.pixel.saturated.any"]))
    good = np.logical_and(np.logical_not(bad), np.logical_not(suspect))

    if False:
        amp30 = np.logical_and(data.xAstrom > 1526, data.xAstrom < 2036)
        amp30 = np.logical_and(amp30, data.yAstrom < 2000)
        good = np.logical_and(good, amp30)

    mag1 = data.getMagsByType(magType1, good)
    mag2 = data.getMagsByType(magType2, good)
    delta = mag1 - mag2

    locus = np.logical_and(np.logical_and(mag1 > magmin, mag1 < maglim), np.abs(delta - meanDelta) < sgVal)

    stellar = data.stellar[good] > 0.99
    nonStellar = np.logical_not(stellar)

    try:
        stats = afwMath.makeStatistics(delta[locus], afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
    except:
        mean, stdev = np.nan, np.nan

    if adjustMean:
        ids = data.ids[good]
        ccds = set([x["ccd"] for x in sum(data.dataSets.values(), [])])
        ccdIds = [butler.mapperInfo.splitId(i)[1] for i in ids]
        for ccd in ccds:
            l = np.logical_and(np.equal(ccdIds, ccd), locus)
            try:
                delta[np.equal(ccdIds, ccd)] += \
                    mean - afwMath.makeStatistics(delta[l], afwMath.MEANCLIP).getValue()
            except Exception, e:
                print "Adjusting mean for CCD %d: %s" % (ccd, e)

    if False:
        axes.plot(mag1[locus], delta[locus], "o", markersize=2*markersize, color="blue")

    axes.plot((magmin, maglim, maglim, magmin, magmin), meanDelta + sgVal*np.array([-1, -1, 1, 1, -1]), "b:")
    axes.plot((0, 30), meanDelta + np.array([0, 0]), "b-")

    color2 = "green"
    axes.plot(mag1[nonStellar], delta[nonStellar], "o", markersize=markersize, markeredgewidth=0, color=color)
    axes.plot(mag1[stellar], delta[stellar], "o", markersize=markersize, markeredgewidth=0, color=color2)
    axes.plot(mag1[suspect], delta[suspect], "o", markersize=markersize, markeredgewidth=0, color="blue")
    axes.plot(mag1[bad], delta[bad], "+", markersize=markersize, markeredgewidth=0, color="red")
    #axes.plot((0, 30), (0, 0), "b-")

    if showMedians:
        binwidth = 1.0
        bins = np.arange(int(magmin), int(max(mag1[stellar])), binwidth)
        vals = np.empty_like(bins)
        for i in range(len(bins) - 1):
            inBin = np.logical_and(mag1 > bins[i], mag1 <= bins[i] + binwidth)
            vals[i] = np.median(delta[np.where(np.logical_and(stellar, inBin))])

        axes.plot(bins + 0.5*binwidth, vals, linestyle="-", marker="o", color="cyan")

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
        objs = np.logical_and(np.logical_and(mag1 > magmin, mag1 < magmax),
                              np.logical_and(delta > dmin, delta < dmax))
        ids = data.ids[good]
        for i in sorted(ids[objs]):
            print i, butler.mapperInfo.splitId(i)

    axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
    axes.set_ylim(meanDelta + (0.2 if ymin is None else ymin), meanDelta + (- 0.8 if ymax is None else ymax))
    axes.set_xlabel(magType1)
    axes.set_ylabel("%s - %s" % (magType1, magType2))
    axes.set_title(re.sub(r"^\+\s*", data.name + " ", title))

    fig.text(0.20, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    flags = data.flags[good]
    x = data.xAstrom[good]
    y = data.yAstrom[good]
    ids = data.ids[good]
    eventHandlers[fig] = EventHandler(axes, mag1, delta, ids, x, y, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotSizes(data, magType="psf",
              xmin=None, xmax=None, ymin=None, ymax=None,
              title="+", markersize=1, fig=None):
    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = np.bitwise_and(data.flags, (flagsDict["flags.pixel.interpolated.center"] |
                                      flagsDict["flags.pixel.edge"]))
    good = np.logical_not(bad)

    fig = getMpFigure(fig)
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    stellar = data.stellar[good] > 0.99
    nonStellar = np.logical_not(stellar)

    mag = data.getMagsByType(magType, good)
    size = np.empty(len(data.shape[0]))

    for i in range(len(size)):
        Ixx = data.shape[0][i]
        Ixy = data.shape[1][i]
        Iyy = data.shape[2][i]
        
        q = afwGeom.ellipses.Quadrupole(Ixx, Iyy, Ixy)
        a = afwGeom.ellipses.Axes(q) # convert to (a, b, theta)

        size[i] = math.sqrt(a.getA()*a.getB())

    color = "red"
    color2 = "green"
    axes.plot(mag[nonStellar], size[nonStellar], "o", markersize=markersize, markeredgewidth=0, color=color)
    axes.plot(mag[stellar], size[stellar], "o", markersize=markersize, markeredgewidth=0, color=color2)
    axes.plot(mag[bad], size[bad], "+", markersize=markersize, markeredgewidth=0, color="red")
    #axes.plot((0, 30), (0, 0), "b-")

    axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
    axes.set_ylim(0.0 if ymin is None else ymin, 10 if ymax is None else ymax)
    axes.set_xlabel(magType)
    axes.set_ylabel(r"$\sqrt{a b}$ (pixels)", fontsize="larger")
    axes.set_title(re.sub(r"^\+\s*", data.name + " ", title))

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmagHistograms(data, magType1="model", magType2="psf", color="red", fig=None):
    """Plot (magType1 - magType2) v. magType2 mags"""

    try:
        data.flags
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = np.bitwise_and(data.flags, (flagsDict["flags.pixel.interpolated.center"] |
                                      flagsDict["flags.pixel.edge"]))
    good = np.logical_not(bad)

    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    mag1 = data.getMagsByType(magType1, good)
    mag2 = data.getMagsByType(magType2, good)
    delta = mag1 - mag2

    sgVal = 0.10
    stellar = np.abs(delta) < sgVal

    nbin = 15
    minmag, maxmag, dmag = 15, 20, 0.5
    binEdge = np.arange(minmag, maxmag - dmag/2, dmag)

    subplots = makeSubplots(fig, nx=1, ny=len(binEdge))

    for edge in binEdge:
        inBin = np.logical_and(np.logical_and(mag1 >= edge, mag1 < edge + dmag), stellar)

        pp = mag2[inBin]        
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

    fig.suptitle("%s - %s: %s" % (magType1, magType2, data.name))

    fig.show()

    return fig

def plotCalibration(data, plotBand=0.05, magType='psf', maglim=20,
                    markersize=1, title="+",
                    frame=None, ctype=None, ds9Size=None, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand
    """
    fig = getMpFigure(fig)

    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    mstars = [m for m in data.matches if
              (m.second and (m.second.getFlagForDetection() & flagsDict["STAR"]))] # data
    realStars = [(m.first.getFlagForDetection() & flagsDict["STAR"]) != 0          # catalogue
                 for m in data.matches if (m.first and \
                                           m.second and (m.second.getFlagForDetection() & flagsDict["STAR"]))]

    if frame is not None:
        kwargs = {}
        if ctype:
            kwargs["ctype"] = ctype
        if ds9Size:
            kwargs["size"] = ds9Size
        for m in mstars:
            s = m.second
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), frame=frame, **kwargs)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars if s.first and s.second])
    ids = np.array([s.second.getId() for s in mstars if s.first and s.second])
    flux = [getFlux(s.second, magType) for s in mstars if s.first and s.second]

    instmag = data.zp - 2.5*np.log10(np.array(flux))
    realStars = np.array([m for m in realStars])

    delta = refmag - instmag

    x = np.array([s.second.getXAstrom() for s in mstars])
    y = np.array([s.second.getYAstrom() for s in mstars])
    flags = np.array([s.second.getFlagForDetection() for s in mstars])
    if False:
        for i in range(len(ids)):
            print "%d4 %.2f %.2f (%.2f, %.2f)" % (ids[i], refmag[i], delta[i], x[i], y[i])
    
    stats = afwMath.makeStatistics(delta[refmag < maglim], afwMath.STDEVCLIP | afwMath.MEANCLIP)
    mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)

    axes.plot(refmag, delta, "m.", markersize=markersize)

    minRef = min(refmag)
    axes.plot((minRef, maglim, maglim, minRef, minRef), 100*np.array([-1, -1, 1, 1, -1]), "g:")

    #axes.plot(refmag[realStars], delta[realStars], "gx")
        
    axes.plot((-100, 100), (0, 0), "g-")
    if plotBand:
        for i in (-plotBand, plotBand):
            axes.plot((-100, 100), i*np.ones(2), "g--")

    axes.set_ylim(-0.6, 0.6)
    #axes.set_ylim(-2.1, 2.1)
    axes.set_xlim(13, 24)
    axes.set_xlabel("Reference")
    axes.set_ylabel("Reference - %s" % magType)
    axes.set_title(data.name)

    fig.text(0.75, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    eventHandlers[fig] = EventHandler(axes, refmag, delta, ids, x, y, flags, [frame,])

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

def plotAstromCalibration(data, plotBand=0.5, magType='psf', maglim=20,
                          markersize=1, title="+",
                          frame=None, ctype=None, ds9Size=None, fig=None):
    """Plot (measured - reference) positions given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand
    """
    fig = getMpFigure(fig)

    if not data.matches:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    mstars = [m for m in data.matches if
              (m.second and (m.second.getFlagForDetection() & flagsDict["STAR"]))] # data
    realStars = [(m.first.getFlagForDetection() & flagsDict["STAR"]) != 0          # catalogue
                 for m in data.matches if (m.first and \
                                           m.second and (m.second.getFlagForDetection() & flagsDict["STAR"]))]

    if frame is not None:
        kwargs = {}
        if ctype:
            kwargs["ctype"] = ctype
        if ds9Size:
            kwargs["size"] = ds9Size
        with ds9.Buffering():
            for m in mstars:
                s = m.second
                ds9.dot("o", s.getXAstrom(), s.getYAstrom(), frame=frame, **kwargs)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars if s.first and s.second])

    refx = np.array([s.first.getXAstrom() for s in mstars if s.first and s.second])
    refy = np.array([s.first.getYAstrom() for s in mstars if s.first and s.second])
    srcx = np.array([s.second.getXAstrom() for s in mstars if s.first and s.second])
    srcy = np.array([s.second.getYAstrom() for s in mstars if s.first and s.second])
        
    ids = np.array([s.second.getId() for s in mstars if s.first and s.second])

    delta = refx - srcx

    flags = np.array([s.second.getFlagForDetection() for s in mstars])
    if False:
        for i in range(len(ids)):
            print "%d4 %.2f %.2f (%.2f, %.2f)" % (ids[i], refmag[i], delta[i], x[i], y[i])
    

    stats = afwMath.makeStatistics(delta[refmag < maglim], afwMath.STDEVCLIP | afwMath.MEANCLIP)
    mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)

    left = np.logical_or(srcx < 512, np.logical_and(srcx > 1024, srcx < 1536))
    right = np.logical_not(left)
    axes.plot(refmag[left], delta[left], "m.", markersize=markersize)
    axes.plot(refmag[right], delta[right], "b.", markersize=markersize)

    minRef = min(refmag)
    axes.plot((minRef, maglim, maglim, minRef, minRef), 100*np.array([-1, -1, 1, 1, -1]), "g:")

    #axes.plot(refmag[realStars], delta[realStars], "gx")
        
    axes.plot((-100, 100), (0, 0), "g-")
    if plotBand:
        for i in (-plotBand, plotBand):
            axes.plot((-100, 100), i*np.ones(2), "g--")

    axes.set_ylim(-1.6, 1.6)
    axes.set_xlim(13, 24)
    axes.set_xlabel("Reference Magnitude")
    axes.set_ylabel(r"$\Delta$ (pixels)")
    axes.set_title(data.name)

    if False:
        fig.text(0.75, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    if False:
        global eventHandlers
        eventHandlers[fig] = EventHandler(axes, refmag, delta, ids, srcx, srcy, flags, [frame,])

    fig.show()

    if False:
        _mstars = [m for m in mstars if math.fabs(-2.5*math.log10(m.first.getPsfFlux()) - 17.5) < 0.5]
        with ds9.Buffering():
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

    with ds9.Buffering():
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
    with ds9.Buffering():
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

def findSource(sourceSet, x, y, radius=2):
    ss = afwDetect.SourceSet()

    s = afwDetect.Source();
    s.setXAstrom(x); s.setYAstrom(y)
    ss.push_back(s)

    matches = [m.first for m in afwDetect.matchXy(sourceSet, ss, radius)]
    if len(matches) == 0:
        raise RuntimeError("Unable to find a source within %g pixels of (%g,%g)" % (radius, x, y))

    if len(matches) == 1:
        return matches[0]
    else:
        return matches

def getRefCatalog(data, dataId):
    """Get the reference catalogue for dataId (a single CCD) as a SourceSet"""

    query = """
        select
            scisql_s2CPolyToBin(llcRa, llcDecl, lrcRa, lrcDecl, urcRa, urcDecl, ulcRa, ulcDecl),
            filterId, visit, raftName, ccdName
        from
           Science_Ccd_Exposure
        where
           visit = %(visit)d and raftName = '%(raft)s' and ccdName = '%(sensor)s'
        into @poly, @filterId, @visit, @raftName, @ccdName;

        select
           refObjectId, sro.ra, sro.decl,
           CASE WHEN @filterId = 0 THEN uMag
                WHEN @filterId = 1 THEN gMag
                WHEN @filterId = 2 THEN rMag
                WHEN @filterId = 3 THEN iMag
                WHEN @filterId = 4 THEN zMag
                WHEN @filterId = 5 THEN yMag
           END as Mag,
           sro.isStar, sro.varClass
        from
           SimRefObject as sro
        where
           scisql_s2PtInCPoly(ra, decl, @poly)
       """ % dataId

    return _getSourceSetFromQuery(data, dataId, query)

def showRefCatalog(data, dataId, wcs=None, frame=0):
    ss = data.getSources(dataId)[0]
    refCat = getRefCatalog(data, dataId)

    if frame is not None:
        if not wcs:
            wcs = data.wcs

        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="+")
        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="x", mask=["PEAKCENTER",])
        showSourceSet(refCat, wcs=wcs, frame=frame, ctype=ds9.GREEN, symb="o", mask=["STAR",])

        showSourceSet(ss, wcs=wcs, frame=frame, ctype=ds9.RED, symb="x")

    return ss, refCat

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def _getSourceSetFromQuery(data, dataId, queryStr):
    """Return a SourceSet given a Data object and query string that returns certain information"""
    
    calexp_md = butler.get('calexp_md', **dataId)
    data.wcs = afwImage.makeWcs(calexp_md)
    data.calib = afwImage.Calib(calexp_md)
    data.zp = data.calib.getMagnitude(1.0)

    ss = afwDetect.SourceSet()

    for Id, ra, dec, mag, isStar, varClass in queryDB(queryStr % dataId):
        s = afwDetect.Source();

        s.setId(Id)
        ra, dec = math.radians(ra), math.radians(dec)
        s.setRa(ra); s.setDec(dec)

        x, y = data.wcs.skyToPixel(ra, dec)
        s.setXAstrom(x), s.setYAstrom(y)

        flux = 10**(-0.4*(mag - data.zp))
        s.setApFlux(flux); s.setPsfFlux(flux); s.setModelFlux(flux)

        flag = flagsDict["BINNED1"]
        if isStar > 0:
            flag |= flagsDict["STAR"]
        if varClass != 0:                   # variable
            flag |= flagsDict["flags.badcentroid"] # XXX
        s.setFlagForDetection(flag)

        ss.push_back(s)

    return ss

def getSpurious(data, dataId):
    """Return a SourceSet of all sources that don't match any in the catalogue for dataId (a single CCD)"""

    return _getSourceSetFromQuery(data, dataId, """
select
   s.sourceId, s.ra, s.decl,
   dnToABMag(s.psfFlux, exp.fluxMag0) as psfMag,
   -1, 0
from
   Source as s join
   Science_Ccd_Exposure as exp on s.scienceCcdExposureId = exp.scienceCcdExposureId join
   RefSrcMatch as rsm on rsm.SourceId = s.SourceId
where
   visit = %(visit)d and raftName = '%(raft)s' and ccdName = '%(sensor)s'
   and rsm.refObjectId is Null
   """ % dataId)   

def getMatched(data, dataId):
    """Return a SourceSet of all the catalogue objects that we detected for dataId (a single CCD)"""

    query = """
select
   sro.refObjectId, sro.ra, sro.decl,
   CASE WHEN exp.filterId = 0 THEN uMag
        WHEN exp.filterId = 1 THEN gMag
	WHEN exp.filterId = 2 THEN rMag
	WHEN exp.filterId = 3 THEN iMag
	WHEN exp.filterId = 4 THEN zMag
	WHEN exp.filterId = 5 THEN yMag
   END as Mag,
   sro.isStar, sro.varClass
from
   Source as s join
   Science_Ccd_Exposure as exp on s.scienceCcdExposureId = exp.scienceCcdExposureId join
   RefSrcMatch as rsm on rsm.SourceId = s.SourceId join
   SimRefObject as sro on sro.refObjectId = rsm.refObjectId
where
   visit = %(visit)d and raftName = '%(raft)s' and ccdName = '%(sensor)s'
   and rsm.refObjectId is not Null
   """

    return _getSourceSetFromQuery(data, dataId, query)
    

def getMissed(data, dataId):
    """Return a SourceSet of all the catalogue objects that we failed to detect for dataId (a single CCD)"""

    query = """
        SELECT
           scisql_s2CPolyToBin(llcRa, llcDecl, lrcRa, lrcDecl, urcRa, urcDecl, ulcRa, ulcDecl),
           scienceCcdExposureId, filterId
        FROM
          Science_Ccd_Exposure
        WHERE
            visit = %(visit)d and raftName = '%(raft)s' and ccdName = '%(sensor)s'
        INTO @poly, @sceId, @filterId;
        
        SELECT
           sro.refObjectId, sro.ra, sro.decl,
           CASE WHEN @filterId = 0 THEN uMag
                WHEN @filterId = 1 THEN gMag
        	WHEN @filterId = 2 THEN rMag
        	WHEN @filterId = 3 THEN iMag
        	WHEN @filterId = 4 THEN zMag
        	WHEN @filterId = 5 THEN yMag
           END as Mag,
           sro.isStar, sro.varClass
        FROM
           SimRefObject AS sro LEFT OUTER JOIN
           (
               SELECT DISTINCT rsm.refObjectId AS refObjectId
               FROM RefSrcMatch AS rsm INNER JOIN
                    Source AS s ON (rsm.sourceId = s.sourceId)
               WHERE s.scienceCcdExposureId = @sceId
           ) AS t ON (sro.refObjectId = t.refObjectId)
        WHERE t.refObjectId IS NULL AND scisql_s2PtInCPoly(sro.ra, sro.decl, @poly) = 1;
        """

    return _getSourceSetFromQuery(data, dataId, query)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def showSourceSet(sourceSet, exp=None, wcs=None, raDec=None, magmin=None, magmax=None, magType="psf",
                  mask=None, symb="+", **kwargs):
    """Show a SourceSet on ds9.

    If mask, it's a set of bitplane names (e.g. INTERP) which must be set to display the source"""
    
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

    if magmin is not None or magmax is not None:
        if not exp:
            raise RuntimeError("You'll have to provide an exposure if you want me to use magnitude limits")
        if magmin is None:
            magmin = -100

    with ds9.Buffering():
        for s in sourceSet:
            if mask:
                if not (s.getFlagForDetection() & mask):
                    continue

            if magmax is not None:
                if exp:
                    mag = exp.getCalib().getMagnitude(getFlux(s, magType))

                    if mag < magmin or mag > magmax:
                        continue

            if raDec:
                ra, dec = s.getRa(), s.getDec()

                x, y = wcs.skyToPixel(ra, dec)
            else:
                x, y = s.getXAstrom(), s.getYAstrom()

            ds9.dot((symb if symb != "id" else ("%d" % s.getId())), x, y, **kwargs)

def showCatalog(data, dataId, calexp=None, refOnly=None, detected=None, spurious=None, ss=None,
                frame=0, **kwargs):
    if not calexp:
        calexp = data.getDataset("calexp", dataId)[0]

    if not ss:
        ss = data.getSources(dataId)[0]

    if not refOnly:
        refOnly = getMissed(data, dataId)
    if not detected:
        detected = getMatched(data, dataId)
    if not spurious:
        spurious = getSpurious(data, dataId)

    showSourceSet(refOnly, calexp, ctype=ds9.BLUE, symb="o", size=3, frame=frame, **kwargs)
    showSourceSet(detected, calexp, ctype=ds9.GREEN, symb="o", size=3, frame=frame, **kwargs)
    showSourceSet(spurious, calexp, ctype=ds9.RED, symb="o", size=3, frame=frame, **kwargs)

    showSourceSet(ss, calexp, ctype=ds9.YELLOW, frame=frame, **kwargs)

def plotCompleteness(data, dataId, refCat=None, calexp=None, matchRadius=2, **kwargs):
    """Plot completeness plots in matplotlib (and maybe ds9)"""

    data.name = butler.mapperInfo.dataIdToTitle(dataId)

    ss = data.getSources(dataId)[0]

    if False:
        if not refCat:
            refCat = getRefCatalog(data, dataId)

        matched, spurious, refOnly = getMatches(ss, refCat, matchRadius)
        detected = matched["src"]
    else:
        refOnly = getMissed(data, dataId)
        detected = getMatched(data, dataId)
        spurious = getSpurious(data, dataId)

    if not calexp:
        calexp = data.getDataset("calexp", dataId)[0]

    orefOnly = refOnly
    blended = afwDetect.SourceSet()
    refOnly = afwDetect.SourceSet()

    for s in orefOnly:
        x, y = int(s.getXAstrom() + 0.5), int(s.getYAstrom() + 0.5)
        try:
            mval = calexp.getMaskedImage().getMask().get(x, y)
        except pexExcept.LsstCppException:
            mval = 0x0

        if mval & afwImage.MaskU.getPlaneBitMask("EDGE"):
            continue

        if mval & afwImage.MaskU.getPlaneBitMask("DETECTED"):
            blended.push_back(s)
            continue

        refOnly.push_back(s)
        
    return (plotCounts(data, detected, refOnly, spurious, blended, ss, calexp,
                       **kwargs), ss, refCat)

def plotCounts(data, matched, refOnly, spurious, blended, ss=None, calexp=None,
               includeSpurious=True, magType="model",
               magmin=None, magmax=None, dmag=0.5, stars=None, galaxies=None,
               log=True, stacked=True, title=None, completeness=False, fig=None, frame=None):

    if magmin is None:
        magmin = 14
    if magmax is None:
        magmax = 25

    if not title:
        title = data.name

    if log:
        stacked = False

    if stars is None and galaxies is None:
        stars, galaxies = False, False

    if stars or galaxies:
        title += " [%s]" % ("stars" if stars else "galaxies")
    if stacked:
        title += " [stacked histograms]"

    bins = np.arange(magmin, magmax + 0.5*dmag, dmag)

    arrays = []                         # our data arrays
    xyPos = []                          # [xy]Astrom for objects
    for s in [matched, blended, refOnly, spurious]:
        ids, flags, stellar, xAstrom, yAstrom, shape, mags = getMagsFromSS(s, data.dataId)

        bad = np.bitwise_and(flags, (flagsDict["flags.pixel.interpolated.center"] |
                                     flagsDict["flags.pixel.edge"]))

        if stars or galaxies:
            isStar = np.bitwise_and(flags, flagsDict["STAR"])
            if stars:
                bad = np.logical_or(bad, np.logical_not(isStar))
            else:
                bad = np.logical_or(bad, isStar)
                 
        mags = data.getMagsByType(magType, np.logical_not(bad), mags)

        xyPos.append(zip(mags, xAstrom[np.logical_not(bad)], yAstrom[np.logical_not(bad)]))
        arrays.append(np.histogram(mags, bins)[0])

    if ss:
        m, x, y = [], [], []
        for s in ss:
            try:
                mag = calexp.getCalib().getMagnitude(getFlux(s, magType))
            except pexExcept.LsstCppException:
                mag = np.nan
                
            m.append(mag)
            x.append(s.getXAstrom())
            y.append(s.getYAstrom())

        xyPos.append(zip(m, x, y))
    #
    # ds9?
    #
    if frame is not None:
        ds9.erase(frame=frame)

        for i, ct in enumerate([ds9.GREEN, ds9.CYAN, ds9.BLUE, ds9.RED, ds9.YELLOW]):
            if i >= len(xyPos):
                break
            
            symb = "+" if ct == ds9.YELLOW else "o"
            with ds9.Buffering():
                for mag, x, y in xyPos[i]:
                    if mag > magmin and mag < magmax:
                        ds9.dot(symb, x, y, size=3, frame=frame, ctype=ct)
    #
    # Time for matplotlib
    #
    fig = getMpFigure(fig)

    plottingArea = (0.1, 0.1, 0.85, 0.80)
    axes = fig.add_axes(plottingArea)

    kwargs = {}
    if stacked:
        alpha = 1.0
    else:
        alpha = 0.7

    left = bins[0:-1] - 0.5*dmag        # left side of bins
    try:
        axes.bar(left, arrays[0], width=dmag, log=log, label="matched", color="green", alpha=alpha)
    except ValueError:
        pass

    if stacked:
        kwargs["bottom"] = arrays[0]
    try:
        axes.bar(left, arrays[1], width=dmag, log=log, label="blended", color="cyan", alpha=alpha, **kwargs)
    except ValueError:
        pass

    if stacked:
        kwargs["bottom"] = arrays[0] + arrays[1]
    try:
        axes.bar(left, arrays[2], width=dmag, log=log, label="refOnly", color="blue", alpha=alpha, **kwargs)
    except ValueError:
        pass

    if includeSpurious:
        axes.bar(left, arrays[3], width=dmag, log=log, label="spurious", alpha=alpha, color="red")

    axes.legend(loc=2)
    axes.set_xlim(magmin - 0.75*dmag, magmax)
    axes.set_ylim(0.5 if log else -1, 1.05*max(arrays[0] + arrays[1] + arrays[2] +
                                               arrays[3]*(1 if includeSpurious else 0)))
    axes.set_xlabel(magType)

    axes.set_title(title)

    if completeness:
        detected = arrays[0] + arrays[1]
        missed = arrays[2]

        detected[detected + missed == 0] = 1
        completeness = 1.0*detected/(detected + missed)

        axes2 = fig.add_axes(plottingArea, frameon=False)
        axes2.yaxis.tick_right()
        axes2.plot(bins[:-1], completeness, color="black")
        axes2.set_xlim(axes.get_xlim())
        axes2.yaxis.set_label_position('right')
        axes2.set_ylabel('Completeness')
        axes2.set_ylim(-0.05, 1.05)
        axes.set_yticks([])

    fig.show()

    return fig

def getMatches(ss, refCat, radius=2):
    """Get the matches between a SourceSet of measured Sources and another, refCat, generated from the DB

    Return three lists:
       matches          As dictionary of two SourceSets (keys ref and src) of Sources from refCat and ss
       spurious     A SourceSet of Sources from ss that don't match
       refOnly          A SourceSet of Sources from refCat that don't match
    """
    matches = afwDetect.matchRaDec(refCat, ss, 2000) # all objects should match

    matchedRefs = {}                    # matched pairs;  reference object is unique
    matchedSrcs = {}                    # list of matched objects from ss
    
    spurious = afwDetect.SourceSet()
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
            # Set the quality flags from the src, and STAR from the ref
            ref.setFlagForDetection( ref.getFlagForDetection() |
                                    (src.getFlagForDetection() & ~flagsDict["STAR"]))
            src.setFlagForDetection( src.getFlagForDetection() |
                                    (ref.getFlagForDetection() & flagsDict["STAR"]))

            matched["ref"].push_back(ref)
            matched["src"].push_back(src)
        else:
            refOnly.push_back(ref)

    sids = set()
    for ref, src, distance in matchedSrcs.values():
        sids.add(src.getId())
        if distance >= radius:
            spurious.push_back(src)
    #
    # It's possible for a src to match to no reference object (think of one in the interior of a triangle)
    #
    for src in ss:
        if src.getId() not in sids:
            spurious.push_back(src)
            
    return matched, spurious, refOnly

def subtractModels(data, dataId, magType="psf", frame=0, subtractedOnly=False, showExposure=True, ccd=None,
                   fitAmplitude=False):
    """Show and return all the images for some objectId and filterId (a maximum of maxDs9 are displayed, starting at ds9 frame0).
If you specify visits, only that set of visits are considered for display
If showExposure is False, show the Image not the Exposure
If subtractedOnly is True, don't show the unsubtracted images
If fitAmplitude, fit the object's amplitudes rather than using the measured fluxes
    """

    if ccd is not None:
        dataId = dataId.copy()
        dataId["ccd"] = ccd

    exp = data.getDataset("calexp", dataId)[0]
    exp.getPsf().setDetector(exp.getDetector())
    psf = exp.getPsf()

    subtracted =  exp.Factory(exp, True)

    for s in data.getSources(dataId)[0]:
        if (s.getFlagForDetection() & flagsDict["flags.pixel.saturated.center"]):
            continue

        flux = np.nan if fitAmplitude else getFlux(s, magType)
        try:
            if s.getApDia() > 0.5:
                measAlg.subtractPsf(psf, subtracted.getMaskedImage(), s.getXAstrom(), s.getYAstrom(), flux)
        except pexExcept.LsstCppException, e:
            pass

    title = butler.mapperInfo.dataIdToTitle([dataId])

    ds9.setMaskTransparency(75)
    if showExposure:
        if not subtractedOnly:
            ds9.mtv(exp, title=title, frame=frame)
        ds9.setMaskTransparency(95)
        ds9.mtv(subtracted, title=title, frame=frame + 1)
    else:
        if not subtractedOnly:
            ds9.mtv(exp.getMaskedImage().getImage(), wcs=ims[visit][0].getWcs(), title=title, frame=frame)
        ds9.setMaskTransparency(95)
        ds9.mtv(subtracted.getMaskedImage().getImage(), title=title, frame=frame + 1)


def plotObject(exp, xc, yc, dir="-", hlen=10, findPeak=0, frame=None, fig=None):
    """Plot a cross section through an object centered at (xc, yc) on the specified figure
    (horizontally if dir=="-", vertically if |"|")

    Search a region of size +-findPeak for the highest pixel
    """
    try:
        image = exp.getMaskedImage().getImage()
    except AttributeError:
        try:
            image = exp.getImage()
        except AttributeError:
            pass

    xc = int(xc + 0.5)
    yc = int(yc + 0.5)

    xc0, yc0, maxVal = xc, yc, image.get(xc, yc)
    for x in range(xc0 - findPeak, xc0 + findPeak + 1):
        for y in range(yc0 - findPeak, yc0 + findPeak + 1):
            if image.get(x, y) > maxVal:
                xc, yc, maxVal = x, y, image.get(x, y)
    
    if dir == "-":
        x0, x1 = xc - hlen, xc + hlen
        y0, y1 = yc, yc
    elif dir == "|":
        x0, x1 = xc, xc
        y0, y1 = yc - hlen, yc + hlen
    else:
        raise RuntimeError("dir should be - or |, not %s" % dir)

    z = np.empty(2*hlen + 1)
    I = np.empty(len(z))
    i = 0
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            I[i] = image.get(x, y)
            z[i] = x if dir == "-" else y
            i += 1

    if frame is not None:
        ds9.pan(xc, yc, frame=frame)
        ds9.dot("+", xc, yc, frame=frame)
        ds9.line([(x0, y0), (x1, y1)], frame=frame)
    
    fig = getMpFigure(fig)
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    axes.plot(z, I, "o", color="b")
    axes.plot(z, 0.5*maxVal + 0*I, "r--")
    axes.set_xlim(min(z) - 1, max(z) + 1)
    axes.set_ylim(-0.01*max(I), 1.05*max(I))
    axes.set_title("%s (%d, %d)" % (butler.mapperInfo.exposureToStr(exp), xc, yc))

    fig.show()

    return fig, z, I

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

    fig.suptitle("PSF eigen components, %s [%s]" % (butler.mapperInfo.exposureToStr(exposure), rerun),
                 fontsize=14) 

    for i in range(1, nfunc):
        axes = subplots.next()
        CS = axes.contour(X, Y, Z[i], 20, colors="red")
        axes.clabel(CS, inline=1, fontsize=10)
        axes.set_xlim(-5, width - 1 + 5)
        axes.set_ylim(-5, height - 1 + 5)

        axes.text(0.1, 0.1, "Component %d" % i, transform=axes.transAxes)

    fig.show()

    return fig

def getFlux(s, magType="psf"):
    """Return the desired type of flux"""
    if magType == "ap":
        return s.getApFlux()
    elif magType == "inst":
        return s.getInstFlux()
    elif magType == "model":
        return s.getModelFlux()
    elif magType == "psf":
        return s.getPsfFlux()
    else:
        raise RuntimeError("Uknown magnitude type %s" % magType)

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

def assembleCcdLsst(dataType, butler, dataId, snap=0, fixAmpLevels=False):
    """Return an Exposure of a raw data frame (or related image stored per-amp such as a flat or bias)

    @param dataType	Desired type (e.g. "raw")
    @param butler       An input butler
    @param snap	        The desired snap, if appropriate
    @param fixAmpLevels       Renormalise each amp image to unit gain
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

    return ipIsr.assembleCcd(ampList, ccd, fixAmpLevels=fixAmpLevels)

def assembleCcdSubaru(dataType, butler, dataId, fixAmpLevels=False):
    """Return an Exposure of a raw data frame (or related image stored per-amp such as a flat or bias)

    @param dataType	Desired type (e.g. "raw")
    @param butler       An input butler
    @param dataId       Dictionary of visit, ccd, etc.
    """
    exp = cameraGeomUtils.trimExposure(butler.get(dataType, **dataId))
    if fixAmpLevels:
        ampWidth = 512
        for x in range(4):
            bbox = afwGeom.BoxI(afwGeom.PointI(ampWidth*x, 0), afwGeom.ExtentI(ampWidth, exp.getHeight()))
            smi = exp.Factory(exp, bbox).getMaskedImage()
            median = afwMath.makeStatistics(smi, afwMath.MEDIAN).getValue()
            if x == 0:
                median0 = median
            else:
                smi += int(median0 - median + 0.5)

    return exp

class EventHandler(object):
    """A class to handle key strokes with matplotlib displays"""
    def __init__(self, axes, xs, ys, ids, x, y, flags, frames=[0]):
        self.axes = axes
        self.xs = xs
        self.ys = ys
        self.ids = ids
        self.x = x
        self.y = y
        self.flags = flags
        self.frames = frames

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
            print ("%.3f %.3f" % (ev.xdata, ev.ydata)), butler.mapperInfo.splitId(objId), \
                  maUtils.explainDetectionFlags(long(self.flags[objId == self.ids][0])),
            sys.stdout.flush()

            if ev.key == "I":
                print ""
            elif ev.key in ("p", "P"):
                x = self.x[objId == self.ids][0]
                y = self.y[objId == self.ids][0]
                for frame in self.frames:
                    ds9.pan(x, y, frame=frame)
                ds9.cmdBuffer.flush()
        else:
            pass

def showPsfResiduals(data, dataId, fluxLim=9e4, sigma=0, frame=0, **kwargs):
    """Suprime cam specific!"""
    dataId = dict(visit = dataId["visit"])

    mos = ds9Utils.Mosaic(gutter=4, background=-10)

    for ccd in (8, 9, 5, 4, 3, 6, 7, 2, 1, 0,):
        dataId["ccd"] = ccd

        calexp = data.getDataset("calexp", dataId)[0]

        ss = [s for s in data.getSources(dataId)[0] if
              s.getPsfFlux() > fluxLim and not (s.getFlagForDetection() &
                                                flagsDict["flags.pixel.saturated.any"]) and
              abs(2.5*math.log10(s.getPsfFlux()/s.getApFlux())) < 0.05]
        im = maUtils.showPsfResiduals(calexp, ss, frame=None, **kwargs)
        if sigma > 0:
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            kernel = afwMath.SeparableKernel(int(5*sigma), int(5*sigma), gaussFunc, gaussFunc)
            cim = im.getImage().Factory(im.getImage(), True)
            afwMath.convolve(cim, im.getImage(), kernel)
            im = cim

        mos.append(im, "%d" % ccd)

    return mos.makeMosaic(frame=frame, title="%s %s" % (data.name.split()[-1], data.name.split()[0][-5:-1]), mode=5)

def showStackedPsfResiduals(data=None, dataId={}, exposure=None, sourceSet=None, sigma=None,
                            magType="psf", magMin=14, magMax=22, dmag=0.5, gutter=4,
                            normalize=True, frame=None):
    """Stack PSF residuals binned by magnitude"""

    if not exposure or not sourceSet:
        if not (data and dataId):
            raise RuntimeError("Please specify Data and dataId or an exposure and sourceSet")
        if not exposure:
            exposure = data.getDataset("calexp", dataId)[0]
                
        if not sourceSet:
            sourceSet = data.getSources(dataId)[0]

    mimIn = exposure.getMaskedImage()
    mimIn = mimIn.Factory(mimIn, True)  # make a copy to subtract from

    Image = mimIn.getImage().Factory    # Image ctor
    
    psf = exposure.getPsf()
    psfWidth, psfHeight = psf.getLocalKernel().getDimensions()
    #
    # Make the images that we'll stack into
    #
    magMax -= 1e-4
    def bin(mag):
        return int((mag - magMin)/dmag)
    nbin = bin(magMax) + 1

    residualImages = []
    for i in range(nbin):
        residualImages.append([Image(psfWidth, psfHeight), Image(psfWidth, psfHeight),
                               magMin + (i + 0.5)*dmag, 0])

    for s in sourceSet:
        x, y = s.getXAstrom(), s.getYAstrom()
        
        try:
            flux = getFlux(s, magType)
            mag = exposure.getCalib().getMagnitude(flux)

            if not (magMin <= mag <= magMax):
                continue

            if (s.getFlagForDetection() & flagsDict["flags.pixel.saturated.any"]) or \
               abs(2.5*math.log10(s.getPsfFlux()/s.getApFlux())) > 0.05:
                continue

            bbox = afwGeom.BoxI(afwGeom.PointI(int(x) - psfWidth//2, int(y) - psfHeight//2),
                                afwGeom.ExtentI(psfWidth, psfHeight))
            dx, dy = x - int(x), y - int(y) # offset from centre of subimage

            residualImages[bin(mag)][0] += afwMath.offsetImage(Image(mimIn.getImage(),
                                                                     bbox, afwImage.PARENT), -dx, -dy)
            flux = np.nan
            chi2 = measAlg.subtractPsf(psf, mimIn, x, y, flux)
        except (pexExcept.LsstCppException, ValueError), e:
            continue

        expIm = Image(mimIn.getImage(), bbox, afwImage.PARENT)
        expIm = afwMath.offsetImage(expIm, -dx, -dy)

        residualImages[bin(mag)][1] += expIm
        residualImages[bin(mag)][3] += 1
        
    objects = ds9Utils.Mosaic(gutter=gutter)
    residuals = ds9Utils.Mosaic(gutter=gutter)

    for obj, res, mag, n in residualImages:
        if normalize and n > 0:
            peak = afwMath.makeStatistics(obj, afwMath.MAX).getValue()
            obj /= peak
            res /= peak

        if sigma > 0:
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            kernel = afwMath.SeparableKernel(int(5*sigma + 1), int(5*sigma + 1), gaussFunc, gaussFunc)
            cres = res.Factory(res, True)
            afwMath.convolve(cres, res, kernel, afwMath.ConvolutionControl(True, True))
            res = cres

        lab = "%.2f %d" % (mag, n) if n else ""
        objects.append(obj, lab)
        residuals.append(res, lab)
        
    title = str(exposure.getDetector().getId())
    mosaics = []
    nx = 2
    mosaics.append(objects.makeMosaic(mode=nx, title=title, frame=frame))
    mosaics.append(residuals.makeMosaic(mode=nx, title=title, frame=None if frame is None else frame+1))
    mosaics.append([(mag, n) for obj, res, mag, n in residualImages])

    return mosaics

def showStackedPsfResidualsCamera(data, dataId, frame=0, overlay=False, normalize=False, **kwargs):
    """
    Show the stackedPsfResiduals laid out in true camera positions (unless overlay is True, in which case the
    chip are simply added

    Position of CCDs in Suprime cam specific!

    """

    subGutter, gutter = 4, 2
    mos = ds9Utils.Mosaic(gutter=gutter, background=0.02)

    dataId = dataId.copy()

    mags, n = None, None

    ccds = (8, 9, 5, 4, 3, 6, 7, 2, 1, 0,)
    for ccd in ccds:
        dataId["ccd"] = ccd

        try:
            object, resid, labels = showStackedPsfResiduals(data, dataId, normalize=not overlay and normalize,
                                                            gutter=subGutter, frame=None, **kwargs)
        except (RuntimeError, TypeError), e:
            print "Failed to return residuals for %s: %s" % (dataId, e)
            continue
            
        if not mags:
            mags  = [m for m, n in labels]
            nstar = np.array([m for m, n in labels])
        else:
            nstar += [n for m, n in labels]

        if overlay:
            if not mos.images:
                mos.append(resid)
            else:
                mos.images[0] += resid
        else:
            mos.append(resid, "%d" % ccd)

    im = mos.makeMosaic(frame=None if overlay else frame,
                        title="%s %s" % (data.name.split()[-1], data.name.split()[0][-5:-1]),
                        mode=1 if overlay else len(ccds)//2)

    if overlay:
        w, h = im.getDimensions()
        npanel = len(mags); nx = 2; ny = npanel//nx
        if ny*nx != npanel:
            ny += 1
            
        w = (w - (nx - 1)*subGutter)//nx
        h = (h - (ny - 1)*subGutter)/ny

        labs = []
        for i, mn in enumerate(zip(mags, nstar)):
            m, n = mn
            ix, iy = i%nx, i//nx
            x, y = ix*(w + subGutter), iy*(h + subGutter)
            bbox = afwGeom.BoxI(afwGeom.PointI(x, y), afwGeom.ExtentI(w, h))
            sim = im.Factory(im, bbox, afwImage.PARENT)

            peak = afwMath.makeStatistics(sim, afwMath.MAX).getValue()
            
            sim /= peak

            labs.append(["%.2f %d" % (m, n), x + 0.5*w, y])

        ds9.mtv(im, frame=frame, title=str(dataId))
        with ds9.Buffering():
            for lab, x, y in labs:
                ds9.dot(lab, x, y, frame=frame)

    return im

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def smooth(x, windowLen, windowType="boxcar"):
    """Adapted from http://www.scipy.org/Cookbook/SignalSmooth"""
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
 
    if x.size < windowLen:
        raise ValueError("Input vector needs to be bigger than window size.")

    if windowLen < 3:
        return x

    if windowType == "boxcar":
        w = np.ones(windowLen,'d')
    elif windowType == "hamming":
        w = np.hamming(windowLen)
    elif windowType == "hanning":
        w = np.hanning(windowLen)
    elif windowType == "bartlett":
        w = np.bartlett(windowLen)
    elif windowType == "blackman":
        w = np.blackman(windowLen)
    else:
        raise ValueError("windowType %s is not one of 'boxcar', 'hanning', 'hamming', 'bartlett', 'blackman'"
                         % windowType)

    s = np.r_[x[windowLen-1:0:-1],x,x[-1:-windowLen:-1]]

    y = np.convolve(w/w.sum(), s, mode='valid')

    return y[windowLen//2:-windowLen//2 + 1]

def plotImageHistogram(calexp, minDN=None, maxDN=None, binwidth=None,
                       bitmask=["DETECTED", "EDGE"], showStats=False,
                       windowLen=0, windowType='hamming', maxBins=1000,
                       title=None, log=False, showUnclipped=False, fig=None):
    """Plot a histogram of image pixels
    @param calexp The Exposure/MaskedImage/Image to process
    @param minDN    The minimum value to plot (default: -maxDN)
    @param maxDN    The maximum value to plot (default: 1000)
    @param binwidth The width of the bins (default: 1)
    @param maxBins  The maximum number of bins in the histogram
    @param bitmask Hex mask for pixels to ignore, or list of bit names (default: ["DETECTED"])
    @param showStats Show median etc. estimated from the Exposure
    @param title    Title of plot
    @param log      Take log of histogram?
    """

    try:
        mi = calexp.getMaskedImage()
    except AttributeError:
        mi = calexp

    try:
        img = mi.getImage()
    except AttributeError:
        img = mi

    try:
        msk = mi.getMask()
    except AttributeError:
        msk = None

    if msk is not None:
        try:
            _bitmask = 0x0
            for n in bitmask:
                _bitmask |= msk.getPlaneBitMask(n)
            bitmask = _bitmask
        except TypeError:
            bitmask = _bitmask

    img = img.getArray()
    if msk is not None:
        msk = msk.getArray()
   
    if maxDN is None:
        maxDN = np.max(img)
    if minDN is None:
        minDN = np.min(img)
    if binwidth is None:
        binwidth = 1

    bins = np.arange(minDN, maxDN, binwidth)
    if len(bins) > maxBins:
        raise RuntimeError("Too many bins: %d > maxBins == %d" % (len(bins), maxBins))

    if msk is None:
        sky = np.histogram(img, bins)[0]
        showUnclipped = False
    else:
        x,y = np.where(np.logical_not(np.bitwise_and(msk, bitmask)))
        sky = np.histogram(img[x,y], bins)[0]
        if showUnclipped:
            unclipped = np.histogram(img, bins)[0]

    fig = getMpFigure(fig)

    plottingArea = (0.1, 0.1, 0.85, 0.80)
    axes = fig.add_axes(plottingArea)

    if showUnclipped:
        axes.bar(bins[0:-1] - 0.5*binwidth, unclipped, width=binwidth, log=log,
                 color="red", linewidth=0, alpha=0.8)
    axes.bar(bins[0:-1] - 0.5*binwidth, sky, width=binwidth, log=log,
             color="green", linewidth=0, alpha=0.8)
    if windowLen > 3:
        axes.bar(bins[0:-1] - 0.5*binwidth, smooth(sky, windowLen, windowType), width=binwidth, log=log,
                 color="blue", linewidth=0, alpha=0.4)
    axes.set_xlim(minDN - 0.75*binwidth, maxDN - 1.25*binwidth)
    axes.set_xlabel("DN")

    axes.set_ylim(axes.get_ylim()[0], 1.05*axes.get_ylim()[1])
    if log:
        ymin, ymax = axes.get_ylim()
        axes.set_ylim(max([0.5, ymin]), ymax)
        
    if title:
        axes.set_title(title)

    if showStats:
        sctrl = afwMath.StatisticsControl()

        vals = ["MEAN", "MEDIAN", "MEANCLIP", "STDEVCLIP",]
        if msk is not None:
            sctrl.setAndMask(bitmask)

        stats = afwMath.makeStatistics(mi,
                                       reduce(lambda x, y:
                                                  x | afwMath.stringToStatisticsProperty(y), vals, 0), sctrl)

        ctypes = ["red", "blue", "green", "yellow", "cyan"]
        for i, v in enumerate(vals):
            val = stats.getValue(afwMath.stringToStatisticsProperty(v))
            axes.plot((val, val), axes.get_ylim(), linestyle="-", color=ctypes[i%len(ctypes)],
                      label="%-9s %.2f" % (v, val))
        i += 1
        axes.plot((0, 0), axes.get_ylim(), linestyle="--", color=ctypes[i%len(ctypes)])
            
        axes.legend(loc=1)

    fig.show()

    return fig
