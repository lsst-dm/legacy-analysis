"""
Some utilities for looking at the outputs of processing runs
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
import lsst.afw.table as afwTable
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
from lsst.meas.algorithms.detection import SourceDetectionTask
import lsst.meas.algorithms.utils as maUtils

try:
    from lsst.obs.lsstSim import LsstSimMapper
    _raw_ = "raw"
    _visit_ = "visit"
except ImportError:
    class LsstSimMapper(object): pass

try:
    from lsst.obs.suprimecam.suprimecamMapper import SuprimecamMapper
    from lsst.obs.suprimecam.suprimecamMapper import SuprimecamMapperMit
    _raw_ = "raw"
    _visit_ = "visit"
except ImportError:
    class SuprimecamMapper(object): pass

try:
    from lsst.obs.sdss.sdssMapper import SdssMapper
    _raw_ = "fpC"
    _visit_ = "run"
except ImportError:
    class SdssMapper(object): pass

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

#
try:
    _prefix_                            # prefix for data type names
except NameError:
    _prefix_ = ""

def dtName(dType, md=False):
    """Get the name of a given data type (e.g. dtName("src"))"""
    dType_base = re.sub(r"_(filename|md)$", "", dType)
    if _prefix_ in ("", "forced",) or dType_base in ("camera", "coaddTempExp", "goodSeeingCoadd",):
        if md:
            dType += "_md"
        return dType

    if _prefix_ == "forced":
        if dType == "src":
            return "forcedsources"
        else:
            return dType

    dType = "%s_%s" % (_prefix_, dType)
    if md:
        dType += "_md"

    return dType
#
# These functions are only useful to ease the transition to the Summer2012 Source schema
#
if True:
    def getFlagForDetection(s, flagName):
        if flagName == "STAR":
            return s.get("classification.extendedness") < 0.5
        elif flagName == "BINNED1":
            return True
        elif flagName == "PEAKCENTER":
            keyName = "centroid.sdss.err"
        elif flagName == "SATUR_CENTER":
            keyName = "flags.pixel.saturated.center"
        else:
            raise RuntimeError("Unknown flag: %s" % (flagName))

        return s.get(keyName)

    def setFlagForDetection(s, flagName, value):
        if flagName == "STAR":
            s.set("classification.extendedness", 0.0 if value else 1.0)
            return
        elif flagName == "BINNED1":
            return
        elif flagName == "PEAKCENTER":
            keyName = "centroid.sdss.err"
        else:
            raise RuntimeError("Unknown flag: %s" % (flagName))

        s.set(keyName, value)
    
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

    def addPairToName(valName, val0, val1):
        """Add a pair of values, val0 and val1, to the valName list"""
        sval1 = str(val1)
        if val0 != val1:
            pre = os.path.commonprefix([str(val0), sval1])
            sval1 = int(sval1[len(pre):])
        valName.append("%s-%s" % (val0, sval1) if val1 != val0 else str(val0))

    valName = []
    val0 = vals[0]; val1 = val0
    for val in vals[1:]:
        if isinstance(val, int) and val == val1 + 1:
            val1 = val
        else:
            addPairToName(valName, val0, val1)
            val0 = val; val1 = val0

    addPairToName(valName, val0, val1)

    return ", ".join(valName)

def getNameOfSRSet(sr, n, omit=[]):
    """Get the name of a set of Sensors or Rafts arranged in an n*n array (with elements of omit missing)"""

    sindex = {}; sname = {}             # generate mapping from Raft/Sensor names to indexes
    k = 0
    for i in range(n):
        for j in range(n):
            s = "%d,%d" % (i, j)
            if s not in omit:
                sindex[s] = k
                sname[str(k)] = s
                k += 1

    for s in ("all",):
        sindex[s] = k
        sname[str(k)] = s
        k += 1

    name = []
    for subsr in getNameOfSet([sindex[s] for s in sr]).split(", "):
        name.append("-".join([sname[i] for i in subsr.split('-')]).replace(",", ""))

    return ", ".join(name)

def makeMapperInfo(mapper):
    """Return an object with extra per-mapper information (e.g. which fields fully specify an exposure)"""

    class MapperInfo(object):
        @staticmethod
        def getColorterm(filterName):
            return None

        def getId(self, src, field="objId"): # can't be static as it calls derived function splitId
            idDict = self.splitId(src.getId(), asDict=True)

            return idDict[field] if field else idDict

        @staticmethod
        def canonicalFiltername(filterName):
            return filterName

        @staticmethod
        def idMask(dataId):
            return 0x0
            
    class LsstSimMapperInfo(MapperInfo):
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
            filters = set()
            sensors = set()
            rafts = set()
            visits = set()
            for dataId in dataIds:
                if dataId.get("sensor") == None:
                    did = dataId.copy(); did["sensor"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **did)).getName())
                    except:
                        filters.add("?")
                    sensors.add("all")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except:
                        filters.add("?")

                    try:
                        sensors.add(dataId["sensor"])
                    except TypeError:
                        for c in dataId["sensor"]:
                            sensors.add(c)

                if dataId.get("raft") == None:
                    did = dataId.copy(); did["raft"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **did)).getName())
                    except:
                        filters.add("?")
                    rafts.add("all")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except:
                        filters.add("?")

                    try:
                        rafts.add(dataId["raft"])
                    except TypeError:
                        for c in dataId["raft"]:
                            rafts.add(c)

                try:
                    visits.add(dataId["visit"])
                except TypeError:
                    for v in dataId["visit"]:
                        visits.add(v)

            sensors = sorted(list(sensors))
            rafts = sorted(list(rafts))
            visits = sorted(list(visits))
            filters = sorted(list(filters))

            if len(visits) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple visits and filters: %s %s" % \
                      (visits, filters)
                visits = visits[0:1]

            title = "%s R%s S%s [%s]" % (getNameOfSet(visits),
                                         getNameOfSRSet(rafts, 5, ['0,0', '4,0', '0,4', '4,4']),
                                         getNameOfSRSet(sensors, 3), ", ".join(filters))
            if rerunName:
                title += " %s" % rerunName

            return title

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
        def splitId(oid, asDict=False):
            """Split an ObjectId into visit, raft, sensor, and objId"""
            objId = int((oid & 0xffff) - 1)     # Should be the same value as was set by apps code
            oid >>= 16
            raftSensorId = oid & 0x1ff
            oid >>= 9
            visit = int(oid)

            raftId, sensorId = int(raftSensorId//10), int(raftSensorId%10)
            raft = "%d,%d" % (raftId//5, raftId%5)
            sensor = "%d,%d" % (sensorId//3, sensorId%3)

            if asDict:
                return dict(visit=visit, raft=raft, sensor=sensor, objId=objId)
            else:
                return visit, raft, sensor, objId

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def splitSdssCcdExposureId(oid, hasFilter=True, asDict=False):
        """Split an ObjectId into run, camcol, [filter], objId

    If hasFilter is True, the ObjectId encodes a filtername
        """
        nbits = 26                      # number of bits reserved for objId
        oid = long(oid)

        omask = 0xffffffffffffffff << nbits
        objId = int(oid & ~omask)     # Should be the same value as was set by apps code
        oid >>= nbits

        field = int(oid % 10000);  oid //= 10000
        camcol = int(oid % 10);    oid //= 10
        filter = int(oid % 10);    oid //= 10
        run = int(oid)

        if hasFilter:
            filterNames = [k for k, v in butler.mapper.filterIdMap.items() if v == filter]
            try:
                filter = filterNames[0]
                assert len(filterNames) == 1
            except IndexError:
                raise RuntimeError("Invalid filter index %d" % filter)

            if asDict:
                return dict(run=run, camcol=camcol, filter=filter, field=field, objId=objId)
            else:
                return run, camcol, filter, field, objId
        else:
            if asDict:
                return dict(run=run, camcol=camcol, field=field, objId=objId)
            else:
                return run, camcol, field, objId

    def splitSdssCoaddId(oid, hasFilter=True, asDict=False):
        """Split an ObjectId into tract, patch, [filter], objId

    If hasFilter is True, the ObjectId encodes a filtername
        """
        nbits = 34                  # number of bits used by patch etc. part of ID
        if hasFilter:
            nbits += 3                  # add 3 bits for filters
        nbits = 64 - nbits          # length
        oid = long(oid)

        omask = 0xffffffffffffffff << nbits
        objId = int(oid & ~omask)     # Should be the same value as was set by apps code
        oid >>= nbits
        if hasFilter:
            filter = int(oid & 0x7)
            oid >>= 3
        patchY = int(oid & 0x1fff)
        oid >>= 13
        patchX = int(oid & 0x1fff)
        oid >>= 13
        tract = int(oid)

        patch = "%d,%d" % (patchX, patchY)

        if hasFilter:
            filterNames = [k for k, v in butler.mapper.filterIdMap.items() if v == filter]
            try:
                filter = filterNames[0]
                assert len(filterNames) == 1
            except IndexError:
                raise RuntimeError("Invalid filter index %d" % filter)

            if asDict:
                return dict(tract=tract, patch=patch, filter=filter, objId=objId)
            else:
                return tract, patch, filter, objId
        else:
            if asDict:
                return dict(tract=tract, patch=patch, objId=objId)
            else:
                return tract, patch, objId


    class SdssMapperInfo(MapperInfo):
        def __init__(self, Mapper):
            SdssMapperInfo.Mapper = Mapper

        @staticmethod
        def getFields(dataType):
            if _prefix_ in ("", "forced",) or dataType in ("coaddTempExp",):
                fields = ["run", "filter", "camcol"]

                if dataType not in ("flat",):
                    fields.append("field")
            elif _prefix_ in ("goodSeeingCoadd",):
                fields = ["patch", "tract", "filter"]
            else:
                raise RuntimeError("I don't know what fields I need to read %s data" % _prefix_)
                pass

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw",)

        @staticmethod
        def photometricTransform(desiredBand, primaryMag, secondaryMag):
            """Return the primary/secondary magnitude transformed into the desiredBand"""
            return SdssMapperInfo._Colorterm.transformMags(desiredBand, primaryMag, secondaryMag)

        @staticmethod
        def dataIdToTitle(dataIds):
            runs = set()
            filters = set()
            camcols = set()
            fields = set()
            for dataId in dataIds:
                if dataId.get("camcol") == None:
                    did = dataId.copy(); did["camcol"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **did)).getName())
                    except:
                        filters.add(dataId.get("filter", "?"))
                    if _prefix_ in ("", "forced"):
                        camcols.add("(all)")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except:
                        filters.add("?")

                    try:
                        camcols.add(dataId["camcol"])
                    except TypeError:
                        for c in dataId["camcol"]:
                            camcols.add(c)

                for k in ["run", "patch", "tract"]:
                    try:
                        fields.add(dataId[k])
                    except KeyError:
                        pass
                    except TypeError:
                        for f in dataId[k]:
                            fields.add(f)

                for k in ["field",]:
                    try:
                        runs.add(dataId[k])
                    except KeyError:
                        pass
                    except TypeError:
                        for f in dataId["field"]:
                            runs.add(f)

            runs = sorted(list(runs))
            fields = sorted(list(fields))
            camcols = sorted(list(camcols))
            filters = sorted(list(filters))

            if len(runs) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple runs and filters: %s %s" % \
                      (runs, filters)
                runs = runs[0:1]

            nameOfFilters = "".join(filters)
            if len(filters) > 1:
                nameOfFilters = "[%s]" % nameOfFilters
            title = "%s %s%s %s" % (getNameOfSet(runs), nameOfFilters, getNameOfSet(camcols),
                                      getNameOfSet(fields))
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

        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            return butler

        @staticmethod
        def splitId(oid, hasFilter=True, asDict=False):
            """Split an ObjectId into run, camcol, [filter], field, objId or tract, patch, [filter], objId

        If hasFilter is True, the ObjectId encodes a filtername
            """

            if _prefix_ in ("goodSeeingCoadd",):
                return splitSdssCoaddId(oid, hasFilter=hasFilter, asDict=asDict)
            else:
                return splitSdssCcdExposureId(oid, hasFilter=hasFilter, asDict=asDict)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    from lsst.meas.photocal.colorterms import Colorterm
    from lsst.obs.suprimecam.colorterms import colortermsData

    class SubaruMapperInfo(MapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.Mapper = Mapper

            SubaruMapperInfo._Colorterm = Colorterm
            SubaruMapperInfo.getColorterm = lambda x, y : Colorterm.getColorterm(y)
            SubaruMapperInfo._Colorterm.setColorterms(colortermsData, "Hamamatsu")

        @staticmethod
        def getFields(dataType):
            if _prefix_ in ("",):
                fields = ["visit", "ccd"]
                if dataType not in ("flat",):
                    fields.append("filter")
            elif _prefix_ in ("stack",):
                fields = ["stack", "patch", "filter"]
            else:
                raise RuntimeError("I don't know what fields I need to read %s data" % _prefix_)

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw",)

        @staticmethod
        def photometricTransform(desiredBand, primaryMag, secondaryMag):
            """Return the primary/secondary magnitude transformed into the desiredBand"""
            return SubaruMapperInfo._Colorterm.transformMags(desiredBand, primaryMag, secondaryMag)

        @staticmethod
        def dataIdToTitle(dataIds):
            if _prefix_ == "stack":
                title = []
                for did in dataIds:
                    title.append("stack %(stack)d patch %(patch)d filter %(filter)s" % did)

                return "[%s]" % "], [".join(title)

            filters = set()
            ccds = set()
            visits = set()
            for dataId in dataIds:
                if dataId.get("ccd") == None:
                    did = dataId.copy(); did["ccd"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **did)).getName())
                    except:
                        filters.add(dataId.get("filter", "?"))
                    ccds.add("(all)")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except:
                        filters.add("?")

                    try:
                        ccds.add(dataId["ccd"])
                    except TypeError:
                        for c in dataId["ccd"]:
                            ccds.add(c)
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
        def idMask(dataId):
            return dataId.get("stack", 0) << 52 # hack hack

        @staticmethod
        def splitId(oid, asDict=False):
            """Split an ObjectId into visit, ccd, and objId"""
            oid = long(oid)
            objId = int(oid & 0xffff)     # Should be the same value as was set by apps code
            oid >>= 32

            if _prefix_ == "stack":
                stack = int(oid >> 20)      # XXX see idMask
                oid &= 0xfffff

                nfilter = len(butler.mapper.filters)
                ifilter = oid % nfilter
                oid //= nfilter
                patch = int(oid)

                filter = [k for k,v  in butler.mapper.filterIdMap.items() if v == ifilter][0]

                if asDict:
                    return dict(stack=stack, patch=patch, filter=filter, objId=objId)
                else:
                    return stack, patch, filter, objId
                
            else:                
                ccd = int(oid % 10)
                oid //= 10
                visit = int(oid)

                if asDict:
                    return dict(visit=visit, ccd=ccd, objId=objId)
                else:
                    return visit, ccd, objId

        @staticmethod
        def canonicalFiltername(filterName):
            mat = re.search(r"W-S-(.)\+", filterName)
            if mat:
                filterName = mat.group(1).lower()

            return filterName

    class SubaruMapperInfoMit(SubaruMapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.__init__(self, None)
            SubaruMapperInfoMit.Mapper = Mapper
            SubaruMapperInfo._Colorterm.setColorterms(colortermsData, "MIT")

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if isinstance(mapper, LsstSimMapper):
        return LsstSimMapperInfo(LsstSimMapper)
    elif isinstance(butler.mapper, SdssMapper):
        return SdssMapperInfo(SdssMapper)
    elif isinstance(butler.mapper, SuprimecamMapper):
        return SubaruMapperInfo(SuprimecamMapper)
    elif isinstance(butler.mapper, SuprimecamMapperMit):
        return SubaruMapperInfoMit(SuprimecamMapperMit)
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
        if kwargs.has_key("cat"):
            self.cat = kwargs["cat"]
            return

        self.setButler(*args, **kwargs)

        self.dataId = []

        self.matched = None
        self.ZP0 = 31
        self.zp = None

        self.name = "??"
        self.astrom = None
            
    def setButler(self, rerun=None, dataRoot=None, registryRoot=None, Mapper=None):
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
        else:
            dataRoot = os.path.join(dataRoot, "rerun", rerun)        

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
                print >> sys.stderr, "I'm unable to find your registry in %s" % registryRoot

        if not Mapper:
            Mapper = getMapper(registryRoot, defaultMapper=LsstSimMapper if False else SdssMapper)
        
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
        for st, v, f, r, s in butler.queryMetadata(_raw_, "skyTile",
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

        dataSets = []
        for did in _dataIdDictOuterProduct(dataId, []):
            if _prefix_ in ("forced", "goodSeeingCoadd", "stack"): # XXXXXXX
                if butler.datasetExists(dtName(dataType), **did):
                    dataSets.append(did)
            else:
                for vals in butler.queryMetadata(_raw_, _visit_, fields, **did):
                    _dataId = dict(zip(fields, vals))
            
                    if butler.datasetExists(dtName(dataType), **_dataId):
                        dataSets.append(_dataId)

        self.dataSets = dataSets

        return self.dataSets

    def getDataset(self, dataType, dataId, ids=True, calibrate=False, setMask=None, fixOrientation=None,
                   fixAmpLevels=False, trim=True):
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

        if trim and dataType in butler.mapperInfo.getTrimmableData():
            return [butler.mapperInfo.assembleCcd(dtName(dataType), inButler,
                                                  dataId, fixAmpLevels=fixAmpLevels)]
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

        dataSets = self.lookupDataByVisit(dataType, dataId)
        self.dataId += dataSets

        data = []
        for did in dataSets:
            dataElem = butler.get(dtName(dataType), **did)

            if dataType == "calexp":
                psf = butler.get(dtName("psf"), **did)
		psf.setDetector(dataElem.getDetector())
                dataElem.setPsf(psf)
                del psf
            elif dataType == "psf":
                calexp = butler.get(dtName("calexp"), **did)
                dataElem.setDetector(calexp.getDetector())
		del calexp
            elif dataType == 'src':
                idMask = butler.mapperInfo.idMask(did)
                if idMask:
                    for s in dataElem:
                        s.set("id", s.getId() | idMask)
                        if s.getParent():
                            s.set("parent", s.getParent() | idMask)
                
            data.append(dataElem)

        return data

    def getSources(self, dataId, setXYfromRaDec=False):
        """Return the list of Sources for the specified dataId
        
        If setXYfromRaDec is true, recalculate (x, y) from (ra, dec)"""

        if setXYfromRaDec:
            calexp_md = butler.get(dtName("calexp", True), **dataId)
            wcs = afwImage.makeWcs(calexp_md)

        sources = []
        dataSets = self.dataSets
        for ss in self.getDataset("src", dataId):
            for s in ss:
                s.setId(butler.mapperInfo.getId(s))
                if setXYfromRaDec:
                    x, y = wcs.skyToPixel(s.getCoord())
                    s.setF("centroid.sdss.x", x)
                    s.setF("centroid.sdss.y", y)

            sources.append(ss)

        if dataSets:
            self.dataSets = dataSets

        return sources

    def getMagsByType(self, magType, good=None, suffix="", magDict=None):
        if magDict:
            mags = magDict[magType]
        else:
            try:
                cat = self.cat
            except:
                cat = self
            if magType == "ap":
                mags = cat.get("apMag" + suffix)
            elif magType == "inst":
                mags = cat.get("instMag" + suffix)
            elif magType == "model":
                mags = cat.get("modelMag" + suffix)
            elif magType == "psf":
                mags = cat.get("psfMag" + suffix)
            else:
                raise RuntimeError("Unknown magnitude type %s" % magType)

        if good is not None:
            mags = mags[good]

        return mags

    def expandDataId(self, dataId):
        """Return a list of fully qualified dataIds, given a possibly-ambiguous one

ccd may be a list"""

        return self.lookupDataByVisit("src", dataId)

    def getMagsByVisit(self, dataId, extraApFlux=0.0, verbose=False):
        """Read the magnitudes for a set of data"""

        dataIds = self.expandDataId(dataId)
        
        catInfo = None
        for did in dataIds:
            if verbose:
                print "Reading %s         \r" % did,
                sys.stdout.flush()

            catInfo = _appendToCatalog(self, did, catInfo, extraApFlux=extraApFlux)
        if verbose:
            print

        cat = catInfo[0] if catInfo else None
            
        if cat is None:
            raise RuntimeError("Failed to read any data for %s" % " ".join([str(d) for d in dataId]))

        try:
            cat.getColumnView()
        except:
            cat = makeContiguous(cat)
        self.cat = cat
        
        self.name = butler.mapperInfo.dataIdToTitle(dataIds)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, dataId, allSources=True,
                        displayType="psf", frame=None, mtv=False, verbose=False):
        """If frame is not None then information about displayType will be plotted on ds9; if
        mtv is true the calexp will be displayed first

        If allSources is True, read the src table and match it to the catalogue
        """
        dataIds = self.expandDataId(dataId)
        if not dataIds:
            raise RuntimeError("%s doesn't expand to any valid dataIds" % dataId)

        self.name = butler.mapperInfo.dataIdToTitle(dataIds)

        if mtv and frame is None:
            frame = 0
        if frame is not None:
            if mtv:
                calexp = self.getDataset("calexp", dataIds[0])[0]
                Id = calexp.getDetector().getId() if calexp.getDetector() else self.name
                ds9.mtv(calexp, title=Id, frame=frame)
            else:
                ds9.erase(frame)
            
        _matched = []
        _zp = []
        for dataId in dataIds:
            if verbose:
                print "Reading %s         \r" % dataId,
                sys.stdout.flush()
                
            matched, zp = self._getCalibObjectsImpl(dataId, displayType, frame, verbose, allSources=allSources)
            
            _matched.append(matched)
            _zp.append(zp)
        if verbose:
            print

        self.matched = sum([len(m) for m in _matched])*[None]
        self.zp = np.empty(len(self.matched))

        start = 0
        for i, m in enumerate(_matched):
            end = start + len(m)
            self.matched[start:end] = m
            self.zp[start:end] = _zp[i]*np.ones(end - start)
            start = end

    def _getCalibObjectsImpl(self, dataId, displayType, frame, verbose, allSources=False):
        calexp_md = butler.get(dtName("calexp", True), **dataId)
        wcs = afwImage.makeWcs(calexp_md)
        imageSize = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
        xy0 = (calexp_md.get("CRVAL1A"), calexp_md.get("CRVAL2A"),)
        filterName = afwImage.Filter(calexp_md).getName()
        calib = afwImage.Calib(calexp_md)

        if not self.astrom:
            self.astrom = measAstrom.Astrometry(measAstrom.Astrometry.ConfigClass())

        if allSources:
            sources = butler.get(dtName("src"), **dataId)
            cat = self.astrom.getReferenceSourcesForWcs(wcs, imageSize, filterName, pixelMargin=50,
                                                        trim=True, allFluxes=True)
            if frame is not None and displayType == "src":
                showSourceSet(sources, xy0=xy0, raDec=False, frame=frame)
                showSourceSet(cat, xy0=xy0, wcs=wcs, raDec=True,  frame=frame, symb="o", ctype=ds9.RED)
            try:
                matched = self.astrom._getMatchList(sources, cat, wcs)
            except Exception, e:
                print "RHL", e

                matchRadius = 2
                x, y = sources.getX(), sources.getY()
                x0, y0 = xy0
                xk = sources.getSchema().find("%s.x" % sources.getTable().getCentroidDefinition()).getKey()
                yk = sources.getSchema().find("%s.y" % sources.getTable().getCentroidDefinition()).getKey()
                for i, s in enumerate(sources):
                    s.setD(xk, x[i] - x0)
                    s.setD(yk, y[i] - y0)

                matched = afwTable.matchRaDec(cat, sources, matchRadius*afwGeom.arcseconds)
        else:
            sources = butler.get(dtName("icSrc"), **dataId)
            matches = butler.get(dtName("icMatch"), **dataId)

            matched = self.astrom.joinMatchListWithCatalog(matches, sources,
                                                           allFluxes=True) # list(ReferenceMatch)

        try:
            zp = calib.getMagnitude(1.0)
        except Exception, e:
            print >> sys.stderr, "Failed to calculate zp: %s" % e
            zp = 0.0

        showDistortion = False
        if displayType == "distortion":
            showDistortion = False
        elif displayType in ("src",):
            pass
        elif displayType in ("psf", "ap"):
            if displayType == "ap":
                displayType = "sinc"
            if matched:
                ref, src, d = matched[0]
                sch = src.getSchema()
                fluxKey = sch.find("flux.%s" % displayType).getKey()
                #
                # Now the reference magnitude, which may need a colour term
                #
                sch = ref.getSchema()

                ct = butler.mapperInfo.getColorterm(filterName)
                if ct:
                    primary, secondary = ct.primary, ct.secondary
                    primaryKey_r = sch.find(primary).getKey()
                    secondaryKey_r = sch.find(secondary).getKey()
                else:
                    fluxKey_r = sch.find("flux").getKey()

        else:
            print >> sys.stderr, "Ignoring unknown displayType %s" % displayType
            frame = None
            
        if False:
            return matched, zp

        if showDistortion:
            calexp = butler.get('calexp', **dataId)
            ccd = cameraGeom.cast_Ccd(calexp.getDetector())
            distortion = ccd.getDistortion()

            import lsst.afw.geom.ellipses as geomEllipses
            quad = geomEllipses.Quadrupole() # used for estimating distortion
            quad.scale(1/math.sqrt(quad.getArea()))

        keepMatched = []
        deltaScale = 1000               # multiplier for mag_ref - mag_src to display on ds9
        with ds9.Buffering():
            for match in matched:
                ref, src, d = match
                x, y = wcs.skyToPixel(ref.getRa(), ref.getDec())
                if x < 0 or y < 0 or x > imageSize[0] or y > imageSize[1]:
                    continue
                
                keepMatched.append(match)

                if frame is not None:
                    if False:
                        ds9.dot("o", x, y, frame=frame, ctype=ds9.GREEN)
                        ds9.dot("+", src.getX(), src.getY(), frame=frame, ctype=ds9.RED)
                    else:
                        if showDistortion:
                            delta = distortion.distort(afwGeom.PointD(x, y), quad, ccd).getArea() - 1
                        else:
                            if ct:
                                p = -2.5*math.log10(ref.get(primaryKey_r))
                                s = -2.5*math.log10(ref.get(secondaryKey_r))

                                refmag = butler.mapperInfo.photometricTransform(filterName, p, s)
                            else:
                                refmag = -2.5*math.log10(ref.get(fluxKey_r))

                            delta = refmag - (zp - 2.5*math.log10(src.get(fluxKey)))

                        size = min([100, max([deltaScale*delta, -100])])
                        ds9.dot("o", src.getX(), src.getY(), size=abs(size),
                                frame=frame, ctype=ds9.CYAN if size > 0 else ds9.MAGENTA)
                        ds9.dot("o", src.getX(), src.getY(), size=deltaScale*0.01, # 10 mmag
                                frame=frame, ctype=ds9.YELLOW)
           
        if verbose > 1:
            print "Kept %d out of %d reference sources for %s" % (len(keepMatched), len(matched), dataId)

        return keepMatched, zp

def _dataIdDictOuterProduct(dataId, expandedDataId=[]):
    """Given a dataId which may contain lists, return a list of dataIds with all lists expanded"""

    for k, v in dataId.items():
        try:
            v[0]                        # can we index it?
            if not isinstance(v, str):  # and it's not a string (damn Guido)
                for x in v:
                    _dataId = dataId.copy(); _dataId[k] = x
                    _dataIdDictOuterProduct(_dataId, expandedDataId)
                return expandedDataId
        except (TypeError, IndexError):
            pass

    expandedDataId.append(dataId)

    return expandedDataId
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getFilters(dataRoot, dataId, filters, verbose=False):
    """Read the data for a number of filters"""

    data = {}
    for f in filters:
        data[f] = Data(dataRoot=dataRoot)
        data[f].getMagsByVisit(dataId.copy(filter=f), verbose=verbose)

    return data

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Convert SourceSets to numpy arrays
#
def _setFlagsFromSource(ssc=None):
    """Return a dictionary of status flags from a Source or table;  if ssc is None, return a
    dictionary of the Source fields used"""
    
    names = dict(
        EDGE =          "flags.pixel.edge",
        INTERP =        "flags.pixel.interpolated.any",
        INTERP_CENTER = "flags.pixel.interpolated.center",
        SATUR =         "flags.pixel.saturated.any",
        SATUR_CENTER =  "flags.pixel.saturated.center",
        )

    if ssc is None:
        return names
    else:
        flags = {}
        for k, n in names.items():
            flags[k] = ssc.get(n)
        return flags

def _appendToCatalog(data, dataId, catInfo=None, scm=None, sourceSet=None, extraApFlux=0.0):
    """Append interesting columns from data's src object (a SourceSet) to a catalogue cat, passed
    in through catInfo == (cat, scm) where scm is the SchemaMapper used to generate cat.

    if data is None, sourceSet must be provided and is assumed to be the sourceSet

    If catInfo is None, a new catalogue (and SchemaMapper) are generated
    """

    if data is None:
        assert sourceSet
    else:
        sourceSet = data.getDataset("src", dataId)[0]
    
    calexp_md = butler.get(dtName("calexp", True), **dataId)
    calib = afwImage.Calib(calexp_md)
    fluxMag0, fluxMag0Err = calib.getFluxMag0()
    if fluxMag0 <= 0.0:
        fluxMag0 = 1e13
        print >> sys.stderr, "Setting fluxMag0 to %g" % fluxMag0
        calib.setFluxMag0(fluxMag0, fluxMag0Err)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if catInfo is not None:
        cat, scm = catInfo
        assert scm
        
        table = cat.table
    else:
        tab = sourceSet.getTable()
        sch = tab.getSchema()
        scm = afwTable.SchemaMapper(sch)
        for f in ["id", "coord", "parent",  # must be first and in this order --- #2154
                  "deblend.nchild",
                  "classification.extendedness",
                  "flags.pixel.edge",
                  "flags.pixel.interpolated.any",
                  "flags.pixel.interpolated.center",
                  "flags.pixel.saturated.any",
                  "flags.pixel.saturated.center",
                  ]:
            scm.addMapping(sch.find(f).getKey())
        for key in [tab.getCentroidKey(), tab.getCentroidErrKey(), tab.getCentroidFlagKey(),
                    tab.getShapeKey(),
                    ]:
            scm.addMapping(key)

        definePhotoSlots = False            # add standard photometric slots?
        if definePhotoSlots:
            for key in [tab.getInstFluxKey(),  tab.getInstFluxErrKey(),  tab.getInstFluxFlagKey(),
                        tab.getModelFluxKey(), tab.getModelFluxErrKey(), tab.getModelFluxFlagKey(),
                        tab.getPsfFluxKey(),   tab.getPsfFluxErrKey(),   tab.getPsfFluxFlagKey(),
                        tab.getApFluxKey(),   tab.getApFluxErrKey(),     tab.getApFluxFlagKey(),
                        ]:
                scm.addMapping(key)

        scm.addOutputField(afwTable.Field["Flag"]("stellar", "Is the source stellar?"))
        scm.addOutputField(afwTable.Field["F"]("apMag", "The magnitude corresponding to apFlux"))
        scm.addOutputField(afwTable.Field["F"]("instMag", "The magnitude corresponding to instFlux"))
        scm.addOutputField(afwTable.Field["F"]("modelMag", "The magnitude corresponding to modelFlux"))
        scm.addOutputField(afwTable.Field["F"]("psfMag", "The magnitude corresponding to psfFlux"))
        scm.addOutputField(afwTable.Field["Flag"]("good", "The object is good"))

        cat = afwTable.SourceCatalog(scm.getOutputSchema())
        table = cat.table

        table.defineCentroid(tab.getCentroidDefinition())

        if definePhotoSlots:
            #
            # Define the slots we were using before
            #
            table.defineApFlux(tab.getApFluxDefinition())    
            table.defineInstFlux(tab.getInstFluxDefinition())    
            table.defineModelFlux(tab.getModelFluxDefinition())    
            table.definePsfFlux(tab.getPsfFluxDefinition())    
    #
    # OK, we have the schema sorted out;  time to read and process
    #
    table.preallocate(len(cat) + len(sourceSet))

    stellarKey = cat.getSchema().find("stellar").getKey()
    apMagKey = cat.getSchema().find("apMag").getKey()
    instMagKey = cat.getSchema().find("instMag").getKey()
    modelMagKey = cat.getSchema().find("modelMag").getKey()
    psfMagKey = cat.getSchema().find("psfMag").getKey()
    goodKey = cat.getSchema().find("good").getKey()

    try:
        throwOnNegative = afwImage.Calib.getThrowOnNegativeFlux()
        afwImage.Calib.setThrowOnNegativeFlux(False)
    except:
        throwOnNegative = None

    try:
        afwTable.SourceRecord.setF
    except AttributeError:
        afwTable.SourceRecord.setD = afwTable.SourceRecord.setF8
        afwTable.SourceRecord.setF = afwTable.SourceRecord.setF4
        afwTable.SourceRecord.setFlag = afwTable.SourceRecord.set_Flag

    try:
        afwTable.SourceRecord.setMag
    except AttributeError:
        typeString = apMagKey.getTypeString()
        print >> sys.stderr, "Checking type of magnitude keys... %s" % typeString
        afwTable.SourceRecord.setMag = \
            afwTable.SourceRecord.setD if typeString == "D" else afwTable.SourceRecord.setF

    try:
        for s in sourceSet:
            cat.append(cat.copyRecord(s, scm))

            cat[-1].setId(s.getId())

            try:
                cat[-1].setMag(apMagKey,    calib.getMagnitude(s.getApFlux() + extraApFlux))
            except Exception, e:
                cat[-1].setMag(apMagKey,  calib.getMagnitude(float("NaN")))
                
            try:
                cat[-1].setMag(instMagKey,  calib.getMagnitude(s.getInstFlux()))
            except Exception, e:
                cat[-1].setMag(instMagKey,  calib.getMagnitude(float("NaN")))

            try:
                cat[-1].setMag(modelMagKey, calib.getMagnitude(s.getModelFlux()))
            except Exception, e:
                cat[-1].setMag(modelMagKey, calib.getMagnitude(float("NaN")))
            
            try:
                cat[-1].setMag(psfMagKey, calib.getMagnitude(s.getPsfFlux()))
            except Exception, e:
                cat[-1].setMag(psfMagKey, calib.getMagnitude(float("NaN")))

            cat[-1].setFlag(goodKey, True)  # for now

            isStar = (s.get("classification.extendedness") < 0.5)
            cat[-1].setFlag(stellarKey, isStar)
    except:
        raise
    finally:
        if throwOnNegative is not None:
            afwImage.Calib.setThrowOnNegativeFlux(throwOnNegative)

    return cat, scm

def zipMatchList(matchList):
    """zip a matchList into a single catalogue

    @param matchList MatchList, as returned by e.g. afwTable.matchRaDec
    """

    records = [matchList[0][i] for i in range(2)]

    requiredFields = ["id", "coord", "parent"] # keys that must be first and in this order --- #2154

    scm = None                          # the SchemaMapper
    keys_2 = []                         # keys from second list that need to be copied over
    for i, rec in enumerate(records):
        suffix = "_%d" % (i + 1)
        
        sch = rec.getSchema()
        if not scm:
            scm = afwTable.SchemaMapper(sch)

            for f in requiredFields: # must be first and in this order --- #2154
                scm.addMapping(sch.find(f).getKey())

        for schEl in [sch.find(n) for n in sch.getNames()]: # just "sch" should work, but there's "getBit" bug...
            inField = schEl.getField()
            if inField.getName() not in requiredFields:
                key = schEl.getKey()
                outputField = type(inField)(inField.getName() + suffix, inField.getDoc())
                if i == 0:
                    scm.addMapping(key, outputField)
                else:
                    keys_2.append((key, scm.addOutputField(outputField),
                                   outputField.getTypeString()))
    #
    # OK, we have the schema sorted out;  time to read and process
    #
    cat = afwTable.SourceCatalog(scm.getOutputSchema())
    cat.table.preallocate(len(matchList))

    for m1, m2, d in matchList:
        cat.append(cat.copyRecord(m1, scm))
        rec = cat[-1]
        for key, okey, typeString in keys_2:
            if typeString == "D":
                rec.setD(okey, m2.getD(key))
            elif typeString == "F":
                rec.setF(okey, m2.getF(key))
            elif typeString == "Flag":
                rec.setFlag(okey, m2.getFlag(key))
            else:
                rec.set(okey, m2.get(key))

    return cat

def getMagsFromSS(ss, dataId, extraApFlux=0.0):
    """Return numpy arrays constructed from SourceSet ss"""

    cat = _appendToCatalog(None, dataId, sourceSet=None, extraApFlux=extraApFlux)

    ids = cat.get("id")
    flags = _setFlagsFromSource(cat)
    stellar = cat.get("classification.extendedness") < 0.5
    apMags = cat.get(apMagKey)
    instMags = cat.get(instMagKey)
    modelMags = cat.get(modelMagKey)
    psfMags = cat.get(psfMagKey)
    x = cat.getX()
    y = cat.getY()
    if False:
        shape = [cat.getIxx(), cat.getIxy(), cat.getIyy()]
    else:
        shape = [np.zeros(len(ids)), np.zeros(len(ids)), np.zeros(len(ids))]

    return ids, flags, stellar, x, y, shape, \
        dict(ap=apMags, inst=instMags, model=modelMags, psf=psfMags)

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
             showMedians=False, sgVal=0.05, meanDelta=0.0, adjustMean=False, parents=False,
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
        data.cat
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)

    if parents:
        good = np.logical_and(good, data.cat.get("parent") == 0)
    else:
        good = np.logical_and(good, data.cat.get("deblend.nchild") == 0)

    mag1 = data.getMagsByType(magType1, good)
    mag2 = data.getMagsByType(magType2, good)
    delta = np.array(mag1 - mag2, dtype='float64') # float64 is needed by makeStatistics --- why?

    locus = np.logical_and(np.logical_and(mag1 > magmin, mag1 < maglim), np.abs(delta - meanDelta) < sgVal)

    stellar = data.cat.get("stellar")
    nonStellar = np.logical_not(stellar)
    stellar = stellar[good]; nonStellar = nonStellar[good]
    ids = data.cat.get("id")[good]

    try:
        stats = afwMath.makeStatistics(delta[locus], afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
    except:
        mean, stdev = np.nan, np.nan

    if adjustMean:
        raise RuntimeError("Fix me")
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
    did = butler.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data.getDataset("calexp_md", did)[0]

    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat
    x = data.cat.getX()[good] + md.get("LTV1")
    y = data.cat.getY()[good] + md.get("LTV2")
    ids = data.cat.get("id")[good]
    eventHandlers[fig] = EventHandler(axes, mag1, delta, ids, x, y, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeSelectCcd(ccd):
    def selectCcd(id):
        return butler.mapperInfo.splitId(id, asDict=True)["ccd"] == ccd

    return lambda id: butler.mapperInfo.splitId(id, asDict=True)["ccd"] == ccd

def plotCM(data1, data2, magType="psf", maglim=20, magmin=14,
           SG="sg", selectObjId=None, showMedians=False, sgVal=0.05, meanDelta=0.0, adjustMean=False,
           xmin=None, xmax=None, ymin=None, ymax=None,
           title="+", markersize=1, color="red", alpha=1.0, frames=[0], fig=None):
    """Plot (data1.magType - data2.magType) v. data1.magType mags (e.g. "psf")

The magnitude limit for the "locus" box used to calculate statistics is maglim, and only objects within
+- sgLim of meanDelta are included in the locus (the axes are also adjusted to include meanDelta in the
visible area).  If adjustMean is True, adjust all the CCDs being plotted to have the same mean within
the "locus" area.

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits (y{min,max} are interpreted relative to meanDelta
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    data = {1 : data1, 2 : data2}
    for i, d in data.items():
        try:
            d.cat
        except AttributeError:
            raise RuntimeError("Please call data%d.getMagsByVisit, then try again" % i)

    matchRadius = 2                     # arcsec
    mat = afwTable.matchRaDec(data[1].cat, data[2].cat, matchRadius*afwGeom.arcseconds)
    matched = Data(cat=zipMatchList(mat))

    bad = None
    for i in range(2):
        suffix = "_%d" % (i + 1)
        for name in ["flags.pixel.edge",
                     #"flags.pixel.interpolated.center",
                     "flags.pixel.saturated.center",]:
            _flg = matched.cat.get(name + suffix)
            if bad is None:
                bad = _flg
            else:
                bad = np.logical_or(bad, _flg)

    good = np.logical_not(bad)

    ids = matched.cat.get("id")
    if selectObjId:
        for i, _id in enumerate(ids):
            if not selectObjId(_id):
                good[i] = False

    ids = ids[good]
    xc = matched.cat.get("centroid.sdss_1.x")[good]
    yc = matched.cat.get("centroid.sdss_1.y")[good]
    stellar = matched.cat.get("stellar_1")[good]

    mag1 = matched.getMagsByType(magType, good, suffix="_1")
    mag2 = matched.getMagsByType(magType, good, suffix="_2")
    delta = mag1 - mag2

    good = good[good]

    nonStellar = np.logical_not(stellar)
    
    try:
        alpha.keys()
    except AttributeError:
        alpha = dict(g=alpha, s=alpha)

    color2 = "green"
    nobj = 0
    if "g" in SG.lower():
        nobj += np.sum(nonStellar)
        axes.plot(delta[nonStellar], mag1[nonStellar], "o",
                  alpha=alpha["g"], markersize=markersize, markeredgewidth=0, color=color)
    if "s" in SG.lower():
        nobj += np.sum(stellar)
        axes.plot(delta[stellar], mag1[stellar], "o",
                  alpha=alpha["s"], markersize=markersize, markeredgewidth=0, color=color2)

    #axes.plot((0, 30), (0, 0), "b-")

    filterNames = [None] + [data[i + 1].dataId[0]["filter"] for i in range(len(data.keys()))]
    filter1, filter2 = filterNames[1], filterNames[2]

    axes.set_xlim(-1 if xmin is None else xmin, 2 if xmax is None else xmax)
    axes.set_ylim(24 if ymax is None else ymax, 14 if ymin is None else ymin)
    axes.set_xlabel("(%s - %s)$_{%s}$" % (filter1, filter2, magType))
    axes.set_ylabel("%s$_{%s}$" % (filter1, magType))

    title += " %d objects" % nobj
    axes.set_title(re.sub(r"^\+\s*", data1.name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat

    did = butler.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data[1].getDataset("calexp_md", did)[0]
    xc = xc[good] + md.get("LTV1")
    yc = yc[good] + md.get("LTV2")
    ids = ids[good]

    if "g" not in SG.lower():
        delta[nonStellar] = -1000
    if "s" not in SG.lower():
        delta[stellar] = -1000

    eventHandlers[fig] = EventHandler(axes, delta, mag1, ids, xc, yc, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotCC(data1, data2, data3, magType="psf", magmax=None, magmin=None,
           SG="sg", selectObjId=None, matchRadius=2, colorCcds=False,
           usePrincipalColor=True, stellarLocusEnds=[], adjustLocus=False, locusLtype="b:", 
           xmin=None, xmax=None, ymin=None, ymax=None,
           title="+", markersize=1, alpha=1.0, color="red", frames=[0], fig=None):
    """Plot (data1.magType - data2.magType) v. (data2.magType - data3.magType) mags (e.g. "psf")

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    data = {1 : data1, 2 : data2, 3 : data3}
    for i, d in data.items():
        try:
            d.cat
        except AttributeError:
            raise RuntimeError("Please call data%d.getMagsByVisit, then try again" % i)

    mat = afwTable.matchRaDec(data[1].cat, data[2].cat, matchRadius*afwGeom.arcseconds)
    mat = afwTable.matchRaDec(zipMatchList(mat), data[3].cat, matchRadius*afwGeom.arcseconds)
    matched = Data(cat=zipMatchList(mat))

    good = None
    for suffix in ["_1_1", "_2_1", "_2"]:
        for name in ["flags.pixel.edge",
                     #"flags.pixel.interpolated.center",
                     "flags.pixel.saturated.center",]:
            _flg = np.logical_not(matched.cat.get(name + suffix))
            if good is None:
                good = _flg
            else:
                good = np.logical_and(good, _flg)

    ids = matched.cat.get("id")
    if selectObjId:
        for i, _id in enumerate(ids):
            if not selectObjId(_id):
                good[i] = False

    ids = ids[good]
    xc = matched.cat.get("centroid.sdss_1_1.x")[good]
    yc = matched.cat.get("centroid.sdss_1_1.y")[good]
    stellar = matched.cat.get("stellar_1_1")[good]

    mag1 = matched.getMagsByType(magType, good, suffix="_1_1")
    mag2 = matched.getMagsByType(magType, good, suffix="_2_1")
    mag3 = matched.getMagsByType(magType, good, suffix="_2")

    good = good[good]

    if magmin is not None:
        good = np.logical_and(good, mag1 > magmin)
    if magmax is not None:
        good = np.logical_and(good, mag1 < magmax)
    
    col12 = mag1 - mag2; col12 = col12[good]
    col23 = mag2 - mag3; col23 = col23[good]

    stellar = stellar[good]
    nonStellar = np.logical_not(stellar)
    
    try:
        alpha.keys()
    except AttributeError:
        alpha = dict(g=alpha, s=alpha)

    filterNames = [None] + [butler.mapperInfo.canonicalFiltername(data[i + 1].dataId[0]["filter"])
                            for i in range(len(data.keys()))]

    color2 = "green"

    nobj = 0
    if "g" in SG.lower():
        nobj += np.sum(nonStellar)

        axes.plot(col12[nonStellar], col23[nonStellar], "o",
                  alpha=alpha["g"], markersize=markersize, markeredgewidth=0, color=color)
    if "s" in SG.lower():
        nobj += np.sum(stellar)

        if colorCcds:
            colors = "rgbcmyk"
            
            xvec = col12[stellar]; yvec = col23[stellar]; 

            for ccd in range(10):
                print ccd, colors[ccd%len(colors)]

            for i, _id in enumerate(ids[good][stellar]):
                ccd = butler.mapperInfo.splitId(_id, asDict=True)["ccd"]
                #if ccd not in (3, 9): continue
                axes.plot(xvec[i], yvec[i], "o", alpha=alpha["g"], markersize=markersize, markeredgewidth=0,
                          color=colors[ccd%len(colors)])

        else:
            axes.plot(col12[stellar], col23[stellar], "o",
                      alpha=alpha["s"], markersize=markersize, markeredgewidth=0, color=color2)
    #axes.plot((0, 30), (0, 0), "b-")
    #
    # Calculate width of blue end of stellar locus
    #
    #stellarLocusEnds = (0.40, 1.00,)    # the blue and red colour defining the straight part of the locus
    #stellarLocusEnds = (0.50, 1.50,)
    if stellarLocusEnds and "s" in SG.lower():
        x, y = col12[stellar], col23[stellar]
        blue = np.logical_and(x > stellarLocusEnds[0], x < stellarLocusEnds[1])
        
        principalColor = ""
        if usePrincipalColor and filterNames[1] == 'g' and filterNames[2] == 'r' and filterNames[3] == 'i':
            pc = -0.227*mag1 + 0.792*mag2 - 0.567*mag3 + 0.050
            principalColor = "w"
            delta = pc[good][stellar]
        elif usePrincipalColor and filterNames[1] == 'r' and filterNames[2] == 'i' and filterNames[3] == 'z':
            pc = -0.270*mag1 + 0.800*mag2 - 0.534*mag3 + 0.054
            principalColor = "y"
            delta = pc[good][stellar]
        else:
            xy = []
            for xx in stellarLocusEnds:
                delta = col23[stellar][abs(col12[stellar] - xx) < 0.1]
                yy = afwMath.makeStatistics(np.array(delta, dtype="float64"), afwMath.MEDIAN).getValue()
                xy.append((xx, yy))

            theta = math.atan2(xy[1][1] - xy[0][1], xy[1][0] - xy[0][0])
            c, s= math.cos(theta), math.sin(theta)
            xp =   x*c + y*s
            yp = - x*s + y*c

            if locusLtype:
                axes.plot([xy[0][0], xy[1][0]], [xy[0][1], xy[1][1]], locusLtype)

            delta = yp

        delta = np.array(delta, dtype="float64")
        try:
            stats = afwMath.makeStatistics(delta[blue], afwMath.STDEVCLIP | afwMath.MEANCLIP)
            mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
        except:
            mean, stdev = float("NaN"), float("NaN")

        #print "%g +- %g" % (mean, stdev)
        fig.text(0.75, 0.85, r"$%s \pm %.3f$" % \
                     ("%s = " % principalColor if principalColor else "", stdev), fontsize="larger")
        
    axes.set_xlim(-1 if xmin is None else xmin, 2 if xmax is None else xmax)
    axes.set_ylim(-1 if ymin is None else ymin, 2 if ymax is None else ymax)

    axes.set_xlabel("(%s - %s)$_{%s}$" % (filterNames[1], filterNames[2], magType))
    axes.set_ylabel("(%s - %s)$_{%s}$" % (filterNames[2], filterNames[3], magType))
    if magmax is not None:
        if magmin is not None:
            title += " [%g < %s < %g]" % (magmin, filterNames[1], magmax)
        else:
            title += " [%s < %g]" % (filterNames[1], magmax)

    title += " %d objects" % nobj
    
    axes.set_title(re.sub(r"^\+\s*", data1.name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat

    did = butler.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data[1].getDataset("calexp_md", did)[0]
    xc = xc[good] + md.get("LTV1")
    yc = yc[good] + md.get("LTV2")
    ids = ids[good]

    if "g" not in SG.lower():
        col12[nonStellar] = -1000
    if "s" not in SG.lower():
        col12[stellar] = -1000
    
    eventHandlers[fig] = EventHandler(axes, col12, col23, ids, xc, yc, flags, frames=frames)

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

    bad = np.bitwise_or(data.flags["INTERP_CENTER"], data.flags["EDGE"])
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

    bad = np.bitwise_and(data.flags, (flagsDict["INTERP_CENTER"] | flagsDict["EDGE"]))
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

def getRefmag(mstars, desiredBand):
    """Get the reference magnitude from a set of match objects"""
    
    try:
        mstars[0]
    except TypeError:
        mstars = [mstars,]

    refmag = np.empty(len(mstars))

    if len(refmag) == 0:
        return refmag

    sch = mstars[0][0].getSchema()

    ct = butler.mapperInfo.getColorterm(desiredBand)
    if False:
        print "RHL Not applying colour terms"
        ct = None

    if ct:
        primaryKey = sch.find(ct.primary).getKey()
        secondaryKey = sch.find(ct.secondary).getKey()

        refMag = -2.5*np.log10(np.array([m[0].get(primaryKey)   for m in mstars]))
        secMag = -2.5*np.log10(np.array([m[0].get(secondaryKey) for m in mstars]))

        refMag = butler.mapperInfo.photometricTransform(desiredBand, refMag, secMag)
    else:
        fluxKey = sch.find("flux").getKey()
        refMag = -2.5*np.log10(np.array([m[0].get(fluxKey) for m in mstars]))

    return refMag

def getCanonicalFilterName(filterName):
    filterId = afwImage.Filter(filterName).getId()

    for name in afwImage.Filter.getNames():
        if afwImage.Filter(name).getId() == filterId:
            return name

    return filterName

def plotCalibration(data, plotBand=0.05, magType='psf', magmin=14, maglim=20, cursor=False,
                    markersize=2, alpha=1.0, title="+", showMedians=False,
                    frame=None, ctype=None, ds9Size=None, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand

If title is provided it's used as a plot title; if it starts + the usual title is prepended
"""
    if not data.matched:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    fig = getMpFigure(fig)

    mstars, zp = [], []
    for m, z in zip(data.matched, data.zp):
        # m[0] == catalogue, m[1] == source
        if m[0].get("stargal") and not m[1].get("flags.pixel.saturated.center"):
            mstars.append(m)
            zp.append(z)

    if frame is not None:
        kwargs = {}
        if ctype:
            kwargs["ctype"] = ctype
        if ds9Size:
            kwargs["size"] = ds9Size
        with ds9.Buffering():
            for i, m in enumerate(mstars):
                s = m[1]
                if zp[i] - 2.5*np.log10(getFlux(s, magType)) < 17.3:
                    ds9.dot("o", s.getX(), s.getY(), frame=frame, **kwargs)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));
    
    refmag = getRefmag(mstars, getCanonicalFilterName(data.dataId[0]["filter"]))

    ids = np.array([s[1].getId() for s in mstars])
    flux = [getFlux(s[1], magType) for s in mstars]

    instmag = zp - 2.5*np.log10(np.array(flux))
    delta = refmag - instmag

    if cursor:
        x = np.empty(len(mstars))
        y = np.empty_like(x)
        flagNames, flags = _setFlagsFromSource(), {}
        for k in flagNames.keys():
            flags[k] = np.empty_like(x)

        for i, s in enumerate(mstars):
            src = s[1]
            x[i], y[i] = src.getX(), src.getY()
            for k, n in flagNames.items():
                flags[k][i] = src.get(n)

        if False:
            for i in range(len(ids)):
                print "%d4 %.2f %.2f (%.2f, %.2f)" % (ids[i], refmag[i], delta[i], x[i], y[i])

    try:
        stats = afwMath.makeStatistics(delta[np.logical_and(refmag > magmin, refmag < maglim)],
                                       afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
    except Exception, e:
        print "Failed to estimate mean: %s" % e
        mean, stdev = float("NaN"), float("NaN")
        
    axes.plot(refmag, delta, "k.", markersize=markersize, alpha=alpha, markeredgewidth=0)

    axes.plot((magmin, maglim, maglim, magmin, magmin), 100*np.array([-1, -1, 1, 1, -1]), "g:")

    axes.plot((-100, 100), (0, 0), "g-")
    if plotBand:
        for i in (-plotBand, plotBand):
            axes.plot((-100, 100), i*np.ones(2), "g--")

    axes.set_ylim(-0.6, 0.6)
    #axes.set_ylim(-2.1, 2.1)
    axes.set_xlim(13, 24)
    axes.set_xlabel("Reference")
    axes.set_ylabel("Reference - %s" % magType)
    axes.set_title(re.sub(r"^\+\s*", data.name + " ", title))

    if showMedians:
        binwidth = 1.0
        bins = np.arange(int(min(refmag)), int(max(refmag)), binwidth)
        vals = np.empty_like(bins)
        for i in range(len(bins) - 1):
            inBin = np.logical_and(refmag > bins[i], refmag <= bins[i] + binwidth)
            vals[i] = np.median(delta[np.where(inBin)])

        axes.plot(bins + 0.5*binwidth, vals, linestyle="-", marker="o", color="cyan")

    fig.text(0.75, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    if cursor:
        global eventHandlers
        eventHandlers[fig] = EventHandler(axes, refmag, delta, ids, x, y, flags, [frame,])

    fig.show()

    return fig

def plotAstromCalibration(data, plotBand=0.5, magType='psf', maglim=20,
                          markersize=1, title="+",
                          frame=None, ctype=None, ds9Size=None, fig=None):
    """Plot (measured - reference) positions given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)
    
If plotBand is provided, draw lines at +- plotBand
    """
    if not data.matched:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    raise RuntimeError("Convert me")

    fig = getMpFigure(fig)

    mstars = [m for m in data.matches if
              m.second and getFlagForDetection(m.second, "STAR")] # data
    realStars = [getFlagForDetection(m.first, "STAR") != 0        # catalogue
                 for m in data.matches if (m.first and \
                                               m.second and getFlagForDetection(m.second, "STAR"))]

    if frame is not None:
        kwargs = {}
        if ctype:
            kwargs["ctype"] = ctype
        if ds9Size:
            kwargs["size"] = ds9Size
        with ds9.Buffering():
            for m in mstars:
                s = m.second
                ds9.dot("o", s.getX(), s.getY(), frame=frame, **kwargs)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    refmag = np.array([-2.5*math.log10(s.first.getPsfFlux()) for s in mstars if s.first and s.second])

    refx = np.array([s.first.getX() for s in mstars if s.first and s.second])
    refy = np.array([s.first.getY() for s in mstars if s.first and s.second])
    srcx = np.array([s.second.getX() for s in mstars if s.first and s.second])
    srcy = np.array([s.second.getY() for s in mstars if s.first and s.second])
        
    ids = np.array([s.second.getId() for s in mstars if s.first and s.second])

    delta = refx - srcx

    flags = np.zeros(len(mstars)) # np.array([s.second.getFlagForDetection() for s in mstars])
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
                print "%.1f %.1f  %.2f %.2f" % (m.second.getX(), m.second.getY(),
                                                -2.5*math.log10(m.first.getPsfFlux()),
                                                -2.5*math.log10(m.first.getPsfFlux()) - 
                                                (data.zp - 2.5*math.log10(m.second.getPsfFlux())))
                ds9.dot("*", m.second.getX(), m.second.getY(), ctype=ds9.MAGENTA, size=10)
                
    return fig

def displayCalibration(data, frame=0):
    """display the calibration objects in Data object data on ds9
Stars are green; galaxies are red based on our processing
    """
    if not data.matched:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    with ds9.Buffering():
        for m in data.matched:
            ref, src = m[0], m[1]
            x1, y1 = src.getX(), src.getY()
            if ref.get("stargal"):      # a star in the reference catalogue
                ptype = "+"
                if src.get("classification.psfstar"):
                    ctype = ds9.GREEN
                else:
                    ctype = ds9.YELLOW
            else:                       # a galaxy in the reference catalogue
                ptype = "o"
                if src.get("classification.psfstar"):
                    ctype = ds9.RED
                else:
                    ctype = ds9.BLUE

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
        for visit, raft, ccd, X, Y, psfMag in queryDB("""
           select
<              visit, raftName, ccdName, XAstrom, YAstrom, dnToABMag(s.psfFlux, exp.fluxMag0) as psfMag
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
                    psfMosaic = maUtils.showPsfMosaic(ims[visit], psf, nx=nx,
                                                      showFWHM=True, frame=None).makeMosaic(mode=nx)
                    sim = im.Factory(im, afwImage.BBox(afwImage.PointI(0, 0), psfMosaic.getWidth(), psfMosaic.getHeight()))
                    sim <<= psfMosaic
                    sim *= 1000
                    del sim

                ds9.mtv(im, wcs=ims[visit].getWcs(), title="%ld %s %s" % (visit, raft, ccd), frame=frame)
                ds9.setMaskTransparency(75)
                ds9.dot("+", X, Y, frame=frame)
                ds9.zoom(4, X, Y, frame=frame)
                nDs9 += 1

            frame += 1

    return ims

def findSource(sourceSet, x, y, radius=2):
    ss = afwDetect.SourceSet()

    s = afwDetect.Source();
    s.setX(x); s.setY(y)
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
    
    calexp_md = butler.get(dtName("calexp", True), **dataId)
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
        s.setX(x), s.setY(y)

        flux = 10**(-0.4*(mag - data.zp))
        s.setApFlux(flux); s.setPsfFlux(flux); s.setModelFlux(flux)

        setFlagForDetection(s, "BINNED1")
        if isStar > 0:
            setFlagForDetection(s, "STAR")
        if varClass != 0:                   # variable
            setFlagForDetection(s, "PEAKCENTER") # XXX

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

def showSourceSet(sourceSet, exp=None, wcs=None, xy0=None, raDec=None, magmin=None, magmax=None, magType="psf",
                  nSource=-1, SG=False, deblend=True, obeyXY0=True,
                  mask=None, symb="+", **kwargs):
    """Show a SourceSet on ds9.

    If nSource > 0, only show the nSource objects with the brightest psf flux
    If mask, it's a set of bitplane names (e.g. INTERP) which must be set to display the source"""
    
    if mask:
        bmask = 0L
        for name in mask:
            bmask |= flagsDict[name]
        mask = bmask

    if exp and not wcs:
        wcs = exp.getWcs()

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

    if nSource > 0:
        fluxInd = sorted(zip([x for x in sourceSet.columns.getPsfFlux()], range(len(sourceSet))),
                          key=lambda x: 0 if np.isnan(x[0]) else -x[0])
        doNotShow = np.ones(len(sourceSet))
        for i in range(nSource):
            doNotShow[fluxInd[i][1]] = 0
    else:
        doNotShow = np.zeros(len(sourceSet))

    # modelMags not available?  Use sourceSet.get("multishapelet.combo.flux") ??
    isStar = (sourceSet.get("classification.extendedness") < 0.5)

    isStar = np.logical_or(isStar, sourceSet.get("flags.pixel.saturated.center"))

    if xy0 is None:
        x0, y0 = exp.getXY0() if exp else [0.0, 0.0]
    else:
        x0, y0 = xy0

    with ds9.Buffering():
        for i, s in enumerate(sourceSet):
            if doNotShow[i]:
                continue

            if mask:
                if not getFlagForDetection(s, mask): # XXX not converted to S2012
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
                x, y = s.getX(), s.getY()

            if obeyXY0:
                x -= x0; y -= y0

            _symb = symb
            if SG:
                kwargs["ctype"] = (ds9.GREEN if isStar[i] else ds9.RED) if s.get("parent") == 0 else \
                    (ds9.CYAN if isStar[i] else ds9.MAGENTA)                    
                _symb = "+" if isStar[i] else "o"

            dx, dy = 0.0, 0.0           # offset of symbol from position
            if symb in ("id", "ID"):
                dx, dy = 0.0, 0.5

                if symb == "id":
                    _symb = butler.mapperInfo.splitId(s.getId(), asDict=True)["objId"]
                else:
                    _symb = "%d" % s.getId()
                if deblend:
                    kwargs["ctype"] = ds9.RED if s.get("parent") == 0 else ds9.MAGENTA                    

            elif symb == "@":
                _symb = "@:%g,%g,%g" % (s.getIxx(), s.getIxy(), s.getIyy())                
            elif deblend:
                kwargs["ctype"] = ds9.RED if s.get("parent") == 0 else ds9.MAGENTA

                pkwargs = kwargs.copy()
                pkwargs["ctype"] = ds9.YELLOW
                pkwargs["size"] = 0.5
                for p in s.getFootprint().getPeaks():
                    ds9.dot("*" if s.get("deblend.deblended-as-psf") else "+", *p.getF(), **pkwargs)

            try:
                ds9.dot(_symb, x + dx, y + dy, **kwargs)
            except:
                pass

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
        x, y = int(s.getX() + 0.5), int(s.getY() + 0.5)
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
        ids, flags, stellar, x, y, shape, mags = getMagsFromSS(s, data.dataId)

        bad = np.bitwise_and(flags, (flagsDict["INTERP_CENTER"] | flagsDict["EDGE"]))

        if stars or galaxies:
            isStar = np.bitwise_and(flags, flagsDict["STAR"])
            if stars:
                bad = np.logical_or(bad, np.logical_not(isStar))
            else:
                bad = np.logical_or(bad, isStar)
                 
        mags = data.getMagsByType(magType, np.logical_not(bad), magDict=mags)

        xyPos.append(zip(mags, x[np.logical_not(bad)], yAstrom[np.logical_not(bad)]))
        arrays.append(np.histogram(mags, bins)[0])

    if ss:
        m, x, y = [], [], []
        for s in ss:
            try:
                mag = calexp.getCalib().getMagnitude(getFlux(s, magType))
            except pexExcept.LsstCppException:
                mag = np.nan
                
            m.append(mag)
            x.append(s.getX())
            y.append(s.getY())

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
            # XXX Need to be converted to S2012 conventions
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
        if getFlagForDetection(s, "SATUR_CENTER"):
            continue

        flux = np.nan if fitAmplitude else getFlux(s, magType)
        try:
            if s.get("classification.extendedness") < 0.5:
                measAlg.subtractPsf(psf, subtracted.getMaskedImage(), s.getX(), s.getY(), flux)
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

def writeRgb(images, rgbFile, min=0, max=50, Q=8, bin=1, scales=[1.0, 1.0, 1.0],
             subtractBkgd=False, boxSize=1024, nsigma=5):
    """Convert the list of images to a true-colour image, and write it to rgbFile (currently .png or .tiff)

Scale the images by scales, if provided

The default stretch is an asinh stretch; Q == 0 corresponds to a linear.  If Q == 0 then min and max are
the range of the stretch; as Q increases they still define the slope of the linear portion of the stretch

If bin is specified, it should be an integer > 1;  the output file will be binned by that factor.
    """
    
    if not afwRgb:
        raise RuntimeError("I was unable to import lsst.afw.extensions.rgb")

    try:
        image = images[:]
    except TypeError:
        images = [images]

    if len(images) == 1:
        grayScale = True
    elif len(images) == 3:
        grayScale = False
    else:
        raise RuntimeError("Please specify one or three images, not %d" % len(images)) 

    for i, image in enumerate(images):
        copied = False                    # have we made a copy?
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
            copied = True

        if scales[i] != 1.0:
            if not copied:
                image = image.Factory(image, True)
                copied = True

            image *= scales[i]

        if subtractBkgd:
            if not copied:
                image = image.Factory(image, True)
                copied = True
 
            sdConfig = SourceDetectionTask.ConfigClass(None)
            sdConfig.reEstimateBackground = True
            sdConfig.thresholdPolarity = "both"
            sdConfig.thresholdType = 'value'
            stats = afwMath.makeStatistics(image, afwMath.MEANCLIP | afwMath.STDEVCLIP)
            image -= stats.getValue(afwMath.MEANCLIP)
            sdConfig.thresholdValue = nsigma*stats.getValue(afwMath.STDEVCLIP)

            sdConfig.background.binSize = boxSize
            try:
                sdConfig.background.isNanSafe = True
            except AttributeError:
                pass
            sdConfig.background.undersampleStyle = "REDUCE_INTERP_ORDER"

            sdTask = SourceDetectionTask(None, config=sdConfig)

            exp = afwImage.makeExposure(afwImage.makeMaskedImage(image))
            sdTask.detectFootprints(exp)
            del exp

        images[i] = image

        if grayScale:
            assert i == 0
            while len(images) < 3:
                images.append(None)
            images[1] = image
            images[2] = image

            break            

    afwRgb.RgbImageF(images[0], images[1], images[2], afwRgb.asinhMappingF(min, max - min, Q)).write(rgbFile)

def grayScale(image, rgbFile, min=0, max=50, Q=8, bin=1):
    """Write a grayscale file (currently png or tiff) to rgbFile

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

    afwRgb.RgbImageF(image, image, image, afwRgb.asinhMappingF(min, max - min, Q)).write(rgbFile)

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
            dist[np.where(np.isnan(dist))] = 1e30
            dmin = min(dist)

            which = np.where(dist == min(dist))
            objId = self.ids[which][0]

            flagsInfo = []
            for k in sorted(self.flags.keys()):
                if self.flags[k][objId == self.ids]:
                    flagsInfo.append(k)
            flagsInfo = ", ".join(flagsInfo)

            print "\r>>>",
            print "%.3f %.3f %s %s%20s\r" % \
                (self.xs[which][0], self.ys[which][0], butler.mapperInfo.splitId(objId), flagsInfo, ""),
            sys.stdout.flush()

            if ev.key == "I":
                print ""
            elif ev.key in ("p", "P"):
                x = self.x[which][0]
                y = self.y[which][0]
                for frame in self.frames:
                    ds9.pan(x, y, frame=frame)
                ds9.cmdBuffer.flush()
        else:
            pass

def showPsfResiduals(data, dataId, fluxLim=9e4, sigma=0, frame=0, **kwargs):
    """Suprime cam specific!"""
    dataId = dict(visit = dataId[_visit_])

    mos = ds9Utils.Mosaic(gutter=4, background=-10)

    for ccd in (8, 9, 5, 4, 3, 6, 7, 2, 1, 0,):
        dataId["ccd"] = ccd

        calexp = data.getDataset("calexp", dataId)[0]

        ss = [s for s in data.getSources(dataId)[0] if
              s.getPsfFlux() > fluxLim and not getFlagForDetection(s, "SATUR") and
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
        x, y = s.getX(), s.getY()
        
        try:
            flux = getFlux(s, magType)
            mag = exposure.getCalib().getMagnitude(flux)

            if not (magMin <= mag <= magMax):
                continue

            if getFlagForDetection(s, "SATUR") or \
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
    @param bitmask Hex mask for pixels to ignore, or list of bit names (default: ["DETECTED", "EDGE"])
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
   
    if maxDN in ("max", None):
        maxDN = np.max(img)

    if minDN == "min":
        minDN = np.min(img)
    elif minDN is None:
        minDN = -maxDN
    if binwidth is None:
        binwidth = 1
    elif binwidth < 0:
        binwidth = (maxDN - minDN)/maxBins

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
def makeContiguous(cat):
    """Make backing storage contiguous, allowing column access"""
    cat = cat.copy(True)
    try:
        cat.getColumnView()
        return cat
    except:
        pass                            # we'll have to do it by hand

    sch = cat.getSchema()
    scm = afwTable.SchemaMapper(sch)
    required = ["id", "coord", "parent",]
    for f in required:
        scm.addMapping(sch.find(f).getKey())

    for f in sch.getNames():
        if f not in required:
            scm.addMapping(sch.find(f).getKey())

    ncat = afwTable.SourceCatalog(scm.getOutputSchema())
    ncat.table.preallocate(len(cat))

    table = cat.table
    ntable = ncat.table
    ntable.defineCentroid(table.getCentroidDefinition())
    if False:
        ntable.defineApFlux(table.getApFluxDefinition())    
        ntable.defineInstFlux(table.getInstFluxDefinition())    
        ntable.defineModelFlux(table.getModelFluxDefinition())    
        ntable.definePsfFlux(table.getPsfFluxDefinition())    
    
    for c in cat:
        ncat.append(ncat.copyRecord(c, scm))

    return ncat
