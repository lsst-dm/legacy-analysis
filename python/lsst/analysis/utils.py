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
import lsst.afw.image.utils as afwImageUtils
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
    import lsst.meas.extensions.multiShapelet.simpleViewer as msViewer
except ImportError:
    msViewer = None
    
try:
    from lsst.obs.lsstSim import LsstSimMapper
    _raw_ = "raw"
    _visit_ = "visit"
except ImportError:
    class LsstSimMapper(object): pass

try:
    from lsst.obs.hsc.hscMapper import HscMapper
    _raw_ = "raw"
    _visit_ = "visit"
except ImportError:
    class HscMapper(object): pass

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
except Exception, e:
    try:
        SdssMapper
    except NameError:
        print >> sys.stderr, "Importing SdssMapper:", e
        class SdssMapper(object): pass

try:
    import lsst.ip.isr as ipIsr
except ImportError:
    ipIsr = None

try:
    warnedHackExtendedness
except NameError:
    warnedHackExtendedness = False      # generate a warning once
    
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
    mpFigures
except NameError:
    mpFigures = {0 : None}              # matplotlib (actually pyplot) figures
    eventHandlers = {}                  # event handlers for matplotlib figures
    eventCallbacks = {}                 # callbacks for event handlers
    eventFrame = 0                      # ds9 frame to use for callbacks

#
try:
    _prefix_                            # prefix for data type names
except NameError:
    _prefix_ = ""

class Prefix(object):
    """An object to use in a with clause to temporarily change _prefix_"""
    def __init__(self, prefix):
        self._saved_prefix_ = _prefix_
        self.prefix = prefix
        
    def __enter__(self):
        global _prefix_
        _prefix_ = self.prefix

    def __exit__(self, *args):
        global _prefix_
        _prefix_ = self._saved_prefix_

def dtName(dType, md=False):
    """Get the name of a given data type (e.g. dtName("src"))"""
    dType_base = re.sub(r"_(filename|md|sub)$", "", dType)
    if _prefix_ in ("", ) or dType_base in ("camera", "coaddTempExp", "goodSeeingCoadd",):
        if md:
            dType += "_md"
        return dType

    dTypeOut = None
    if _prefix_ == "" or dType_base in ("camera", "coaddTempExp", "goodSeeingCoadd",):
        dTypeOut = dType
    elif _prefix_ == "forced":
        if dType == "src":
            dTypeOut = "forcedsources"
        else:
            dTypeOut = dType

    if not dTypeOut:
        dTypeOut = "%s_%s" % (_prefix_, dType)

    if md:
        dTypeOut += "_md"

    return dTypeOut
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
    
def findMapper(root, mapperFile):
    """Search root, looking for mapperFile; follow _parent links if they exist"""
    root0 = root
    while root:
        mapper_py = os.path.join(root, mapperFile)
        if os.path.exists(mapper_py):
            return mapper_py

        parent = os.path.join(root, "_parent")
        if os.path.exists(parent):
            root = parent
        else:
            raise RuntimeError("Unable to find file %s starting at %s" % (mapperFile, root0))

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
    mapper_py = findMapper(registryRoot, mapperFile)
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

    def addPairToName(valName, val0, val1, stride=1):
        """Add a pair of values, val0 and val1, to the valName list"""
        if isinstance(val0, str) and isinstance(val1, str):
            if val0 != val1:
                pre = os.path.commonprefix([val0, val1])
                sval0 = val0[len(pre):]
                sval1 = val1[len(pre):]
        else:
            sval0 = str(val0)
            sval1 = str(val1)
        if sval1 == sval0:
            dvn = str(val0)
        else:
            dvn = "%s-%s" % (sval0, sval1)
            if stride > 1:
                dvn += ":%d" % stride
        valName.append(dvn)
    #
    # Find the minimum spacing between values and interpret it as a stride
    #
    if len(vals) <= 1 or not isinstance(vals[0], int):
        stride = 1
    else:
        stride = vals[1] - vals[0]
        oval = vals[1]
    for val in vals[2:]:
        if val - oval < stride:
            stride = val - oval
        if stride == 1:
            break
        oval = val

    valName = []
    val0 = vals[0]; val1 = val0
    for val in vals[1:]:
        if isinstance(val, int) and val == val1 + stride:
            val1 = val
        else:
            addPairToName(valName, val0, val1, stride=stride)
            val0 = val; val1 = val0

    addPairToName(valName, val0, val1, stride)

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

def idsToDataIds(data, ids):
    """Convert an array of 64-bit IDs to a list of dataIds
    (e.g. [{'ccd': 43, 'visit': 902040}, {'ccd': 1, 'visit': 902040}, ...])
    """
    idDict = data.mapperInfo.splitId(ids, True); del(idDict["objId"])
    return map(dict, set(zip(*[zip(len(v)*[k], v) for k, v in idDict.items()])))

def makeMapperInfo(butler):
    """Return an object with extra per-mapper information (e.g. which fields fully specify an exposure)"""

    mapper = butler.mapper
    
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
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]

            filters = set()
            sensors = set()
            rafts = set()
            visits = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("sensor") == None:
                    dataId["sensor"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
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
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]

            runs = set()
            filters = set()
            camcols = set()
            fields = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("camcol") == None:
                    dataId["camcol"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
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
                        runs.add(dataId[k])
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
                            fields.add(f)

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
            if _prefix_ in ("", "forced"):
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
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]


            if _prefix_ == "stack":
                title = []
                for did in dataIds:
                    title.append("stack %(stack)d patch %(patch)d filter %(filter)s" % did)

                return "[%s]" % "], [".join(title)

            filters = set()
            ccds = set()
            visits = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("ccd") == None:
                    dataId["ccd"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
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
            return 0x0

        @staticmethod
        def splitId(oid, asDict=False):
            """Split an ObjectId into visit, ccd, and objId.
            See obs/subaru/python/lsst/obs/suprimecam/suprimecamMapper.py"""
            oid = np.array(oid, dtype='int64')
            objId = np.bitwise_and(oid, 0xffff) # Should be the same value as was set by apps code
            oid = np.right_shift(oid, 22).astype('int32')

            if _prefix_ == "stack":
                print "Warning: not vectorized"
                nfilter = len(butler.mapper.filters)
                nPatches = 1000000L

                ifilter = oid % nfilter
                oid //= nfilter

                patch = oid % nPatches
                oid //= nPatches

                stack = int(oid)

                filter = [k for k,v  in butler.mapper.filterIdMap.items() if v == ifilter][0]

                if asDict:
                    return dict(stack=stack, patch=patch, filter=filter, objId=objId)
                else:
                    return stack, patch, filter, objId
                
            else:
                oid = np.right_shift(oid, 10).astype('int32')
                ccd = oid % 10
                oid //= 10
                visit = oid

                if visit.size == 1:     # sqlite doesn't like numpy types
                    visit = int(visit)
                    ccd = int(ccd)
                    objId = int(objId)

                if asDict:
                    return dict(visit=visit, ccd=ccd, objId=objId)
                else:
                    return visit, ccd, objId

        @staticmethod
        def canonicalFiltername(filterName):
            mat = re.search(r"W-J-(.)", filterName)
            if mat:
                return mat.group(1)

            mat = re.search(r"W-S-(.)\+", filterName)
            if mat:
                return mat.group(1).lower()

            return filterName

    class SubaruMapperInfoMit(SubaruMapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.__init__(self, None)
            SubaruMapperInfoMit.Mapper = Mapper
            SubaruMapperInfo._Colorterm.setColorterms(colortermsData, "MIT")

    from lsst.obs.hsc.colorterms import colortermsData

    class HscMapperInfo(SubaruMapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.__init__(self, None)
            HscMapperInfo.Mapper = Mapper

            HscMapperInfo._Colorterm.setColorterms(colortermsData, "Hamamatsu")
            
        @staticmethod
        def exposureToStr(exposure):
            try:
                ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getSerial()
                visit = re.sub(r"^HSC", "", exposure.getMetadata().get("FRAMEID"))
            except AttributeError:
                return "??"

            return "%s %s" % (visit, ccdId)

        @staticmethod
        def splitId(oid, asDict=False):
            """Split an ObjectId (maybe an numpy array) into visit, ccd, and objId.
            See obs/subaru/python/lsst/obs/hscSim/hscMapper.py"""
            oid = np.array(oid, dtype='int64')
            objId = np.bitwise_and(oid, 0xffffffff) # Should be the same value as was set by apps code
            oid = np.right_shift(oid, 32).astype('int32')

            if _prefix_ == "stack":
                raise RuntimeError("Please teach HscMapperInfo how to process splitId on a %s" % _prefix_)
            else:
                ccd = (oid % 200).astype('int32')
                oid //= 200
                visit = oid.astype('int32')

                if visit.size == 1:     # sqlite doesn't like numpy types
                    visit = int(visit)
                    ccd = int(ccd)
                    objId = int(objId)

                if asDict:
                    return dict(visit=visit, ccd=ccd, objId=objId)
                else:
                    return visit, ccd, objId

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if isinstance(mapper, LsstSimMapper):
        return LsstSimMapperInfo(LsstSimMapper)
    elif isinstance(butler.mapper, SdssMapper):
        return SdssMapperInfo(SdssMapper)
    elif isinstance(butler.mapper, SuprimecamMapper):
        return SubaruMapperInfo(SuprimecamMapper)
    elif isinstance(butler.mapper, SuprimecamMapperMit):
        return SubaruMapperInfoMit(SuprimecamMapperMit)
    elif isinstance(butler.mapper, HscMapper):
        return HscMapperInfo(HscMapper)
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

    def copy(self, **kwargs):
        """Return a copy, setting any fields specified by kwargs"""
        new = DataId(self)
        for k, v in kwargs.items():
            new[k] = v

        return new

def fixTypesForButler(did):
    for k, v in did.items():
        if isinstance(v, np.int64):
            did[k] = long(v)

    return did

class Data(object):
    def __init__(self, *args, **kwargs):
        if kwargs.has_key("cat"):
            self.cat = kwargs["cat"]
            return

        self._dataRoot = None
        self._rerun = None
        self.butler = None
        self.setButler(*args, **kwargs)
        if self.butler:
            mat = re.search(r"/rerun/(.*)", self.butler.mapper.root)
            self._rerun = mat.group(1)

        self.dataId = []
        self.matches = None

        self.matchedCalibs = None
        self.ZP0 = 31
        self.zp = []

        self.name = "??"
        self.astrom = None
            
    def setButler(self, dataRoot=None, rerun=None):
        if (self._dataRoot and dataRoot != self._DataRoot) or (self._rerun and self._rerun != rerun):
            self.butler = None
            self._dataRoot = dataRoot
            self._rerun = rerun

        if self.butler is not None:
            return

        if dataRoot is None:
            dataRoot = self._dataRoot
        if rerun is None:
            rerun = self._rerun

        if rerun is None and dataRoot is None:
            raise RuntimeError("Please specify a rerun or root")

        if rerun is not None:
            dataRoot = os.path.join(dataRoot, "rerun", rerun)
        if dataRoot:
            dataRoot = os.path.normpath(dataRoot)

        self.butler = dafPersist.Butler(dataRoot)
        self.mapperInfo = makeMapperInfo(self.butler)

        self._dataRoot = dataRoot
        self._rerun = rerun

        if rerun:
            rerunName = rerun
        else:
            dr0 = dataRoot
            if os.path.islink(dr0):
                bdr0 = os.readlink(dr0)

            rerunName = os.path.basename(dr0)

    def get(self, *args, **kwargs):
        return self.butler.get(*args, **kwargs)

    def lookupDataBySkytile(self, dataType):
        """N.b. not converted to use dataId --- Lsst specific"""
        dataSets = {}
        for st, v, f, r, s in self.butler.queryMetadata(_raw_, "skyTile",
                                                   ["skyTile", "visit", "filter", "raft", "sensor"]):
            if self.butler.datasetExists(dataType, visit=v, filter=f, raft=r, sensor=s):
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
        fields = self.mapperInfo.getFields(dataType)

        dataSets = []
        for did in _dataIdDictOuterProduct(dataId, []):
            if _prefix_ in ("forced", "goodSeeingCoadd", "stack"): # XXXXXXX
                if self.butler.datasetExists(dtName(dataType), **did):
                    dataSets.append(did)
            else:
                for vals in self.butler.queryMetadata(_raw_, _visit_, fields, **fixTypesForButler(did)):
                    _dataId = dict(zip(fields, vals))
            
                    if self.butler.datasetExists(dtName(dataType), **_dataId):
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

        if trim and dataType in self.mapperInfo.getTrimmableData():
            return [self.mapperInfo.assembleCcd(dtName(dataType), self.butler,
                                                  dataId, fixAmpLevels=fixAmpLevels)]
        elif dataType in ("eimage",):
            sdataId = dataId.copy(); sdataId["snap"] = dataId.get("snap", 0)
            raw_filename = self.butler.get('raw_filename', channel='0,0', **sdataId)[0]
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
            dataElem = self.butler.get(dtName(dataType), **did)

            if dataType == "psf":
                calexp = self.butler.get(dtName("calexp"), **did)
                dataElem.setDetector(calexp.getDetector())
		del calexp
            elif dataType == 'src':
                cat = dataElem
                sch = cat.getSchema()
                idKey = sch.find("id").getKey()
                parentKey = sch.find("parent").getKey()
                extendednessKey = sch.find("classification.extendedness").getKey()

                idMask = self.mapperInfo.idMask(did)
                if idMask:
                    for s in cat:
                        s.set(idKey, s.getId() | idMask)
                        if s.getParent():
                            s.set(parentKey, s.getParent() | idMask)
                
                if _prefix_ == "forced":
                    poiKey = sch.find("parentObjectId").getKey()
                    for s in cat:
                        s.set(parentKey, s.get(poiKey))

                if True:
                    modelFluxName, corrFac = "flux.kron", 10**(0.4*0.047)
                    extendedness = np.where(cat.get(modelFluxName)*corrFac/cat.getPsfFlux() > 1.05, 1.0, 0.0)
                    extendedness = np.where(np.isfinite(cat.get(modelFluxName)),
                                            extendedness, cat.get(extendednessKey))
                    global warnedHackExtendedness
                    if not warnedHackExtendedness:
                        print >> sys.stderr, "Hacking the extendeness value to use %.2f*%s" % \
                            (corrFac, modelFluxName)
                        warnedHackExtendedness = True
                    for s, e in zip(cat, extendedness):
                        s.set(extendednessKey, e)
                    
            data.append(dataElem)

        return data

    def getSources(self, dataId, setXYfromRaDec=False):
        """Return the list of Sources for the specified dataId
        
        If setXYfromRaDec is true, recalculate (x, y) from (ra, dec)"""

        if setXYfromRaDec:
            calexp_md = self.butler.get(dtName("calexp", True), **dataId)
            wcs = afwImage.makeWcs(calexp_md)

        sources = []
        dataSets = self.dataSets
        for ss in self.getDataset("src", dataId):
            for s in ss:
                s.setId(self.mapperInfo.getId(s))
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
            elif magType == "kron":
                mags = cat.get("kronMag" + suffix)
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

    def getMagsByVisit(self, dataId={}, extraApFlux=0.0, verbose=False):
        """Read the magnitudes for a set of data"""

        dataIds = self.expandDataId(dataId)
        
        catInfo = None
        for did in dataIds:
            if verbose:
                print "Reading %s         \r" % did,
                sys.stdout.flush()

            try:
                catInfo = _appendToCatalog(self, did, catInfo, extraApFlux=extraApFlux)
            except Exception, e:
                print e
                import pdb; pdb.set_trace() 
        if verbose:
            print

        cat = catInfo[0] if catInfo else None
            
        if cat is None:
            raise RuntimeError("Failed to read any data for %s" % " ".join([str(d) for d in dataId]))

        try:
            cat.getColumnView()
        except:
            cat = cat.copy(True)
        self.cat = cat
        
        self.name = self.mapperInfo.dataIdToTitle(dataIds, self._rerun)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def getCalibObjects(self, dataId={}, allSources=True,
                        displayType="psf", frame=None, mtv=False, erase=True, maglim=None, showGalaxies=False,
                        verbose=False):
        """If frame is not None then information about displayType will be plotted on ds9; if
        mtv is true the calexp will be displayed first

        If allSources is True, read the src table and match it to the catalogue
        """
        dataIds = self.expandDataId(dataId)
        if not dataIds:
            raise RuntimeError("%s doesn't expand to any valid dataIds" % dataId)

        self.name = self.mapperInfo.dataIdToTitle(dataIds, self._rerun)

        if mtv and frame is None:
            frame = 0
            
        _matched = []
        _zp = []
        for dataId in dataIds:
            if verbose:
                print "Reading %s         \r" % dataId,
                sys.stdout.flush()
                
            if frame is not None:
                if mtv:
                    calexp = self.getDataset("calexp", dataId)[0]
                    Id = calexp.getDetector().getId() if calexp.getDetector() else self.name
                    ds9.mtv(calexp, title=Id, frame=frame)
                elif erase:
                    ds9.erase(frame)
                else:
                    pass 
                    
            matched, zp, frame = self._getCalibObjectsImpl(dataId, displayType, frame, verbose=verbose,
                                                           maglim=maglim, showGalaxies=showGalaxies,
                                                           allSources=allSources)

            if frame is not None and len(dataIds) > 1:
                while True:
                    try:
                        res = raw_input("%s: [ chnpq] " % str(dataId)[1:-1]).strip()
                    except EOFError:
                        res = ""
                        
                    if res == "c":
                        frame = None
                    elif res == "h":
                        print """' ' \\n c[ontinue without plotting] h[elp] n[ext] p[db] q[uit reading data]

Plotted symbols are:
   if displayType is "src":
      detected sources: x green
      catalog sources:  o red
   else:
      plot a 0.01 error at each matched point
          In red if bad else yellow; symbol is * if stellar else o
      plot error of specified type in cyan if error > 0 else magenta

      if displayType is "jacobian":
         The error is the (Jacobian - 1)
      elif displayType in ("ap", "psf", "sinc")
         error is (ref - mag) of specified type"""
                        continue
                    elif res == "p":
                        import pdb; pdb.set_trace()
                        continue
                    elif res == "q":
                        break
                    elif not res or res == "n":
                        break
                    else:
                        print >> sys.stderr, "Unexpected response: %s" % res
                        continue

                    break

                if res == "q":
                    break                
            
            _matched.append(matched)
            _zp.append(zp)
        if verbose:
            print

        matchedCalibs_old, zp_old = self.matchedCalibs, self.zp
        self.matchedCalibs = sum([len(m) for m in _matched])*[None]
        self.zp = np.empty(len(zp_old) + len(self.matchedCalibs))

        self.zp[0:len(zp_old)] = zp_old

        zstart = len(zp_old)
        start = 0
        for i, m in enumerate(_matched):
            end = start + len(m)
            zend = zstart + len(m)
            self.matchedCalibs[start:end] = m
            self.zp[zstart:zend] = _zp[i]*np.ones(end - start)
            start, zstart = end, zend

        self.matchedCalibs = zipMatchList(self.matchedCalibs, swap=True)
        if matchedCalibs_old:
            if True:                                               # work around a bug
                matchedCalibs_new = afwTable.SourceCatalog(matchedCalibs_old.getSchema())
                matchedCalibs_new.reserve(len(matchedCalibs_old) + len(self.matchedCalibs))
                
                for s in matchedCalibs_old:
                    matchedCalibs_new.append(s)

                for s in self.matchedCalibs:
                    matchedCalibs_new.append(s)

                self.matchedCalibs = matchedCalibs_new.copy(deep=True) # make contiguous
            else:
                matchedCalibs_old.extend(self.matchedCalibs, afwTable.SchemaMapper(matchedCalibs_old.getSchema()))
                self.matchedCalibs = matchedCalibs_old.copy(deep=True) # make contiguous

    def _getCalibObjectsImpl(self, dataId, displayType, frame, maglim=None,
                             showGalaxies=False, verbose=False, allSources=False):
        calexp_md = self.butler.get(dtName("calexp", True), **dataId)
        wcs = afwImage.makeWcs(calexp_md)
        imageSize = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
        xy0 = (calexp_md.get("CRVAL1A"), calexp_md.get("CRVAL2A"),)
        filterName = afwImage.Filter(calexp_md).getName()
            
        calib = afwImage.Calib(calexp_md)

        if not self.astrom:
            astromConfig = measAstrom.Astrometry.ConfigClass()
            astromConfig.catalogMatchDist = 10
            self.astrom = measAstrom.Astrometry(astromConfig)

        ct = self.mapperInfo.getColorterm(filterName)
        primaryFilterName = ct.primary if ct else filterName

        if allSources:
            sources = self.butler.get(dtName("src"), **dataId)
            cat = self.astrom.getReferenceSourcesForWcs(wcs, imageSize, primaryFilterName, pixelMargin=50,
                                                        trim=True, allFluxes=True)
            if frame is not None and displayType == "src":
                showSourceSet(sources, xy0=xy0, raDec=False, frame=frame)
                showSourceSet(cat, xy0=xy0, wcs=wcs, raDec=True,  frame=frame, symb="o", ctype=ds9.RED)
                displayType = None

            try:
                matched = self.astrom._getMatchList(sources, cat, wcs)
            except Exception, e:
                print "RHL", e

                matchRadius = 3
                x, y = sources.getX(), sources.getY()
                x0, y0 = xy0
                xk = sources.getSchema().find("%s.x" % sources.getTable().getCentroidDefinition()).getKey()
                yk = sources.getSchema().find("%s.y" % sources.getTable().getCentroidDefinition()).getKey()
                for i, s in enumerate(sources):
                    s.setD(xk, x[i] - x0)
                    s.setD(yk, y[i] - y0)

                matched = afwTable.matchRaDec(cat, sources, matchRadius*afwGeom.arcseconds)
        else:
            sources = self.butler.get(dtName("icSrc"), **dataId)
            matches = self.butler.get(dtName("icMatch"), **dataId)

            matched = self.astrom.joinMatchListWithCatalog(matches, sources,
                                                           allFluxes=True) # list(ReferenceMatch)

        try:
            zp = calib.getMagnitude(1.0)
        except Exception, e:
            print >> sys.stderr, "Failed to calculate zp: %s" % e
            zp = 0.0

        showJacobian = False
        if displayType == "jacobian":
            showJacobian = True
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

                ct = self.mapperInfo.getColorterm(filterName)
                if ct:
                    primary, secondary = ct.primary, ct.secondary
                    primaryKey_r = sch.find(primary).getKey()
                    secondaryKey_r = sch.find(secondary).getKey()
                else:
                    fluxKey_r = sch.find("flux").getKey()

        else:
            if displayType:
                print >> sys.stderr, "Ignoring unknown displayType %s" % displayType
            frame = None
            
        if showJacobian:
            calexp = self.butler.get('calexp', **dataId)
            ccd = cameraGeom.cast_Ccd(calexp.getDetector())
            distortion = ccd.getDistortion()

            import lsst.afw.geom.ellipses as geomEllipses
            quad = geomEllipses.Quadrupole() # used for estimating jacobian
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

                if frame is not None and displayType:
                    if showJacobian:
                        delta = distortion.distort(afwGeom.PointD(x, y), quad, ccd).getArea() - 1
                    else:
                        if ct:
                            p = -2.5*math.log10(ref.get(primaryKey_r))
                            s = -2.5*math.log10(ref.get(secondaryKey_r))

                            refmag = self.mapperInfo.photometricTransform(filterName, p, s)
                        else:
                            refmag = -2.5*math.log10(ref.get(fluxKey_r))

                        try:
                            delta = refmag - (zp - 2.5*math.log10(src.get(fluxKey)))
                        except:
                            continue

                        if maglim is not None:
                            if refmag > maglim:
                                continue
                        
                    size = min([100, max([deltaScale*delta, -100])])
                    bad = src.get("flags.pixel.interpolated.center") | src.get("flags.pixel.edge")
                    stellar = src.get("classification.extendedness") < 0.5
                    
                    if not showJacobian and not (stellar or showGalaxies):
                        continue

                    ds9.dot("*" if stellar else "o", src.getX(), src.getY(), size=deltaScale*0.01, # 10 mmag
                            frame=frame, ctype=ds9.RED if bad else ds9.YELLOW)

                    if not showJacobian and bad:
                        continue

                    ds9.dot("o", src.getX(), src.getY(), size=abs(size),
                            frame=frame, ctype=ds9.CYAN if size > 0 else ds9.MAGENTA)
           
        if verbose > 1:
            print "Kept %d out of %d reference sources for %s" % (len(keepMatched), len(matched), dataId)

        return keepMatched, zp, frame

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
        BAD =           "flags.pixel.bad",
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

    if _prefix_ == "forced":            # hack hack --- hscMosaic
        calib = data.butler.get("wcs", **dataId).getCalib()
    else:
        calexp_md = data.butler.get(dtName("calexp", True), **dataId)
        calib = afwImage.Calib(calexp_md)

    fluxMag0, fluxMag0Err = calib.getFluxMag0()
    if fluxMag0 <= 0.0:
        fluxMag0 = 1e13
        print >> sys.stderr, "Setting fluxMag0 to %g" % fluxMag0
        calib.setFluxMag0(fluxMag0, fluxMag0Err)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    magTypes = ["ap", "inst", "model", "psf", "kron",]

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
                  "multishapelet.combo.devfrac",
                  "flux.kron.radius",
                  "flags.pixel.bad",
                  "flags.pixel.edge",
                  "flags.pixel.interpolated.any",
                  "flags.pixel.interpolated.center",
                  "flags.pixel.saturated.any",
                  "flags.pixel.saturated.center",
                  "jacobian",
                  ]:
            try:
                scm.addMapping(sch.find(f).getKey())
            except KeyError, e:
                print >> sys.stderr, e
        for key in [tab.getCentroidKey(), tab.getCentroidErrKey(), tab.getCentroidFlagKey(),
                    tab.getShapeKey(),
                    ]:
            try:
                scm.addMapping(key)
            except Exception, e:
                #print >> sys.stderr, e
                pass

        definePhotoSlots = not False            # add standard photometric slots?
        if definePhotoSlots:
            for key in [tab.getInstFluxKey(),  tab.getInstFluxErrKey(),  tab.getInstFluxFlagKey(),
                        tab.getModelFluxKey(), tab.getModelFluxErrKey(), tab.getModelFluxFlagKey(),
                        tab.getPsfFluxKey(),   tab.getPsfFluxErrKey(),   tab.getPsfFluxFlagKey(),
                        tab.getApFluxKey(),    tab.getApFluxErrKey(),     tab.getApFluxFlagKey(),
                        ]:
                scm.addMapping(key)

        scm.addOutputField(afwTable.Field["Flag"]("stellar", "Is the source stellar?"))
        for x in magTypes:
            scm.addOutputField(afwTable.Field["F"]("%sMag" % x,
                                                   "The magnitude corresponding to %sFlux" % x))
            scm.addOutputField(afwTable.Field["F"]("%sMagErr" % x,
                                                   "The error in the magnitude corresponding to %sFlux" % x))

        scm.addOutputField(afwTable.Field["Flag"]("good", "The object is good"))

        cat = afwTable.SourceCatalog(scm.getOutputSchema())
        table = cat.table

        try:
            table.defineCentroid(tab.getCentroidDefinition())
        except:
            pass                        # there may not be a centroid slot defined

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
    size = len(sourceSet)
    if len(cat) == 0 or size > table.getBufferSize():
        table.preallocate(size if size > len(cat) else len(cat)) # add room for new sources

    stellarKey = cat.getSchema().find("stellar").getKey()
    goodKey = cat.getSchema().find("good").getKey()

    try:
        afwTable.SourceRecord.setF
    except AttributeError:
        afwTable.SourceRecord.setD = afwTable.SourceRecord.setF8
        afwTable.SourceRecord.setF = afwTable.SourceRecord.setF4
        afwTable.SourceRecord.setFlag = afwTable.SourceRecord.set_Flag

    tab = sourceSet.table

    with afwImageUtils.CalibNoThrow():
        oldLen = len(cat)

        cat.extend(sourceSet, True, scm)

        if False and sourceSet.get("deblend.nchild") > 0:
            pass

        # ensure that cat is contiguous
        try:
            cat.getColumnView()
        except:
            cat = cat.copy(True)

        cat.get("id")[oldLen:] = sourceSet.get("id")
        if False:                           # can't set bool columns (because they are bitmaps)
            cat.get("good")[oldLen:] = True # for now
            cat.get("stellar")[oldLen:] = s.get("classification.extendedness") < 0.5
        else:
            isStar = (sourceSet.get("classification.extendedness") < 0.5)

            i = oldLen
            for s in sourceSet:
                cat[i].setFlag(goodKey, True)  # for now

                cat[i].setFlag(stellarKey, isStar[i - oldLen])

                i += 1

        for x in magTypes:
            try:
                fluxFieldName =    tab.__getattribute__("get%sFluxDefinition"    % x.title())()
            except AttributeError:          # not available as a slot
                fluxFieldName = "flux.%s" % x

            try:
                tab.getSchema().find(fluxFieldName) # was the flux measured?
            except KeyError, e:
                if catInfo is None:     # we're reading the first CCD
                    print >> sys.stderr, e                
                continue

            flux = sourceSet.get(fluxFieldName)
            fluxErr = sourceSet.get("%s.err" % fluxFieldName)
            if x == "ap":
                flux += extraApFlux

            flux[np.logical_not(np.isfinite(flux))] = np.nan
            fluxErr[fluxErr > 1e38] = np.nan

            mag_magErr = calib.getMagnitude(flux, fluxErr)
            cat.get("%sMag" % x)[oldLen:] = mag_magErr[0]
            cat.get("%sMagErr" % x)[oldLen:] = mag_magErr[1]

    return cat, scm

def zipMatchList(matchList, suffixes=None, swap=False, verbose=True):
    """zip a matchList into a single catalogue

    @param matchList MatchList, as returned by e.g. afwTable.matchRaDec
    @param swap      Interchange [0] and [1]
    """

    records = [matchList[0][i] for i in range(2)]
    if suffixes is None:
        suffixes = ["_1", "_2"]
    if swap:
        records = list(reversed(records))
        suffixes = list(reversed(suffixes))

    requiredFields = ["id", "coord", "parent"] # keys that must be first and in this order --- #2154

    scm = None                          # the SchemaMapper
    keys_2 = []                         # keys from second list that need to be copied over
    for i, rec in enumerate(records):
        suffix = suffixes[i]
        
        sch = rec.getSchema()
        if not scm:
            scm = afwTable.SchemaMapper(sch)

            for f in requiredFields: # must be first and in this order --- #2154
                scm.addMapping(sch.find(f).getKey())

        for schEl in sch:
            inField = schEl.getField()
            name = inField.getName()
            if name == "id" or name not in requiredFields:
                key = schEl.getKey()
                try:
                    outputField = inField.copyRenamed(name + suffix)
                except pexExcept.LsstCppException, e:
                    print >> sys.stderr, "Mapping %s: %s" % (inField.getName(), e)
                    continue

                if name == "id":
                    scm.addOutputField(outputField)
                else:
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

    id_1Key = scm.getOutputSchema().find("id_1").getKey()
    id_2Key = scm.getOutputSchema().find("id_2").getKey()

    for m1, m2, d in matchList:
        if swap:
            m1, m2 = m2, m1
        cat.append(cat.copyRecord(m1, scm))
        rec = cat[-1]

        rec.setL(id_1Key, m1.getId())
        rec.setL(id_2Key, m2.getId())
        
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
    apMags = cat.get(magKeys["ap"])
    instMags = cat.get(magKeys["inst"])
    modelMags = cat.get(magKeys["model"])
    psfMags = cat.get(magKeys["psf"])
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

def overlayCcds(camera, axes, names=False):
    """Overlay camera's CCDs' outlines and serials (or names) using maplotlib axes"""
    for raft in camera:
        raft = cameraGeom.cast_Raft(raft)
        for ccd in raft:
            ccd = cameraGeom.cast_Ccd(ccd)
            ccd.setTrimmed(True)
            
            width, height = ccd.getAllPixels(True).getDimensions()

            corners = ((0.0,0.0), (0.0, height), (width, height), (width, 0.0), (0.0, 0.0))
            for (x0, y0), (x1, y1) in zip(corners[0:4],corners[1:5]):
                if x0 == x1 and y0 != y1:
                    yList = np.linspace(y0, y1, num=2)
                    xList = [x0] * len(yList)
                elif y0 == y1 and x0 != x1:
                    xList = np.linspace(x0, x1, num=2)
                    yList = [y0] * len(xList)
                else:
                    raise RuntimeError("Should never get here")

                xOriginal = []; yOriginal = []
                for x, y in zip(xList, yList):
                    position = ccd.getPositionFromPixel(afwGeom.Point2D(x,y)) # focal plane position

                    xOriginal.append(position.getMm().getX())
                    yOriginal.append(position.getMm().getY())

                axes.plot(xOriginal, yOriginal, 'k-')

            x,y = ccd.getPositionFromPixel(afwGeom.Point2D(width/2, height/2)).getMm()
            cid = ccd.getId()
            if names:
                axes.text(x, y, cid.getName(), ha='center', rotation=90 if height > width else 0,
                        fontsize="smaller")
            else:
                axes.text(x, y, cid.getSerial(), ha='center', va='center')

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmag(data, magType1="model", magType2="psf", maglim=20, magmin=14,
             showMedians=False, sgVal=0.05, criticalFracDeV=0.5,
             selectObjId=None, colorCcds=False, showCamera=False,
             meanDelta=0.0, adjustMean=False, parents=False,
             xmin=None, xmax=None, ymin=None, ymax=None, overlay=True,
             title="+", markersize=1, SG="sg",  color="red", color2="green", frames=[0], wcss=[], fig=None):
    """Plot (magType1 - magType2) v. magType1 mags (e.g. "model" and "psf")

If selectObjId is provided, it's a function that returns True or False for each object. E.g.
    sel = makeSelectCcd(ccd=2)
    plotDmag(..., selectObjId=makeSelectCcd(2), ...)

The magnitude limit for the "locus" box used to calculate statistics is maglim, and only objects within
+- sgLim of meanDelta are included in the locus (the axes are also adjusted to include meanDelta in the
visible area).  If adjustMean is True, adjust all the CCDs being plotted to have the same mean within
the "locus" area.

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If overlay is True, plot the CCDs' outlines and IDs

If non-None, [xy]{min,max} are used to set the plot limits (y{min,max} are interpreted relative to meanDelta
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call data.getMagsByVisit, then try again")

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)

    if parents:
        good = np.logical_and(good, data.cat.get("parent") == 0)
    else:
        good = np.logical_and(good, data.cat.get("deblend.nchild") == 0)

    if selectObjId:
        good = np.logical_and(good, selectObjId(data, data.cat.get("id")))

    mag1 = data.getMagsByType(magType1, good)
    mag2 = data.getMagsByType(magType2, good)
    delta = np.array(mag1 - mag2, dtype='float64') # float64 is needed by makeStatistics --- why?

    locus = np.logical_and(np.logical_and(mag1 > magmin, mag1 < maglim), np.abs(delta - meanDelta) < sgVal)

    stellar = data.cat.get("stellar")
    nonStellar = np.logical_not(stellar)
    stellar = stellar[good]; nonStellar = nonStellar[good]
    ids = data.cat.get("id")[good]

    if "s" not in SG.lower():
        mean = None
    else:
        try:
            stats = afwMath.makeStatistics(delta[locus], afwMath.STDEVCLIP | afwMath.MEANCLIP)
            mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
        except:
            mean, stdev = np.nan, np.nan

    if False:
        axes.plot(mag1[locus], delta[locus], "o", markersize=2*markersize, color="blue")

    if "d" in SG.lower() or "e" in SG.lower():           # split by deV/exp?
        SG += "g"

    plotKW = dict(markersize=markersize, markeredgewidth=0, zorder=1)
    nobj = 0
    if "g" in SG.lower():
        if "d" in SG.lower() or "e" in SG.lower():           # split by deV/exp
            try:
                fracDeV = data.cat.get("multishapelet.combo.devfrac")[good]
            except KeyError:
                print >> sys.stderr, "No frac_deV is available"
                fracDeV = np.ones_like(nonStellar)
                
            deV = fracDeV > criticalFracDeV # it's a deV

            plotKW["alpha"] = 0.5
            if "d" in SG.lower():
                tmp = np.logical_and(deV, nonStellar)
                nobj += tmp.sum()
                axes.plot(mag1[tmp], delta[tmp], "o", color="red", **plotKW)
            if "e" in SG.lower():
                tmp = np.logical_and(np.logical_not(deV), nonStellar)
                nobj += tmp.sum()
                axes.plot(mag1[tmp], delta[tmp], "o", color="blue", **plotKW)
        else:
            nobj += nonStellar.sum()
            axes.plot(mag1[nonStellar], delta[nonStellar], "o", color=color, **plotKW)
    if "s" in SG.lower():
        nobj += stellar.sum()

        if showCamera:
            #
            # Get focal-plane positions ("mm" -- for HSC they are actually in units of 15 micron pixels)
            #
            inRange = np.logical_and(stellar, np.logical_and(mag1 > magmin, mag1 < maglim))

            camera = data.butler.get("camera")
            xc = data.cat.getX()[good][inRange]
            yc = data.cat.getY()[good][inRange]

            xmm = np.empty_like(xc); ymm = np.empty_like(yc)
            for i, _id in enumerate(ids[inRange]):
                ccdId = data.mapperInfo.splitId(_id, True)["ccd"]
                ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(ccdId))
                xmm[i], ymm[i] = ccd.getPositionFromPixel(afwGeom.PointD(xc[i], yc[i])).getMm()


            sc = axes.scatter(xmm, ymm, c=delta[inRange], norm=pyplot.Normalize(ymin, ymax),
                              cmap=pyplot.cm.rainbow, marker="o", s=10*markersize,
                              edgecolors="none")
            fig.colorbar(sc)
            axes.set_aspect('equal')

            axes.set_xlabel("X (mm)")
            axes.set_ylabel("Y (mm)")

            title += " %s - %s" % (magType1, magType2)

            if overlay:
                overlayCcds(camera, axes) 
        elif colorCcds:
            ccdIds = data.mapperInfo.splitId(ids[stellar], True)["ccd"]
            sc = axes.scatter(mag1[stellar], delta[stellar], c=ccdIds,
                              cmap=pyplot.cm.rainbow, marker="o", s=5*markersize,
                              edgecolors="none")
            fig.colorbar(sc)
        else:
            axes.plot(mag1[stellar], delta[stellar], "o", markersize=markersize, markeredgewidth=0,
                      color=color2, zorder=1)

    if not showCamera:
        axes.axhline(meanDelta, linestyle=":", color="black")
        if mean is not None:
            axes.plot((magmin, maglim, maglim, magmin, magmin),
                      meanDelta + sgVal*np.array([-1, -1, 1, 1, -1]), linestyle=":", color="black")

        if showMedians:
            if sum(stellar) == 0:
                print >> sys.stderr, "Without stars I can't show the stellar sequence"
                showMedians = False
        if showMedians:
            binwidth = 1.0
            bins = np.arange(np.floor(magmin), np.floor(max(mag1[stellar])), binwidth)
            vals = np.empty_like(bins)
            err = np.empty_like(bins)
            for i in range(len(bins) - 1):
                inBin = np.logical_and(mag1 > bins[i], mag1 <= bins[i] + binwidth)
                tmp = delta[np.where(np.logical_and(stellar if 's' in SG else nonStellar, inBin))]
                tmp.sort()

                if len(tmp) == 0:
                    median, iqr = np.nan, np.nan
                else:
                    median = tmp[int(0.5*len(tmp) + 0.5)] if len(tmp) > 1 else tmp[0]
                    iqr = (tmp[int(0.75*len(tmp) + 0.5)] - tmp[int(0.25*len(tmp) + 0.5)]) \
                        if len(tmp) > 2 else np.nan

                vals[i] = median
                err[i] = 0.741*iqr

            axes.errorbar(bins + 0.5*binwidth, vals, yerr=err, zorder=3,
                          linestyle="-", marker="o", color="black")

        axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
        axes.set_ylim(meanDelta + (0.2 if ymin is None else ymin),
                      meanDelta + (- 0.8 if ymax is None else ymax))
        axes.set_xlabel(magType1)
        axes.set_ylabel("%s - %s" % (magType1, magType2))

    title += " %d objects" % nobj
    name = data.name if selectObjId is None else \
        data.mapperInfo.dataIdToTitle(idsToDataIds(data, ids), data._rerun)
    axes.set_title(re.sub(r"^\+\s*", name + " ", title))

    if mean is not None:
        fig.text(0.20, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    if "g" not in SG.lower():
        mag1[nonStellar] = -1000        # we don't want to pick these objects
    if "s" not in SG.lower():
        mag1[stellar] = -1000
        if "d" in SG.lower():
            if "e" in SG.lower():       # show deV and /exp
                pass
            else:
                mag1[np.logical_not(deV)] = -1000
        elif "e" in SG.lower():         # show only exp
            mag1[deV] = -1000

    did = data.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data.getDataset("calexp_md", did)[0]

    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat
    try:
        x = data.cat.getX()[good] + md.get("LTV1")
        y = data.cat.getY()[good] + md.get("LTV2")
    except pexExcept.LsstCppException, e:
        if not re.search(r"pex::exceptions::LogicErrorException", e.message.getType()):
            raise e
        x = np.zeros_like(mag1)
        y = np.zeros_like(mag1)
        
    ids = data.cat.get("id")[good]
    eventHandlers[fig] = EventHandler(data, axes, mag1, delta, ids, x, y, flags, frames=frames, wcss=wcss)

    fig.show()

    return fig

def plotClassification(data, magType1="model", maglim=20, magmin=14,
             showMedians=False, parents=False, criticalFracDeV=0.5,
             xmin=None, xmax=None, ymin=None, ymax=None,
             title="+", markersize=1, SG="sg",  color="red", color2="green", frames=[0], fig=None):
    """Plot objects' classification v. magType1

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)

    if parents:
        good = np.logical_and(good, data.cat.get("parent") == 0)
    else:
        good = np.logical_and(good, data.cat.get("deblend.nchild") == 0)

    mag1 = data.getMagsByType(magType1, good)
    fracDeV = data.cat.get("multishapelet.combo.devfrac")[good]
    deV = fracDeV > criticalFracDeV

    nonStellar = np.logical_not(data.cat.get("stellar"))[good]
    ids = data.cat.get("id")[good]

    plotKW = dict(markersize=markersize, markeredgewidth=0, zorder=1)
    if "d" in SG.lower() or "e" in SG.lower():           # split by deV/exp
        plotKW["alpha"] = 0.5
        if "d" in SG.lower():
            tmp = np.logical_and(deV, nonStellar)
            axes.plot(mag1[tmp], fracDeV[tmp], "o", color="red", **plotKW)
        if "e" in SG.lower():
            tmp = np.logical_and(np.logical_not(deV), nonStellar)
            axes.plot(mag1[tmp], fracDeV[tmp], "o", color="blue", **plotKW)
    else:
        axes.plot(mag1[nonStellar], fracDeV[nonStellar], "o", color=color, **plotKW)

    axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
    axes.set_ylim(-0.1 if ymin is None else ymin, 1.1 if ymax is None else ymax)
    axes.set_xlabel(magType1)
    axes.set_ylabel("$frac_{deV}$")
    axes.set_title(re.sub(r"^\+\s*", data.name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    did = data.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data.getDataset("calexp_md", did)[0]

    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat
    x = data.cat.getX()[good] + md.get("LTV1")
    y = data.cat.getY()[good] + md.get("LTV2")
    ids = data.cat.get("id")[good]
    eventHandlers[fig] = EventHandler(data, axes, mag1, fracDeV, ids, x, y, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotRadii(data, magType1="model", maglim=20, magmin=14,
              showMedians=False, parents=False, criticalFracDeV=0.5,
              xmin=None, xmax=None, ymin=None, ymax=None,
              title="+", markersize=1, SG="sg",  color="red", color2="green", frames=[0], fig=None):
    """Plot objects' radii v. magType1

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits
    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    if "d" in SG.lower() or "e" in SG.lower():           # split by deV/exp?
        SG += "g"

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)

    if parents:
        good = np.logical_and(good, data.cat.get("parent") == 0)
    else:
        good = np.logical_and(good, data.cat.get("deblend.nchild") == 0)

    mag1 = data.getMagsByType(magType1, good)
    radii = data.cat.get("flux.kron.radius")[good]
    fracDeV = data.cat.get("multishapelet.combo.devfrac")[good]
    deV = fracDeV > criticalFracDeV

    stellar = data.cat.get("stellar")[good]
    nonStellar = np.logical_not(stellar)
    ids = data.cat.get("id")[good]

    plotKW = dict(markersize=markersize, markeredgewidth=0, zorder=1)
    if "g" in SG.lower():
        if "d" in SG.lower() or "e" in SG.lower():           # split by deV/exp
            plotKW["alpha"] = 0.5
            if "d" in SG.lower():
                tmp = np.logical_and(deV, nonStellar)
                axes.plot(mag1[tmp], radii[tmp], "o", color="red", **plotKW)
            if "e" in SG.lower():
                tmp = np.logical_and(np.logical_not(deV), nonStellar)
                axes.plot(mag1[tmp], radii[tmp], "o", color="blue", **plotKW)
        else:
            axes.plot(mag1[nonStellar], radii[nonStellar], "o", color=color, **plotKW)
    if "s" in SG.lower():
        axes.plot(mag1[stellar], radii[stellar], "o", markersize=markersize, markeredgewidth=0,
                  color=color2, zorder=1)

    axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
    axes.set_ylim(-0.1 if ymin is None else ymin, 1.1 if ymax is None else ymax)
    axes.set_xlabel(magType1)
    axes.set_ylabel("$R_{%s}$" % (magType1.title()))
    axes.set_title(re.sub(r"^\+\s*", data.name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    did = data.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data.getDataset("calexp_md", did)[0]

    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat
    x = data.cat.getX()[good] + md.get("LTV1")
    y = data.cat.getY()[good] + md.get("LTV2")
    ids = data.cat.get("id")[good]
    eventHandlers[fig] = EventHandler(data, axes, mag1, radii, ids, x, y, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeSelectCcd(ccds, include=True):
    def selectCcd(data, id, ccds=ccds, include=include):
        try:
            ccds[0]
        except TypeError:
            ccds = [ccds]

        if ccds[0] is None:
            return True

        if not isinstance(id, dict):
            id = data.mapperInfo.splitId(id, asDict=True)
        ccd = id["ccd"]
        matches = reduce(lambda m, c: np.logical_or(m, ccd == c), ccds, False)

        return matches if include else np.logical_not(matches)

    return selectCcd

def makeSelectVisit(visits=None, ccds=None, include=True):
    try:
        visits[0]
    except TypeError:
        visits = [visits]

    if ccds is None:
        selectCcd = None
    else:
        selectCcd = makeSelectCcd(ccds)

    def selectVisit(data, id, visits=visits, include=include):
        if visits[0] is None and selectCcd is None:
            return True

        if not isinstance(id, dict):
            id = data.mapperInfo.splitId(id, asDict=True)

        if visits[0] is None:
            matches = np.ones_like(id)
        else:
            visit = id["visit"]
            matches = reduce(lambda m, v: np.logical_or(m, visit == v), visits, False)

        if selectCcd:
            matches = np.logical_and(matches, selectCcd(data, id))

        return matches if include else np.logical_not(matches)

    return selectVisit

def plotCM(data, select1, select2, magType="psf", maglim=20, magmin=14,
           SG="sg", showMedians=False,
           xmin=None, xmax=None, ymin=None, ymax=None,
           title="+", markersize=1, color="red", alpha=1.0, frames=[0], verbose=False, fig=None):
    """Plot (data[select1].magType - data[select2].magType) v. data1.magType mags (e.g. "psf")
where data[selectN] means, "the objects in data for which selectN(ids) returns True".  E.g.

    plotCM(data, makeSelectVisit(902030), makeSelectVisit(902032)...)

If a select "function" is actually an int, it is interpreted as makeSelectVisit(selectN) -- i.e.
as a visit IDs

The magnitude limits for the box used to calculate statistics are magmin..maglim

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If non-None, [xy]{min,max} are used to set the plot limits
    """

    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call data.getMagsByVisit, then try again")

    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))
    #
    selections = {1: select1,
                  2: select2,
                  }

    allIds = data.cat.get("id")
    for k, v in selections.items():
        if isinstance(v, int):
            selections[k] = makeSelectVisit(v)

        selections[k] = selections[k](data, allIds)

    cats = dict([(k, data.cat[v]) for k,v in selections.items()])
    #
    # See if we've done this match before
    #
    if data.matches:
        if set(selections.keys()) != set(data.matches[0].keys()):
            data.matches = None
        else:
            dselections = data.matches[0]
            for k in selections.keys():
                if not (selections[k] == dselections[k]).all():
                    data.matches = None
                    break

    if data.matches:
        matched = data.matches[1]
    else:
        if verbose:
            print "Matching selected sources"
        matchRadius = 2                     # arcsec
        mat = afwTable.matchRaDec(cats[1], cats[2], matchRadius*afwGeom.arcseconds)
        matched = Data(cat=zipMatchList(mat))
        
        data.matches = selections, matched

    bad = None
    for i in range(2):
        suffix = "_%d" % (i + 1)
        for name in ["flags.pixel.edge",
                     "flags.pixel.bad",
                     #"flags.pixel.interpolated.center",
                     "flags.pixel.saturated.center",]:
            _flg = matched.cat.get(name + suffix)
            if bad is None:
                bad = _flg
            else:
                bad = np.logical_or(bad, _flg)

    good = np.logical_not(bad)
    stellar = matched.cat.get("stellar_1")
    if SG.lower() == 's':
        good = np.logical_and(good, stellar)
    elif SG.lower() == 'g':
        good = np.logical_and(good, np.logical_not(stellar))

    ids = matched.cat.get("id")[good]
    xc = matched.cat.get("centroid.sdss_1.x")[good]
    yc = matched.cat.get("centroid.sdss_1.y")[good]
    stellar = stellar[good]

    mag1 = matched.getMagsByType(magType, good, suffix="_1")
    mag2 = matched.getMagsByType(magType, good, suffix="_2")
    delta = np.array(mag1 - mag2, dtype='float64') # float64 is needed by makeStatistics --- why?

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
                  alpha=alpha["g"], markersize=markersize, markeredgewidth=0, color=color, zorder=-1)
    if "s" in SG.lower():
        nobj += np.sum(stellar)
        axes.plot(delta[stellar], mag1[stellar], "o",
                  alpha=alpha["s"], markersize=markersize, markeredgewidth=0, color=color2, zorder=-1)

    if showMedians:
        binwidth = 0.5
        bins = np.arange(np.floor(magmin), np.floor(np.nanmax(mag1)), binwidth)
        vals = np.empty_like(bins)
        err = np.empty_like(bins)
        for i in range(len(bins) - 1):
            inBin = np.logical_and(mag1 > bins[i], mag1 <= bins[i] + binwidth)
            if SG.lower() == 's':
                tmp = delta[np.where(np.logical_and(stellar, inBin))]
            else:
                tmp = delta[inBin]
            tmp.sort()

            if len(tmp) == 0:
                median, iqr = np.nan, np.nan
            else:
                median = tmp[int(0.5*len(tmp) + 0.5)] if len(tmp) > 1 else tmp[0]
                iqr = (tmp[int(0.75*len(tmp) + 0.5)] - tmp[int(0.25*len(tmp) + 0.5)]) \
                    if len(tmp) > 2 else np.nan

            vals[i] = median
            err[i] = 0.741*iqr

        if True:
            axes.errorbar(vals, bins + 0.5*binwidth, xerr=err, zorder=3,
                          linestyle="-", marker="o", color="black")
        else:
            axes.plot(vals, bins + 0.5*binwidth, zorder=3, marker="o", color="black")
            for pm in (-1, 1):
                axes.plot(vals + pm*err, bins + 0.5*binwidth, zorder=3, linestyle="-", color="blue")

        axes.axvline(0.0, color="blue", ls=":")

        try:
            inBin = np.logical_and(mag1 > magmin, mag1 <= maglim)
            stats = afwMath.makeStatistics(delta[inBin], afwMath.STDEVCLIP | afwMath.MEANCLIP)
            mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
            print "%-5s %5.2f%% +- %5.2f%%" % (magType, 100*mean, 100*stdev)
        except Exception, e:
            print "Failed to estimate mean: %s" % e
            mean, stdev = float("NaN"), float("NaN")

        for m in (magmin, maglim,):
            axes.axhline(m, color="black", ls=":")

        fig.text(0.20, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")
        
    filterNames = [None]
    for k in sorted(cats.keys()):
        dataId = idsToDataIds(data, allIds[selections[k]][0:2])[0]
        for k, v in dataId.items():
            if isinstance(v, np.int32):
                dataId[k] = int(v)
        filterNames.append(afwImage.Filter(data.butler.get(dtName("calexp", True), **dataId)).getName())
    filter1, filter2 = filterNames[1], filterNames[2]

    axes.set_xlim(-1 if xmin is None else xmin, 2 if xmax is None else xmax)
    axes.set_ylim(24 if ymax is None else ymax, 14 if ymin is None else ymin)
    axes.set_xlabel("(%s - %s)$_{%s}$" % (filter1, filter2, magType))
    axes.set_ylabel("%s$_{%s}$" % (filter1, magType))

    title += " %d objects" % nobj
    name = ", ".join([data.mapperInfo.dataIdToTitle(
                idsToDataIds(data, allIds[selections[i]])) for i in (1, 2)])

    axes.set_title(re.sub(r"^\+\s*", name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat

    did = data.mapperInfo.splitId(ids[0], asDict=True); del did["objId"]
    md = data.getDataset("calexp_md", did)[0]
    xc = xc[good] + md.get("LTV1")
    yc = yc[good] + md.get("LTV2")
    ids = ids[good]

    if "g" not in SG.lower():
        delta[nonStellar] = -1000
    if "s" not in SG.lower():
        delta[stellar] = -1000

    eventHandlers[fig] = EventHandler(data, axes, delta, mag1, ids, xc, yc, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def robustStd(y):
    """Return a robust estimate of y's standard deviation based on the IQR"""
    
    y = sorted(y)
    n = len(y)
    return 0.741*(y[int(0.75*n)] - y[int(0.25*n)])

def plotStellarLocus(axes, mags, k1, k2, k3, stellar, filterNames, stellarLocusEnds=[], locusLtype=False,
                     usePrincipalColor=False, plotRaDec=False):
    """
    Calculate width of blue end of stellar locus
    """
    col12 = mags[k1] - mags[k2]
    col23 = mags[k2] - mags[k3]

    #stellarLocusEnds = (0.40, 1.00,)    # the blue and red colour defining the straight part of the locus
    #stellarLocusEnds = (0.50, 1.50,)

    # Find the stellar locus
    xy = []
    try:
        stellarLocusEnds[0][0]  # both x and y are provided
        stellarLocusEnds[0][0]  # both x and y are provided
        xy = stellarLocusEnds
        stellarLocusEnds = (xy[0][0], xy[1][0])
    except TypeError:
        dxx = 0.1*abs(stellarLocusEnds[1] - stellarLocusEnds[0])
        for xx in stellarLocusEnds:
            delta = col23[stellar][abs(col12[stellar] - xx) < dxx]
            try:
                yy = afwMath.makeStatistics(np.array(delta, dtype="float64"), afwMath.MEDIAN).getValue()
            except pexExcept.LsstCppException:
                yy = np.nan

            if False and not plotRaDec:
                axes.axvline(xx, color="blue", ls=":")

            xy.append((xx, yy))

    if False:
        print "Stellar locus: [(%.3f, %.3f), (%.3f, %.3f)]" % (xy[0][0], xy[0][1], xy[1][0], xy[1][1])

    principalColor = ""
    if usePrincipalColor and filterNames[k1] == 'g' and filterNames[k2] == 'r' and filterNames[k3] == 'i':
        pc = -0.227*mags[k1] + 0.792*mags[k2] - 0.567*mags[k3] + 0.050
        principalColor = "w"
        delta = pc[stellar]
    elif usePrincipalColor and filterNames[k1] == 'r' and filterNames[k2] == 'i' and filterNames[k3] == 'z':
        pc = -0.270*mags[k1] + 0.800*mags[k2] - 0.534*mags[k3] + 0.054
        principalColor = "y"
        delta = pc[stellar]
    else:
        theta = math.atan2(xy[1][1] - xy[0][1], xy[1][0] - xy[0][0])
        c, s = math.cos(theta), math.sin(theta)

        x, y = col12[stellar], col23[stellar]
        xp =   x*c + y*s
        yp = - x*s + y*c

        delta = yp

    if locusLtype and not plotRaDec:
        axes.plot([xy[0][0], xy[1][0]], [xy[0][1], xy[1][1]], locusLtype)


    delta = np.array(delta, dtype="float64")
    blue = np.logical_and(col12 > stellarLocusEnds[0], col12 < stellarLocusEnds[1])[stellar]
    try:
        stats = afwMath.makeStatistics(delta[blue], afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
    except:
        mean, stdev = float("NaN"), float("NaN")

    #print "%g +- %g" % (mean, stdev)
    axes.text(xy[0][0], 0.5*(xy[0][1] + xy[1][1]) + 0.25, r"$%s \pm %.3f$" % \
                  ("%s = " % principalColor if principalColor else "", stdev), fontsize="larger")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotCC(data, select1, select2, select3, magType="psf", SG="sg", select4=None,
           verbose=False, fig=None, *args, **kwargs):
    """
Plot (data[select1].magType - data[select2].magType) v. (data[select1].magType - data[select3].magType)
where data[selectN] means, "the objects in data for which selectN(ids) returns True".  E.g.
    plotCC(data, makeSelectVisit(902030), makeSelectVisit(902032), makeSelectVisit(902034), "psf", ...)

If selectN is an int it's interpreted as makeSelectVisit(selectN), so this is equivalent to
    plotCC(data, 902030, 902032, 902034, "psf", ...)

(let's write data[selectN].magType as magN).  If provided, idN specifies the data from which selection
function should to be used for magnitude limits (e.g. idN == 3 => visit 902034 in the example)

This can be used to plot 3-colour diagrams or to compare 3 epochs.

If select4 is provided, compare (mag1 - mag2) with (mag3 - mag4)
If both select1 and select3 are the same band and select2 and select4 are the same (different) band, plot
   (mag1 - mag2) - (mag3 - mag4) against mag1

If title is provided it's used as a plot title; if it starts + the usual title is prepended and if it's
None no title is provided; if it starts "T:" the label will be written at the top of the plot

showXlabel may be "bottom", "top", or None;  showYlabel may be "left", "right", or None.

If non-None, [xy]{min,max} are used to set the plot limits
    """

    if isinstance(magType, str):
        magType = [magType,]
    if isinstance(SG, str):
        SG = [SG,]

    fig = getMpFigure(fig)

    selections = {1: select1,
                  2: select2,
                  3: select3
                  }
    if select4 is not None:
        selections[4] = select4

    allIds = data.cat.get("id")
    for k, v in selections.items():
        if isinstance(v, int):
            selections[k] = makeSelectVisit(v)

        selections[k] = selections[k](data, allIds)

    dataKeys = sorted(selections.keys())

    cats = dict([(k, data.cat[v]) for k,v in selections.items()])

    empty = []
    for i, d in cats.items():
        if len(d) == 0:
            empty.append(str(i))

    if empty:
        raise RuntimeError("No objects were found with selection%s %s" % ("s" if len(empty) > 1 else "",
                                                                          ", ".join(empty)))
    #
    # See if we've done this match before
    #
    if data.matches:
        if set(selections.keys()) != set(data.matches[0].keys()):
            data.matches = None
        else:
            dselections = data.matches[0]
            for k in selections.keys():
                if not (selections[k] == dselections[k]).all():
                    data.matches = None
                    break

    if data.matches:
        matched = data.matches[1]
    else:
        if verbose:
            print "Matching selected sources"
        matchRadius = 2                     # arcsec
        mat = afwTable.matchRaDec(cats[1], cats[2], matchRadius*afwGeom.arcseconds)
        cat12=zipMatchList(mat)
        mat = afwTable.matchRaDec(cat12, cats[3], matchRadius*afwGeom.arcseconds)
        cat123 = zipMatchList(mat)
        if select4 is None:
            matched = Data(cat=cat123)
        else:
            mat = afwTable.matchRaDec(cat123, cats[4], matchRadius*afwGeom.arcseconds)
            matched = Data(cat=zipMatchList(mat))
        
        data.matches = selections, matched

    filterNames = [None]
    visitNames = [None]
    for k in sorted(cats.keys()):
        dataId = idsToDataIds(data, allIds[selections[k]][0:2])[0]
        for k, v in dataId.items():
            if isinstance(v, np.int32):
                dataId[k] = int(v)
        visitNames.append(dataId["visit"])
        filterNames.append(afwImage.Filter(data.butler.get(dtName("calexp", True), **dataId)).getName())

    datasetName = ", ".join([data.mapperInfo.dataIdToTitle(
                idsToDataIds(data, allIds[selections[i]])) for i in dataKeys])

    subplots = makeSubplots(fig, nx=len(magType), ny=len(SG))
    j = -1
    for _magType in magType:
        for _sg in SG:
            axes = subplots.next(); j += 1
            
            if len(SG) == 1 and len(magType) == 1:
                showXlabel = "bottom"
                showYlabel = "left"
            else:
                showXlabel = \
                    "top"    if False and j//len(SG) == 0      else \
                    "bottom" if j//len(SG) == len(magType) - 1 else None
                showYlabel = \
                    "left"  if j//len(SG) == len(magType)//2 and _sg == SG[0] else \
                    "right" if j//len(SG) == len(magType)//2 and _sg == SG[len(SG)-1] else None

            _, matched = _plotCCImpl(data, matched, dataKeys, _magType, filterNames, visitNames, _sg,
                                     verbose=verbose, fig=axes,
                                     title="T:+" if j == 0 else None, datasetName=datasetName,
                                     showXlabel=showXlabel, showYlabel=showYlabel,
                                     *args, **kwargs
                                     )
                                     
    fig.show()

    return fig

def _plotCCImpl(data, matched, dataKeys, magType, filterNames, visitNames, SG, fig=None, show=True,
                magmax=None, magmin=None, 
                verbose=False,
                idN=None, idColorN=None, selectObjId=None, matchRadius=2, plotRaDec=False,
                showStatistics=False, show_r_xy=True, colorCcds=False, colorVisits=False,
                usePrincipalColor=True, stellarLocusEnds=[],
                adjustLocus=False, locusLtype="b:",
                xmin=None, xmax=None, ymin=None, ymax=None,
                showXlabel="bottom", showYlabel="left",
                datasetName="", title="+",
                markersize=1, alpha=1.0, color="red", frames=[0], wcss=[]):
    """
    \param idN   The (1-based) index into the selection functions to use for magnitude limits etc.
    """

    subplot = isinstance(fig, pyplot.Axes)
    if subplot:
        axes = fig
        fig = axes.figure
    else:
        fig = getMpFigure(fig)
        
        axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    if idN is None:
        idN = dataKeys[0]
    if idColorN is None:
        idColorN = idN

    suffixes = ["_2" + i*"_1" for i in range(len(dataKeys) - 1, -1, -1)]
    suffixes[0] = suffixes[0][2:]       # first name has no "_2".  Grrr
    suffixes = dict(zip(dataKeys, suffixes))

    if len(dataKeys) == 3:
        k1, k2, k3 = dataKeys
    else:
        k1, k2, k3, k4 = dataKeys

    # are we looking at e.g. (g-r)_1 - (g-r)_2?  If so, plot against magnitude
    deltaColor = len(dataKeys) > 3 and \
        filterNames[k1] != filterNames[k2] and \
        (filterNames[k1], filterNames[k2]) == (filterNames[k3], filterNames[k4])

    mag1 = matched.getMagsByType(magType, suffix=suffixes[1])
    good = True

    for suffix in suffixes.values():
        for name in ["flags.pixel.edge",
                     "flags.pixel.bad",
                     #"flags.pixel.interpolated.center",
                     "flags.pixel.saturated.center",]:
            good = np.logical_and(good, np.logical_not(matched.cat.get(name + suffix)))

    good = np.logical_and(good, sum([matched.cat.get("deblend.nchild%s" % s) for s in suffixes.values()]) == 0)
    good = np.logical_and(good, matched.cat.get("parent") == 0)

    centroidStr = "centroid.sdss%s" % suffixes[idN]
    xc = matched.cat.get("%s.x" % centroidStr)[good]
    yc = matched.cat.get("%s.y" % centroidStr)[good]
    stellar = matched.cat.get("stellar%s" % suffixes[idN])[good]

    mags, ids = {}, {}
    for k, s in suffixes.items():
        mags[k] = matched.getMagsByType(magType, good, suffix=s)
        ids[k]  = matched.cat.get("id%s" %s )[good]

    good = good[good]

    if not deltaColor:
        if magmin is not None:
            good = np.logical_and(good, mags[idN] > magmin)
        if magmax is not None:
            good = np.logical_and(good, mags[idN] < magmax)

    for k, s in suffixes.items():
        mags[k] = mags[k][good]

    stellar = stellar[good]
    nonStellar = np.logical_not(stellar)

    for k in ids.keys():
        ids[k] = ids[k][good]

    if plotRaDec:
        ra = np.degrees(matched.cat.get("coord.ra"))[good]
        dec = np.degrees(matched.cat.get("coord.dec"))[good]

        xvec, yvec = ra, dec
    else:
        if len(dataKeys) == 3:
            xvec = mags[k1] - mags[k2]
            yvec = mags[k2] - mags[k3]
        else:
            xvec = mags[k1] - mags[k2]
            yvec = mags[k3] - mags[k4]

            if deltaColor:
                show_r_xy = False
                yvec -= xvec
                xvec = mags[k1]
    try:
        alpha.keys()
    except AttributeError:
        alpha = dict(g=alpha, s=alpha)

    idsForColor = ids[idColorN]

    splitId = data.mapperInfo.splitId(idsForColor, asDict=True)
    ccds = splitId["ccd"]
    visits = splitId["visit"]
    #
    # Are we dealing with multi-band or multi-epoch data?
    #
    multiEpoch = False if len(set(filterNames)) == len(filterNames) else True

    nobj = 0
    for c, l, ptype, markersize, color in [("g", nonStellar, "h",     markersize, color),
                                           ("s", stellar,    "*", 1.5*markersize, "green",)]:
        if c not in SG.lower():
            continue

        if colorCcds or colorVisits:
            colors = ["red", "green", "blue", "cyan", "magenta", "yellow", "black", "orange", "brown", "gray"]
            
            if colorCcds:
                for i, ccd in enumerate(sorted(set(ccds))):
                    color = colors[i%len(colors)]
                    tmp = (ccd == ccds)

                    axes.plot(xvec[np.logical_and(l, tmp)],
                              yvec[np.logical_and(l, tmp)],
                              ptype, alpha=alpha[c], markersize=markersize, markeredgewidth=0,
                              markeredgecolor=color, color=color, label=str(ccd))
            else:
                for i, visit in enumerate(sorted(set(visits))):
                    color = colors[i%len(colors)]
                    tmp = (visit == visits)

                    axes.plot(xvec[np.logical_and(l, tmp)],
                              yvec[np.logical_and(l, tmp)],
                              ptype, alpha=alpha[c], markersize=markersize, markeredgewidth=0,
                              color=color, label=str(visit))
            if nobj == 0:
                axes.legend(loc="upper left", numpoints=1,
                            ncol=2, columnspacing=0,
                            markerscale=2,
                            borderpad=0.1, labelspacing=0, handletextpad=0, handlelength=1, borderaxespad=0,
                            )
        else:
            axes.plot(xvec[l], yvec[l], ptype,
                      alpha=alpha[c], markersize=markersize, markeredgewidth=0, color=color)

        nobj += np.sum(l)

    #axes.plot((0, 30), (0, 0), "b-")

    if stellarLocusEnds and "s" in SG.lower():
        plotStellarLocus(axes, mags, k1, k2, k3, stellar, filterNames,
                         stellarLocusEnds, locusLtype, usePrincipalColor, plotRaDec)
        
    if plotRaDec:
        if xmin is not None and xmax is not None:
            axes.set_xlim(xmin, xmax)
        if ymin is not None and ymax is not None:
            axes.set_ylim(ymin, ymax)
    elif multiEpoch:
        axes.set_xlim(-0.3 if xmin is None else xmin, 0.3 if xmax is None else xmax)
        axes.set_ylim(-0.3 if ymin is None else ymin, 0.3 if ymax is None else ymax)
    else:
        axes.set_xlim(-1 if xmin is None else xmin, 2 if xmax is None else xmax)
        axes.set_ylim(-1 if ymin is None else ymin, 2 if ymax is None else ymax)

    if xmax is not None and xmin is not None and ymax is not None and ymin is not None and \
            abs(xmax - xmin) == abs(ymax - ymin):
        axes.set_aspect('equal')
        
    if showStatistics:
        mean, stdev = [], []
        delta = []
        for v in (xvec, yvec,):
            l = True
            if "g" not in SG.lower():
                l = stellar
            if "s" not in SG.lower():
                l = nonStellar

            if deltaColor:
                l = np.logical_and(l, np.logical_and(xvec >= magmin, xvec < magmax))

            delta.append(np.array(v[l], dtype="float64"))
            try:
                stats = afwMath.makeStatistics(delta[-1], afwMath.STDEVCLIP | afwMath.MEANCLIP)
                mean.append(stats.getValue(afwMath.MEANCLIP))
                stdev.append(stats.getValue(afwMath.STDEVCLIP))
            except:
                mean.append(np.nan); stdev.append(np.nan)
        #
        # now the covariance
        #
        if deltaColor:
            if mean is not None:
                for m in (magmin, magmax,):
                    axes.axvline(m, linestyle=":", color="black")

            msg = r"$%.3f \pm %.3f$" % (mean[1], stdev[1])
        else:
            delta = (delta[0] - mean[0])*(delta[1] - mean[1])
            try:
                stdev.append(afwMath.makeStatistics(delta, afwMath.MEANCLIP).getValue())
            except:
                stdev.append(np.nan)

            ax = afwGeom.ellipses.Axes(afwGeom.ellipses.Quadrupole(stdev[0]**2, stdev[1]**2, stdev[2]))
            theta, A, B = ax.getTheta(), ax.getA(), ax.getB()
            if show_r_xy:
                r_xy = stdev[2]/(stdev[0]*stdev[1]) # Pearson's correlation coefficient

            if True:
                from matplotlib.patches import Ellipse
                sig1, sig2 = 1.5151729039613389, 2.4859755240637766 # 68.3% and 95.4% contours
                for nSig in [sig1, sig2]:
                    ell = Ellipse(xy=mean, width=2*nSig*A, height=2*nSig*B, angle=np.rad2deg(theta), zorder=10)
                    axes.add_artist(ell)
                    ell.set_clip_box(axes.bbox)
                    if not True:
                        ell.set_fill(False)
                    else:
                        ell.set_alpha(0.2)
                        ell.set_facecolor("black")
            else:
                _a = np.linspace(0, 2*np.pi, 100)
                ct, st = np.cos(theta), np.sin(theta)
                _x, _y = A*np.cos(_a), B*np.sin(_a)
                axes.plot(mean[0] + _x*ct - _y*st, mean[1] + _x*st + _y*ct, "k")
            #
            # Done with covariance
            #
            msg = "$%.3f, %.3f$\n$\\pm %.3f, %.3f$" % (mean[0], mean[1], stdev[0], stdev[1])
            axes.plot(mean[0], mean[1], "k+", markersize=10)

        axes.text(0.5, 0.70 if show_r_xy else 0.85, msg,
                  fontsize="larger", ha="center", transform = axes.transAxes)

        axes.axvline(0, color="black", ls=":")
        axes.axhline(0, color="black", ls=":")
    else:
        show_r_xy = False

    if plotRaDec:
        xlabel = r"$\alpha$"
        ylabel = r"$\delta$"
    else:
        if len(dataKeys) == 3:
            lk1, lk2, lk3, lk4 = k1, k2, k2, k3
        else:
            lk1, lk2, lk3, lk4 = k1, k2, k3, k4

        if multiEpoch:
            if deltaColor:
                xlabel = "%s (%s)" % (filterNames[lk1], visitNames[lk1])
                ylabel = r"$\Delta %s - %s$  [(%s-%s) - (%s-%s)]" % (filterNames[lk1], filterNames[lk2],
                                                                     visitNames[lk1], visitNames[lk2],
                                                                     visitNames[lk3], visitNames[lk4])
            else:
                xlabel = ("%s" % filterNames[lk1]) if filterNames[lk1] == filterNames[lk2] else \
                    ("%s - %s" % (filterNames[lk1], filterNames[lk2]))
                ylabel = ("%s" % filterNames[lk3]) if filterNames[lk3] == filterNames[lk4] else \
                    ("%s - %s" % (filterNames[lk3], filterNames[lk4]))

                xlabel += " [%s - %s]" % (visitNames[lk1], visitNames[lk2])
                ylabel += " [%s - %s]" % (visitNames[lk3], visitNames[lk4])

        else:
            xlabel = "%s - %s" % (filterNames[lk1], filterNames[lk2])
            ylabel = "%s - %s" % (filterNames[lk3], filterNames[lk4])
        
    axes.text(0.05, 0.1, magType + (r"  $r_{xy}=%.2f$" % r_xy if show_r_xy else ""),
              #fontsize="larger",
              ha="left", transform = axes.transAxes)
    if idN != idColorN:
        axes.text(0.80, 0.1, "legend: %s" % idColorN,
                  fontsize="larger", ha="center", transform = axes.transAxes)

    if showXlabel is not None:
        axes.set_xlabel(xlabel)
        axes.xaxis.set_label_position(showXlabel)
    if showYlabel is not None:
        axes.set_ylabel(ylabel)
        axes.yaxis.set_label_position(showYlabel)

    if title is not None:
        if re.search(r"^T:", title):
            title = title[2:]
            titlePos = "figure"
        else:
            titlePos = "axes"
            
        if not deltaColor and magmax is not None:
            if magmin is not None:
                title += " [%g < %s < %g]" % (magmin,
                                              data.mapperInfo.canonicalFiltername(filterNames[idN]), magmax)
            else:
                title += " [%s < %g]" % (filterNames[idN], magmax)

        title += " %d objects" % nobj
        if data._rerun:
            title += " rerun %s" % data._rerun

        title = re.sub(r"^\+\s*", datasetName + "\n", title)
        if titlePos == "axes":
            axes.set_title(title)
        else:
            fig.suptitle(title)
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    global eventHandlers
    flags = {}
    if False:
        for k, v in data.flags.items():
            flags[k] = data.flags[k][good] # needs to be converted to use data.cat

    if len(xc):
        xc = xc[good]
        yc = yc[good]

    canonicalIds = ids[idN]
    if len(canonicalIds):
        did = data.mapperInfo.splitId(canonicalIds[0], asDict=True); del did["objId"]
        md = data.getDataset("calexp_md", did)[0]
        xc += md.get("LTV1")
        yc += md.get("LTV2")

    if "g" not in SG.lower():
        xvec[nonStellar] = -1000
    if "s" not in SG.lower():
        xvec[stellar] = -1000
    
    eventHandlers[fig] = EventHandler(data, axes, xvec, yvec, canonicalIds, xc, yc, flags,
                                      frames=frames, wcss=wcss, selectWcs=idN-1)

    if not subplot:                     # they didn't pass in an Axes object
        show = False

    if show:
        fig.show()

    return fig, matched

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotSizes(data, magType="psf", selectObjId=None,
              xmin=None, xmax=None, ymin=None, ymax=None,
              title="+", markersize=1, fig=None, frames=[0]):
    """
If selectObjId is provided, it's a function that returns True or False for each object. E.g.
    sel = makeSelectCcd(ccd=2)
    plotSizes(..., selectObjId=utils.makeSelectCcd(2), ...)
"""
    
    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)
    if selectObjId:
        ids = data.cat.get("id")
        good = np.logical_and(good, selectObjId(data, data.cat.get("id")))
        ids = ids[good]
        
    fig = getMpFigure(fig)
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    stellar = data.cat.get("classification.extendedness")[good] < 0.5
    nonStellar = np.logical_not(stellar)

    mag = data.getMagsByType(magType, good)

    Ixx = data.cat.get("shape.sdss.xx")[good]
    Ixy = data.cat.get("shape.sdss.xy")[good]
    Iyy = data.cat.get("shape.sdss.yy")[good]

    size = np.empty(len(Ixx))
    for i in range(len(size)):
        q = afwGeom.ellipses.Quadrupole(Ixx[i], Iyy[i], Ixy[i])
        a = afwGeom.ellipses.Axes(q) # convert to (a, b, theta)

        size[i] = math.sqrt(a.getA()*a.getB())

    color = "red"
    color2 = "green"
    axes.plot(mag[nonStellar], size[nonStellar], "o", markersize=markersize, markeredgewidth=0, color=color)
    axes.plot(mag[stellar], size[stellar], "o", markersize=markersize, markeredgewidth=0, color=color2)
    #axes.plot(mag[bad], size[bad], "+", markersize=markersize, markeredgewidth=0, color="red")
    #axes.plot((0, 30), (0, 0), "b-")

    axes.set_xlim(14 if xmin is None else xmin, 26 if xmax is None else xmax)
    axes.set_ylim(0.0 if ymin is None else ymin, 10 if ymax is None else ymax)
    axes.set_xlabel(magType)
    axes.set_ylabel(r"$\sqrt{a b}$ (pixels)", fontsize="larger")
    name = data.name if selectObjId is None else \
        data.mapperInfo.dataIdToTitle(idsToDataIds(data, ids), data._rerun)
    axes.set_title(re.sub(r"^\+\s*", name + " ", title))

    global eventHandlers
    flags = {}
    try:
        x = data.cat.getX()[good] # + md.get("LTV1")
        y = data.cat.getY()[good] # + md.get("LTV2")
    except pexExcept.LsstCppException, e:
        if not re.search(r"pex::exceptions::LogicErrorException", e.message.getType()):
            raise e
        x = np.zeros_like(mag)
        y = np.zeros_like(mag)
        
    eventHandlers[fig] = EventHandler(data, axes, mag, size, ids, x, y, flags, frames=frames)

    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotDmagHistograms(data, magType1="model", magType2="psf",
                       nbin=15, magmin=15, magmax=20, dmag=0.5,
                       xmin=-0.1, xmax=None, color="red", fig=None):
    """Plot (magType1 - magType2) v. magType2 mags"""

    try:
        data.cat
    except AttributeError:
        raise RuntimeError("Please call da.getMagsByVisit, then try again")

    bad = reduce(lambda x, y: np.logical_or(x, data.cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  #"flags.pixel.interpolated.center",
                  "flags.pixel.saturated.center",],
                 False)
    good = np.logical_not(bad)

    fig = getMpFigure(fig)

    mag1 = data.getMagsByType(magType1, good)
    mag2 = data.getMagsByType(magType2, good)
    delta = np.array(mag1 - mag2, dtype='float64') # float64 is needed by makeStatistics --- why?

    stellar = data.cat.get("classification.extendedness")[good] < 0.5
    #
    # Get a characteristic (we hope) Calib object
    #
    if False:
        calexp_md = data.butler.get(dtName("calexp", True), data.dataSets[0])
        calib = afwImage.Calib(calexp_md)
        psf = None                      # can't unpersist the PSF without the pixels; #2777
    else:
        calexp = data.butler.get(dtName("calexp"), data.dataSets[0])
        calib = calexp.getCalib()
        psf = calexp.getPsf()
        alpha = psf.computeShape().getDeterminantRadius() # psf ~ N(0, alpha^2)
        print "# alpha = %.3f" % (alpha)

    if xmax is None:
        xmax = -xmin
    binEdge = np.arange(magmin, magmax - dmag/2, dmag)

    subplots = makeSubplots(fig, nx=1, ny=len(binEdge))
    fig.subplots_adjust(wspace=0.001, top=0.875)

    for edge in binEdge:
        inBin = np.logical_and(np.logical_and(mag1 >= edge, mag1 < edge + dmag), stellar)
        inBin = np.logical_and(np.logical_and(delta >= xmin, delta < xmax), inBin)

        pp = mag2[inBin]        
        dd = delta[inBin]

        if len(dd) == 0:
            continue

        stats = afwMath.makeStatistics(dd, afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)
        
        axes = subplots.next()
        axes.xaxis.set_major_formatter(pyplot.NullFormatter())
        axes.yaxis.set_major_formatter(pyplot.NullFormatter())

        hist, bins = np.histogram(dd, bins=nbin); hist = np.array(hist, dtype='d')
        hist /= max(hist)
        width =  bins[1] - bins[0]
        center = (bins[:-1] + bins[1:])/2
        axes.bar(center, hist, align='center', width=width, linewidth=0, color=color)

        axes.axvline(0.0, color="blue", ls=":")
        axes.axvline(mean, color="green", ls=":")

        axes.set_xlim(xmin, xmax)
        axes.set_ylim(0, 1.05)

        meanMag = edge + 0.5*dmag
        axes.set_xlabel("%.2f" % meanMag)
        axes.set_title("%.3f\n$\pm$%.3f" % (mean, stdev), fontsize="smaller") 
        #
        # Put out information for intensity-dependent PSFs
        #
        print "%.4f %6.2e %.4f %.4f" % (meanMag, calib.getFlux(meanMag), mean, stdev)
        
    fig.suptitle("%s - %s: %s" % (magType1, magType2, data.name))

    fig.show()

    return fig

def getRefmag(data, matchedCat, desiredBand, referenceSuffix=""):
    """Get the reference magnitude from a catalog objects"""

    if len(matchedCat) == 0:
        return np.zeros(0)

    sch = matchedCat.getSchema()

    ct = data.mapperInfo.getColorterm(desiredBand)
    if False:
        print "RHL Not applying colour terms"
        ct = None

    if ct:
        primaryKey = sch.find(ct.primary + referenceSuffix).getKey()
        secondaryKey = sch.find(ct.secondary + referenceSuffix).getKey()

        refMag = -2.5*np.log10(matchedCat.get(primaryKey))
        secMag = -2.5*np.log10(matchedCat.get(secondaryKey))

        refMag = data.mapperInfo.photometricTransform(desiredBand, refMag, secMag)
    else:
        fluxKey = sch.find("flux" + referenceSuffix).getKey()
        refMag = -2.5*np.log10(matchedCat.get(fluxKey))

    return refMag

def getCanonicalFilterName(filterName):
    filterId = afwImage.Filter(filterName).getId()

    for name in afwImage.Filter.getNames():
        if afwImage.Filter(name).getId() == filterId:
            return name

    return filterName

def plotCalibration(data, plotBand=0.05, magType='psf', magmin=14, maglim=20, eventHandler=False,
                    selectObjId=None, showCamera=False, correctRadial=False, plotRadial=False,
                    showJacobian=False, correctJacobian=False, useDistortion=False, showZP=False,
                    xmin=None, xmax=None, ymin=None, ymax=None, meanDelta=0.0,
                    markersize=2, alpha=1.0, title="+", showMedians=False, overlay=True,
                    mtv=True, frame=None, ctype=None, ds9Size=None, fig=None):
    """Plot (instrumental - reference) v. reference magnitudes given a Data object.

The Data must have been initialised with a call to the getCalibObjects method

If selectObjId is provided, it's a function that returns True or False for each object. E.g.
    plotCalibration(..., selectObjId=makeSelectCcd(2), ...)

Use the matplotlib Figure object fig; if none is provided it'll be created and saved in data (and if it's true, a new one will be created)

If showZP is True, the data plotted is the per-CCD zeropoint

if eventHandler is true, make the points in the plot live
    
If plotBand is provided, draw lines at +- plotBand

If title is provided it's used as a plot title; if it starts + the usual title is prepended

If overlay is true, overlay the CCDs' boundaries and IDs

If frame is not None plot the objects in ds9 (if mtv is True as well, display the calexp first)
"""
    if not data.matchedCalibs:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    sourceSuffix = "_2"
    referenceSuffix="_1"
    good = np.logical_and(data.matchedCalibs["stargal" + referenceSuffix],
                          np.logical_not(data.matchedCalibs["flags.pixel.saturated.center" + sourceSuffix]))
    good = np.logical_and(good,
                          np.logical_not(data.matchedCalibs["flags.pixel.edge" + sourceSuffix]))
    if selectObjId:
        good = np.logical_and(good, selectObjId(data, data.matchedCalibs.get("id")))
    magTypeName = "sinc" if magType == "ap" else "deconvolvedPsf" if magType == "inst" else magType
    flux = data.matchedCalibs.get("flux.%s%s" % (magTypeName, sourceSuffix))[good]
    
    ids = data.matchedCalibs.get("id")[good]
    xc = data.matchedCalibs.get("centroid.sdss%s.x" % sourceSuffix)[good]
    yc = data.matchedCalibs.get("centroid.sdss%s.y" % sourceSuffix)[good]
    #
    # Get focal-plane positions ("mm" -- for HSC they are actually in units of 15 micron pixels)
    #
    camera = data.butler.get("camera")
    xmm = np.empty_like(xc); ymm = np.empty_like(yc)
    for i, _id in enumerate(ids):
        ccdId = data.mapperInfo.splitId(_id, True)["ccd"]
        ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(ccdId))
        xmm[i], ymm[i] = ccd.getPositionFromPixel(afwGeom.PointD(xc[i], yc[i])).getMm()

    if correctJacobian or showJacobian:
        if useDistortion:
            title += " (Jacobian corr from Distort)"
            jacobian = np.empty_like(xmm)
            dist = camera.getDistortion()
            for i in range(len(jacobian)):
                jacobian[i] = dist.computeQuadrupoleTransform(afwGeom.PointD(xmm[i], ymm[i]), False).computeDeterminant()
        else:
            title += " (Jacobian corr)"
            jacobian = data.matchedCalibs.get("jacobian%s" % sourceSuffix)[good]

        if correctJacobian:
            perChipJacobian = True      # set the average Jacobian to 1.0 for every CCD
            
        if perChipJacobian:
            ccdId = data.mapperInfo.splitId(ids, True)["ccd"]
            for _id in set(ccdId):
                thisCcd = np.where(ccdId == _id)
                jacobian[thisCcd] /= np.median(jacobian[thisCcd])

        flux *= jacobian

    zp = data.zp[good]
    instmag = zp - 2.5*np.log10(flux)

    if frame is not None:
        if mtv:
            ccdId = data.mapperInfo.splitId(ids[0], True)
            ds9.mtv(data.butler.get(dtName("calexp"), visit=ccdId["visit"], ccd=ccdId["ccd"]), frame=frame)

        kwargs = {}
        if ctype:
            kwargs["ctype"] = ctype
        if ds9Size:
            kwargs["size"] = ds9Size
        with ds9.Buffering():
            for i in range(len(xc)):
                if instmag[i] < 21:        # XXXXXXXX
                    ds9.dot("o", xc[i], yc[i], frame=frame, **kwargs)

    fig = getMpFigure(fig)
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    did = data.mapperInfo.splitId(ids[0], True)
    refmag = getRefmag(data, data.matchedCalibs,
                       getCanonicalFilterName(data.butler.queryMetadata('raw', 'filter',
                                                                        visit=did['visit'])[0]),
                       referenceSuffix=referenceSuffix)[good]

    if correctRadial:               # apply a radial correction
        if correctJacobian:
            if useDistortion:
                x0, y0 = 3000, 0               # center of correction pattern
                r = np.hypot(xmm - x0, ymm - y0)

                instmag += -0.009 - 0.045*(r/15000)**2
            else:
                x0, y0 = 3000, 0               # center of correction pattern
                r = np.hypot(xmm - x0, ymm - y0)

                instmag +=  0.005 - 0.050*(r/15000)**2
        else:
            x0, y0 = -2000, 0               # center of correction pattern
            r = np.hypot(xmm - x0, ymm - y0)

            instmag += -0.010 + 0.065*(r/15000)**2
        title += " (radial correction)"

    if showJacobian:
        delta = jacobian
        title = "+ Jacobian"
    elif showZP:
        delta = zp
        title += " Zero Points"
        if correctJacobian:
            delta += 2.5*np.log10(jacobian)

    else:
        delta = refmag - instmag

    try:
        stats = afwMath.makeStatistics(delta[np.logical_and(refmag >= magmin, refmag <= maglim)] - meanDelta,
                                       afwMath.STDEVCLIP | afwMath.MEANCLIP)
        mean, stdev = stats.getValue(afwMath.MEANCLIP) + meanDelta, stats.getValue(afwMath.STDEVCLIP)
    except Exception, e:
        print "Failed to estimate mean: %s" % e
        mean, stdev = float("NaN"), float("NaN")
        
    if showCamera or plotRadial:
        inRange = np.logical_and(refmag >= magmin, refmag <= maglim)

        if plotRadial:
            r = np.hypot(xmm, ymm)
            axes.plot(r[inRange], delta[inRange], "k.", markersize=markersize, alpha=alpha, markeredgewidth=0)

            axes.set_ylim(-0.6 if ymin is None else ymin,
                           0.6 if ymax is None else -ymin if ymax is None else ymax)
            axes.axhline(0.0, linestyle=":", color="cyan")
            axes.axhline(mean, linestyle=":", color="red")
            axes.set_xlabel("Radius (mm)")
            axes.set_ylabel("Reference - %s" % magType)
        else:
            sc = axes.scatter(xmm[inRange], ymm[inRange], c=delta[inRange], norm=pyplot.Normalize(ymin, ymax),
                              cmap=pyplot.cm.rainbow, marker="o", s=10*markersize,
                              edgecolors="none", alpha=alpha)
            fig.colorbar(sc)
            axes.set_aspect('equal')
            
            axes.set_xlabel("X (mm)")
            axes.set_ylabel("Y (mm)")
            
            title += " %s" % magType

            if overlay:
                overlayCcds(camera, axes) 
    else:
        axes.plot(refmag, delta, "k.", markersize=markersize, alpha=alpha, markeredgewidth=0)

        for m in (magmin, maglim,):
            axes.axvline(m, color="blue", ls=":")

        axes.axhline(mean, color="green", linestyle="-")
        if plotBand:
            for i in (-plotBand, plotBand):
                axes.axhline(mean + i, color="green", linestyle="--")

        axes.set_xlim(magmin - 0.5 if xmin is None else xmin,
                      26  if xmax is None else xmax)
        axes.set_ylim(-0.6         if ymin is None else ymin,
                      0.6 if ymax is None else -ymin if ymax is None else ymax)
        axes.set_xlabel("Reference")
        axes.set_ylabel("Reference - %s" % magType)

        if showMedians:
            binwidth = 1.0
            bins = np.arange(int(min(refmag)), int(max(refmag)), binwidth)
            vals = np.empty_like(bins)
            for i in range(len(bins) - 1):
                inBin = np.logical_and(refmag > bins[i], refmag <= bins[i] + binwidth)
                vals[i] = np.median(delta[np.where(inBin)])

            axes.plot(bins + 0.5*binwidth, vals, linestyle="-", marker="o", color="cyan")

    fig.text(0.2, 0.85, r"$%.3f \pm %.3f$" % (mean, stdev), fontsize="larger")

    name = data.name if selectObjId is None else \
        data.mapperInfo.dataIdToTitle(idsToDataIds(data, ids), data._rerun)

    axes.set_title(re.sub(r"^\+\s*", name + " ", title))
    #
    # Make "i" print the object's ID, p pan ds9, etc.
    #
    if eventHandler:
        global eventHandlers
        flags ={}
        eventHandlers[fig] = EventHandler(data, axes,
                                          xmm if showCamera else refmag,
                                          ymm if showCamera else delta, ids, xc, yc, flags, [frame,])

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
    if not data.matchedCalibs:
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

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

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
        eventHandlers[fig] = EventHandler(data, axes, refmag, delta, ids, srcx, srcy, flags, [frame,])

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
    if not data.matchedCalibs:
        raise RuntimeError("You need to call the Data.getCalibObjects method before calling this routine")

    raise RuntimeError("Convert me to zipped matches")

    with ds9.Buffering():
        for m in data.matchedCalibs:
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

def plotZeroPoints(data, selectObjId=None,
                   correctJacobian=False,
                   ymin=None, ymax=None, title="+", markersize=1, fig=None):
    """Plot the per-CCD zeropoints (only for visit == visit if non-None)

If selectObjId is provided, it's a function that returns True or False for each CCD. E.g.
    plotZeroPoints(..., selectObjId=makeSelectVisit(visit=902040), ...)   

    """
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    dataSets = data.dataSets
    if selectObjId:
        dataSets = [did for did in dataSets if selectObjId(data, did)]
    #
    # Get focal-plane positions ("mm" -- for HSC they are actually in units of 15 micron pixels)
    #
    camera = data.butler.get("camera")
    dist = camera.getDistortion()

    xmm = np.empty(len(dataSets)); ymm = np.empty_like(xmm)
    zp = np.empty_like(xmm); jacobian = np.empty_like(xmm)
    for i, did in enumerate(dataSets):
        calexp_md = data.butler.get(dtName("calexp", True), **did)
        w, h = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")

        ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(did["ccd"]))
        xmm[i], ymm[i] = ccd.getPositionFromPixel(afwGeom.PointD(0.5*w, 0.5*h)).getMm()

        calib = afwImage.Calib(calexp_md)
        zp[i] = calib.getMagnitude(1.0)

        jacobian[i] = dist.computeQuadrupoleTransform(afwGeom.PointD(xmm[i], ymm[i]),
                                                      False).computeDeterminant()

    if correctJacobian:
        zp += 2.5*np.log10(jacobian)
        title += " (Jacobian corr from Distort)"

    if len(xmm) == 0:
        raise RuntimeError("I'm afraid that I have no points to plot")

    sc = axes.scatter(xmm, ymm, c=zp, norm=pyplot.Normalize(ymin, ymax),
                      cmap=pyplot.cm.rainbow, marker="o", s=10*markersize, edgecolors="none")
    fig.colorbar(sc)
    axes.set_aspect('equal')

    name = data.name if selectObjId is None else data.mapperInfo.dataIdToTitle(dataSets, data._rerun)
    axes.set_title(re.sub(r"^\+\s*", "Zeropoints " + name + " ", title))
    
    #
    # Make "z" print the z-axis (zp)
    #
    global eventHandlers
    eventHandlers[fig] = EventHandler(data, axes, xmm, ymm, zp)
    fig.show()

    return fig

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
    
    calexp_md = data.butler.get(dtName("calexp", True), **dataId)
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

def drawKron(s, xy0=(0, 0), **kwargs):
    x, y = s.getCentroid()
    x -= xy0[0]; y -= xy0[1]

    if s.get("flux.kron.flags"): # failed; draw an \odot
        kw = kwargs.copy()
        kw["size"] = 5; kw["ctype"] = ds9.GREEN
        ds9.dot('x', x, y, **kw)
        ds9.dot('o', x, y, **kw)
        return

    shape = s.getShape()
    R_K = s.get("flux.kron.radius")
    if False:
        if math.isnan(R_K):
            import pdb; pdb.set_trace() 
        assert not math.isnan(R_K)

    nRadiusForFlux = 2.5
    rDet = shape.getDeterminantRadius()
    for r, ct in [(R_K, ds9.BLUE), (R_K*nRadiusForFlux, ds9.MAGENTA),
                  (s.get("flux.kron.radiusForRadius"), ds9.GREEN)]:
        shape = shape.clone()
        shape.scale(r/shape.getDeterminantRadius())

        try:
            ds9.dot(shape, x, y, ctype=ct, silent=True, **kwargs) # requires a ds9.py >= 2013-01-15
        except:
            pass

def kronEventCallback(key, source, im, frame):
    """Callback for event handlers; e.g. utils.EventHandler.callbacks['d'] = kronEventCallback"""
    drawKron(source, im.getXY0(), frame=frame)

def showSourceSet(sourceSet, exp=None, wcs=None, xy0=None, raDec=None, magmin=None, magmax=None, magType="psf",
                  nSource=-1, SG=False, deblend=True, obeyXY0=True, mapperInfo=None,
                  mask=None, symb="+", **kwargs):
    """Show a SourceSet on ds9.

    If nSource > 0, only show the nSource objects with the brightest psf flux
    If mask, it's a set of bitplane names (e.g. INTERP) which must be set to display the source"""
    
    try:
        sourceSet.getColumnView()
    except:
        sourceSet = sourceSet.copy(True)

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
    try:
        isStar = (sourceSet.get("classification.extendedness") < 0.5)
        isStar = np.logical_or(isStar, sourceSet.get("flags.pixel.saturated.center"))
    except KeyError:
        isStar = True

    try:
        sourceSet.get("parent")
    except:
        deblend = False

    if xy0 is None:
        x0, y0 = exp.getXY0() if exp else [0.0, 0.0]
    else:
        x0, y0 = xy0

    mat = re.search(r"^R(\d+(\.\d+)?)", symb)
    if mat:
        symb = "o"
        kwargs["size"] = float(mat.group(1))
    
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

                _id = s.getId()
                if symb == "id":
                    if mapperInfo:
                        _symb = mapperInfo.splitId(_id, asDict=True)["objId"]
                    else:
                        _id = _id & 0xffff # guess wildly
                else:
                    _symb = "%d" % _id
                if deblend:
                    kwargs["ctype"] = ds9.RED if s.get("parent") == 0 else ds9.MAGENTA                    

            elif symb == "@":
                _symb = s.getShape()    # requires a ds9.py >= 2013-01-15
            elif symb.lower() == "kron":
                if not deblend and s.get("deblend.nchild") > 0:
                    continue

                drawKron(s, xy0=(x0, y0) if obeyXY0 else (0, 0), **kwargs)
                continue
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

    data.name = data.mapperInfo.dataIdToTitle(dataId, self._rerun)

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

def subtractModels(data, dataId, magType="model", frame=0, subtractedOnly=False, showExposure=True, ccd=None,
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
    smi = subtracted.getMaskedImage()
    si = smi.getImage()

    for s in data.getSources(dataId)[0]:
        if getFlagForDetection(s, "SATUR_CENTER"):
            continue

        if s.get("deblend.nchild") > 0: # subtract the children, not the parent
            continue

        if magType == "model":
            try:
                comboImage, expImage, devImage = msViewer.makeGalaxyImages(s, imageBBox=exp.getBBox())
            except ValueError, e:
                continue
            
            try:
                if True:
                    modelImage = comboImage
                elif True:
                    modelImage = devImage
                else:
                    modelImage = expImage

                sub = si.Factory(si, modelImage.getBBox(afwImage.PARENT))
                sub -= modelImage
                del sub
            except pexExcept.LsstCppException, e:
                #print e
                pass

            continue

        flux = np.nan if fitAmplitude else getFlux(s, magType)
        try:
            if s.get("classification.extendedness") < 0.5:
                measAlg.subtractPsf(psf, smi, s.getX(), s.getY(), flux)
        except pexExcept.LsstCppException, e:
            pass

    title = data.mapperInfo.dataIdToTitle([dataId], self._rerun)

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

    si -= exp.getMaskedImage().getImage()
    si *= -1

    return subtracted    

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
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    axes.plot(z, I, "o", color="b")
    axes.plot(z, 0.5*maxVal + 0*I, "r--")
    axes.set_xlim(min(z) - 1, max(z) + 1)
    axes.set_ylim(-0.01*max(I), 1.05*max(I))
    axes.set_title("%s (%d, %d)" % (data.mapperInfo.exposureToStr(exp), xc, yc))

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

    fig.suptitle("PSF eigen components, %s [%s]" % (data.mapperInfo.exposureToStr(exposure), rerun),
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
    elif magType == "kron":
        return s.get("flux.kron")
    elif magType == "model":
        return s.getModelFlux()
    elif magType == "psf":
        return s.getPsfFlux()
    else:
        raise RuntimeError("Unknown magnitude type %s" % magType)

def subtractBackground(image, nsigma=5, boxSize=1024, method="fit"):
    if method == "min":
        arr = image.getArray().flatten()
        arr.sort()
        bkgd = float(arr[len(arr)//100])
        image -= bkgd
    elif method == "plane":
        arr = image[0:1000, 0:1000].getArray().flatten()
        arr.sort()
        valLL = float(arr[len(arr)//100])

        arr = image[-1000:-1, -1000:-1].getArray().flatten()
        arr.sort()
        valUR = float(arr[len(arr)//100])

        width, height = image.getDimensions()
        xx = np.arange(width)
        for y in range(height):
            image.getArray()[y, xx] -= valLL + 0.5*(y/(height - 1.0) + xx/(width - 1.0))*(valUR - valLL)

        bkgd = None
    elif method == "fit":
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

        try:
            image.getImage
        except AttributeError:
            image = afwImage.makeMaskedImage(image)
            
        exp = afwImage.makeExposure(image)
        bkgd = sdTask.detectFootprints(exp, sigma=1).background
        del exp
    else:
        raise RuntimeError("I don't know how to subtract a background using method %s" % method)

    return bkgd

def writeRgb(images, rgbFile, min=0, max=50, Q=8, bin=1, scales=[1.0, 1.0, 1.0],
             fixSaturation=False, saveFmt=None,
             subtractBkgd=False, boxSize=1024, nsigma=5):
    """Convert the list of images to a true-colour image, and write it to rgbFile (currently .png or .tiff)
The order is R, G, B (e.g. i, r, g)

Scale the images by scales, if provided

The default stretch is an asinh stretch; Q == 0 corresponds to a linear.  If Q == 0 then min and max are
the range of the stretch; as Q increases they still define the slope of the linear portion of the stretch

If bin is specified, it should be an integer > 1;  the output file will be binned by that factor.

If saveFmt is non-None it's taken to be a string with a single %s which will be replaced by "R", "G", and "B"
when writing the final images to disk
    """
    
    if not afwRgb:
        raise RuntimeError("I was unable to import lsst.afw.extensions.rgb")

    R, G, B = 0, 1, 2

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
        #
        # Handle Exposures
        #
        try:
            image = image.getMaskedImage()    # maybe it's an Exposure
        except AttributeError:
            pass

        images[i] = image

    if fixSaturation:
        bctrl = afwMath.BackgroundControl(image.getWidth()//boxSize, image.getHeight()//boxSize)

        brightLevel = 130
        for i, image in enumerate(images):
            image = image.clone()
            images[i] = image

            backobj = afwMath.makeBackground(image.getImage(), bctrl)
            bkgd = backobj.getImageF(afwMath.Interpolate.AKIMA_SPLINE, afwMath.REDUCE_INTERP_ORDER)

            image.getImage()[:] -= bkgd

            img = image.getImage().getArray()
            msk = image.getMask().getArray()
            SAT = image.getMask().getPlaneBitMask(["SAT", "INTRP"])

            msk[np.where(img < brightLevel)] &= ~SAT

            image += bkgd

        saturValue = np.nan if True else 2000
        afwRgb.replaceSaturatedPixels(images[R], images[G], images[B], 2, saturValue)

    for i, image in enumerate(images):
        copied = False                    # have we made a copy?

        if fixSaturation:
            copied = True
        #
        # Handle MaskedImages
        #
        try:
            image = image.getImage()          # maybe it's (now) a MaskedImage
        except AttributeError:
            pass

        if bin > 1:
            image = afwMath.binImage(image, bin)
            copied = True

        if scales[i] != 1.0:
            if not copied:
                image = image.clone()
                copied = True

            image *= float(scales[i])

        if subtractBkgd:
            if not copied:
                image = image.clone()
                copied = True

            if isinstance(subtractBkgd, bool):
                subtractBkgd = "fit"
            subtractBackground(image, nsigma, boxSize, subtractBkgd)

        images[i] = image

        if grayScale:
            assert i == 0
            while len(images) < 3:
                images.append(None)
            images[1] = image
            images[2] = image

            break            

    afwRgb.RgbImageF(images[R], images[G], images[B], afwRgb.asinhMappingF(min, max - min, Q)).write(rgbFile)

    if saveFmt:
        images[R].writeFits(saveFmt % "R")
        images[G].writeFits(saveFmt % "G")
        images[B].writeFits(saveFmt % "B")

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
    def __init__(self, data, axes, xs, ys, ids, x=None, y=None, flags={}, frames=[0], wcss=[], selectWcs=0):
        self.data = data
        self.axes = axes
        self.xs = xs
        self.ys = ys
        self.ids = ids
        self.selectWcs = selectWcs                  # which of frames[], wcss[] the ids correspond to
        self.x = x
        self.y = y
        self.flags = flags

        if frames is None:
            frames = []
        else:
            try:
                frames[0]
            except TypeError:
                frames = [frames]
            
        if not wcss:
            wcss = [None]*len(frames)

        self.frames = frames
        self.wcss = wcss

        self.cid = self.axes.figure.canvas.mpl_connect('key_press_event', self)

    def __call__(self, ev):
        if ev.inaxes != self.axes:
            return
        
        if not ev.key:
            return

        if ev.key.lower() == 'h':
            print """
Options:
   d, D    Display the object in ds9, frame utils.eventFrame.  If D, show the deblended child
   h, H    Print this message
   p, P    Pan to the object in ds9 in all the frames specified to the event handler; 'P' also prints a newline
   i, I    Print info about the object; 'I' also prints a newline
   z, Z    Print raw "ids" array about the object; 'Z' also prints a newline

If utils.eventCallbacks[ev.key] is defined it'll be called with arguments:
   ev.key, source, maskedImage, frame
(see utils.kronCallback for an example)
""",
            return
        elif ev.key.lower() in ("dipz"):
            dist = np.hypot(self.xs - ev.xdata, self.ys - ev.ydata)
            dist[np.where(np.isnan(dist))] = 1e30
            dmin = min(dist)

            which = np.where(dist == min(dist))
            objId = self.ids[which][0]

            if ev.key in "zZ":
                print "\r>>>",
                print "%.3f %.3f %g%20s\r" % (self.xs[which][0], self.ys[which][0], objId, ""),
                sys.stdout.flush()
                if ev.key.isupper():
                    print ""
                return

            flagsInfo = []
            for k in sorted(self.flags.keys()):
                if self.flags[k][objId == self.ids]:
                    flagsInfo.append(k)
            flagsInfo = ", ".join(flagsInfo)

            print "\r>>>",
            print "%.3f %.3f %s %s (%6.1f, %6.1f)%20s\r" % \
                (self.xs[which][0], self.ys[which][0], self.data.mapperInfo.splitId(objId), flagsInfo,
                 self.x[which][0], self.y[which][0], ""),
            sys.stdout.flush()

            if ev.key in ('d', 'D'):
                dataId = self.data.mapperInfo.splitId(objId, asDict=True)
                title = re.sub(r"[' {}]", "", str(dataId)).replace(",", " ")
                oid = dataId["objId"]; del dataId["objId"]
                ss = self.data.butler.get(dtName("src"), **dataId)
                i = int(np.where(objId == ss.get("id"))[0][0])
                s = ss[i]

                bbox = s.getFootprint().getBBox()
                grow = 0.5
                bbox.grow(afwGeom.ExtentI(int(grow*bbox.getWidth()), int(grow*bbox.getHeight())))

                md = self.data.butler.get(dtName("calexp_md"), **dataId)
                bbox.clip(afwGeom.BoxI(afwGeom.PointI(0, 0),
                                       afwGeom.ExtentI(md.get("NAXIS1"), md.get("NAXIS2"))))

                if ev.key == 'd':
                    calexp = self.data.butler.get(dtName("calexp_sub"), bbox=bbox, **dataId).getMaskedImage()
                else:
                    import deblender
                    try:
                        calexp = deblender.footprintToImage(s.getFootprint())
                    except:
                        calexp = self.data.butler.get(dtName("calexp"), **dataId).getMaskedImage()
                        calexp = deblender.footprintToImage(s.getFootprint(), calexp)

                ds9.mtv(calexp, title=title, frame=eventFrame)
                ds9.pan(s.getX() - calexp.getX0(), s.getY() - calexp.getY0(), frame=eventFrame)

                callback = eventCallbacks.get(ev.key)
                if callback:
                    callback(ev.key, s, calexp, frame=eventFrame)
            elif ev.key == "I":
                print ""
            elif ev.key in ("p", "P"):
                x = self.x[which][0]
                y = self.y[which][0]

                for frame, wcs in zip(self.frames, self.wcss):
                    if wcs:
                        raDec = self.wcss[self.selectWcs].pixelToSky(x, y)
                        ds9.pan(*wcs.skyToPixel(raDec[0], raDec[1]), frame=frame)
                    else:
                        ds9.pan(x, y, frame=frame)
                ds9.cmdBuffer.flush()
        else:
            pass

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class PsfModelImage(cameraGeomUtils.GetCcdImage):
    """A class to return an Image of a given Ccd based on its cameraGeometry"""
    
    def __init__(self, butler, bin=10, background=np.nan, verbose=False, stampSize=0, *args, **kwargs):
        """Initialise
        gravity  If the image returned by the butler is trimmed (e.g. some of the SuprimeCam CCDs)
                 Specify how to fit the image into the available space; N => align top, W => align left
        background  The value of any pixels that lie outside the CCDs
        """
        super(PsfModelImage, self).__init__(*args)
        self.butler = butler
        self.verbose = verbose
        self.kwargs = kwargs
        self.bin = bin
        self.stampSize = stampSize
        self.imageIsBinned = True       # i.e. images returned by getImage are already binned

        self.isRaw = False
        
        self.gravity = None
        self.background = background

    def getImage(self, ccd, amp=None, imageFactory=None):
        """Return an image of the specified amp in the specified ccd"""

        ccdNum = ccd.getId().getSerial()

        try:
            if self.kwargs.get("ccd") is not None and self.kwargs.get("ccd") != ccdNum:
                raise RuntimeError

            dataId = self.kwargs
            if dataId.has_key("ccd"):
                dataId = self.kwargs.copy()
                del dataId["ccd"]

            bbox = afwGeom.Box2I(afwGeom.Point2I(), afwGeom.Extent2I(1,1))
            psf = self.butler.get("calexp_sub", bbox=bbox,
                                  ccd=ccd.getId().getSerial(), **self.kwargs).getPsf()

            w, h =  ccd.getAllPixels(True).getDimensions()
            binnedW, binnedH = w//self.bin, h/self.bin
            if self.stampSize:
                psfW, psfH = self.stampSize, self.stampSize
            else:
                psfW, psfH = psf.computeImage(afwGeom.PointD(0, 0)).getDimensions()

            nx = binnedW/psfW
            while nx*psfW > binnedW:
                if nx == 1:
                    break
                nx -= 1
            mos = maUtils.showPsfMosaic((w, h), psf, nx=nx,
                                        showCenter=False, showEllipticity=False, showFwhm=False,
                                        stampSize=self.stampSize)
            mos.gutter = 0
            im = mos.makeMosaic(mode=nx)
        except Exception, e:
            if self.verbose and str(e):
                print e

            im = afwImage.ImageF(*ccd.getAllPixels(self.isTrimmed).getDimensions())
            if self.bin:
                im = afwMath.binImage(im, self.bin)

        return im

        
def showPsfModelImage(butler, visit, bin=20, stampSize=0, frame=0, verbose=False, title=None):
    """Show a mosaic of PSF mosaics"""
    cameraImage = cameraGeomUtils.showCamera(butler.get("camera"),
                                    PsfModelImage(butler, visit=visit, verbose=verbose,
                                                  stampSize=stampSize, bin=bin),
                                    frame=frame, bin=bin, title=title)
    if frame is not None:
        xy0 = afwGeom.ExtentI(cameraImage.getXY0())

        if butler.get("camera").getId().getName() == "HSC":
            ds9.dot("o", 0/bin - xy0[0], -136/bin - xy0[1], size=18280/bin,
                    frame=frame, ctype=ds9.RED)

    return cameraImage

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ResidualImage(cameraGeomUtils.GetCcdImage):
    """A class to return an Image of a given Ccd based on its cameraGeometry"""
    
    def __init__(self, butler, magMin=None, magLim=30, magType="psf", apCorr=1.0, showPsfStars=False,
                 showPsfModels=True,
                 bin=10, stampSize=0, background=np.nan, sigma=0, verbose=False, *args, **kwargs):
        """Initialise
        gravity  If the image returned by the butler is trimmed (e.g. some of the SuprimeCam CCDs)
                 Specify how to fit the image into the available space; N => align top, W => align left
        background  The value of any pixels that lie outside the CCDs
        """
        super(ResidualImage, self).__init__(*args)
        self.butler = butler
        self.verbose = verbose
        self.kwargs = kwargs
        self.sigma = sigma
        self.bin = bin
        self.magMin = magMin
        self.magLim = magLim
        self.magType = magType
        self.imageIsBinned = True       # i.e. images returned by getImage are already binned
        self.apCorr = apCorr
        self.showPsfStars = showPsfStars
        self.showPsfModels = showPsfModels
        self.stampSize = stampSize
        self.isRaw = False
        
        self.gravity = None
        self.background = background

    def getImage(self, ccd, amp=None, imageFactory=None):
        """Return an image of the specified amp in the specified ccd"""

        ccdNum = ccd.getId().getSerial()

        try:
            if self.kwargs.get("ccd") is not None and not ccdNum in self.kwargs.get("ccd"):
                raise RuntimeError

            dataId = self.kwargs
            if dataId.has_key("ccd"):
                dataId = self.kwargs.copy()
                del dataId["ccd"]

            calexp = self.butler.get("calexp", ccd=ccdNum, immediate=True, **dataId)
        except Exception, e:
            if self.verbose and str(e):
                print e

            ccdImage = afwImage.ImageF(*ccd.getAllPixels(self.isTrimmed).getDimensions())
            if self.bin:
                ccdImage = afwMath.binImage(ccdImage, self.bin)

            return ccdImage

        ss = self.butler.get("src",  ccd=ccdNum, **dataId)
        
        calib = calexp.getCalib()
        if calib.getFluxMag0()[0] > 0.0:
            with afwImageUtils.CalibNoThrow():
                psfMag = np.array(calib.getMagnitude(ss.getPsfFlux()))
        else:
            print >> sys.stderr, "CCD %d has fluxMag0 <= 0.0; showing all stars" % ccdNum
            psfMag = np.zeros(len(ss)) + 0.9999*self.magLim
        psfWidth, psfHeight = calexp.getPsf().getLocalKernel().getDimensions()

        good = np.logical_not(ss.get("flags.pixel.saturated.any"))
        if self.magLim is not None:
            good = np.logical_and(psfMag < self.magLim, good)
        if self.magMin is not None:
            good = np.logical_and(good, psfMag > self.magMin)
        isStar = ss.get("classification.extendedness") < 0.5
        good = np.logical_and(good, isStar)
        good = np.logical_and(good, ss.get("deblend.nchild") == 0)
        if self.showPsfStars:
            good = np.logical_and(good, ss.get("calib.psf.used"))

        if self.showPsfModels:
            stampSize = self.stampSize

            w, h = calexp.getWidth()//self.bin, calexp.getHeight()//self.bin
            nx, ny = w//(stampSize if stampSize else psfWidth), h//(stampSize if stampSize else psfHeight)
            if nx == 1:
                nx += 1
            return maUtils.showPsfMosaic(calexp, nx=nx, ny=ny,
                                         showCenter=True, showEllipticity=False, showFwhm=True,
                                         stampSize=stampSize, frame=None, title=None).makeMosaic(mode=nx)
        else:
            im = maUtils.showPsfResiduals(calexp, ss[good], scale=self.bin,
                                          magType=self.magType,
                                          #apCorr=self.apCorr,
                                          frame=None).getImage()

        if self.sigma > 0:
            gaussFunc = afwMath.GaussianFunction1D(self.sigma)
            kernel = afwMath.SeparableKernel(int(5*self.sigma), int(5*self.sigma), gaussFunc, gaussFunc)
            cim = im.clone()
            afwMath.convolve(cim, im, kernel)
            im = cim

        # showPsfResiduals pads the image to allow psfs that hang over the edge
        im = im[afwGeom.BoxI(afwGeom.PointI(psfWidth//2, psfHeight//2),
                             afwGeom.ExtentI(calexp.getWidth()//self.bin, calexp.getHeight()//self.bin))]

        return im
    

def showPsfResiduals(data, dataId, magMin=None, magLim=23, magType="psf", showPsfStars=False,
                     showPsfModels=True, stampSize=0, radius=17920,
                     onlySerials=[], nJob=0, bin=10, sigma=0, frame=0, verbose=False):
    """Show the residuals across the camera resulting from subtracting PSFs from an image (the background
    is scaled down by bin)

    If showPsfModels is True, show the PSF model instead
    """
    mos = cameraGeomUtils.showCamera(data.butler.get("camera"),
                                      ResidualImage(data.butler, bin=bin, sigma=sigma, verbose=verbose,
                                                    magMin=magMin, magLim=magLim, magType=magType,
                                                    showPsfStars=showPsfStars,
                                                    showPsfModels=showPsfModels, stampSize=stampSize,
                                                    **dataId),
                                      onlySerials=onlySerials, nJob=nJob, bin=bin, frame=frame,
                                      title=data.mapperInfo.dataIdToTitle([dataId]))
    if frame is not None and radius:
        ds9.dot("o", -mos.getX0(), -mos.getY0(), size=radius/float(bin), ctype=ds9.RED, frame=frame)
        
    return mos

def showStackedPsfResiduals(data=None, dataId={}, exposure=None, sourceSet=None, sigma=None,
                            magType="psf", magMin=14, magMax=22, dmag=0.5, gutter=4, title="+",
                            referencePsf=0, normalizeReferenceFlux=False,
                            normalize=True, frame0=None):
    """Stack PSF residuals binned by magnitude"""

    if not exposure or not sourceSet:
        if not (data and dataId):
            raise RuntimeError("Please specify Data and dataId or an exposure and sourceSet")
        if exposure:
            exposures = [exposure]
        else:
            exposures = data.getDataset("calexp", dataId)
                
        if sourceSet:
            sourceSets = [sourceSet]
        else:
            sourceSets = data.getSources(dataId) # should make a single sourceSet here

    exposure = exposures[0]
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
        residualImages.append([afwImage.vectorImageF(), afwImage.vectorImageF(), magMin + (i + 0.5)*dmag, 0])

    for i, sourceSet in enumerate(sourceSets):
        exposure = exposures[i]
        mimIn = exposure.getMaskedImage()[:]
        psf = exposure.getPsf()

        for s in sourceSet:
            x, y = s.getX(), s.getY()

            try:
                flux = getFlux(s, magType)
                mag = exposure.getCalib().getMagnitude(flux)

                if not (magMin <= mag <= magMax):
                    continue

                if s.get("flags.pixel.saturated.center") or \
                        abs(2.5*math.log10(s.getPsfFlux()/s.getApFlux())) > 0.05:
                    continue

                bbox = afwGeom.BoxI(afwGeom.PointI(int(x) - psfWidth//2, int(y) - psfHeight//2),
                                    afwGeom.ExtentI(psfWidth, psfHeight))
                dx, dy = x - int(x), y - int(y) # offset from centre of subimage

                residualImages[bin(mag)][0].push_back(afwMath.offsetImage(Image(mimIn.getImage(),
                                                                                bbox, afwImage.PARENT),
                                                                          -dx, -dy))
                flux = np.nan           # i.e. fit the amplitude
                chi2 = measAlg.subtractPsf(psf, mimIn, x, y, flux)
            except (pexExcept.LsstCppException, ValueError), e:
                continue

            expIm = Image(mimIn.getImage(), bbox, afwImage.PARENT)
            expIm = afwMath.offsetImage(expIm, -dx, -dy)

            residualImages[bin(mag)][1].push_back(expIm)
            residualImages[bin(mag)][3] += 1
        
    objects = ds9Utils.Mosaic(gutter=gutter)
    residuals = ds9Utils.Mosaic(gutter=gutter)

    sctrl = afwMath.StatisticsControl()
    for i, vals in enumerate(residualImages):
        obj, res, mag, n = vals
        try:
            obj = afwMath.statisticsStack(obj, afwMath.MEANCLIP, sctrl)
            res = afwMath.statisticsStack(res, afwMath.MEANCLIP, sctrl)
        except:
            obj, res = Image(psfWidth, psfHeight), Image(psfWidth, psfHeight)

        if normalize and n > 0:
            peak = afwMath.makeStatistics(obj, afwMath.MAX).getValue()
            obj /= peak
            res /= peak

        residualImages[i] = (obj, res, mag, n,)

        if sigma > 0:
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            kernel = afwMath.SeparableKernel(int(5*sigma + 1), int(5*sigma + 1), gaussFunc, gaussFunc)
            cres = res.Factory(res, True)
            afwMath.convolve(cres, res, kernel, afwMath.ConvolutionControl(True, True))
            res = cres

        lab = "%.2f %d" % (mag, n) if n else ""
        objects.append(obj, lab)
        residuals.append(res, lab)
    #
    # Now subtract the brightest star from all the others
    #
    starResiduals = ds9Utils.Mosaic(gutter=gutter)
    referencePsfNo = referencePsf
    referencePsf = objects.images[referencePsf].clone()
    if normalizeReferenceFlux:
        referencePsf /= float(referencePsf.getArray().sum())

    for i, val in enumerate(residualImages):
        obj, _, mag, n = val
        res = obj.clone()

        if res.getArray().sum() > 0:
            if i == referencePsfNo:
                X, Y = np.meshgrid(range(res.getWidth()), range(res.getHeight()))
                res.getArray()[:] = 2e-3*np.where((X + Y)%2 == 1, 1, -1)
            else:
                if normalizeReferenceFlux:
                    res.getArray()[:] -= referencePsf.getArray()*res.getArray().sum()
                else:
                    res -= referencePsf

        lab = "%.2f %d" % (mag, n) if n else ""
        starResiduals.append(res, lab)

    title = re.sub(r"^\+\s*", " " + str(exposure.getDetector().getId()), title)
    mosaics = []
    nx = int(1/dmag + 0.5)
    mosaics.append(residuals.makeMosaic(    mode=nx, title=title, frame=frame0))
    mosaics.append(starResiduals.makeMosaic(mode=nx, title=title, frame=None if frame0 is None else frame0+1))
    mosaics.append(objects.makeMosaic(      mode=nx, title=title, frame=None if frame0 is None else frame0+2))
    mosaics.append([(mag, n) for obj, res, mag, n in residualImages])

    return mosaics

def showStackedPsfResidualsCamera(data, dataId, frame=0, overlay=False, normalize=False, **kwargs):
    """
    Show the stackedPsfResiduals laid out in true camera positions (unless overlay is True, in which case the
    chip are simply added

    Position of CCDs in Suprime cam specific! (See rewrite of showPsfResiduals)
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
                                                            gutter=subGutter, frame0=None, **kwargs)
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


def comparePsfModels(da1, da2, dataId, nx=3, stampSize=0):
    """Compare the PSF models from da1 and da2, returning a mosaic of the difference"""
    dataIds = da1.expandDataId(dataId)

    n = 0
    ims = {}
    onePixel = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.PointI(1, 1))
    dims = da1.butler.get("calexp", **dataIds[0]).getDimensions()
    for did in dataIds:
        for da in (da1, da2):
            try:
                psf = da.butler.get("calexp_sub", bbox=onePixel, **did).getPsf()
            except Exception, e:
                print "comparePsfModels reading %s: %s", (did, e)
                ims[da] = None          # continue the outer loop
                continue
            try:
                ims[da] = maUtils.showPsfMosaic(dims, psf, stampSize=stampSize, nx=nx).makeMosaic(mode=nx)
            except Exception, e:
                print "comparePsfModels, %s: %s", (did, e)
                ims[da] = None          # continue the outer loop
                continue

        if ims[da] == None:
            continue
            
        if n == 0:
            res = ims[da1].clone(); res[:] = 0
            im1 = ims[da1].clone()
            im2 = ims[da2].clone()
        else:
            im1 += ims[da1]
            im2 += ims[da2]

        ims[da1] -= ims[da2]
        res += ims[da1]
        n += 1

    res /= n
    im1 /= n
    im2 /= n

    return res, [im1, im2]

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
            bitmask = msk.getPlaneBitMask(bitmask)
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
            axes.axvline(val, linestyle="-", color=ctypes[i%len(ctypes)],
                      label="%-9s %.2f" % (v, val))
        i += 1
        axes.axvline(0, linestyle="--", color=ctypes[i%len(ctypes)])
            
        axes.legend(loc=1)

    fig.show()

    return fig

def plotFrames(das, fig=None):
    fig = getMpFigure(fig)
    
    plottingArea = (0.1, 0.1, 0.85, 0.80)
    axes = fig.add_axes(plottingArea)

    for k, da in das.items():
        for did in da.dataId:
            calexp_md = data.butler.get(dtName("calexp", True), **did)
            w, h = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
            wcs = afwImage.makeWcs(calexp_md)
            xvec, yvec = [], []
            for x, y in [(0, 0), (w, 0), (w, h), (0, h),]:
                pos = wcs.pixelToSky(x, y)
                xvec.append(pos[0].asDegrees())
                yvec.append(pos[1].asDegrees())
            #print k, did, w, h, xvec, yvec
            xvec.append(xvec[0]); yvec.append(yvec[0])

            axes.plot(xvec, yvec)

            p = pyplot.Polygon(zip(xvec, yvec), alpha=0.2 )
            axes.add_artist(p)

    axes.set_xlabel("RA")
    axes.set_ylabel("Dec")
    
    fig.show()

    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def findMissedSatur(data, maglim=15, deltaMin=-0.4, deltaMax=-0.06, frame=None, **dataId):
    ids = data.cat.get('id')
    stellar = data.cat.get('stellar')
    psf = data.cat.get('psfMag')
    ap = data.cat.get('apMag')

    sids = data.mapperInfo.splitId(ids, True)
    dataIdKeys = [k for k in sids.keys() if k != "objId"]

    missed = reduce(np.logical_and, [np.logical_not(stellar),
                                     psf < maglim,
                                     ap - psf >= deltaMin,
                                     ap - psf < deltaMax,
                                     data.cat.get("flags.pixel.bad") == False,
                                     data.cat.get("flags.pixel.saturated.center") == False],
                    )
    #
    # Possible restriction to some dataId
    #
    for k, v in dataId.items():
        try:
            inV = sids[k] == v[0]

            for vv in v[1:]:
                inV = np.logical_or(inV, sids[k] == vv)

            missed = np.logical_and(missed, inV)
        except TypeError:
            missed = np.logical_and(missed, sids[k] == v)
            
    #
    # OK, we've got our bad objects
    #
    x = (data.cat.getX()[missed] + 0.5).astype(int)
    y = (data.cat.getY()[missed] + 0.5).astype(int)
    nchild = data.cat.get("deblend.nchild")[missed]

    psf = psf[missed]
    ap =  ap[missed]

    sids = data.mapperInfo.splitId(ids[missed], True)
    if len(psf) == 1:                   # => the elements of dids are ints, not arrays
        for k, v in sids.items():
            sids[k] = [v]
            
    ids = [dict(zip(sids.keys(), _)) for _ in zip(*sids.values())]

    del sids["objId"]
    for s in sorted(set([tuple(zip(sids.keys(), _)) for _ in zip(*sids.values())])):
        did = dict([(k, int(v)) for k, v in s]) # n.b. sqlite3 doesn't like numpy integral types

        raw = cameraGeomUtils.trimExposure(data.butler.get('raw', **did),
                                           subtractBias=True).getMaskedImage().getImage()
        if frame is not None:
            ds9.mtv(raw, title=" ".join(["%s:%s" % kv for kv in did.items()]), frame=frame)
            pass

        for i, _id in enumerate(ids):
            if sum([_id[k] != v for k, v in did.items()]) == 0:
                if nchild[i] > 0:
                    continue
                print "%-50s (%6.1f, %6.1f) %6.2f %6.2f %g" % (_id, x[i], y[i], psf[i], ap[i] - psf[i],
                                                               int(raw[x[i], y[i]]))
                if frame is not None:
                    ds9.dot("+", x[i], y[i], size=4, frame=frame, ctype=ds9.RED)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plotImageCrossSection(im, dir='-', topBottom=False, fig=1, ymin=None, ymax=None, title=None):
    w, h = im.getDimensions()
    try:
        im = im.getMaskedImage()
    except AttributeError:
        pass

    try:
        im = im.getImage()
    except AttributeError:
        pass

    vals2 = None
    if dir == '-':
        if topBottom:
            vals1 = np.empty(w); vals2 = np.empty(w)
            r10, r11 = int(0.1*h), int(0.2*h)
            r20, r21 = int(0.8*h), int(0.9*h)
            for x in range(w):

                vals1[x] = afwMath.makeStatistics(im[x:x+1, r10:r11], afwMath.MEANCLIP).getValue()
                vals2[x] = afwMath.makeStatistics(im[x:x+1, r20:r21], afwMath.MEANCLIP).getValue()

            label1 = "Rows %d..%d" % (r10, r11)
            label2 = "Rows %d..%d" % (r20, r21)
        else:
            vals1 = np.empty(w)
            for x in range(w):
                vals1[x] = afwMath.makeStatistics(im[x:x+1, :], afwMath.MEANCLIP).getValue()
    elif dir == '|':
        vals1 = np.empty(h)
        for y in range(h):
            vals1[y] = afwMath.makeStatistics(im[:, y:y+1], afwMath.MEANCLIP).getValue()
    else:
        raise RuntimeError("Direction %s is not yet implemented" % dir)

    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))
    axes.plot(vals1, label=label1)
    if vals2.all() is not None:
        vals2 *= vals1[100]/vals2[100]
        axes.plot(vals2, color='red', label="%s scaled to bottom" % label2)
        axes.legend(loc="upper left")
        
    axes.set_xlim(-0.02*len(vals1), 1.02*len(vals1))
    axes.set_ylim(ymin, ymax)
    axes.set_xlabel("column" if dir == '-' else "row")
    axes.set_ylabel("DN")
    if title:
        axes.set_title(title)

    fig.show()

    return fig

    
def standardEllipticity(ss, calexp, magLim=100, frame=0):
    n = 0
    sIxx, sIyy, sIxy = 0.0, 0.0, 0.0
    for s in ss[ss.get("calib.psf.used")]:
        Q = s.getShape()                # a afwGeom.ellipses.Quadrupole

        sIxx += Q.getIxx()
        sIyy += Q.getIyy()
        sIxy += Q.getIxy()
        n += 1
        
    sIxx /= n
    sIyy /= n
    sIxy /= n

    Qstar = afwGeom.ellipses.Quadrupole(sIxx, sIyy, sIxy)
    print "Qstar =", afwGeom.ellipses.Axes(Qstar)

    calib = calexp.getCalib()
    if calib and magLim is not None:
        with afwImageUtils.CalibNoThrow():
            modelMag = calib.getMagnitude(ss.getModelFlux())
    else:
        modelMag = np.ones_like(ss.getModelFlux())
        magLim = modelMag + 1

    rms_standard = 0.7/0.200
    Izz_standard = math.pow(rms_standard, 2)

    scale = 20                           # scale for plotting ellipses
    ds9.erase(frame=frame)
    with ds9.Buffering():
        for s in ss[np.less(modelMag, magLim)]:
            Q = s.getShape()                # a afwGeom.ellipses.Quadrupole

            if \
                    Q.getIxx() < Qstar.getIxx() or \
                    Q.getIyy() < Qstar.getIyy():
                Ixx, Iyy, Ixy = 0, 0, 0
            else:
                Ixx = Q.getIxx() - Qstar.getIxx()
                Iyy = Q.getIyy() - Qstar.getIyy()
                Ixy = Q.getIxy() - Qstar.getIxy()

            Qstandard = afwGeom.ellipses.Quadrupole(Ixx + Izz_standard, Iyy + Izz_standard, Ixy)

            #Q.scale(scale/Q.getDeterminantRadius())
            Qstandard.scale(scale/Qstandard.getDeterminantRadius())

            ds9.dot(Q, *s.getCentroid(), ctype=ds9.RED, frame=frame)
            ds9.dot(Qstandard, *s.getCentroid(), ctype=ds9.GREEN, frame=frame)

        
        
def showApcorr(data,  ymin=None, ymax=None, markersize=1, fig=None, **dataId):
    """Show the aperture corrections for a given detector"""
    
    fig = getMpFigure(fig)

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    ss = data.butler.get("src", **dataId)
    xmm = ss.get("focalplane.x")
    ymm = ss.get("focalplane.y")
    apcorr = ss.get("correctfluxes.apcorr")
        
    sc = axes.scatter(xmm, ymm, c=apcorr, norm=pyplot.Normalize(ymin, ymax),
                      cmap=pyplot.cm.rainbow, marker="o", s=10*markersize,
                      edgecolors="none")   
    fig.colorbar(sc)
    axes.set_aspect('equal')

    psfStar = ss.get("calib.psf.used")
    axes.plot(xmm[psfStar], ymm[psfStar], "*", color="black", markerfacecolor="none")
    
    axes.set_xlabel("X (mm)")
    axes.set_ylabel("Y (mm)")
    axes.set_title("Aperture correction %s" % data.mapperInfo.dataIdToTitle([dataId]))
    
    fig.show()
    return fig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def medianFilterImage(img, nx, ny=None, verbose=False):
    if nx%1 == 0:
        nx += 1
    if ny:
        if ny%1 == 0:
            ny += 1
    else:
        ny = nx

    w, h = img.getDimensions()

    mimg = img.clone()

    imga = img.getArray()
    mimga = mimg.getArray()
    sim = afwImage.ImageF(nx, ny); sima = sim.getArray() # permits me to mtv it

    for y in range(nx//2, h - nx//2):
        if verbose:
            print "%d\r" % y,; sys.stdout.flush()
        for x in range(ny//2, w - ny//2):
            sima[:] = imga[y - ny//2:y + ny//2, x - nx//2:x + nx//2]
            mimga[y, x] = np.median(sima[np.where(np.isfinite(sima))])

            if False:
                ds9.mtv(sim)

    return mimg

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def addBackgroundCallback(im, ccd=None, butler=None, imageSource=None):
    """A callback function that adds the background back into a calexp"""

    assert butler

    if not ccd:
        ccd = cameraGeom.cast_Ccd(im.getDetector())

    bkgd = butler.get("calexpBackground", ccd=ccd.getId().getSerial(), **imageSource.kwargs)

    im += bkgd.getImage()

    return im

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# makeVignettingImage getVignetting correctVignettingCallback
# are now in lsst/obs/subaru/ccdTesting.py
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def assembleTiles(images):
    """Assemble a list of tiles according to their XY0 values"""
    bigBBox = afwGeom.BoxI()

    for im in images:
        bigBBox.include(im.getBBox(afwImage.PARENT))

    bigIm = afwImage.MaskedImageF(bigBBox)
    for im in images:
        if True:
            sub = bigIm.Factory(bigIm, im.getBBox(afwImage.PARENT), afwImage.PARENT)
            sub <<= im.getMaskedImage()
            del sub
        else:
            bigIm[im.getBBox(afwImage.PARENT)] = im

    return bigIm

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def readPatches(butler, tract=0, patches={}, pps=None):
    if pps is None:
        pps=["%d,%d" % (i , j) for i in range(1, 8) for j in range(1, 8)]

    for pp in pps:
        for f in "gri":
            patches.update([("%s-%s" % (pp, f), butler.get("deepCoadd", filter='HSC-%s' % f.upper(),
                                                          tract=tract, patch=pp))])

    pps = set([k[0:3] for k in patches.keys() if re.search(r"^.,.-.$", k)])
    patches.update([(f, assembleTiles([patches["%s-%s" % (pp, f)] for pp in pps])) for f in "gri"])

    return patches

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def patchNucleus(images, pp, threshold=2000, name=""):
    for f in "gri":
        mi = images["%s-%s" % (pp, f)].getMaskedImage()
        SAT = mi.getMask().getPlaneBitMask(["SAT", "INTRP"])
        fs = afwDetect.FootprintSet(mi.getMask(), afwDetect.Threshold(SAT, afwDetect.Threshold.BITMASK))
        for foot in fs.getFootprints():
            if foot.getNpix() > 15000 and len(foot.getSpans()) > 100:
                print "Adjusting Footprint at", foot.getCentroid(), f, name
                xv, yv = [], []
                for s in foot.getSpans():
                    y = s.getY()
                    for x in range(s.getX0(), s.getX1() + 1):
                        xv.append(x)
                        yv.append(y)

                xv = np.array(xv); yv = np.array(yv)
                yv -= mi.getY0();  xv -= mi.getX0()

                tooLow = np.where(mi.getImage().getArray()[yv, xv] < threshold)[0]
                mi.getMask().getArray()[yv[tooLow], xv[tooLow]] &= ~SAT

        #ds9.mtv(mi, title=f); import pdb; pdb.set_trace() #return

def readM31Images(butler, images={}, pps=None):
    if pps is None:
        pps=["%d,%d" % (i , j) for i in range(1, 8) for j in range(1, 8)]

    for pp in pps:
        for f in "gri":
            if False:                   # old and tested
                images.update([("%s-%s" % (pp, f), afwImage.ExposureF("HSC-%s/0/%s.fits" % (f.upper(), pp)))])
            else:                       # better to use the butler
                images.update([("%s-%s" % (pp, f), butler.get("deepCoadd", filter='HSC-%s' % f.upper(),
                                                              tract=0, patch=pp))])


    # Patch saturation for nuclei of M31/M32
    for pp, threshold in (["4,2", 1300],
                          ["4,4", 1500]): # 1000 OK, 2000 bad
        if pp in pps:
            patchNucleus(images, pp, threshold, name=pp)

    pps = set([k[0:3] for k in images.keys() if re.search(r"^.,.-.$", k)])
    images.update([(f, assembleTiles([images["%s-%s" % (pp, f)] for pp in pps])) for f in "gri"])

    return images

def hackM31Sky(fileFmt="sat-16-%s.fits", outfile="fixed.png", *args, **kwargs):
    """Hack up the M31 sky levels"""
    images = {}
    images.update([(f, afwImage.MaskedImageF(fileFmt % f)) for f in "RGB"])

    sim=afwImage.MaskedImageF(5,5)

    im = images["B"]
    sim[0:3,0:3]=10; sim[3:,0:3]=35; sim[0:3,3:]=24; sim[3:,3:]=35
    sim[0:3,0:3] -= 0
    sim[0:3,0:3] -= 2
    sim[3:,3:] -= 2
    sim[0, 0] += 2
    sim[4, 0] -= 13
    sim[0,1] += 2
    sim[1,2] += 2
    sim[0,3] -= 2
    sim[4,0:4] += 10
    sim[0,3] -= 2
    sim[4,1] -= 6
    sim[1,4] -= 3
    sim[1,1] += 2
    sim[4,1] += 2
    backobj = afwMath.BackgroundMI(im.getBBox(), sim); im -= backobj.getImageF("LINEAR")

    im = images["G"]
    sim[0:3,0:3]=8; sim[3:,0:3]=8; sim[0:3,3:]=14; sim[3:,3:]=8
    sim[0:3,0:3] -= 2
    sim[0:3,0:3] -= 2
    sim[3:,3:] -= 2
    sim[4, 0] -= 8
    sim[0,2] += 2
    sim[3,3] += 3
    sim[3,2] += 10
    sim[4,1] += 4
    sim[0,4] += 5
    sim[4,1] -= 3
    sim[1,4] -= 2
    sim[1,1] += 2
    sim[3,0] += 5
    sim[4,0] += 7
    backobj = afwMath.BackgroundMI(im.getBBox(), sim); im -= backobj.getImageF("LINEAR")

    im = images["R"]
    sim[0:3,0:3]=7; sim[3:,0:3]=17; sim[0:3,3:]=1; sim[3:,3:]=5;
    sim[:, 0] += 8
    sim[3, 0] += 10
    sim[3, 1] += 3
    sim[4, 0] -= 5
    sim[0:3,0:3] -= 2
    sim[0:2,0:3] -= 2
    sim[2:,3:] += 2
    sim[0,1] += 2
    sim[1,1] += 2
    sim[1,1] += 2
    sim[4,1] += 10
    sim[4,1] -= 6
    sim[1,4] -= 3
    sim[1,1] += 2
    sim[3,3] += 2
    sim[1,1] -= 2
    sim[2,1] += 5
    sim[2,0] += 3

    sim[2,4] += 5
    sim[2,3] += 7
    sim[3,3] += 7
    sim[4,2] += 5
    backobj = afwMath.BackgroundMI(im.getBBox(), sim); im -= backobj.getImageF("LINEAR")

    writeRgb([images["R"], images["G"], images["B"]], outfile, *args, **kwargs)
    
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
#import lsst.meas.algorithms.defects as defects

def trimRemoveCrCallback(im, ccd=None, butler=None, imageSource=None, subtractBackground=False):
    if hasattr(im, "convertF"):
        im = im.convertF()
        
    im = cameraGeomUtils.trimRawCallback(im, ccd, butler, imageSource=imageSource, subtractBackground=True)
    if True:                            # XXX
        return im

    if False:
        defects = afwDetect.FootprintSet(im, afwDetect.Threshold(10, afwDetect.Threshold.STDEV, False), "SAT")
        del defects

    if subtractBackground:
        bctrl = afwMath.BackgroundControl(2, 5)
        for a in ccd:
            aim = im[a.getAllPixels()]
            backobj = afwMath.makeBackground(aim, bctrl)
            bkgd = backobj.getImageF(afwMath.Interpolate.AKIMA_SPLINE, afwMath.REDUCE_INTERP_ORDER)

            aim[:] -= bkgd

    FWHM = 5                   # pixels
    psf = measAlg.DoubleGaussianPsf(29, 29, FWHM/(2*math.sqrt(2*math.log(2))))

    stats = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    background = stats.getValue(afwMath.MEANCLIP)

    if hasattr(im, "getImage"):
        mi = im
    else:
        mi = afwImage.makeMaskedImage(im)
        mi.getVariance()[:] = stats.getValue(afwMath.STDEVCLIP)**2

    crConfig = measAlg.FindCosmicRaysConfig()
    crConfig.nCrPixelMax = 100000
    crConfig.nCrPixelMax *= 10
    crs = measAlg.findCosmicRays(mi, psf, background, pexConfig.makePolicy(crConfig))
    if not False:
        print "CCD %03d %5d CRs" % (ccd.getId().getSerial(), len(crs))

    return im

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def analyzeForced(dataId, sourceSet_forced, sourceSet_sfm, magType="psf",
                  showAstrometricResiduals=False, matchRadius=1,
                  SG="sg", maglim=None, xmin=None, xmax=None, ymin=None, ymax=None,
                  frames=[0], markersize=3, title="+", fig=None):

    sourceSet_forced = _appendToCatalog(None, dataId, sourceSet=sourceSet_forced)[0]
    sourceSet_sfm    = _appendToCatalog(None, dataId, sourceSet=sourceSet_sfm)[0]
    
    matched = afwTable.matchRaDec(sourceSet_forced, sourceSet_sfm, matchRadius*afwGeom.arcseconds)

    if showAstrometricResiduals:
        calexp_md = butler.get("calexp_md", **dataId)
        wcs = afwImage.makeWcs(calexp_md)
        
        #wcs = improveWcs(wcs, matched, sipOrder=5)

    matched = zipMatchList(matched, verbose=False)
    
    mag_f = matched.get("%sMag_1" % magType)
    mag_s = matched.get("%sMag_2" % magType)

    stellar = matched.get("stellar_1")
    x = matched.get("centroid.sdss_2.x")
    y = matched.get("centroid.sdss_2.y")
    ids = matched.get("id")

    if showAstrometricResiduals:
        x2 = np.empty_like(x)
        y2 = np.empty_like(y)
        for i, m in enumerate(matched):
            x2[i], y2[i] = wcs.skyToPixel(m.getCoord())

        xvec = np.hypot(x - x2, y - y2)
        if xmin is None: xmin = -0.002
        if xmax is None: xmax = 1
        xlab = "position error (pixels)"
    else:
        xvec = mag_f + 0.0              # adding 0.0 makes xvec writable
        if xmin is None: xmin = 14
        if xmax is None: xmax = 23
        xlab = magType

    yvec = mag_f - mag_s
    if ymin is None: ymin = -0.01
    if ymax is None: ymax =  0.10
        
    if maglim is not None:
        good = mag_f < maglim

        xvec = xvec[good]
        yvec = yvec[good]

        x = x[good]
        y = y[good]
        ids = ids[good]
        stellar = stellar[good]

    nonStellar = np.logical_not(stellar)
    #
    # Time to plot
    #
    fig = getMpFigure(fig)    
    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    if "s" in SG.lower():
        axes.plot(xvec[stellar], yvec[stellar], "g.", markersize=markersize, markeredgewidth=0)
    if "g" in SG.lower():
        axes.plot(xvec[nonStellar], yvec[nonStellar], "r.", markersize=markersize, markeredgewidth=0)

    if True:
        axes.set_xlim(xmin, xmax)
        axes.set_ylim(ymin, ymax)
    axes.set_xlabel(xlab)
    axes.set_ylabel("(forced - sfm) %s" % (magType))

    if maglim is not None:
        title += "%s < %g" % (magType, maglim)
        
    axes.set_title(re.sub(r"^\+\s*", str(dataId) + " ", title))    
    #
    # Make plot live
    #
    global eventHandlers
    flags = {}

    if "g" not in SG.lower():
        xvec[nonStellar] = -1000
    if "s" not in SG.lower():
        xvec[stellar] = -1000
    
    eventHandlers[fig] = EventHandler(axes, xvec, yvec, ids, x, y, flags, frames=frames)

    fig.show()

    return matched

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

from lsst.pipe.tasks.astrometry import AstrometryTask
import lsst.meas.astrom.sip as astromSip

def improveWcs(wcs, matched, sipOrder=3):
    rmv = afwTable.ReferenceMatchVector()
    rmv.reserve(len(matched))
    for m in matched:
        rm = afwTable.ReferenceMatch(m.first, m.second, m.distance)
        rmv.push_back(rm)
    matched = rmv

    medToOneDRms = 1/(math.sqrt(2*math.log(2))) # convert median(hypot(dx, dy)) to a 1-D rms for a Gaussian

    sipOrderMax = sipOrder
    sipOrder = 1
    iter = 0
    bestWcs, bestRms = wcs, np.nan
    while sipOrder <= sipOrderMax:
        try:
            sipObject = astromSip.CreateWcsWithSip(matched, wcs, sipOrder)
            wcs = sipObject.getNewWcs()

            sipRms = sipObject.getScatterInPixels()*medToOneDRms

        except pexExcept.LsstCppException, e:
            print >> sys.stderr, ('Failed to calculate distortion terms. Error: %s' % e)
            sipRms = None
            #break

        dx = np.empty(len(matched))
        dy = np.empty_like(dx)
        for i, m in enumerate(matched):
            dx[i], dy[i] = wcs.skyToPixel(m.first.getCoord())
            dx[i] -= m.first.get("centroid.sdss.x")
            dy[i] -= m.first.get("centroid.sdss.y")

        dz = np.hypot(dx, dy)
        i = np.arange(len(dz))
        
        nsigma = 4
        rms = np.median(dz)*medToOneDRms
        good = dz < nsigma*rms
        bad = np.logical_not(good)

        if sipRms is None:
            sipRms = rms

        if not np.isfinite(bestRms) or sipRms < bestRms:
            bestWcs, bestRms = wcs, sipRms

        print 'Sip iteration %i: order %d  %i objects match. rms scatter is %g pixels' % \
            (iter, sipOrder, len(matched), sipRms)

        if False:
            fig = getMpFigure(2)
            axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));
            axes.plot(i[good], dz[good], "g.")
            axes.plot(i[bad], dz[bad], "r.")
            fig.show()
            import pdb; pdb.set_trace() 
        sipOrder += 1
        
        rmv = afwTable.ReferenceMatchVector()
        rmv.reserve(int(sum(good)))
        for i, m in enumerate(matched):
            if good[i]:
                rmv.push_back(m)
        matched = rmv

    print "bestRms = %g" % bestRms
    return bestWcs

