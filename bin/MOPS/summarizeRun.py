#!/usr/bin/env python
"""
Given a DC2 runId, summarize quantities of interest from the database

Run as summarizeRun.py runId

"""
import string, sys, os

def execQuery(query):

    # set up mysql access
    mySqlHost = 'lsst10.ncsa.uiuc.edu'
    mySqlUser = 'test'
    mySqlPasswd = 'globular.test'

    mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

    print os.popen(mySqlCmd).read()
    

# get input args
runId = sys.argv[1]

print 'Run: %s' % (runId)
print '--------------------------------------------------------\n'

# Summarize visits that were processed


mySqlQuery = \
           'USE %s;\
           SELECT visitId, ra, decl, filterId, dateObs, expTime, airmass\
           FROM Visit v, Raw_FPA_Exposure e\
           WHERE v.exposureId = e.rawFPAExposureId;' % (runId)

execQuery(mySqlQuery)

# Get number of DIASources matching existing Objects and number that are new

print 'Number of DIASources matching existing Objects; number that are new'
print '------------------------------------------------------------------\n'

mySqlQuery = \
           'USE %s;\
           SELECT COUNT(DISTINCT first) Num_Matching\
           FROM DiaSourceToObjectMatches;\
           SELECT COUNT(*) Num_New FROM NewObjectIdPairs;' % (runId)

execQuery(mySqlQuery)

# Get average match multiplicity

print 'Average and max number of sources matched to each object'
print '--------------------------------------------------------\n'

mySqlQuery = \
           'USE %s; \
           SELECT AVG(a.nmatch) Avg_match_multiplicity, MAX(a.nmatch) Max_match_multiplicity\
           FROM\
           (SELECT count(*) nmatch, second FROM DiaSourceToObjectMatches GROUP BY second, visitId) a;' % (runId)

execQuery(mySqlQuery)

# Get number of MopsPreds that were detected

print 'MopsPreds that were detected'
print '----------------------------\n'

mySqlQuery = \
           'USE %s; \
           SELECT COUNT(DISTINCT first,visitId) Num_Detected \
           FROM MopsPredToDiaSourceMatches;' % (runId)

execQuery(mySqlQuery)

mySqlQuery = 'use %s; select DISTINCT first MopsId,visitId from MopsPredToDiaSourceMatches;' % (runId)

execQuery(mySqlQuery)

# Get number of DIASources for each Visit and CCD

print 'Number of DIASources for each Visit and CCD'
print '-------------------------------------------\n'


mySqlQuery =\
           'USE %s;\
           SELECT visitId, url, count(*) nSources\
           FROM Visit v, Raw_CCD_Exposure e, DIASource s\
           WHERE v.exposureId = e.rawFPAExposureId\
           AND e.rawCCDExposureId = s.ccdExposureId\
           GROUP BY e.rawCCDExposureId;' % (runId)


execQuery(mySqlQuery)
