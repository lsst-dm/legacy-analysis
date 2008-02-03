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

print 'Run: %s\n' % (runId)

# Summarize visits that were processed


mySqlQuery = 'use %s; select count(*) Num_Visits from Raw_FPA_Exposure; select count(*) Num_CCD_Visits from Raw_CCD_Exposure; select ROUND(rawFPAExposureId/2) visit, ra, decl, filterId, dateObs, expTime, airmass from Raw_FPA_Exposure;' % (runId)

execQuery(mySqlQuery)

# Get number of DIASources matching existing Objects and number that are new

print 'Number of DIASources matching exiting Objects; number that are new\n'

mySqlQuery = 'use %s; select count(*) Num_Matching from DiaSourceToObjectMatches; select count(*) Num_New from NewObjectIdPairs;' % (runId)

execQuery(mySqlQuery)

# Get average match multiplicity

print 'Average and max number of sources matched to each object\n'

mySqlQuery = 'use %s; select AVG(a.nmatch) Avg_match_multiplicity, MAX(a.nmatch) Max_match_multiplicity from (select count(*) nmatch, second from DiaSourceToObjectMatches group by second, visitId) a;' % (runId)

execQuery(mySqlQuery)

# Get number of MopsPreds that were detected

print 'Number of MopsPreds that were detected\n'

mySqlQuery = 'use %s; select count(*) Num_Detected from MopsPredToDiaSourceMatches;' % (runId)

execQuery(mySqlQuery)

# Get number of DIASources for each Visit and CCD

print 'Number of DIASources for each Visit and CCD\n'

mySqlQuery = 'use %s; select ROUND(jj.rawFPAExposureId/2) visit, jj.url, count(*) nSources from DIASource s inner join (select f.rawFPAExposureId, c.rawCCDExposureId, c.url from Raw_FPA_Exposure f inner join Raw_CCD_Exposure c on f.rawFPAExposureId=c.rawFPAExposureId) jj on s.ccdExposureId = jj.rawCCDExposureId group by jj.url;;' % (runId)

execQuery(mySqlQuery)
