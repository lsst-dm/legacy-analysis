#!/usr/bin/env python
"""
Given a DC2 runId, summarize quantities of interest from the database

Run as summarizeRun.py runId

"""
import string, sys, os

# set up mysql access
mySqlHost = 'lsst10.ncsa.uiuc.edu'
mySqlUser = 'test'
mySqlPasswd = 'globular.test'

# get input args
runId = sys.argv[1]

print 'Run: %s\n' % (runId)

# Summarize visits that were processed


mySqlQuery = 'use %s; select count(*) Num_Visits from Raw_FPA_Exposure; select count(*) Num_CCD_Visits from Raw_CCD_Exposure; select ROUND(rawFPAExposureId/2) visit, ra, decl, filterId, dateObs, expTime, airmass from Raw_FPA_Exposure;' % (runId)

mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

print os.popen(mySqlCmd).read()

# Get number of DIASources matching existing Objects and number that are new

print 'Number of DIASources matching exiting Objects; number that are new\n'

mySqlQuery = 'use %s; select count(*) Num_Matching from DiaSourceToObjectMatches; select count(*) Num_New from NewObjectIdPairs;' % (runId)

mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

print os.popen(mySqlCmd).read()

# Get number of MopsPreds that were detected

print 'Number of MopsPreds that were detected\n'

mySqlQuery = 'use %s; select count(*) Num_Detected from MopsPredToDiaSourceMatches;' % (runId)

mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

print os.popen(mySqlCmd).read()

# Get number of DIASources for each Visit and CCD

print 'Number of DIASources for each Visit and CCD\n'

mySqlQuery = 'use %s; select ROUND(jj.rawFPAExposureId/2) visit, jj.url, count(*) nSources from DIASource s inner join (select f.rawFPAExposureId, c.rawCCDExposureId, c.url from Raw_FPA_Exposure f inner join Raw_CCD_Exposure c on f.rawFPAExposureId=c.rawFPAExposureId) jj on s.ccdExposureId = jj.rawCCDExposureId group by jj.url;;' % (runId)

mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

print os.popen(mySqlCmd).read()






