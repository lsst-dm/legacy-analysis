#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

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
