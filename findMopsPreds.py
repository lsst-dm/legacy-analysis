#!/usr/bin/env python
"""
Given a DC2 runId and a visitId, checks which MopsPreds are to be found in the
images that are present for that runId/visitId, and prints their coordinates
in the images.

Run as findMopsPreds.py runId visitId

NOTE: currently works with DC2 schema prior to 1/30/08 only!
"""
import string, sys, os

# image directory setup

imageDirBase = '/lsst/DC2root'
## imageDirBase = './test'

# set up mysql access
mySqlHost = 'lsst10.ncsa.uiuc.edu'
mySqlUser = 'test'
mySqlPasswd = 'globular.test'

# get input args
runId = sys.argv[1]
visitId = sys.argv[2]

print 'Run: %s Visit: %s\n' % (runId, visitId)

# get list of image files

imageDir = '%s/%s/ipd/output/%s' % (imageDirBase, runId, visitId)
imageCcdDirs = os.listdir(imageDir)

imageList = []
for dir in imageCcdDirs:
    # need to split ccd from path of form /lsst/DC2root/rlp0083/ipd/output/704893/003
    dirComps = string.split(dir, '/')
    nComps = len(dirComps)
    ccd = dirComps[nComps-1]
    imageName = '%s/%s/%sp_%s_diff_img.fits' % (imageDir,dir,visitId,ccd)
##     print imageName
    imageList.append(imageName)

# get list of MopsPreds

mySqlQuery = 'use %s; select orbit_id, ra_deg, dec_deg from MopsPreds_visit%s;' % (runId, visitId)

mySqlCmd = 'mysql -h %s -u %s -p%s -e \'%s\'' % (mySqlHost, mySqlUser, mySqlPasswd, mySqlQuery)

## print mySqlCmd

for foo in os.popen(mySqlCmd).readlines():
    (orbit_id, ra, dec) = string.split(foo)
    if ra != 'ra_deg':
##         print ra, dec
        for image in imageList:
            sky2xyCmd = 'sky2xy %s %s %s' % (image, ra, dec)
            result = os.popen(sky2xyCmd).read()  # want form that gives all as one line
            if string.find(result,'off') != -1:
                # not found
##                 print result
                pass
            else:
                # found
                print orbit_id, image, ra, dec, result




