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

#
# Tools useful for SuprimeCam (e.g. they may hardwire cameraGeom) but useful at least as
# a starting point for other cameras
#
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms.utils as maUtils
import lsst.analysis as analysis

def showPsf(da, distort=False, stampSize=19, nx=7, ny=None, gutterLevel=0.05):
    mos = displayUtils.Mosaic(background=gutterLevel, gutter=1)

    if ny is None:
        ny = 2*nx

    dataId = da.dataId.copy()
    for ccd in (8, 9, 5, 4, 3,  6, 7, 2, 1, 0):
        dataId.update(ccd=ccd)
        calexp = da.getDataset("calexp", dataId)[0]
        if False:
            print calexp.getDetector().getId()
            
        mos.append(maUtils.showPsfMosaic(calexp, stampSize=stampSize,
                                         nx=nx, ny=ny, distort=distort).makeMosaic(mode=nx))

    return mos.makeMosaic(mode=5)       # Camera is 5x2
