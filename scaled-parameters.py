import pyfits
import numpy as np
import sys
import math
import os
import re
from scipy import interpolate
import scipy.signal
from optparse import OptionParser
from kapteyn import wcs
from pylab import *

import actutils

# TODO apply image amplitude and local distance cuts!

# note can use wcs.coordmap() to generate the coordinate map input to
# scipy.interpolate if we want to do interpolation properly (see Kapteyn manual)

class TelLookupTable(object):
    """ 
    Read and interpolate telescope-wise lookup tables
    """

    _valueDict = dict()
    _proj = dict()
    values = dict()
    telID = 0

    def __init__(self, lookupName, telID, 
                 lookupDir= None,useSmooth=True):
        """
        initialize a lookup table

        TODO: eventually these will be 4D data cubes probably! (for zenith
        and azimuth dependence)

        `lookupName`: name of lookup table to load (e.g. HIL_TEL_WIDTH)
        `telID`: telescope ID number
        """
        if lookupDir==None:
            lookupDir = os.environ["HOME"]+"/Analysis/FITSEventLists/Lookups"

        self.values = ["mean","stddev","count"]    
        self._valueDict = dict()
        self.telID = telID

        print "Loading Lookups for CT",tel

        smooth=""
        if useSmooth: 
            smooth="-gauss"

        for what in self.values:
            fname = "CT%03d-%s-lookup-%s%s.fits" % (tel,lookupName,what,smooth)    
            print "\t",fname
            hdu = pyfits.open(lookupDir+"/"+fname)[0]
            edgesAndHist = actutils.histFromFITS( hdu )
            self._valueDict[what] = edgesAndHist
            self._proj[what] = wcs.Projection( hdu.header )

#            xed,yed,val = self._valueDict[what]
#            figure()
#            pcolor( xed,yed,val.transpose() ) 
#            title("CT%d: %s" % (self.telID, what))
#            show()           



        #self.extrapolateLookups()

    def extrapolateLookups(self, minCounts=10):
        """
        Make the lookup tables nicer by extrapolating unfilled values
        
        Arguments:
        - `self`:
        """

        xed,yed,counts = self._valueDict["count"]

        for what in self.values:
            if what=='count': continue
            print "Extrapolating: ", what
            xed,yed,val = self._valueDict[what]

            # todo: do this as an array op, not for-loop

            for ii in xrange(counts.shape[0]):
                lastgood = (0,0)
                for jj in xrange(counts.shape[1]):
                    if counts[ii,jj] > minCounts:
                        lastgood=(ii,jj)
                    else:
                        val[ii,jj] = val[lastgood]
                        counts[ii,jj]=minCounts+1

            for ii in xrange(counts.shape[0]):
                lastgood = (0,0)
                for jj in xrange(counts.shape[1]-1,0,-1):
                    if counts[ii,jj] > minCounts:
                        lastgood=(ii,jj)
                    else:
                        val[ii,jj] = val[lastgood]
                        counts[ii,jj]=minCounts+1



    def getValue(self, coord, what="mean" ):
        """
        
        Arguments:
        - `coord`: coordinate or array of coordinates  in logSize,Dist space
        - `impactDist`: single value or array of values

        returns value in the lookup table
        """
        
        world = np.array( coord )
        bin = np.trunc(self._proj[what].topixel( world ) - 1.0).astype(int)
        xed,yed,val = self._valueDict[what]
        shape = np.array(val.shape)
        
        # outliers:
        bin[bin>shape] = shape-1
        v = val[bin]

        #print "coord=",coord,"bin=",bin, " value=",v

        return v
        
# ===========================================================================

def calcMeanReducedScaledValue( tels, coords, vals, lookup):
    """
    
    `tels`: array of telescope-ids
    `vals`: array of values corresponding to each tel
    `lookup`: the big lookup table structure
    """

    mrsv = 0.0
    ntels = vals.shape[0]

    if ntels==0:
        return -10000

    for itel in xrange( ntels ):
        vMean = lookup[tels[itel]].getValue( coords[itel], "mean" )
        vSigma = lookup[tels[itel]].getValue( coords[itel], "stddev" )
        mrsv += (vals[itel] - vMean)/vSigma

    return mrsv/float(ntels)
    
                        


if __name__ == '__main__':
    
    ineventlistfile = sys.argv.pop(1)
    evfile = pyfits.open(ineventlistfile)
    events = evfile["EVENTS"]
    telarray = evfile["TELARRAY"]

    tposx = telarray.data.field("POSX")
    tposy = telarray.data.field("POSY")
    telid = np.array(events.header['TELLIST'].split(",")).astype(int)

    # load the lookups:
    lookupName = "HIL_TEL_WIDTH"

    telLookup = dict()

    for tel in telid:
        telLookup[tel] = TelLookupTable(lookupName,tel)

    # THE FOLLOWING IS SIMILAR TO generate-lookup-tables (should combine them)
         
    telMask   = events.data.field("TELMASK") 
    allValues = events.data.field(lookupName) 
    allSizes = events.data.field("HIL_TEL_SIZE") 
    allCoreX = events.data.field("COREX") 
    allCoreY = events.data.field("COREY") 
    allMSV = events.data.field("HIL_MSW")
    cogx = events.data.field("HIL_TEL_COGX")
    cogy = events.data.field("HIL_TEL_COGY")
    localDistance = sqrt( cogx**2 +cogy**2)

    # impacts distances need to be calculated for each telescope (the
    # impact distance stored is relative to the array center)
    allImpacts = np.zeros( allValues.shape )
    for ii in range(allValues.shape[1]):
        print "t",ii, " at ", tposx[ii],tposy[ii]
        allImpacts[:,ii] = np.sqrt( (allCoreX-tposx[ii])**2 +
                                    (allCoreY-tposy[ii])**2 )

    nevents,ntels = allValues.shape
    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = allValues > -100 
    telMask *= valueMask  # mask off bad values
 
    localDistMask = localDistance < 0.025
#    telMask *= localDistMask  # mask off bad values

    # now, for each event, we want to calculate the MRSW/MRSL value,
    # which is just 1/Ntelsinevent sum( (V[tel] -
    # Vmean[tel])/Vsigma[tel] )

    print "Calculating MRSV from the lookups...."

    mrsv = np.zeros( allValues.shape[0] )

    for evnum in xrange(allValues.shape[0]):
        vals = allValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        lsizes = np.log10(allSizes[evnum][ telMask[evnum] ])
        limpacts = np.log10(allImpacts[evnum][ telMask[evnum] ])
        
        if 0:
            print "-----------------"
            print "VALS:",vals
            print "TELS:",tels
            print "LSIZ:",lsizes
            print "IMPT:",limpacts
            print "ZIP:", zip(lsizes,limpacts)
        mrsv[evnum] = calcMeanReducedScaledValue( tels,zip(lsizes,limpacts), 
                                                  vals, lookup=telLookup )
    
    mrsv[ np.isnan(mrsv) ] = -10000
    print "Done"
    

    hist(mrsv,   range=[-2,2], bins=100, histtype="step", label="calculated MRSV")
    hist(allMSV, range=[-2,2], bins=100, histtype="step", label="Eventlist MRSV")
    legend()

    figure()
    scatter( allMSV, mrsv )
    xlim([-3,3])
    ylim([-3,3])

    figure()
    semilogy()
    hist( mrsv-allMSV, range=[-5,5],bins=50, histtype='step', label='Residuals' )
    legend()
