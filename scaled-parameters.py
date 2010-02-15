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


    def __init__(self, lookupName, telID, 
                 lookupDir= None):
        """
        initialize a lookup table

        TODO: eventually these will be 4D data cubes probably! (for zenith
        and azimuth dependence)

        `lookupName`: name of lookup table to load (e.g. HIL_TEL_WIDTH)
        `telID`: telescope ID number
        """
        if lookupDir==None:
            lookupDir = os.environ["HOME"]+"/Analysis/FITSEventLists/Lookups"

        values = ["mean","stddev","count"]    
        self._valueDict = dict()

        print "Loading Lookups for CT",tel

        for val in values:
            fname = "CT%03d-%s-lookup-%s.fits" % (tel,lookupName,val)    
            print "\t",fname
            hdu = pyfits.open(lookupDir+"/"+fname)[0]
            edgesAndHist = actutils.histFromFITS( hdu )
            self._valueDict[val] = edgesAndHist
            self._proj[val] = wcs.Projection( hdu.header )

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
        v = val[bin[0],bin[1]]

        #print "coord=",coord,"bin=",bin, " value=",v

        return v
        
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
    allMSW = events.data.field("HIL_MSW")
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
        impacts = allImpacts[evnum][ telMask[evnum] ]
        
        mrsv[evnum] = calcMeanReducedScaledValue( tels,zip(lsizes,impacts), 
                                                  vals, lookup=telLookup )
    
    mrsv[ np.isnan(mrsv) ] = -10000
    print "Done"
    

    hist(mrsv, range=[-2,2], bins=50,histtype="step", label="calculated MRSV")
    hist(allMSW, range=[-2,2], bins=50,histtype="step", label="Eventlist MRSV")
    legend()

    figure()
    scatter( allMSW, mrsv )
    xlim([-3,3])
    ylim([-3,3])
