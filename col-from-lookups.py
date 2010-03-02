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
                 lookupDir= None,useSmooth=False):
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
        bin[bin>=shape] = shape-1
        bin[bin<0] = 0
        v = val[bin[0],bin[1]] # correct order?

        #print "coord=",coord,"bin=",bin, " value=",v

        return v

    def display(self, what="mean"):
        """
    
        Arguments:
        - `self`:
        - `what`:
        """

        xed,yed,val = self._valueDict[what]
        figure()
        pcolormesh( xed, yed, val.transpose())
        title("CT%d: %s" % (self.telID, what))
        xlabel("log10(SIZE)")
        ylabel("IMPACT (m)")
        colorbar()
        show()           
        
        
# ===========================================================================

def calcMeanReducedScaledValue( tels, coords, vals, lookupDict):
    """
    Return the mean-reduced-scaled value using the lookup tables. The
    MRSV is defined as \sum_i (x_i - \bar{x}) / \sigma^i_x
    
    This is used e.g. for calculating mean-reduced-scaled length and width values

    `tels`: array of telescope-ids
    `coords`: array of the dependent coordinates (e.g. [log(size), log(impact)])
    `vals`: array of values corresponding to each tel
    `lookup`: dictionary containing lookup table for each telescope
    for the given parameter
    """

    debug=0
    mrsv = 0.0
    ntels = vals.shape[0]

    if ntels==0:
        return -10000

    if debug:
        print "========================",tels

    for itel in xrange( ntels ):
        vMean = lookupDict[tels[itel]].getValue( coords[itel], "mean" )
        vSigma = lookupDict[tels[itel]].getValue( coords[itel], "stddev" )
        mrsv += (vals[itel] - vMean)/vSigma

        if debug:
            print "CT",tels[itel],"coord=",coords[itel],
            print " val=",vals[itel],"mean=",vMean,"sig=",vSigma

    if debug:
        print "-------------------------"
        print "     MRSV=",mrsv
    return mrsv/float(ntels), 0.0
    

def calcWeightedAverage( tels, coords, vals, lookupDict):
    """ 
    Returns the weighted average over telescopes (and the weighted
    standard deviation) of the parameter. This is used, e.g., for
    energy reconstruction.
    
   `tels`: array of telescope-ids
    `coords`: array of the dependent coordinates (e.g. [log(size), log(impact)])
    `vals`: array of values corresponding to each tel
    `lookup`: dictionary containing lookup table for each telescope
    for the given parameter
 
    """
    ntels = vals.shape[0]
    sumweight = 0

    vMean = np.zeros( ntels )
    vSigma = np.zeros( ntels )

    for itel in xrange(ntels):
        vMean[itel] = lookupDict[tels[itel]].getValue( coords[itel], "mean" )
        vSigma[itel] =lookupDict[tels[itel]].getValue( coords[itel], "stddev" )

    wmean, wsum = np.average( vMean, weights = 1.0/vSigma**2, returned=True)
    sumsqr = np.sum( vMean**2 )
    vMean /= wsum
    wstddev = math.sqrt( sumsqr / (wsum-np.sum(vMean)**2)  )

    return wmean, wstddev


def testValue(value,error,trueValue):
    """
    display some tests

    Arguments:
    - `value`: array of values calculated 
    - `error`: array of errors
    - `trueValue`: array true value for comparison
    """
    
    from pylab import *

    percentError = (trueValue-value)/trueValue * 100.0
    
    figure()
    hist( percentError, range=[-100,100], bins=50 )

    figure()
    H,x,y = histogram2d( trueValue,value,
                         range=[[-4,4],[-4,4]], bins=[200,200] )
    pcolormesh( x,y, H)
    plot(x,x)

    figure()
    H2,x,y = histogram2d( trueValue,(trueValue-value)/trueValue,
                         bins=[200,200] ) #range=[[-4,4],[-4,4]]
    pcolormesh( x,y, H2)
    print x,y
    

if __name__ == '__main__':

    from pylab import *

    debug=1
    paramType = "msw"

    # open the input eventlist

    ineventlistfile = sys.argv.pop(1)
    evfile = pyfits.open(ineventlistfile)
    events = evfile["EVENTS"]
    telarray = evfile["TELARRAY"]

    # load telescope information
    tposx = telarray.data.field("POSX")
    tposy = telarray.data.field("POSY")
    telid = np.array(events.header['TELLIST'].split(",")).astype(int)

    lookupName = ""
    outputName = ""
    valueScale = 1.0

    if (paramType == "energy" ):
        lookupName = "MC_ENERGY"
        outputName = "ENERGY"
        reductionFunction=calcWeightedAverage
    elif (paramType == "msl"):
        lookupName = "HIL_TEL_LENGTH"
        outputName = "HIL_MSL"
        reductionFunction=calcMeanReducedScaledValue
        valueScale=1000.0
    elif (paramType == "msw") :
        lookupName = "HIL_TEL_WIDTH"
        outputName = "HIL_MSW"
        reductionFunction=calcMeanReducedScaledValue
    else:
        print "Unknown type:",paramType
        sys.exit(1)

    # load the lookups:
    telLookup = dict()

    for tel in telid:
        telLookup[tel] = TelLookupTable(lookupName,tel)
        telLookup[tel].display()

    # THE FOLLOWING IS SIMILAR TO generate-lookup-tables (should combine them)
         
    telMask   = events.data.field("TELMASK") 
    telValues = events.data.field(lookupName) # values of the requested parameter

    try:
        telSizes = events.data.field("HIL_TEL_SIZE") 
    except KeyError:
        telSizes = events.data.field("INT_TEL_SIZE") 


    if (telValues.ndim == 1):
        # this is not a telescope-wise parameter, like ENERGY, so need
        # to make it one:
        tmp = np.ones_like(telSizes)
        for ii in xrange(tmp.shape[1]):
            tmp[:,ii] *= telValues
        telValues = tmp


    coreX = events.data.field("COREX") 
    coreY = events.data.field("COREY") 

    # impacts distances need to be calculated for each telescope (the
    # impact distance stored is relative to the array center)
    telImpacts = np.zeros_like( telValues )
    for ii in range(telValues.shape[1]):
        print "CT%03d"%ii, " at ", tposx[ii],tposy[ii]
        telImpacts[:,ii] = np.sqrt( (coreX-tposx[ii])**2 +
                                    (coreY-tposy[ii])**2 )

    nevents,ntels = telValues.shape
    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = telValues > -100 
    telMask *= valueMask  # mask off bad values

    # now, for each event, we want to calculate the MRSW/MRSL value,
    # which is just 1/Ntelsinevent sum( (V[tel] -
    # Vmean[tel])/Vsigma[tel] )

    print "Calculating", outputName,"from the",lookupName,
    print "lookups for %d events" % nevents,
    print "and %d telescopes" % ntels,"using",reductionFunction.func_name,"..."

    value  = np.zeros( nevents ) # the reduced value
    error  = np.zeros( nevents ) # the error on the reduced value

    for evnum in xrange(telValues.shape[0]):
        vals = telValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        lsizes = np.log10(telSizes[evnum][ telMask[evnum] ])
        impacts = telImpacts[evnum][ telMask[evnum] ]
        coords = zip(lsizes,impacts)

        # call the appropriate reduction function (defined earlier)
        value[evnum],error[evnum] = reductionFunction( tels,
                                                       coords=coords, 
                                                       vals=vals, 
                                                       lookupDict=telLookup )
        
    value[ np.isnan(value) ] = -10000
    print "Done."
    

    gmask = value > -1000
    gmask *= value < 10000
#    testValue( value[gmask], error[gmask], 
#               np.log10(events.data.field("MC_ENERGY")[gmask]) )

    testValue( value[gmask], error[gmask], 
               (events.data.field("HIL_MSW")[gmask]*valueScale) )
