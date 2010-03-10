import pyfits
import numpy as np
import sys
import math
import os
import re
from scipy import interpolate

from optparse import OptionParser
from kapteyn import wcs

import actutils
from fitshistogram import Histogram



# TODO apply image amplitude and local distance cuts!

# note can use wcs.coordmap() to generate the coordinate map input to
# scipy.interpolate if we want to do interpolation properly (see Kapteyn manual)


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()



class TelLookupTable(object):
    """ 
    Read and interpolate telescope-wise lookup tables
    """

    _valueDict = dict()
    values = dict()
    telID = 0

    def __init__(self, lookupName, telID, valueScale=1.0,
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
        self._valueScale = valueScale

        print ("Loading Lookups for CT%03d"%tel).center(79,"-")

        smooth=""
        if useSmooth: 
            smooth="-gauss"

        for what in self.values:
            fname = "CT%03d-%s-lookup-%s%s.fits" % (tel,lookupName,what,smooth)    
            print "\t",fname
            hdu = pyfits.open(lookupDir+"/"+fname)[0]
            hist = Histogram(initFromFITS=hdu)
            self._valueDict[what] = hist

    

        self.extrapolateLookups()
#        self.smoothLookups(3)

    def extrapolateLookups(self, minCounts=10):
        """
        Make the lookup tables nicer by extrapolating unfilled values
        
        Arguments:
        - `self`:
        """

        for what in self.values:
            counts = self._valueDict["count"].hist.copy()

            if what=='count': continue
            print "\tExtrapolating:", what
            val = self._valueDict[what].hist

            # todo: do this as an array op, not for-loop

            for ii in xrange(counts.shape[0]):
                lastgood = (-1,-1)
                for jj in xrange(counts.shape[1]):
                    if counts[ii,jj] > minCounts and np.isfinite(val[ii,jj]):
                        lastgood=(ii,jj)
                    else:
                        if (lastgood != (-1,-1)):
                            val[ii,jj] = val[lastgood]
                            counts[ii,jj]=minCounts+1

            for ii in xrange(counts.shape[0]):
                lastgood = (-1,-1)
                for jj in xrange(counts.shape[1]-1,0,-1):
                    if counts[ii,jj] > minCounts and np.isfinite(val[ii,jj]):
                        lastgood=(ii,jj)
                    else:
                        if (lastgood != (-1,-1)):
                            val[ii,jj] = val[lastgood]
                            counts[ii,jj]=minCounts+1

            for jj in xrange(counts.shape[1]):
                lastgood = (-1,-1)
                for ii in xrange(counts.shape[0]):
                    if counts[ii,jj] > minCounts and np.isfinite(val[ii,jj]):
                        lastgood=(ii,jj)
                    else:
                        if (lastgood != (-1,-1)):
                            val[ii,jj] = val[lastgood]
                            counts[ii,jj]=minCounts+1

            for jj in xrange(counts.shape[1]):
                lastgood = (-1,-1)
                for ii in xrange(counts.shape[1]-1,0,-1):
                    if counts[ii,jj] > minCounts and np.isfinite(val[ii,jj]):
                        lastgood=(ii,jj)
                    else:
                        if (lastgood != (-1,-1)):
                            val[ii,jj] = val[lastgood]
                            counts[ii,jj]=minCounts+1


    def smoothLookups(self, pixels):
        """
        smooth the lookup tables with a gaussian filter
        """
        for what in self.values:
            print "\tSmoothing:",what
            im = self._valueDict[what].hist 
            gauss = gauss_kern(pixels)
            self._valueDict[what].hist = signal.convolve(im,gauss, mode='same')


    def getValue(self, coord, what="mean" ):
        """
        
        Arguments:
        - `coord`: coordinate or array of coordinates  in logSize,Dist space
        - `what`: what type of value (mean, stddev)

        returns value in the lookup table
        """

        return (self._valueDict[what].getValue( coord,outlierValue=-10000 )
                / self._valueScale)


    def display(self, what="mean"):
        """
    
        Arguments:
        - `self`:
        - `what`:
        """

        val = self._valueDict[what].hist

        if (len(val.shape)==2):
            xed,yed =  self._valueDict[what].binLowerEdges
            figure()
            pcolormesh( xed, yed, val.transpose())
            title("CT%d: %s" % (self.telID, what))
            xlabel("log10(SIZE)")
            ylabel("IMPACT (m)")
            colorbar()
            show()                   
        
# ===========================================================================

def calcMeanReducedScaledValue( tels, coords, vals, lookupDict,debug=0):
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

    mrsv = 0.0
    ntels = vals.shape[0]
    EPSILON = 1e-12

    if ntels==0:
        return -10000

    if debug:
        print "========================",tels

    for itel in xrange( ntels ):
        vMean = lookupDict[tels[itel]].getValue( coords[itel], "mean" )
        vSigma = lookupDict[tels[itel]].getValue( coords[itel], "stddev" )

        if (np.isfinite(vMean) == False or np.isfinite(vSigma)==False
            or np.abs(vSigma <EPSILON)):
            ntels -= 1 # skip telescopes that don't have a good value
            continue

        mrsv += (vals[itel] - vMean)/vSigma

        if debug:
            print "CT",tels[itel],"coord=",coords[itel],
            print " val=",vals[itel],"mean=",vMean,"sig=",vSigma

    if debug:
        print "-------------------------"
        print "     MRSV=",mrsv

    if (ntels >0):
        return mrsv/float(ntels), 0.0
    else:
        return (-100000,-100000)
    

def calcWeightedAverage( tels, coords, vals, lookupDict,debug=0):
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

    EPSILON = 1e-12
    ntels = vals.shape[0]

    wsum = 0 # weighted sum
    sumw = 0  # sum of weights

    for itel in xrange(ntels):
        vMean  = lookupDict[tels[itel]].getValue( coords[itel], "mean" )
        vSigma = lookupDict[tels[itel]].getValue( coords[itel], "stddev" )

        if (np.isfinite(vMean) == False or np.isfinite(vSigma)==False
            or np.abs(vSigma <EPSILON)):
            ntels -= 1 # skip telescopes that don't have a good value
            continue
        
        weight = 1/vSigma**2
        wsum += vMean*weight
        sumw += weight

    if (math.fabs(sumw) > EPSILON):
        wmean = wsum/sumw
        wstddev = math.sqrt( sumw )
        return wmean, wstddev
    else:
        return (-100000,-100000)


def testValue(value,error,trueValue):
    """
    display some tests

    Arguments:
    - `value`: array of values calculated 
    - `error`: array of errors
    - `trueValue`: array true value for comparison
    """
    
    import pylab

    percentError = (trueValue-value)/trueValue 
    
    pylab.figure()
    pylab.hist( percentError, range=[-5,5], bins=50 )

    pylab.figure()
    H,x,y = np.histogram2d( trueValue,value,
                            range=[[-4,4],[-4,4]], bins=[100,100] )
    pylab.pcolormesh( x,y, H)
    pylab.plot(x,x)

#    figure()
#    scatter( trueValue, percentError ) #range=[[-4,4],[-4,4]]
    print x,y



if __name__ == '__main__':

    from pylab import *

    parser = OptionParser()
    parser.add_option("-t","--type", dest="paramType", 
                      help="column to generate (energy, msw, msl")
    (opts,args) = parser.parse_args()

    debug=0
    paramType = opts.paramType

    # open the input eventlist

    ineventlistfile = args.pop(0)
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
        outputName = "ENERGY2"
        reductionFunction=calcWeightedAverage
    elif (paramType == "msl"):
        lookupName = "HIL_TEL_LENGTH"
        outputName = "HIL_MSL2"
        reductionFunction=calcMeanReducedScaledValue
        valueScale=1000.0
    elif (paramType == "msw") :
        lookupName = "HIL_TEL_WIDTH"
        outputName = "HIL_MSW2"
        reductionFunction=calcMeanReducedScaledValue
        valueScale=1000.0
    else:
        print "Unknown type:",paramType
        sys.exit(1)

    # load the lookups:
    telLookup = dict()

    for tel in telid:
        telLookup[tel] = TelLookupTable(lookupName,tel, valueScale=valueScale)
        #telLookup[tel].display()

    telLookup[1].display("mean")
    telLookup[1].display("stddev")
    telLookup[1].display("count")


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

    valueMask = (telValues > -100) * isfinite(telValues)
    if ((valueMask==False).any()):
        print "WARNING: %d values were undefined?" % sum(valueMask==False)
#    telMask *= valueMask  # mask off bad values

    # now, for each event, we want to calculate the MRSW/MRSL value,
    # which is just 1/Ntelsinevent sum( (V[tel] -
    # Vmean[tel])/Vsigma[tel] )

    print "="*70
    print "Calculating", outputName,"from the",lookupName,
    print "lookups for %d events" % nevents
    print "and %d telescopes" % ntels,"using",reductionFunction.func_name,"..."
    print "="*70

    value  = np.zeros( nevents ) # the reduced value
    error  = np.zeros( nevents ) # the error on the reduced value

    evdebug=0

    for evnum in xrange(telValues.shape[0]):
        vals = telValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        lsizes = np.log10(telSizes[evnum][ telMask[evnum] ])
        impacts = telImpacts[evnum][ telMask[evnum] ]
        coords = zip(lsizes,impacts)

        if(evnum==421): evdebug=1
        if(evnum==422): evdebug=0

        # call the appropriate reduction function (defined earlier)
        value[evnum],error[evnum] = reductionFunction( tels,
                                                       coords=coords, 
                                                       vals=vals, 
                                                       lookupDict=telLookup,
                                                       debug=evdebug)
        if (evnum%1000==0):
            print "    event",evnum,"of",telValues.shape[0]
        
    value[ np.isnan(value) ] = -10000
    print "Done."
    

    gmask = value > -1000
    gmask *= value < 1000
    
    if (paramType == "energy"):
        testValue( value[gmask], error[gmask], 
                   log10(events.data.field("ENERGY")[gmask]) )
    elif(paramType=="msw"):
        testValue( value[gmask], error[gmask], 
                   (events.data.field("HIL_MSW")[gmask]) )


    if (paramType == "energy"):
        value = 10**value # not log scale
    
    outputCol = pyfits.Column( name=outputName, format='D', array=value )
    cols = events.columns + outputCol
    outputTable = pyfits.new_table( cols )
    outputTable.writeto("testtable.fits", clobber=True)
    
