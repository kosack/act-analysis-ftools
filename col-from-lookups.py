import pyfits
import numpy as np
import sys
import math
import os
import re
from scipy import interpolate,signal

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
    
    by default, this loads `<TYPE>-<VALUE>-lookup-<param><tag>.fits`
    where tag is generally blank, and <TYPE> is either the telescope
    id or the telescope type (class,subclass)
    
    `tag`: extra text to add (like "-gauss" for externally smoothed lookups)
    """

    _valueDict = dict()
    values = dict()
    telID = 0

    def __init__(self, lookupName, telID, valueScale=1.0,
                 lookupDir= None,tag="", byTelType=False):
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

        for what in self.values:
            fname = "CT%03d-%s-lookup-%s%s.fits" % (tel,lookupName,what,tag)    
            if byTelType:
                fname = "TYPE%02d_%02d-%s-lookup-%s%s.fits" % (telID[0],
                                                               telID[1],
                                                               lookupName,what,tag)

            print "\t",fname
            hdu = pyfits.open(lookupDir+"/"+fname)[0]
            hist = Histogram(initFromFITS=hdu)
            self._valueDict[what] = hist

            

        self.extrapolateLookups()
#        self.smoothLookups(1) # note smoothing seems to cause an energy bias

    def extrapolateLookups(self, minCounts=5):
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
            # can do it with numpy's apply_along_axis() function!

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

        hdu = self._valueDict[what].asFITS()
        actutils.displayFITS( hdu.header, hdu.data ) 

        # if (len(val.shape)==2):
        #     xed,yed =  self._valueDict[what].binLowerEdges
        #     figure()
        #     pcolormesh( xed, yed, val.transpose())
        #     title("CT%d: %s" % (self.telID, what))
        #     xlabel("log10(SIZE)")
        #     ylabel("IMPACT (m)")
        #     colorbar()
        #     show()                   
        
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

    # TODO: modify this to do the whole thing at once? E.g., calculate
    # the vMean and vSigma for *all* events first, and then use them
    # to get the weighted sum

    EPSILON = 1e-12
    ntels = vals.shape[0]

    if ntels == 0:   
        return (-100000,-100000)

    if debug:
        print "========================",tels

    vMean = array([lookupDict[tels[itel]].getValue( coords[itel], "mean" )[0] 
                   for itel in xrange(ntels)])

    vSigma = array([lookupDict[tels[itel]].getValue( coords[itel], "stddev" )[0] 
                    for itel in xrange(ntels)])

    if debug:
        for t,c,m,s in zip(tels,coords,vMean,vSigma):
            print "CT%03d  coord="%t,c,"mean=",m,"sig=",s

    wmean,sumOfWeights = np.average( vMean, weights=1.0/vSigma, returned=True)
    wstddev = 1.0/sumOfWeights

    if (math.fabs(sumOfWeights) > EPSILON):
        if debug:
            print "-------------------------"
            print "     WeightedAvg=",wmean
            print "     WeightSigma=",wstddev
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

#    percentError = (trueValue-value)/trueValue 
#    pylab.figure()
#    pylab.hist( percentError, range=[-5,5], bins=50 )

    H = Histogram( range=[[-1,2],[-1,2]], bins=[70,70],axisNames=["log10(true)",
                                                                  "log10(reco)"])
    H.fill( (log10(trueValue),log10(value)) )
    hdu=H.asFITS()
    actutils.displayFITS( hdu.header, hdu.data )
    title("Reconstructed vs Simulated energy")
    return H

                                 
#    figure()
#    scatter( trueValue, percentError ) #range=[[-4,4],[-4,4]]
#    print x,y



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
    tel2type,type2tel = actutils.getTelTypeMap( telarray )

    telImpacts,telSizes,telid,telMask=actutils.loadLookupTableColumns( events, telarray )

    lookupName = ""
    outputName = ""
    valueScale = 1.0
    computeFromValue = True

    if (paramType == "energy" ):
        lookupName = "MC_ENERGY"
        outputName = "ENERGY"
        reductionFunction=calcWeightedAverage
        computeFromValue = False
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
        ttype = tel2type[tel]
        telLookup[tel] = TelLookupTable(lookupName,ttype, 
                                        valueScale=valueScale, byTelType=False)
        #telLookup[tel].display()

    telLookup[1].display("mean")
    telLookup[1].display("stddev")
#    telLookup[1].display("count")


    # THE FOLLOWING IS SIMILAR TO generate-lookup-tables (should combine them)
         
    if (computeFromValue==True):
        telValues = events.data.field(lookupName) # values of the
                                                  # requested
                                                  # parameter
        telMask *= telValues > -1000  # mask off some bad values
    else:
        # for e.g. energy, there's no "value" to compute from
        telValues = np.zeros_like(telSizes) # dummy

    nevents,ntels = telSizes.shape
    print "="*70
    print "Calculating", outputName,"from the",lookupName,
    print "lookups for %d events" % nevents
    print "and %d telescopes" % ntels,"using",reductionFunction.func_name,"..."
    print "="*70

    value  = np.zeros( nevents ) # the reduced value
    error  = np.zeros( nevents ) # the error on the reduced value

    evdebug=0


    for evnum in xrange(nevents):
        vals = telValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        sizes = telSizes[evnum][ telMask[evnum] ]
        impacts = telImpacts[evnum][ telMask[evnum] ]
        coords = zip(sizes,impacts)

        if(evnum==4): evdebug=1
        if(evnum==10): evdebug=0

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
#    gmask *= error < 0.5

    if (paramType == "energy"):
        value = 10**value # want it to be not log scale
    
    if (paramType == "energy"):
        testh = testValue( value[gmask], error[gmask], 
                           events.data.field("MC_ENERGY")[gmask] )
    elif(paramType=="msw"):
        testValue( value[gmask], error[gmask], 
                   (events.data.field("HIL_MSW")[gmask]) )

        
    
    # check if column already exists, if so rename it
    coldic = dict()
    for ii in range(len(events.columns)):
        coldic[events.columns[ii].name] = ii

    if coldic.has_key(outputName):
        icol = coldic[outputName]
        events.columns[icol].name=outputName+"_ORIG"

    outputCol = pyfits.Column( name=outputName, format='D', array=value )
    cols = events.columns + outputCol
    outputTable = pyfits.new_table( cols )
    outputTable.writeto("testtable.fits", clobber=True)
    
