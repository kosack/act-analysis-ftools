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

        return v/1000.0

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
        ylabel("log10(IMPACT)")
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
    return mrsv/float(ntels)
    

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



if __name__ == '__main__':

    debug=0
    
    ineventlistfile = sys.argv.pop(1)
    evfile = pyfits.open(ineventlistfile)
    events = evfile["EVENTS"]
    telarray = evfile["TELARRAY"]

    tposx = telarray.data.field("POSX")
    tposy = telarray.data.field("POSY")
    telid = np.array(events.header['TELLIST'].split(",")).astype(int)

    # load the lookups:
#    lookupName = "HIL_TEL_WIDTH"
    lookupName = "MC_ENERGY"

    telLookup = dict()

    for tel in telid:
        telLookup[tel] = TelLookupTable(lookupName,tel)

    # THE FOLLOWING IS SIMILAR TO generate-lookup-tables (should combine them)
         
    telMask   = events.data.field("TELMASK") 
    allValues = events.data.field(lookupName) 


    try:
        allSizes = events.data.field("HIL_TEL_SIZE") 
    except KeyError:
        allSizes = events.data.field("INT_TEL_SIZE") 


    if (allValues.ndim == 1):
        # this is not a telescope-wise parameter, like ENERGY, so need
        # to make it one:
        tmp = np.ones_like(allSizes)
        for ii in xrange(tmp.shape[1]):
            tmp[:,ii] *= allValues
        allValues = tmp


    allCoreX = events.data.field("COREX") 
    allCoreY = events.data.field("COREY") 
    allMSV = events.data.field("HIL_MSW")

    # impacts distances need to be calculated for each telescope (the
    # impact distance stored is relative to the array center)
    allImpacts = np.zeros_like( allValues )
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
 
#    localDistMask = localDistance < 0.025
#    telMask *= localDistMask  # mask off bad values

    # now, for each event, we want to calculate the MRSW/MRSL value,
    # which is just 1/Ntelsinevent sum( (V[tel] -
    # Vmean[tel])/Vsigma[tel] )

    print "Calculating MRSV from the lookups for %d events " % nevents,
    print "and %d telescopes..." % ntels

    mrsv  = np.zeros( nevents )
    energy = np.zeros( nevents )
    sig_e  = np.zeros( nevents )

    for evnum in xrange(allValues.shape[0]):
        vals = allValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        lsizes = np.log10(allSizes[evnum][ telMask[evnum] ])
        limpacts = np.log10(allImpacts[evnum][ telMask[evnum] ])
        
        coords = zip(lsizes,limpacts)

        mrsv[evnum] = calcMeanReducedScaledValue( tels,coords=coords, 
                                                  vals=vals, lookupDict=telLookup )
        
        energy[evnum], sig_e[evnum] = calcWeightedAverage(tels,coords=coords, 
                                                          vals=vals, 
                                                          lookupDict=telLookup )
        
    mrsv[ np.isnan(mrsv) ] = -10000
    print "Done"
    
    

    if (debug):
        from pylab import *

        hist(mrsv, range=[-5,5], bins=100, histtype="step", label="calculated MRSV")
        hist(allMSV,range=[-5,5], bins=100, histtype="step", label="Eventlist MRSV")
        legend()

        figure()
        h,x,y = histogram2d( mrsv, allMSV, range=[[-2,5],[-2,5]], bins=[100,100]  )
        pcolormesh( x,y,h.transpose() )

        figure()
        semilogy()
        hist( mrsv-allMSV,range=[-8,8],bins=50, histtype='step', label='Residuals' )
        legend()

        # more tests: a cut on mrsv:
        cut =  (mrsv > -2.0) * (mrsv < 0.5)
        cut2 = (allMSV > -2.0) * (allMSV<0.5)
        X = events.data.field("DETX")
        Y = events.data.field("DETY")
        posraw,x,y = histogram2d( X,Y,range=[[-4,4],[-4,4]], bins=[200,200] )
        poscut,x,y = histogram2d( X[cut],Y[cut],
                                  range=[[-4,4],[-4,4]], bins=[200,200])
        poscut2,x,y = histogram2d( X[cut2],Y[cut2], 
                                   range=[[-4,4],[-4,4]], bins=[200,200])
        figure()
        pcolor( x,y, posraw )
        title("Raw")
        colorbar()
        figure()
        pcolor( x,y,poscut)
        title ("With Cut" )
        colorbar()
        figure()
        pcolor( x,y,poscut2)
        title ("With Cut (orig MRSV)" )
        colorbar()

