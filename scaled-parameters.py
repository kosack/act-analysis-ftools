import pyfits
import numpy as np
import sys
import math
import re
from scipy import interpolate
import scipy.signal
from optparse import OptionParser
from kapteyn import wcs

import actutils


# TODO: make LookupTable a class, and we can then load from it:
#    loadFromFile()
#    interpolate( tels, vals )

def loadLookups(name,tel):
    """
    
    TODO: eventually these will be 4D data cubes probably! (for zenith
    and azimuth dependence)

    Arguments:
    - `name`:
    - `tel`:
    """

    lookupDir = "/home/kosack/Analysis/FITSEventLists/Lookups"

    values = ["mean","stddev","count"]    

    valueDict = dict()

    print "Loading Lookups for CT",tel

    for val in values:
        fname = "CT%03d-%s-lookup-%s.fits" % (tel,lookupName,val)    
        print "\t",fname
        hdu = pyfits.open(lookupDir+"/"+fname)[0]
        edgesAndHist = loadFITSHistogram( hdu )
        valueDict[val] = edgesAndHist

    return valueDict
        
        

def loadFITSHistogram(hdu):
    """
    returns xed,yed,Hist (as in np.histogram2d)

    Arguments:
    - `hdu`: hdu containing the histogram
    """
    
    hist = hdu.data
    bins = hdu.data.shape
    proj = wcs.Projection ( hdu.header ) # for going between world and pix coords
    
    xed,tmp = proj.toworld( (np.arange(bins[0])+0.5, np.zeros(bins[1])+0.5) )
    tmp,yed = proj.toworld( (np.zeros(bins[0])+0.5,np.arange(bins[1])+0.5) )
    
    return xed,yed, hist.transpose() # transposed according to
                                     # np.histogram2d convention


def calcMeanReducedScaledValue( tels, vals, lookup):
    """
    
    `tels`: array of telescope-ids
    `vals`: array of values corresponding to each tel
    `lookup`: the big lookup table structure
    """
    
    mrsv = np.average(vals)

    return mrsv
    
                        


if __name__ == '__main__':
    
    ineventlistfile = sys.argv.pop(1)
    evfile = pyfits.open(ineventlistfile)
    events = evfile["EVENTS"]

    telid = np.array(events.header['TELLIST'].split(",")).astype(int)

    # load the lookups:
    lookupName = "HIL_TEL_WIDTH"


    telLookup = dict()

    for tel in telid:
        telLookup[tel] = loadLookups(lookupName,tel)

    # now we have for example (note H[y,x] is the correct order):
    xed, yed, H = telLookup[2]["stddev"]

    # THE FOLLOWING IS SIMILAR TO GENERATE LOOKP TABLES (should combine them)
         
    telMask   = events.data.field("TELMASK") 
    allValues = events.data.field(lookupName) 
    allSizes = events.data.field("HIL_TEL_SIZE") 
    allCoreX = events.data.field("COREX") 
    allCoreY = events.data.field("COREY") 
    allImpacts = np.sqrt( allCoreX**2 + allCoreY**2 )
    nevents,ntels = allValues.shape
    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = allValues > -100 
    telMask *= valueMask  # mask off bad values
 

    # now, for each event, we want to calculate the MRSW/MRSL value, which is just
    # 1/Ntelsinevent sum( (V[tel] - Vmean[tel])/Vsigma[tel] )

    for evnum in xrange(allValues.shape[0]):
        vals = allValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]

        mrsv = calcMeanReducedScaledValue( tels,vals, lookup=telLookup )
        print zip(tels,vals),"--> ",mrsv
    
    
