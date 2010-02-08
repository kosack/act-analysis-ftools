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


class TelLookupTable(object):
    """ 
    Read and interpolate telescope-wise lookup tables
    """

    _valueDict = dict()
    _proj = dict()


    def __init__(self, lookupName, telID, 
                 lookupDir= "/home/kosack/Analysis/FITSEventLists/Lookups"):
        """
        initialize a lookup table

        TODO: eventually these will be 4D data cubes probably! (for zenith
        and azimuth dependence)

        `lookupName`: name of lookup table to load (e.g. HIL_TEL_WIDTH)
        `telID`: telescope ID number
        """

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

    def getValue(self, coord, tel, what="mean" ):
        """
        
        Arguments:
        - `tel`: telescope id
        - `coord`: coordinate or array of coordinates  in logSize,Dist space
        - `impactDist`: single value or array of values

        returns value in the lookup table
        """
        
        world = np.array( coord )
        bin = np.trunc(self._proj[what].topixel( world ) - 1.0).astype(int)
        xed,yed,val = self._valueDict[what]
        
        print "bin=",bin, " value=",va;[bin]


        



def calcMeanReducedScaledValue( tels, coords, vals, lookup):
    """
    
    `tels`: array of telescope-ids
    `vals`: array of values corresponding to each tel
    `lookup`: the big lookup table structure
    """
    
    print "TELS:",tels
    print "COORDS:",coords
    print "VALS:",vals

    for ii in xrange( vals.shape[0] ):
        print "\tTEL",tels[ii]
        print "\tVAL",vals[ii]
        print "\tCOORD", coords[ii]
#        lval = lookup.getValue( 

#    mrsv = np.average(vals)

    return 1.0
    
                        


if __name__ == '__main__':
    
    ineventlistfile = sys.argv.pop(1)
    evfile = pyfits.open(ineventlistfile)
    events = evfile["EVENTS"]

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
        print "--------------------------------"
        vals = allValues[evnum][ telMask[evnum] ]
        tels = telid[ telMask[evnum] ]
        sizes = allSizes[evnum][ telMask[evnum] ]
        impacts = allImpacts[evnum][ telMask[evnum] ]
        
        mrsv = calcMeanReducedScaledValue( tels,zip(sizes,impacts), 
                                           vals, lookup=telLookup )
        print zip(tels,vals,zip(sizes,impacts)),"--> ",mrsv
    
    
