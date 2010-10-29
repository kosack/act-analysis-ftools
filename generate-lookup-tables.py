import pyfits
import numpy as np
from optparse import OptionParser
import sys
import math
from scipy import interpolate
from scipy import spatial
import scipy.signal
import re
from pylab import *
import os.path

import actutils
from fitshistogram import Histogram
 
# Note this assumes that the local-distance and min SIZE cuts have
# already been appluied before generating the eventlist (no
# localDistance cut is made here)

# old histrange =[[0.5,6],[0,2000.0]],

def generateTelLookupTables(events,varName="HIL_TEL_WIDTH",
                            bins = [100,100], 
                            histrange =[[0.5,6.0],[0,2000.0]],
                            debug=False, namebase=None, 
                            valueScale=1.0, useLogScale=False,
                            includeTelTypeInFilename=True):
    """
    generates lookup table for the given variable (average and sigma
    as a function of logSIZE and IMPACT DISTANCE)

    Arguments:
    - `events`: event-list FITS HDU 
    - `varName`: variable to histogram (HIL_TEL_WIDTH or HIL_TEL_LENGTH)
    - `bins`: number of bins for logSIZE,DISTANCE
    - `histrange`: logSIZE and DISTANCE ranges

    - `valueScale`: scale factor to multiply the value by, if
      requested (to keep it in a nice range)

    - `useLogScale`: apply log10 function to the variable being
      histogrammed (e.g. for MC_ENERGY)
    """

    evfile = pyfits.open(infile)
    events = evfile["EVENTS"]
    telarray = evfile["TELARRAY"]

    # these are the data fields as NxM arrays (where N is number of events,
    # and M is number of telescopes):    
    telImpacts,telSizes,telid,telMask=actutils.loadLookupTableColumns( events, telarray )

    print "-------------------------------------"
    print "generate-lookup-tables: "
    print "DIRECTORY:",os.path.dirname(infile)
    print "EVENTLIST:",os.path.basename(infile)

    nevents,ntels = telSizes.shape

    telValues = events.data.field(varName) 

    if (useLogScale):
        telValues = np.log10(telValues)

    # TEST:::::::::::::::::::::::::
#    telValues /= telSizes

    if (telValues.ndim == 1):
        # this is not a telescope-wise parameter, like ENERGY, so need
        # to make it one:
        tmp = np.ones_like(telSizes)
        for ii in xrange(tmp.shape[1]):
            tmp[:,ii] *= telValues
        telValues = tmp


    # exclude some bad values
    telMask *= telValues > -1000
    
    if (debug):
        figure( figsize=(15,10))
        subplot(1,1,1)

    tel2type,type2tel = actutils.getTelTypeMap( telarray )

    # now, generate the lookup table for each telescope separately:
    for itel in range(ntels):

        goodEvents = telMask[:,itel]  
        value = telValues[:,itel][goodEvents]  * valueScale # scale to mrad
        impact = telImpacts[:,itel][goodEvents]
        size = telSizes[:,itel][goodEvents]

        names = ["LSIZE","DIST"]
        sumHist    = Histogram( range=histrange, bins=bins,name="SUM", 
                                axisNames=names )
        sumSqrHist = Histogram( range=histrange, bins=bins,name="SUMSQR",
                                axisNames=names )
        countHist  = Histogram( range=histrange, bins=bins,name="COUNT",
                                axisNames=names )

        sumHist.fill(    (size,impact), weights=value, normed=False )
        sumSqrHist.fill( (size,impact), weights=value**2, normed=False )
        countHist.fill(  (size,impact), normed=False )
        
        # write it out as a FITS file with 2 image HDUs VALUE and SIGMA

        tag = "CT%03d" % telid[itel]
        if (includeTelTypeInFilename):
            tag += "-TYPE%02d_%02d" % tel2type[telid[itel]]
        
        if namebase:
            filename = "%s-%s-%s-lookup" % (namebase,tag, varName)
        else:
            filename = "CT%03d-%s-lookup" % (telid[itel], varName)

        print "CT",telid[itel],varName,"Nevents=",len(value),
        print "outliers=",len(value)-sum(countHist.hist)
#        print "out:",filename

        sumHist.asFITS().writeto( filename+"-sum.fits", clobber=True )
        sumSqrHist.asFITS().writeto( filename+"-sum2.fits", clobber=True )
        countHist.asFITS().writeto( filename+"-count.fits", clobber=True )

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.set_usage( "generate-lookup-table.py [options] EVENTFILE")
    parser.add_option( "-d","--debug", dest="debug", action="store_true",
                       default=False, help="display debug plots interactively")    
    parser.add_option( "-o","--output", dest="output", 
                       default=None, 
                       help="output name base (rest of name will be appended)")
    (opts, args) = parser.parse_args()

    if len(args)==0:
        parser.error("Please specify a filename")

    infile = args.pop(0)

    if (opts.output==None):
        match = re.search( r"_([\d]+)_", infile )
        if (match):
            runid = match.group(1)
            namebase = "run_%s" % runid
        else:
            namebase=None
    else:
        namebase = opts.output
    

    generateTelLookupTables( infile, varName="HIL_TEL_WIDTH", 
                             debug=opts.debug,namebase=namebase,valueScale=1000.0 )

    generateTelLookupTables( infile, varName="HIL_TEL_LENGTH" , 
                             debug=opts.debug, namebase=namebase,valueScale=1000.0 )

    generateTelLookupTables( infile, varName="MC_ENERGY", 
                             debug=opts.debug, namebase=namebase,
                             valueScale=1.0, useLogScale=True)


 
