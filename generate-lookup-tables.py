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

import actutils
from fitshistogram import Histogram
 
# Note this assumes that the local-distance and min SIZE cuts have
# already been appluied before generating the eventlist (no
# localDistance cut is made here)


def generateTelLookupTables(events,varName="HIL_TEL_WIDTH",
                            bins = [100,100], histrange =[[0.5,6],[0,2000.0]],
                            debug=False, namebase=None, 
                            valueScale=1.0, useLogScale=False):
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

    print "-------------------------------------"
    print infile

    tposx = telarray.data.field("POSX")
    tposy = telarray.data.field("POSY")
    telid = np.array(events.header['TELLIST'].split(",")).astype(int)
    
    # these are the data fields as NxM arrays (where N is number of events,
    # and M is number of telescopes):
    telMask   = events.data.field("TELMASK") 
    telSizes = events.data.field("HIL_TEL_SIZE") 
    telCoreX = events.data.field("COREX") 
    telCoreY = events.data.field("COREY") 
    mccorex = events.data.field("MC_COREX") 
    mccorey = events.data.field("MC_COREY") 
    nevents,ntels = telSizes.shape

    telValues = events.data.field(varName) 

    if (useLogScale):
        telValues = np.log10(telValues)

    if (telValues.ndim == 1):
        # this is not a telescope-wise parameter, like ENERGY, so need
        # to make it one:
        tmp = np.ones_like(telSizes)
        for ii in xrange(tmp.shape[1]):
            tmp[:,ii] *= telValues
        telValues = tmp



    cogx = events.data.field("HIL_TEL_COGX")
    cogy = events.data.field("HIL_TEL_COGY")

    # impacts distances need to be calculated for each telescope (the
    # global impact distance stored is relative to the array center)
    telImpacts = np.zeros_like( telValues )

    for itel in range(telValues.shape[1]):
        nev = telImpacts.shape[0]
        telImpacts[:,itel] = np.sqrt( (telCoreX-tposx[itel])**2 +
                                      (telCoreY-tposy[itel])**2 )


    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = telValues > -100 
    telMask *= valueMask  # mask off bad values
    
    if (debug):
        figure( figsize=(15,10))
        subplot(1,1,1)

    # now, generate the lookup table for each telescope separately:
    for itel in range(ntels):

        goodEvents = telMask[:,itel]  
        value = telValues[:,itel][goodEvents]  * valueScale # scale to mrad
        impact = telImpacts[:,itel][goodEvents]
        logsiz = np.log10(telSizes[:,itel][goodEvents])

        sumHist    = Histogram( range=histrange, bins=bins,name="SUM" )
        sumSqrHist = Histogram( range=histrange, bins=bins,name="SUMSQR" )
        countHist  = Histogram( range=histrange, bins=bins,name="COUNT" )

        sumHist.fill(    (logsiz,impact), weights=value, normed=False )
        sumSqrHist.fill( (logsiz,impact), weights=value**2, normed=False )
        countHist.fill(  (logsiz,impact), normed=False )
        
        # write it out as a FITS file with 2 image HDUs VALUE and SIGMA

        if namebase:
            filename = "%s-CT%03d-%s-lookup" % (namebase,telid[itel], varName)
        else:
            filename = "CT%03d-%s-lookup" % (telid[itel], varName)

        print "CT",telid[itel],varName,"Nevents=",len(value),
        print "outliers=",len(value)-sum(countHist.hist),
        print "out:",filename

        hdu1=sumHist.asFITS()
        hdu2=sumSqrHist.asFITS()
        hdu3=countHist.asFITS()

        # hdu1=actutils.histToFITS( sumHist, bins=bins,
        #                           histrange=histrange,name="SUM",
        #                           valueScale=valueScale )
        # hdu2=actutils.histToFITS( sumSqrHist, bins=bins,
        #                           histrange=histrange, name="SUMSQR",
        #                           valueScale=valueScale )
        # hdu3=actutils.histToFITS( countHist, bins=bins,
        #                           histrange=histrange, name="COUNTS" )


        hdu1.writeto( filename+"-sum.fits", clobber=True )
        hdu2.writeto( filename+"-sum2.fits", clobber=True )
        hdu3.writeto( filename+"-count.fits", clobber=True )
        

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


 
