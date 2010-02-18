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
 
# Note this assumes that the local-distance and min SIZE cuts have
# already been appluied before generating the eventlist (no
# localDistance cut is made here)


def generateTelLookupTables(events,varName="HIL_TEL_WIDTH",
                            bins = [60,60], histrange =[[0.5,5],[0,5.0]],
                            debug=False, namebase=None, 
                            singleFile=False):
    """
    generates lookup table for the given variable (average and sigma
    as a function of logSIZE and IMPACT DISTANCE)

    Arguments:
    - `events`: event-list FITS HDU 
    - `varName`: variable to histogram (HIL_TEL_WIDTH or HIL_TEL_LENGTH)
    - `bins`: number of bins for logSIZE,DISTANCE
    - `histrange`: logSIZE and DISTANCE ranges
    - `singleFile`: write to a single FITS file  with multiple HDUs
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
    allValues = events.data.field(varName) 
    allSizes = events.data.field("HIL_TEL_SIZE") 
    allCoreX = events.data.field("COREX") 
    allCoreY = events.data.field("COREY") 
    mccorex = events.data.field("MC_COREX") 
    mccorey = events.data.field("MC_COREY") 

    cogx = events.data.field("HIL_TEL_COGX")
    cogy = events.data.field("HIL_TEL_COGY")

    # impacts distances need to be calculated for each telescope (the
    # global impact distance stored is relative to the array center)
    allImpacts = zeros( allValues.shape )
    allImpactstest = zeros( allValues.shape )
    for itel in range(allValues.shape[1]):
        nev = allImpacts.shape[0]
        allImpacts[:,itel] = np.sqrt( (allCoreX-tposx[itel])**2 +
                                      (allCoreY-tposy[itel])**2 )


    nevents,ntels = allValues.shape

    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = allValues > -100 
    telMask *= valueMask  # mask off bad values
#    telMask *= localDistMask  # mask off bad values
    
    if (debug):
        figure( figsize=(15,10))
        subplot(1,1,1)

    # now, generate the lookup table for each telescope separately:
    for itel in range(ntels):

        goodEvents = telMask[:,itel]  
        value = allValues[:,itel][goodEvents] * 1000.0 # scale to mrad
        logimpact = np.log10(allImpacts[:,itel][goodEvents])
        logsiz = np.log10(allSizes[:,itel][goodEvents])

        sumHist,edX,edY = np.histogram2d( logsiz,logimpact, 
                                           weights=value,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        sumSqrHist,edX,edY = np.histogram2d( logsiz,logimpact, 
                                           weights=(value)**2,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        countHist,edX,edY = np.histogram2d( logsiz,logimpact, 
                                            weights=None,
                                            range=histrange, 
                                            bins=bins,
                                            normed=False)


        
        # write it out as a FITS file with 2 image HDUs VALUE and SIGMA

        if namebase:
            filename = "%s-CT%03d-%s-lookup" % (namebase,telid[itel], varName)
        else:
            filename = "CT%03d-%s-lookup" % (telid[itel], varName)

        print "CT",telid[itel],varName,"Nevents=",len(value),
        print "outliers=",len(value)-sum(countHist),
        print "out:",filename

        if (singleFile) :
            meanHist = sumHist/(countHist.astype(float)+1e-10)
            sigmaHist = (sumHist**2 - sumSqrHist)/(countHist.astype(float)+1e-10)
            hdu1=histToFITS( meanHist, bins=bins,
                             histrange=histrange,name="MEAN" )
            hdu2=histToFITS( sigmaHist, bins=bins,
                             histrange=histrange, name="STDDEV")
            hdu3=histToFITS( countHist, bins=bins,
                             histrange=histrange, name="COUNTS" )
            hdulist = pyfits.HDUList()
            hdulist.append( pyfits.PrimaryHDU() )
            hdulist.append( hdu1 )
            hdulist.append( hdu2 )
            hdulist.append( hdu3 )
            hdulist.writeto( filename+".fits", clobber=True )
        else:
            hdu1=actutils.histToFITS( sumHist, bins=bins,
                                      histrange=histrange,name="MEAN" )
            hdu2=actutils.histToFITS( sumSqrHist, bins=bins,
                                      histrange=histrange, name="STDDEV")
            hdu3=actutils.histToFITS( countHist, bins=bins,
                                      histrange=histrange, name="COUNTS" )
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
                             debug=opts.debug,namebase=namebase )
    generateTelLookupTables( infile, varName="HIL_TEL_LENGTH" , 
                             debug=opts.debug, namebase=namebase )

    


 
