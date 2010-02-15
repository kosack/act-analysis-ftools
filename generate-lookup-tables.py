import pyfits
import numpy as np
from optparse import OptionParser
import sys
import math
from scipy import interpolate
import scipy.signal
import re
from pylab import *

import actutils
 
# TODO: apply local-distance cut! (maybe make a selection tool for
# telescope cuts: telselect. Want localdist<0.525 (but better as a
# fraciton or something to scale with bigger telescopes)


def generateTelLookupTables(events,varName="HIL_TEL_WIDTH",
                            bins = [60,60], histrange =[[0,7],[0,1500]],
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

    tposx = telarray.data.field("POSX")
    tposy = telarray.data.field("POSY")
    print "TPOSX",tposx
    print "TPOSY",tposy
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
    localDistance = sqrt( cogx**2 +cogy**2)
    localDistMask = localDistance < 0.025


    # impacts distances need to be calculated for each telescope (the
    # impact distance stored is relative to the array center)
    allImpacts = zeros( allValues.shape )
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
    telMask *= localDistMask  # mask off bad values
    
    if (debug):
        figure( figsize=(15,10))
        subplot(1,1,1)

    # now, generate the lookup table for each telescope separately:
    for itel in range(ntels):

        goodEvents = telMask[:,itel]  
        value = allValues[:,itel][goodEvents]
        impact = allImpacts[:,itel][goodEvents]
        logsiz = np.log10(allSizes[:,itel][goodEvents])

        sumHist,edX,edY = np.histogram2d( logsiz,impact, 
                                           weights=value,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        sumSqrHist,edX,edY = np.histogram2d( logsiz,impact, 
                                           weights=value**2,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        countHist,edX,edY = np.histogram2d( logsiz,impact, 
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

    


 
