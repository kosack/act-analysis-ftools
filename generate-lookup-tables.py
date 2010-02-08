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
 

# input is simulated event list (should apply some basic cuts to avoid
# crazy values)



def showIt(x,y,hist):
    pcolor( x,y, hist.transpose() )
    xlabel("log10(SIZE)")
    ylabel("Impact distance")
    colorbar()


def gaussKern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def generateTelLookupTables(events,varName="HIL_TEL_WIDTH",
                            bins = [20,20], histrange =[[0,7],[0,1500]],
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

    # impacts distances need to be calculated for each telescope (the
    # impact distance stored is relative to the array center)
    allImpacts = zeros( allValues.shape )
    for ii in range(allValues.shape[1]):
        print "t",ii
        allImpacts[:,ii] = np.sqrt( (allCoreX-tposx[ii])**2 +
                                    (allCoreY-tposy[ii])**2 )


    nevents,ntels = allValues.shape

    # we want to do some basic cuts on width, removing ones with bad
    # values (why do bad value widths still exist for triggered
    # events?)
    valueMask = allValues > -100 
    telMask *= valueMask  # mask off bad values
    
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

        meanHist = sumHist/(countHist.astype(float)+1e-10)
        sigmaHist = (sumHist**2 - sumSqrHist)/(countHist.astype(float)+1e-10)
        
        if (debug):
            subplot(4,ntels,ntels*0 + itel + 1)
            showIt( edX, edY, meanHist )
            title ("CT%d %s" % (telid[itel],varName))

            subplot(4,ntels,ntels*1 + itel + 1)
            showIt( edX, edY, sigmaHist )
            title ("CT %d Sigma"  % telid[itel])

            subplot(4,ntels,ntels*2 + itel + 1)
            showIt( edX, edY, meanHist )
            title ("CT%d %s" % (telid[itel],varName))



        # write it out as a FITS file with 2 image HDUs VALUE and SIGMA

        if namebase:
            filename = "%s-CT%03d-%s-lookup" % (namebase,telid[itel], varName)
        else:
            filename = "CT%03d-%s-lookup" % (telid[itel], varName)

        print "CT",telid[itel],varName,"Nevents=",len(value),
        print "outliers=",len(value)-sum(countHist),
        print "out:",filename

        
        if (singleFile) :
        
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
            hdu1.writeto( filename+"-sum.fits" )
            hdu2.writeto( filename+"-sum2.fits" )
            hdu3.writeto( filename+"-count.fits" )


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

    


 
