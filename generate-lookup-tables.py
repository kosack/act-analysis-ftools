import pyfits
import numpy as np
from optparse import OptionParser
import sys
import math
from scipy import interpolate
import scipy.signal
import re
from pylab import *



# input is simulated event list (should apply some basic cuts to avoid
# crazy values)



def showIt(x,y,hist):
    pcolor( x,y, hist.transpose() )
    xlabel("log10(SIZE)")
    ylabel("Impact distance")
    colorbar()


def smoothLookup(values,counts, extend=False):
    """
    Arguments:
    - `value`:
    - `counts`:
    
    """
    kern = gaussKern( 5,5 )
    values = values.copy()
    cx = kern.shape[0]/2.0
    cy = kern.shape[1]/2.0
    kern /= kern[ cx,cy ] # correct normalization 

    smoothed = scipy.signal.convolve( values, kern, mode='same' ).transpose()
    ecounts = counts.transpose()

    if extend:
        # extend to the right
        for ii in arange(1,smoothed.shape[0]):
            badpix = ecounts[ii]<10
            smoothed[ii][badpix] = smoothed[ii-1][badpix] 

    return smoothed.transpose()

    
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
                            bins = [80,80], histrange =[[0,7],[0,1500]],
                            debug=False, namebase=None):
    """
    generates lookup table for the given variable (average and sigma
    as a function of logSIZE and IMPACT DISTANCE)

    Arguments:
    - `events`: event-list FITS HDU 
    - `varName`: variable to histogram (HIL_TEL_WIDTH or HIL_TEL_LENGTH)
    - `bins`: number of bins for logSIZE,DISTANCE
    - `histrange`: logSIZE and DISTANCE ranges
    """

    telid = np.array(events.header['TELLIST'].split(",")).astype(int)
    
    # these are the data fields as NxM arrays (where N is number of events,
    # and M is number of telescopes):
    telMask   = events.data.field("TELMASK") 
    allValues = events.data.field(varName) 
    allSizes = events.data.field("HIL_TEL_SIZE") 
    allCoreX = events.data.field("COREX") 
    allCoreY = events.data.field("COREY") 
    allImpacts = np.sqrt( allCoreX**2 + allCoreY**2 )
    mccorex = events.data.field("MC_COREX") 
    mccorey = events.data.field("MC_COREY") 

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
        impact = allImpacts[goodEvents]
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
        
        meanHistSmooth = smoothLookup(meanHist, countHist)
        sigmaHistSmooth = smoothLookup(sigmaHist, countHist)

        if (debug):
            subplot(4,ntels,ntels*0 + itel + 1)
            showIt( edX, edY, meanHist )
            title ("CT%d %s" % (telid[itel],varName))

            subplot(4,ntels,ntels*1 + itel + 1)
            showIt( edX, edY, sigmaHist )
            title ("CT %d Sigma"  % telid[itel])

            subplot(4,ntels,ntels*2 + itel + 1)
            showIt( edX, edY, meanHistSmooth )
            title ("CT%d %s" % (telid[itel],varName))

            subplot(4,ntels,ntels*3 + itel + 1)
            showIt( edX, edY, sigmaHistSmooth )
            title ("CT %d Sigma (smooth)"  % telid[itel])

        print "CT",telid[itel],varName,"Nevents=",len(value),
        print "outliers=",len(value)-sum(countHist)

        # write it out as a FITS file with 2 images VALUE and SIGMA

        if namebase:
            filename = "%s-CT%03d-%s-lookup.fits" % (namebase,telid[itel], varName)
        else:
            filename = "CT%03d-%s-lookup.fits" % (telid[itel], varName)
        print "    --> ",filename
        
        



if __name__ == '__main__':
    
    parser = OptionParser()
    parser.set_usage( "generate-lookup-table.py [options] EVENTFILE")
    parser.add_option( "-d","--debug", dest="debug", action="store_true",
                       default=False, help="display debug plots interactively")    
    (opts, args) = parser.parse_args()

    if len(args)>0:
        infile = args.pop(0)
    else:
        parser.error("Please specify a filename")

    # lookup table definition:
    evfile = pyfits.open(infile)
    events = evfile["EVENTS"]

    match = re.search( r"_([\d]+)_", infile )
    if (match):
        runid = match.group(1)
        namebase = "run%s" % runid
    else:
        namebase=None
    


    generateTelLookupTables( events, varName="HIL_TEL_WIDTH", 
                             debug=False,namebase=namebase )
    generateTelLookupTables( events, varName="HIL_TEL_LENGTH" , 
                             debug=False, namebase=namebase )

    


 
