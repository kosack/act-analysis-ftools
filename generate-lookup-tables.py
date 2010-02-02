import pyfits
import numpy as np
from optparse import OptionParser
import sys
import math
import scipy.interpolate
import scipy.signal

from pylab import *


Debug=False

# input is simulated event list (should apply some basic cuts to avoid
# crazy values)



def showit(x,y,hist):
    figure()
    pcolor( x,y, hist.transpose() )


def smoothLookup(values,counts, extend=False):
    """
    Arguments:
    - `value`:
    - `counts`:
    
    """

    print "Smoothing"

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
    print "Done"

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


    


if __name__ == '__main__':
    
    parser = OptionParser()

    (opts, args) = parser.parse_args()

    infile = args.pop(0)

    evfile = pyfits.open(infile)
    events = evfile["EVENTS"]

    telid = np.array(events.header['TELLIST'].split(",")).astype(int)
    trig   = events.data.field("TELMASK") 
    awid = events.data.field("HIL_TEL_WIDTH") 
    asiz = events.data.field("HIL_TEL_SIZE") 
    acorex = events.data.field("COREX") 
    acorey = events.data.field("COREY") 
    aimpact = np.sqrt( acorex**2 + acorey**2 )

    mccorex = events.data.field("MC_COREX") 
    mccorey = events.data.field("MC_COREY") 

    nevents,ntels = awid.shape

    widmask = awid > -100 # remove bad widths (why are they there?)
    mask = widmask * trig

    bins = [80,80]
    histrange =[[0,7],[0,1500]]
    

    if (Debug):
        figure( figsize=(15,10))
        subplot(1,1,1)

    telwidmean = dict()
    telwidsigma = dict()

    for itel in range(ntels):

        trigmask = mask[:,itel]  # cuts out non-triggered telescopes
        wid = awid[:,itel][trigmask]
        impact = aimpact[trigmask]
        logsiz = np.log10(asiz[:,itel][trigmask])

        widsum,edX,edY = np.histogram2d( logsiz,impact, 
                                           weights=wid,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        wid2sum,edX,edY = np.histogram2d( logsiz,impact, 
                                           weights=wid**2,
                                           range=histrange, 
                                           bins=bins,
                                           normed=False)

        count,edX,edY = np.histogram2d( logsiz,impact, 
                                        weights=None,
                                        range=histrange, 
                                        bins=bins,
                                        normed=False)

        meanwid = widsum/(count.astype(float)+1e-10)
        sigmawid = (widsum**2 - wid2sum)/(count.astype(float)+1e-10)

        telwidmean[telid[itel]] = meanwid
        telwidsigma[telid[itel]] = sigmawid

        #meanwid = fixLookup( value=meanwid, count=count, range=histrange)
        #sigmawid = fixLookup( value=sigmawid, count=count, range=histrange)

        if (Debug):
            subplot(2,ntels,itel+1)
            pcolor( edX, edY, meanwid.transpose() )
            colorbar()
            title ("CT%d Width" % telid[itel])
            xlabel("log10(SIZE)")
            ylabel("Impact distance")

            subplot(2,ntels,ntels*1+itel+1)
            pcolor( edX, edY, sigmawid.transpose() )
            colorbar()
            title ("CT %d Sigma"  % telid[itel])
            xlabel("log10(SIZE)")
            ylabel("Impact distance")

        print "CT",telid[itel]," width:",len(wid),sum(widsum), sum(count)

 
    ii = np.arange(edX.shape[0])
    deltaSize = ((edX[ii]-edX[ii-1])/2.0)[1:]
    csize = edX[0:-1]+deltaSize # center of size bin
    deltaDist = ((edY[ii]-edY[ii-1])/2.0)[1:]
    cdist = edY[0:-1]+deltaDist # center of impact distance bin

    msize,mdist = np.meshgrid( csize,cdist )


    smooth = smoothLookup(telwidsigma[telid[0]], count)

    valueinterp = scipy.interpolate.interp2d(x=msize.flatten(),y=mdist.flatten(),  
                                             z=smooth,
                                             kind='linear')
    

    interpSizes = arange(0,7,0.1)
    interpDists = arange(0,1470,0.1)
    interpVal = valueinterp(interpSizes,interpDists)

    figure( figsize=(15,5))
    subplot(1,3,1)
    pcolor( edX, edY,telwidsigma[telid[0]].transpose() )
    colorbar()
    subplot(1,3,2)
    pcolor( edX, edY,smooth.transpose() )
    colorbar()
    subplot(1,3,3)
    pcolor( interpSizes,interpDists, interpVal.transpose() )
    colorbar()
