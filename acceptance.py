import pyfits
import numpy as np

from pylab import *
#from astLib import astPlots
from astLib import astWCS
from astLib import astCoords

def makeRadialProfile(events,bins=10,range=[0,7]):
    """ 
    Generates an radial profile from the events given in Detector
    coordinates (which are assumed to have (0,0) as the origin)
    
    Arguments:
    - `events`: hdu of cut and masked eventlist (uses header and event info)
    """
    
    # build 1D histogram in detector coordinates (pointing dir is 0,0,
    # and distances are in degrees from center) 

    # Note: we could also use the obspos + RA/DEC coords and calculate
    # the angular distance using WCS information, but this is easier
    # if DETX, and DETY exist

    runhdr = events.header
    obspos = array([runhdr.get("RA_PNT"), runhdr.get("DEC_PNT")])
    
    X = events.data.field("DETX")
    Y = events.data.field("DETY")
    D = sqrt(X**2 + Y**2)

    th2hist,ed = histogram( D*D,bins=bins, range=range,normed=False, new=True)
    print "  HIST: ",th2hist
    th2hist = th2hist.astype(np.float64)

    nevents = sum(th2hist)
    areaperbin = math.pi*(ed[1]) # ed[1]=r^2 (each bin has equal 2D area)
    th2hist /= areaperbin      # now in units of counts/deg^2 (density)

    # TODO: really should correct by exclusion area

    print "edges in t^2: ", ed
    print "edges in   t: ", sqrt(ed)
    print "total events in list: ",len(X)
    print "total events in hist: ",nevents
    print " area per radial bin: ", areaperbin,"deg^2"
    print "SC HIST: ",th2hist

#    scatter( ed, h )
#    show()


    return th2hist,ed

def makeAcceptanceMapFromEvents(events,imagehdu, debug=False):
    """
    Returns an acceptance map
    
    Arguments:
    - `events`: hdu containing *excluded* eventlist
    - `imagehdu`: hdu for the source image 
    """
    runhdr = events.header
    obspos = array([runhdr.get("RA_PNT"), runhdr.get("DEC_PNT")])
    
    profile,edges = makeRadialProfile( events, range=[0,4.0**2])

    # TODO: need to correct the profile for exclusion regions! Use
    # make-flat-eventlist.py and pass events through the region
    # filter? For now it's just ignored...

    wcs = astWCS.WCS( imagehdu.header, mode='pyfits' )
    binarea = wcs.getXPixelSizeDeg() * wcs.getYPixelSizeDeg()
    print "2D BIN AREA:",binarea, "deg^2"

    outhdu = imagehdu.copy()
    acc = outhdu.data
    nx  = acc.shape[0]
    ny  = acc.shape[1]

    # hack to get in list format for wcs routines
    ia,ja = np.mgrid[0:nx,0:ny]+0.5
    z=(ia+ja*1j).flatten() 

    # transform pixel grid to RA/Dec coordinates and calculate angular
    # separation
    radec=array(wcs.pix2wcs( z.real,z.imag ))

    # have to do a loop here since calcAngSepDeg only takes a
    # scalar.. could use map() maybe, but it's reasonably fast.  Note
    # that calcAngSepDeg only works for RA/Dec coordinates and assumes
    # a tangent projection (see astLib documentation). SHould be done
    # properly in 3D!

    dists2 = zeros(radec.shape[0])
    for ii in range(radec.shape[0]):
        dists2[ii] = astCoords.calcAngSepDeg( radec[ii][0], radec[ii][1],
                                              obspos[0], obspos[1] )**2
    
    print dists2.shape,ia.shape
    dists2.shape = ia.shape

    # interpolate the squared distances to each bin from the radial
    # profile. 

    acc = np.interp( dists2,  edges[:edges.shape[0]-1], profile ) 

    # [counts/area] * area = [counts]"
    acc *= binarea

    # want integrals to match: sum(acc) = sum(profile)
    print "INTEGRALS: profile=",sum(profile), " acc=", sum(acc)
    

    if (debug):
        testx=dists2[dists2.shape[0]/2] #arange(0,5.5,0.1)**2
        testa = np.interp( testx, edges[:edges.shape[0]-1], profile )
        reala = acc[dists2.shape[0]/2]/binarea
        #print zip(testx,testa)
        scatter( edges,profile,label='Radial Acc' )
        scatter( testx, testa, color='r', s=1, label='Interpolated' )
        scatter( testx, reala, color='g', s=1, label='2D Interpolated' )
        xlabel("$\Theta^2$")
        ylabel("$Counts/deg^2$")
        legend()
        show()


    return acc


    
if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-o","--output", dest="output", help="output filename")
    (opts, args) = parser.parse_args()

    evfile = args.pop(0)
    imfile = args.pop(0)
    
    print "EVENTS: ",evfile
    print "IMAGE : ",imfile
    print "OUTPUT: ", opts.output
    

    evhdu = pyfits.open(evfile)['EVENTS']
    imhdu = pyfits.open(imfile)[0]

    A = makeAcceptanceMapFromEvents( evhdu, imhdu )

    if (opts.output):
        pyfits.writeto( opts.output, header=imhdu.header, data=A,clobber=True )
        

