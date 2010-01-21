import pyfits
import numpy as np

from pylab import *

from kapteyn import wcs

def makeRadialProfile(events,bins=14,range=[0,10]):
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

def angSepDeg(lambda0, beta0, lambda1, beta1):
    """ calculate angular separation between two ra/dec coordinates
    using the Vincenti formula (accurate for large and small
    separations).  lambda = longitudinal coordinates (e.g. RA or gall)
    and beta=latitudinal coordinate (e.g. dec or galb). These must be
    in degrees already.
    
    EXAMPLE: 
        ra = arange(100)*0.1
        dec = ones(100)*-20.0
        center_radec = (120.0,15.0) # distance to a single point (could be array)
        dists = angSepDeg( ra, dec, center_radec[0], center_radec[1] )

    Arguments:
    - `lambda0,beta0`: first ra/dec coordinate in deg (or arrays)
    - `lambda1, beta1`: second ra/dec coordinates in deg (or arrays)
    """
    
    print "DEBUG: ",lambda0, lambda1

    l0 = np.radians( np.array(lambda0) )
    l1 = np.radians( np.array(lambda1) )
    b0 = np.radians( np.array(beta0) )
    b1 = np.radians( np.array(beta1) )

    cosb0 = np.cos(b0)
    cosb1 = np.cos(b1)
    sinb0 = np.sin(b0)
    sinb1 = np.sin(b1)
    cosdl = np.cos(l1-l0)
    
    numer = np.sqrt( (cosb1*np.sin(l1-l0))**2 +
                     (cosb0*sinb1 - sinb0*cosb1*cosdl)**2 )

    denom = sinb0*sinb1 + cosb0*cosb1*cosdl
    
    return np.degrees(np.arctan2( numer, denom ))
    

def makeDistanceMap(imagehdu,pos):
    """ 
    Generate a map of angular distances to a given observation
    position.
    
    Arguments:
    - `imagehdu`: input image HDU (data + header)
    - `pos`:

    returns: data array of distances from each bin center to the position
    """
    
    proj = wcs.Projection( imagehdu.header )
    binarea = proj.cdelt[0] * proj.cdelt[1]
    (nx,ny) = imagehdu.data.shape # CHECK THAT THIS IS NOT BACKWARD!
    
    ia,ja = np.mgrid[0:nx,0:ny]+0.5
    (ra,dec) = proj.toworld( (ia.flatten(), ja.flatten()) )

    distmap = angSepDeg( ra, dec,  pos[0], pos[1] )
    distmap.shape = ia.shape

    return distmap
    


def makeAcceptanceMapFromEvents(events,imagehdu, obspos, rmax=3.0,debug=True):
    """
    Returns an acceptance map
    
    Arguments:
    - `events`: hdu containing *excluded* eventlist
    - `imagehdu`: hdu for the source image 
    - `obspos`: observation position in RA/Dec (position about which
      to make the acceptance)
    - `rmax`: maximum radius (acceptance is set to 0 outside this)
    """
    runhdr = events.header
     
    profile,edges = makeRadialProfile( events, range=[0,6.0**2])

    proj = wcs.Projection( imagehdu.header )
    binarea = abs(proj.cdelt[0] * proj.cdelt[1])

    print "2D BIN AREA:",binarea, "deg^2"

    dists2 = makeDistanceMap( imagehdu, obspos )**2
    
    # interpolate the squared distances to each bin from the radial
    # profile. 

    acc = np.interp( dists2,  edges[:edges.shape[0]-1], profile ) 

    # [counts/area] * area = [counts]"
    acc *= binarea

    # apply cutoff if requested:
    acc[dists2>rmax**2] = 0.0

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
    parser.add_option( "-r","--rmax", dest="rmax", default="5.0",
                       help="maximum radius")
    parser.add_option( "-e","--exflat", dest="exflat", 
                       help="excluded flat eventlist (for exclusion region correction)")
    
    (opts, args) = parser.parse_args()

    evfile = args.pop(0)
    imfile = args.pop(0)
    
    print "EVENTS: ",evfile
    print "IMAGE : ",imfile
    print "OUTPUT: ", opts.output
    

    evhdu = pyfits.open(evfile)['EVENTS']
    imhdu = pyfits.open(imfile)[0]

    try:
        runhdr = evhdu.header
        obspos = array([runhdr["RA_PNT"], runhdr["DEC_PNT"]])
    except KeyError:
        runhdr = imhdu.header
        obspos = array([runhdr["RA_PNT"], runhdr["DEC_PNT"]])        

    rmax = float(opts.rmax)

    A = makeAcceptanceMapFromEvents( evhdu, imhdu, obspos, rmax=rmax )

    if (opts.exflat):
        # apply  correction for  exclusion regions  from  the excluded
        # flat map (which  samples the image in a  "flat" way, but not
        # covering exclusion regions. 
        exfile = pyfits.open( opts.exflat)['EVENTS']
        E = makeAcceptanceMapFromEvents( exfile, imhdu, obspos, rmax=rmax )
        # need to scale it?
        A /= E+1e-30

    if (opts.output):
        pyfits.writeto( opts.output, header=imhdu.header, data=A,clobber=True )
        

