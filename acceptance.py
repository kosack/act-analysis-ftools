import pyfits
import numpy as np
from pylab import *
from kapteyn import wcs

import actutils


# TODO: do the exclusion correction correctly! Not sure the flatlist is the best idea

    


def makeAcceptanceMapFromEvents(events,imagehdu, obspos, rmax=3.0,debug=False):
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
     
    profile,edges = actutils.makeRadialProfile( events, range=[0,6.0**2])

    proj = wcs.Projection( imagehdu.header )
    binarea = abs(proj.cdelt[0] * proj.cdelt[1])

    print "2D BIN AREA:",binarea, "deg^2"

    dists2 = actutils.makeDistanceMap( imagehdu, obspos )**2
    
    # interpolate the squared distances to each bin from the radial
    # profile. 

    acc = np.interp( dists2,  edges[:edges.shape[0]-1], profile ) 

    # [counts/area] * area = [counts]"
    acc *= binarea

    # apply cutoff if requested:
    acc[dists2>rmax**2] = 0.0

    # want integrals to match: sum(acc) = sum(profile)
    print "INTEGRALS: profile=",sum(profile), " acc=", sum(acc)
    
    print "debug:",debug
    if debug==True:
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
    parser.add_option( "-d","--debug", dest="debug", action="store_true",
                       default=False, help="display debug plots interactively")
    
    (opts, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("incorrect number of arguments")

    evfile = args.pop(0)
    imfile = args.pop(0)
    
    print "EVENTS: ",evfile
    print "IMAGE : ",imfile
    print "OUTPUT: ", opts.output
    print "DEBUG: ",opts.debug

    evhdu = pyfits.open(evfile)['EVENTS']
    imhdu = pyfits.open(imfile)[0]

    try:
        runhdr = evhdu.header
        obspos = array([runhdr["RA_PNT"], runhdr["DEC_PNT"]])
    except KeyError:
        runhdr = imhdu.header
        obspos = array([runhdr["RA_PNT"], runhdr["DEC_PNT"]])        

    rmax = float(opts.rmax)



    A = makeAcceptanceMapFromEvents( evhdu, imhdu, obspos, rmax=rmax, 
                                     debug=int(opts.debug) )

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
        

