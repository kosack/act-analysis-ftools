import pyfits
import numpy as np
import math
from astLib import astWCS
from astLib import astCoords

from numpy.fft import fft2, ifft2
from scipy import signal

def calcRadiiFromArea(onRadiusDeg, minGapDeg=0.1, areaFactor=7.0):
    """ 
    Returns good inner and outter radius (r1,r2) given on-region
    radius
    
    Arguments:
    - `onRadiusDeg`: radius of ON region
    - `minGapDeg`: minimum gap between onRadius and r1
    - `areaFactor`: how many times the area of ON should the OFF be?
    """
    
    r1 = onRadiusDeg + minGapDeg
    r2 = math.sqrt(areaFactor*onRadiusDeg**2 + r1**2)

    return (r1,r2)


def ringArea( radii ):
    return math.pi*(radii[1]**2-radii[0]**2)
    

def makeRingMap(hdu, radii):
    """ Makes a ring map that can be convolved with another image
    
    Arguments:
    - `hdu`: input image hdu (WCS, etc info is used to construct the ring)
    - `radii`: radii in degrees (r1,r2) of the ring. 
    """

    nx,ny = hdu.data.shape 
    wcs = astWCS.WCS( hdu.header, mode='pyfits' )
    center = wcs.getCentreWCSCoords()
    ia,ja = np.mgrid[0:nx,0:ny].astype(np.double)+0.5 # grid with half-bin shift
    z=(ia+ja*1j).flatten() # wcs needs a list not a grid for whatever reason
    radec = np.array(wcs.pix2wcs( z.real, z.imag )) 

    dists = np.zeros(radec.shape[0])
    for ii in range(radec.shape[0]):
        dists[ii] = astCoords.calcAngSepDeg( radec[ii][0], radec[ii][1],
                                             center[0], center[1] )
    dists.shape = ia.shape # back to grid
 
    ring = np.ones(dists.shape)
    ring[dists<radii[0]] = 0.0
    ring[dists>radii[1]] = 0.0

    print np.sum(ring)
    return ring

    pass


if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-o","--output", dest="output", help="output filename")
    parser.add_option( "-m","--onradius", dest="onradius", help="ON region radius")
    parser.add_option( "-f","--areafactor", dest="areafact", 
                       help="Area factor", default=7.0)
    parser.add_option( "-g","--gap", dest="gap", default=0.1,
                       help="Gap between r_ON and r1 in deg")
    parser.add_option( "-O","--make-on", dest="makeon", action="store_true",
                       default=False,
                       help="instead of a ring, make an ON region (tophat)")

    (opts, args) = parser.parse_args()
    
    if (opts.onradius):
        r_on = float(opts.onradius)
    else:
        r_on = 0.1
    
    

    imfile = args.pop(0)
    radii = calcRadiiFromArea( r_on, areaFactor=float(opts.areafact), 
                               minGapDeg=float(opts.gap))

    if (opts.makeon):
        radii = (0,r_on)

    area_ring = ringArea( radii )
    area_on = math.pi*(r_on**2)

    print "  INPUT IMAGE: ",imfile
    print "       OUTPUT: ", opts.output
    
    print "        RADII: ", radii
    print "    RING AREA: ", area_ring, "deg^2 (before exclusions)"
    print "      ON AREA: ", area_on, "deg^2 (before exclusions)"
    print "  AREA FACTOR: ", area_ring/area_on

    imhdu = pyfits.open(imfile)[0]

    ring = makeRingMap(imhdu, radii)

    #ring /= area_ring   

    if (opts.output):
        pyfits.writeto( opts.output, header=imhdu.header,
                        data=ring,clobber=True )



