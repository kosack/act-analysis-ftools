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


def convolveImages(image1,image2,exclusionMap=None):
    """
    Convolves image1 with image2, which should be 2D arrays (the data
    part of a FITS file), producing a new image, taking into account
    an exclusion map if provided
    
    Arguments:
    - `image1`: NxM array input image
    - `image2`: NxM array convolution map
    - `exclusionMap`:  NxM array of boolean values for exclusion
    """
    
    if (image1.shape != image2.shape):
        raise Exception("Image shapes don't agree")

    if (exclusionMap == None):
        excl=np.zeros( image1.shape )
    else:
        excl = exclusionMap

    # note that the ring is centered at (nx/2, ny/2),
    
    # can probably do this in one big vector op with mgrid, but for
    # now let's try it with for-loops (slow but easier)

    nx,ny = image1.shape
    cx,cy = nx/2,ny/2
    conv = np.zeros( image1.shape )

    print "Convolving... ",image1.shape, cx,cy
    for ii in xrange(nx):
        for jj in xrange(ny):
            shifted = image2[ii-cx:ii-cx+nx,jj-cy:jj-cy+ny]
            conv[ii,jj] = np.sum(image1 * shifted )
            if (conv[ii,jj] > 0): print ii,jj,conv[ii,jj]

    return conv


if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-o","--output", dest="output", help="output filename")
    parser.add_option( "-e","--exclmap", dest="exclusion", 
                       help="exclusion map (FITS image)")
    parser.add_option( "-a","--accmap", dest="accmap", 
                       help="acceptance map (FITS image)")
    parser.add_option( "-m","--onradius", dest="onradius", help="ON region radius")
    parser.add_option( "-f","--areafactor", dest="areafact", help="Area factor")

    (opts, args) = parser.parse_args()
    
    if (opts.onradius):
        r_on = float(opts.onradius)
    else:
        r_on = 0.1
    

    imfile = args.pop(0)
    radii = calcRadiiFromArea( r_on, areaFactor=7)
    area_ring = ringArea( radii )
    area_on = math.pi*(r_on**2)

    print "  INPUT IMAGE: ",imfile
    print "EXCLUSION MAP: ", opts.exclusion
    print "       OUTPUT: ", opts.output
    
    print "        RADII: ", radii
    print "    RING AREA: ", area_ring, "deg^2 (before exclusions)"
    print "      ON AREA: ", area_on, "deg^2 (before exclusions)"
    print "  AREA FACTOR: ", area_ring/area_on

    imhdu = pyfits.open(imfile)[0]

    excl = None
    if (opts.exclusion):
        excl =  pyfits.open(opts.exclusion)[0]

    acc = pyfits.open( opts.accmap)[0]

    ring = makeRingMap(imhdu, radii)
    ring /= area_ring


#    conv = convolveImages( imhdu.data, ring, exclusionMap=excl )
    print "Convolving ON..."
    conv = signal.fftconvolve( imhdu.data, ring, mode='same' )
    print "Convolving ACC..."
    aconv = signal.fftconvolve( acc.data, ring, mode='same' )
    print "Convolving EXCL..."
    econv = signal.fftconvolve( excl.data, ring, mode='same' )
                              
    
    if (opts.output):
        pyfits.writeto( "ring_"+opts.output, header=imhdu.header,
                        data=ring,clobber=True )
        pyfits.writeto( "acc_"+opts.output, header=imhdu.header, 
                        data=aconv,clobber=True )
        pyfits.writeto( opts.output, header=imhdu.header, 
                        data=conv,clobber=True )
        pyfits.writeto( "ringexcl_"+opts.output, header=imhdu.header, 
                        data=econv,clobber=True )  
        pyfits.writeto( "ringbg_"+opts.output, header=imhdu.header, 
                        data=conv/econv,clobber=True )  



