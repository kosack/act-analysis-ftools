#!/usr/bin/python2.4

import pyfits
import numpy as np

from kapteyn import wcs


# def dumpImagePixels(image, oversample=4):
#     wcs = astWCS.WCS( image.header, mode='pyfits' )
#     ia = arange(image.data.shape[0])
#     ja = arange(image.data.shape[0])

#     print wcs.pix2wcs( zip(ia,ja) )



if __name__ == "__main__":

    
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-v","--verbose", dest="verbose", help="more output")
    parser.add_option( "-s","--oversample", dest="oversample", 
                       help="Number of times to oversample each pixel", default=1)
    parser.set_usage("make-flat-eventlist.py [options] <input fits image> "
                     +"<output fits evlist>")

    (opts, args) = parser.parse_args()

    filename = args.pop(0)
    outfile = args.pop(0)

    print " INPUT: ",filename
    print "OUTPUT: ",outfile

    fits = pyfits.open( filename )
    image = fits[0]
    shape = np.array(image.data.shape)
    
    oversample = int(opts.oversample)

    print "     SHAPE: ",shape
    print "OVERSAMPLE: ",oversample

    nx = (shape[0])*oversample
    ny = (shape[1])*oversample

    proj = wcs.Projection( image.header )
    ipix = np.mgrid[0:nx,0:ny].astype(float)+0.5 # grid with half-bin shift
    ipix /= float(oversample) # put in correct range

    pix = np.array( proj.toworld( (ipix[0].flatten(), ipix[1].flatten() ) ) )
    detx = (ipix[0].flatten() - shape[0]/2.0) * proj.cdelt[0]
    dety = (ipix[1].flatten() - shape[1]/2.0) * proj.cdelt[1]

    vals = image.data.flatten()/(oversample**2)

    # now write out an eventlist with columns RA,DEC,DETX,DETY

    names = np.array( ['RA', 'DEC'] )
    c1 = pyfits.Column( name='RA', format='D', array = pix[0] )
    c2 = pyfits.Column( name='DEC', format='D', array = pix[1] )
    c3 = pyfits.Column( name='DETX', format='D', array = detx )
    c4 = pyfits.Column( name='DETY', format='D', array = dety )
    c5 = pyfits.Column( name="VALUE", format='D', array=vals )

    evlist = pyfits.new_table([c1,c2,c3,c4,c5])
    evlist.name = "EVENTS"

    print evlist.header.ascardlist()
    
    evlist.writeto(outfile, clobber=True)
    
    
