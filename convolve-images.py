import pyfits
import numpy as np
import math
from astLib import astWCS
from astLib import astCoords

from numpy.fft import fft2, ifft2
from scipy import signal

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-o","--output", dest="output", help="output filename")
    parser.add_option( "-f","--nofft", dest="nofft", 
                       help="don't use FFT to convolve (slower)")

    (opts, args) = parser.parse_args()
    
    if len(args) < 2:
        print "convolve-images.py [--output <output>] <image1> <image2>"
        exit(1)

    imfile1 = args.pop(0)
    imfile2 = args.pop(0)

    imhdu1 = pyfits.open(imfile1)[0]
    imhdu2 = pyfits.open(imfile2)[0]

    convolve = signal.fftconvolve
    if (opts.nofft):
        convolve = signal.convolve2d

    print "CONVOLVING:",imfile1,"with",imfile2
    conv = convolve( imhdu1.data, imhdu2.data, mode='same' )
                             
    if (opts.output):
        pyfits.writeto( opts.output, header=imhdu1.header,
                        data=conv,clobber=True )


