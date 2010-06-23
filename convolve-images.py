import pyfits
import numpy as np
import math

from numpy.fft import fft2, ifft2
from scipy import signal

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-o","--output", dest="output", help="output filename")
    parser.add_option( "-f","--nofft", dest="nofft", action="store_true",
                       default=False, help="don't use FFT to convolve (slower)")
    parser.add_option( "-0","--round-to-zero", dest="zero", action="store_true",
                       default=False, help="round small values to 0")    

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


    imhdu1.data[np.isnan(imhdu1.data)] = 0
    imhdu2.data[np.isnan(imhdu2.data)] = 0

    print "CONVOLVING:",imfile1,"with",imfile2
    conv = convolve( imhdu1.data, imhdu2.data, mode='same' )

    if (opts.zero==True):
        conv[abs(conv)<1e-10] = 0.0
                             
    if (opts.output):
        pyfits.writeto( opts.output, header=imhdu1.header,
                        data=conv,clobber=True )



