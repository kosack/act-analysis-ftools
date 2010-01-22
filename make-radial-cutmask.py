import actutils
import pyfits
from optparse import OptionParser
from kapteyn import wcs
import numpy as np

# this program takes an input image and creates an output image of the
# same geometry that has a value of 1.0 inside the given radius and
# 0.0 outside. The radius is calculated from the pointing direction f
# the run (from the RA_PNT and DEC_PNT header keys) or can be
# specified in RA/Dec coordinates via the --center option. If no
# center is given and no pointing direction is found in the map, the
# map center is used as the reference position.

parser = OptionParser(usage="%prog [options] R_max input-image output-image")
parser.add_option( "-c","--center", dest="center", 
                   help="image center in degrees",
                   metavar="RA,DEC")

(opts, args) = parser.parse_args()

if len(args) != 3:
    parser.error("incorrect number of arguments")

rmax = float(args.pop(0))
inimg = args.pop(0)
outimg = args.pop(0)

hdu = pyfits.open(inimg)[0]

if (opts.center):
    center = np.array(options.center.split(",")).astype(float)
else:
    center=None

mask = actutilsmakeRadialFOVMask( hdu, rmax, center=center )

pyfits.writeto( outimg, header=hdu.header, data=mask, clobber=True )
    
