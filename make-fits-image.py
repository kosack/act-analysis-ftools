#!/usr/bin/python

import pyfits
import numpy as np
import sys

import actutils
from optparse import OptionParser

from kapteyn import wcs

#
# TODO:
#    - add a makeCube() function, with logEnergy as the 3rd dimension
#





center=[83.633333,22.014444]
geom = [200,200]
FOV= [5.0,5.0]

parser = OptionParser()
parser.add_option( "-f","--fov", dest="fov", help="field of view",
                   metavar="X,Y")
parser.add_option( "-g","--geom", dest="geom", help="image geometry",
                   metavar="NX,NY")
parser.add_option( "-c","--center", dest="center", help="image center degrees",
                   metavar="RA,Dec")

parser.add_option( "-o","--output", dest="output", help="output filename")
parser.add_option( "-p","--projection", dest="proj", help="projection",
                   default="CAR")
parser.add_option( "-v","--verbose", dest="verbose", help="more output")
parser.add_option( "-r","--rmax", dest="rmax", help="maximum radius in deg")
parser.add_option( "-s","--system", dest="system", 
                   help="coordinate sys of output (e.g. 'fk5' or 'galactic')",
                   default="fk5");

parser.set_usage("makefits.py [options] eventsfile.fits")


(options, args) = parser.parse_args()

if (options.fov):
    FOV = np.array(options.fov.split(",")).astype(float)

if (options.geom):
    geom = np.array(options.geom.split(",")).astype(float)

if (options.center):
    center = np.array(options.center.split(",")).astype(float)

if (options.rmax):
    rmax = float(options.rmax)


if len(sys.argv)>1:
    inputfile = args.pop()
else:
    print "Please specify at least an input file. See --help";
    sys.exit(1)

if options.verbose:
    print "CENTER:",center
    print "  GEOM:",geom
    print "   FOV:",FOV

# generate blank output image:
hdu = actutils.makeImageHDU(center=center, geom=geom, 
                            FOV=FOV,projection=options.proj, 
                            system=options.system)

# get events
ff=pyfits.open(inputfile)
events = ff['EVENTS']
ra = events.data.field("RA").astype(float)
dec = events.data.field("DEC").astype(float)
ehdr = events.header

# update header with event information:
actutils.copyHeaders( ehdr, hdu.header )

# make count map:

newdata = actutils.makeCountMap( hdu, ra,dec, verbose=options.verbose)

# make FOV mask:
if (options.rmax):
    mask = actutils.makeRadialFOVMask( hdu, rmax )
    newdata *= mask

if (options.output):
    pyfits.writeto(options.output, header=hdu.header, 
                   data=newdata.transpose(), clobber=True )








