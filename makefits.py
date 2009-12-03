#!/usr/bin/python2.4

import pyfits
from numpy import *

from pylab import *
from astLib import astPlots
from astLib import astWCS

def makeFITS(centerRADec,geom=(300,300), FOV=(1.0,1.0), projection="CAR",
             output=None):
    """
    Generate a blank FITS image. Returns a HDU with the given geometry
    (so far in RA/Dec only). 

    Arguments:
    - `centerRADec`: (lambda,beta) center of image
    - `geom`: (N,M) size in pixels of the image
    - `FOV_deg`: (dx,dy) size of FOV in degrees
    - `output`: output filename, or None for no output
    """
  
    centerRADec=array(centerRADec)
    geom = array(geom)
    FOV = array(FOV)
    delta = FOV/geom
    
    hdu = pyfits.NP_pyfits.PrimaryHDU()
    hdu.data = zeros( geom )
    
    hdu.update_header()
    hdu.header.add_comment("Generated by makefits.py")
    hdu.header.update( "CTYPE1", "RA---"+projection )
    hdu.header.update( "CTYPE2", "DEC--"+projection )
    hdu.header.update( "RADESYS", "fk5" )
    hdu.header.update( "EQUINOX", 2000 )
    
    hdu.header.update( "CRVAL1", centerRADec[0] )
    hdu.header.update( "CRVAL2", centerRADec[1] )
    
    hdu.header.update( "CRPIX1", geom[0]/2 )
    hdu.header.update( "CRPIX2", geom[1]/2 )
    
    hdu.header.update( "CDELT1", delta[0] )
    hdu.header.update( "CDELT2", delta[1] )
    
    hdu.header.update( "BSCALE", 1.0 )
    hdu.header.update( "BZERO", 0.0 )
    
    print hdu.header.items();
    
    if output:
        pyfits.writeto(output, header=hdu.header, data=hdu.data )
        
    return hdu

def makeCountMap(hdu, ra,dec, output=None):
    """
    Arguments:
    - `hdu`: hdu containing data array and header with proper WCS info 
    - `ra`: array of RA coordinates
    - `dec`: array of Dec coordinates

    todo: add events coordsystem as an input (defaults to RA/Dec)
    """

    wcs = astWCS.WCS( hdu.header, mode='pyfits' )
    coords = array(wcs.wcs2pix( ra,dec ))

    bins = array(hdu.data.shape)
    pixsize = array( (wcs.getXPixelSizeDeg(),wcs.getYPixelSizeDeg()))
    
    
    imrange = zip( (0,0), bins) # TODO: probably need some half bins here to make it
                                                # right?  check
                                                # histogram code
    print "range: ",imrange
    
    
    H, xedges, yedges = histogram2d( coords[:,0],coords[:,1],
                                     bins=bins, range=imrange )
    
    print "   H: ", H.shape
    print "data: ", hdu.data.shape

    if output:
        print "Writing count map:",output
        pyfits.writeto(output, header=hdu.header, data=H, clobber=True )

    return H


def makeAcceptanceFromEvents(events,image,bins=20,range=[0,7]):
    """ 
    Generates an acceptance map from the events of a single run
    (doesn't use lookups, just a simple accepance-from-data
    
    Arguments:
    - `eventshdu`: hdu of cut and masked eventlist (uses header and event info)
    - `imagehdu`: hdu of output image
    """
    
    # build 1D histogram in detector coordinates (pointing dir is 0,0,
    # and distances are in degrees from center)

    runhdr = events.header
    obspos = array([runhdr.get("RA_PNT"), runhdr.get("DEC_PNT")])
    
    X = events.data.field("DETX")
    Y = events.data.field("DETY")
    D = sqrt(X**2 + Y**2)

    subplot(211)
    h,ed = histogram( D*D,bins=bins, range=range)
    scatter( ed,h )
    subplot(212)
    scatter(X,Y,s=1)
    show()

    return D


if __name__ == "__main__":

    from optparse import OptionParser

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

    parser.add_option( "-d","--display", action="store_true",
                       dest="display", help="show the image")

    parser.add_option( "-o","--output", dest="output", help="output filename")
    parser.add_option( "-p","--projection", dest="proj", help="projection",
                       default="CAR")
    parser.set_usage("makefits.py [options] eventsfile.fits")


    (options, args) = parser.parse_args()

    if (options.fov):
        FOV = array(options.fov.split(",")).astype(float)

    if (options.geom):
        geom = array(options.geom.split(",")).astype(float)

    if (options.center):
        center = array(options.center.split(",")).astype(float)

    print "ARGV:",args

    if len(sys.argv)>1:
        inputfile = args.pop()
    else:
        print "Please specify at least an input file. See --help";
        sys.exit(1)

    print "CENTER:",center
    print "  GEOM:",geom
    print "   FOV:",FOV

    # generate blank output image:
    hdu = makeFITS( centerRADec=center, geom=geom, FOV=FOV,projection=options.proj)
    wcs = astWCS.WCS( hdu.header, mode='pyfits' )

    # get events
    ff=pyfits.open(inputfile)
    events = ff['EVENTS']
    ra = events.data.field("RA").astype(float)
    dec = events.data.field("DEC").astype(float)

    # make count map:

    newdata = makeCountMap( hdu, ra,dec, output=options.output )

    #display it
    if options.display:
        img = astPlots.ImagePlot( newdata, wcs, 
                                  cutLevels=["relative", 99.5], 
                                  colorMapName="jet"   )
        img.draw()
        title("Crab Nebula")
        show()
