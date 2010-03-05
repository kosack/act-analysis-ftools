import pyfits
import numpy as np
import math

from kapteyn import wcs,maputils

# various utility functions 

# A note from the MatPlotLib manual for the pcolor() function:
#-------------
# Grid Orientation
#     The orientation follows the Matlab(TM) convention: an
#     array C with shape (nrows, ncolumns) is plotted with
#     the column number as X and the row number as Y, increasing
#     up; hence it is plotted the way the array would be printed,
#     except that the Y axis is reversed.  That is, C is taken
#     as C(y,x).
#     Similarly for meshgrid:
#         x = arange(5)
#         y = arange(3)
#         X, Y = meshgrid(x,y)
#     is equivalent to
#         X = array([[0, 1, 2, 3, 4],
#                   [0, 1, 2, 3, 4],
#                   [0, 1, 2, 3, 4]])
#         Y = array([[0, 0, 0, 0, 0],
#                   [1, 1, 1, 1, 1],
#                   [2, 2, 2, 2, 2]])
#     so if you have
#         C = rand( len(x), len(y))
#     then you need
#         pcolor(X, Y, transpose(C))
#     or
#         pcolor(transpose(C))
# -----------
# mgrid works the opposite of meshgrid()! 

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
    obspos = np.array([runhdr.get("RA_PNT"), runhdr.get("DEC_PNT")])
    
    X = events.data.field("DETX")
    Y = events.data.field("DETY")
    D = np.sqrt(X**2 + Y**2)

    th2hist,ed = np.histogram( D*D,bins=bins, range=range,normed=False, new=True)
    th2hist = th2hist.astype(np.float64)

    nevents = sum(th2hist)
    areaperbin = math.pi*(ed[1]) # ed[1]=r^2 (each bin has equal 2D area)
    th2hist /= areaperbin      # now in units of counts/deg^2 (density)

    print "edges in t^2: ", ed
    print "edges in   t: ", np.sqrt(ed)
    print "total events in list: ",len(X)
    print "total events in hist: ",nevents
    print " area per radial bin: ", areaperbin,"deg^2"
    print "SC HIST: ",th2hist

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
    

def makeDistanceMap(imagehdu,posRADec=None):
    """ 
    Generate a map of angular distances to a given observation
    position.
    
    Arguments:
    - `imagehdu`: input image HDU (data + header)

    - `pos`: ra/dec position of observation position (will be
      transformed properly to map coordinates)

    returns: data array of distances from each bin center to the position
    """

    # note in Kapteyn/FITS, the center of the first bin is (1,1) in
    # pixel coordinates.

    proj = wcs.Projection( imagehdu.header )
#    tran = wcs.Transformation( wcs.fk5, proj.skysys )
#    pos = tran(posRADec) # transform to map's coordinates

    if (posRADec==None):
        # get the position from the map center:
        pos= proj.toworld(np.array(proj.naxis)/2.0+1.0)
    else:
        pos = posRADec

    print "DISTANCE FROM POSITION: ", pos
    

    (nx,ny) = np.array(imagehdu.data.transpose()).shape # needs the
                                                        # transpose of
                                                        # FITS data!
    
    ia,ja = makePixCoordGrid(nx,ny)

    (ra,dec) = proj.toworld( (ia.flatten(), ja.flatten()) )

    distmap = angSepDeg( ra, dec,  pos[0], pos[1] )
    distmap.shape = ia.shape

    return distmap


def makePixCoordGrid( nx, ny ):
    """
    returns a grid of pixel centers in FITS pixel coordinates, given
    nx,ny bins . The grid returned is (X,Y), where X=projected RA
    coordinate, and Y is projected Dec coordinate (for J2000 system,
    or l,b for galactic)

    EXAMPLE: ia,ja = makePixCoordGrid( 10, 50 )

    """
    return np.mgrid[0:nx,0:ny]+1.0   # shift +1 to fits standard (1,1)
                                     # is center of first bin


def makeDataCubeHDU( centerRADec, erange=(0.1,100), nlogebins=10, **kwargs ):
    """
    Generate a blank FITS datacube, with energy as the 3rd
    dimension. Takes the same key arguments as makeImageHDU, plus
    energy range info.

    """
    imhdu = makeImageHDU( centerRADec, *kwargs )
    imhdu.header.update("CTYPE3", "logE")
    imhdu.header.update("CUNIT3", "log10(TeV)")
    imhdu.header.update("CRPIX3", 0 )
    imhdu.header.update("CRVAL3", log10(erange(0)))
    
    erange = np.array(erange)
    logerange = np.log10(erange)
    
    dlogE = (logerange[1]-logerange[0])/float(nlogebins)
    imhdu.header.update("CDELT3", dlogE )
    
    return imhdu
    
    
    

def makeImageHDU(centerRADec,geom=(300,300), FOV=(1.0,1.0), projection="CAR",
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
  
    centerRADec = np.array(centerRADec)
    geom = np.array(geom)
    FOV = np.array(FOV)
    delta = FOV/geom
    
    hdu = pyfits.NP_pyfits.PrimaryHDU()
    hdu.data = np.zeros( geom ).transpose()
    
    hdu.update_header()

    # Set up a cartesian projection that is basically flat at the
    # camera pointing position (so little distortion). This projection
    # isn't so good if you want to stitch together maps to make a
    # scan, but it minimizes distortion in the field of view, and
    # makes things like the acceptance and ring maps circular. Note
    # that RA/Dec lines will not be parallel in this projection in
    # general.

    hdu.header.add_comment("Generated by actutils.py makeImageHDU")
    hdu.header.update( "CTYPE1", "RA---"+projection )
    hdu.header.update( "CTYPE2", "DEC--"+projection )
    hdu.header.update( "CUNIT1", "degrees")
    hdu.header.update( "CUNIT2", "degrees")
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


    print "== NEW IMAGE HDU: ==========="
    print hdu.header
    print "============================="

    # NOTE: the following could be used for a projection with pixels
    # square in RA/Dec (as in the HESS software), but generally that's
    # not what we want since at high declination, it gives distorted
    # images (and the acceptance will be an oval). However, since the
    # tools take into account the projection correctly, it should
    # still work

    #    hdu.header.update( "CRVAL1", 0.0 )
    #    hdu.header.update( "CRVAL2", 0.0 )
    #    hdu.header.update( "LONPOLE", 0.0 )
    #    hdu.header.update( "CRPIX1", (geom[0]+1)/2 - centerRADec[0]/delta[0] )
    #    hdu.header.update( "CRPIX2", (geom[1]+1)/2 - centerRADec[1]/delta[1] )

    if output:
        pyfits.writeto(output, header=hdu.header, data=hdu.data )
        
    return hdu

def makeCountMap(hdu, lam,bet, output=None, insystem=wcs.fk5, verbose=False):
    """
    Arguments:
    - `hdu`: hdu containing data array and header with proper WCS info 
    - `ra`: array of RA coordinates
    - `dec`: array of Dec coordinates

    todo: add events coordsystem as an input (defaults to RA/Dec)
    
    """

    proj = wcs.Projection( hdu.header )
    proj.skyout = insystem  # allow galactic or other coordinates
    coords = proj.topixel( zip(lam,bet) ) # convert to pixel coordinates
    coords = np.array(coords)

    bins = np.array(hdu.data.transpose().shape)
    
    lowerLeftEdge = (0.5,0.5) # lower-left edge of the histogram (note
                          # that 1.0,1.0 is the center of the
                          # lower-left bin in pixel coordinates

    upperRightEdge = bins +0.5  # the upper-right 
    imrange = zip( lowerLeftEdge, upperRightEdge) 

    if verbose:
        print "bins:", bins
        print "range: ",imrange
    
    
    H, xedges, yedges = np.histogram2d( coords[:,0],coords[:,1],
                                        bins=bins, range=imrange )

    # now H is in standard numpy coordinates. To write it to a FITS
    # image, need to take the transpose (H.T)

    if verbose:
        print "   H: ", H.shape
        print "data: ", hdu.data.shape
        
    if output:
        print "Writing count map:",output
        pyfits.writeto(output, header=hdu.header, data=H.transpose(), clobber=True )

    return H

def copyHeaders(ihdr,ohdr, keys=None):
    """
    copies useful info form eventlist to image header

    Arguments:
    - `ihdu`: input header
    - `ohdu`: output header
    - `keys`: a list of keys to copy. If None, a default list will be used
    """
    
    if keys==None:
        keys = [ "RUN_ID", "DATE_OBS", "TIME_OBS",
                 "DATE_END", "TIME_END", "TSTART", "TSTOP",
                 "MJDREFI","MJDREFF", "TIMEUNIT", "TIMESYS", "TIMEREF",
                 "TASSIGN", "TELAPSE", "ONTIME", "LIVETIME", "DEADC",
                 "OBJECT", "RA_OBJ", "DEC_OBJ", "RA_PNT", "DEC_PNT",
                 "ALT_PNT", "AZ_PNT", "RADESYS", "OBS_MODE", "TELLIST"]
             
    for key in keys:
        try:
            ohdr.update( key, ihdr[key] )
        except KeyError:
            pass

def makeRadialFOVMask(imagehdu,radius,centerWorld=None):
    """ Generates a image mask that is 1.0 inside the given radius, and 0.0 outside
    
    Arguments:
    - `imagehdu`:
    - `centerWorld`: observation position (if None, then taken from
      RA_PNT/DEC_PNT or map center)
    - `radius`: in degrees
    """

    if centerWorld==None:
        # try to use the pointing position:
        try:
            centerWorld = (imagehdu.header["RA_PNT"],imagehdu.header["DEC_PNT"])
        except KeyError:
            # use center of map in RA/Dec
            print "No RA_PNT keys found in header: using map center"
            proj = wcs.Projection( imagehdu.header )
            tran = wcs.Transformation( proj.skysys, wcs.fk5 )
            centerPix = np.array( proj.naxis )/2.0 + 0.5
            centerWorld = tran(proj.toworld( centerPix ))

    dists = makeDistanceMap( imagehdu, centerWorld )
    mask = np.zeros_like( dists )
    mask[dists<radius] = 1.0
    return mask

def histToFITS(histdata, bins, histrange, name="", valueScale=None):
    """
    turn the histogram output from numpy.histogram2d into a FITS image extension

    Arguments:
    - `histdata`: the data output from numpy.histogram2d (not transposed)
    - `bins`: bin ranges given to numpy.histogram2d
    - `range`: range that was given to histogram2d
    """
    
    fullrange =  np.array( (histrange[0][1]- histrange[0][0], 
                            histrange[1][1]- histrange[1][0]) )
    nbins = np.array(bins)
    delta = fullrange/nbins.astype(float)


    # the lower-left edge of the first bin is (Xed[0],Yed[0]), which
    # is (0.5,0.5) in FITS pixel coordinates (the center of the bin is
    # at (1.0,1.0)):
    bin0pix = np.array((0.5,0.5)) # reference coordinate is lower left
                                  # corner of bin
    bin0coord = np.array( (histrange[0][0], 
                           histrange[1][0]) ) # world coord of
                                              # lower-left corner is
                                              # the lower-part of the
                                              # range we specified

    # now note that this defines the first pixel in FITS coordinates
    # with the center (1.0,1.0). in integer python coordinates it is [0,0]

    # to transform a world value, you need to subtract 1.0 and round
    # down to get the bin number:
    #       ibin = round( Vpix-1.0 )
    # To get the value of ibin, you need to go the other way:
    #       Vworld[ibin] = proj.toworld( ibin+0.5 )

    #  note that histogram2d returns the histogram with X vertical and
    #  Y horizontal (the reason is that for a general DD histogram,
    #  the bin order makes more sense that way), so H[j,i] is the
    #  value for bin [i,j], where i is in the x-direction, and j in
    #  the y-direction for this reason, we transpose here, since FITS
    #  uses the opposite convention:
    ohdu = pyfits.ImageHDU(data=histdata.transpose())
    ohdu.name = name
    ohdu.header.update( "CTYPE1", "LSIZ", "log10(SIZE)" );
    ohdu.header.update( "CUNIT1", "PE");
    ohdu.header.update( "CTYPE2", "LDIS", "log10(impact_param)" );
    ohdu.header.update( "CUNIT2", "m");
    ohdu.header.update( "CDELT1", delta[0] )
    ohdu.header.update( "CDELT2", delta[1] )
    ohdu.header.update( "CRVAL1", bin0coord[0] )
    ohdu.header.update( "CRVAL2", bin0coord[1] )
    ohdu.header.update( "CRPIX1", bin0pix[0] )
    ohdu.header.update( "CRPIX2", bin0pix[1] )

    if (valueScale):
        ohdu.header.update( "BSCALE", 1.0/float(valueScale) )


    return ohdu

def histFromFITS(hdu):
    """
    returns xed,yed,Hist (as in np.histogram2d)

    Arguments:
    - `hdu`: hdu containing the histogram
    """
    
    hist = hdu.data.transpose()
    bins = hist.shape
    proj = wcs.Projection ( hdu.header ) # for going between world and pix coords
    
    # bin edges (not centers) start at 0.5 in pixel coordinates, since
    # centers are at 1.0
    xed,tmp = proj.toworld( (np.arange(bins[0]) + 0.5, np.zeros(bins[0])  + 0.5) )
    tmp,yed = proj.toworld( (np.zeros(bins[1])  + 0.5, np.arange(bins[1]) + 0.5) )
    
    return xed,yed,hist # same format as numpy.histogram2d()
        

def displayFITS(header, data):
    import pylab

    f = pylab.maputils.FITSimage( externalheader=header, externaldata=data)
    frame = pylab.figure()
    img = f.Annotatedimage(frame)
    colorbar = img.Colorbar()
    img.Image()
    img.Graticule()

    img.plot()
    img.interact_imagecolors()
    img.interact_toolbarinfo()
    img.interact_writepos()

    pylab.plt.title( fname )
    pylab.plt.show()


