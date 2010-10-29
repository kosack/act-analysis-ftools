import pyfits
import numpy as np
import math

from kapteyn import wcs,maputils
from matplotlib import pyplot as plt

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

def makeRadialProfile(events,bins=14,range=[0,10], 
                      verbose=False, offset=[0,0],
                      squaredBins=False):
    """ 
    Generates an radial profile from the events given in Detector
    coordinates (which are assumed to have (0,0) as the origin)
    
    Arguments:
    - `events`: hdu of cut and masked eventlist (uses header and event info)
    """
    
    # build 1D histogram in detector coordinates (pointing dir is 0,0,
    # and distances are in degrees from center) 

    X = events.data.field("DETX") + offset[0]
    Y = events.data.field("DETY") + offset[1]
    if (squaredBins):
        D = X**2 + Y**2
    else:
        D = np.sqrt(X**2 + Y**2)

    th2hist,ed = np.histogram( D*D,bins=bins, range=range,normed=False, new=True)
    th2hist = th2hist.astype(np.float64)

    nevents = sum(th2hist)
    areaperbin = math.pi*(ed[1]) # ed[1]=r^2 (each bin has equal 2D area)
    th2hist /= areaperbin      # now in units of counts/deg^2 (density)

    if verbose:
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

    if (posRADec==None):
        # get the position from the map center:
        pos= proj.toworld(np.array(proj.naxis)/2.0+1.0)
    else:
        tran = wcs.Transformation( wcs.fk5, proj.skysys )
        pos = tran(posRADec) # transform to map's coordinates


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
    
    
    

def makeImageHDU(center,geom=(300,300), FOV=(1.0,1.0), projection="CAR",
                 output=None, system="fk5"):
    """
    Generate a blank FITS image. Returns a HDU with the given geometry
    (so far in RA/Dec only). 

    Arguments:
    - `center`: (lambda,beta) center of image
    - `geom`: (N,M) size in pixels of the image
    - `FOV_deg`: (dx,dy) size of FOV in degrees
    - `output`: output filename, or None for no output
    """
  
    center = np.array(center)
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

    if (system == "fk5"):
        hdu.header.update( "CTYPE1", "RA---"+projection )
        hdu.header.update( "CTYPE2", "DEC--"+projection )
    elif (system == "galactic"):
        hdu.header.update( "CTYPE1", "GLON-"+projection )
        hdu.header.update( "CTYPE2", "GLAT-"+projection )
    else:
        raise ValueError("Unknown system: %s" %system)

    hdu.header.update( "CUNIT1", "degrees")
    hdu.header.update( "CUNIT2", "degrees")
    hdu.header.update( "RADESYS", "fk5" )
    hdu.header.update( "EQUINOX", 2000 )
    
    hdu.header.update( "CRVAL1", center[0] )
    hdu.header.update( "CRVAL2", center[1] )
    hdu.header.update( "CRPIX1", geom[0]/2 )
    hdu.header.update( "CRPIX2", geom[1]/2 )

    hdu.header.update( "CDELT1", delta[0] )
    hdu.header.update( "CDELT2", delta[1] )
    
#    hdu.header.update( "BSCALE", 1.0 )
#    hdu.header.update( "BZERO", 0.0 )


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

def makeCountMap(hdu, lam,bet, output=None, verbose=False):
    """
    Arguments:
    - `hdu`: hdu containing data array and header with proper WCS info 
    - `ra`: array of RA coordinates
    - `dec`: array of Dec coordinates

    todo: add events coordsystem as an input (defaults to RA/Dec)
    
    """

    proj = wcs.Projection( hdu.header )

    tran = wcs.Transformation( wcs.fk5, proj.skysys )
    print "Transforming from system: ",wcs.fk5,"to",proj.skysys
    X,Y = tran((lam,bet)) # transform to target coordinate system

    coords = proj.topixel( zip(X,Y) ) # convert to pixel coordinates
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




def displayFITS(header, data):
    """
    Display a FITS image using Kapteyn
    """
    import pylab

    f = maputils.FITSimage( externalheader=header, externaldata=data)
    fig = plt.figure()
    frame = fig.add_subplot(1,1,1)
    img = f.Annotatedimage(frame)
    colorbar = img.Colorbar()
    img.Image()
    img.Graticule()

    img.plot()
    img.interact_imagecolors()
    img.interact_toolbarinfo()
    img.interact_writepos()

    if header.has_key("EXTNAME"):
        plt.title( header["EXTNAME"] )
    plt.show()

def displayFITSHDU( hdu ):        
    displayFITS(hdu.header, hdu.data)

def getTelTypeMap( telarray_hdu ):
    """
    returns two dictionaries: telId2Type,type2TelId mapping telescope
    id number to telescope type (class), and vice-versa, given the
    TELARRAY hdu of an eventlist.

    Telescope types are a tuple consisting of the (type,subtype)
        
    """

    telid    = telarray_hdu.data.field("TELID")
    tclass    = telarray_hdu.data.field("CLASS")
    tsubclass = telarray_hdu.data.field("SUBCLASS")
    
    telid_to_type = dict()
    type_to_telid = dict()

    for tid,tcl,tsub in zip(telid,tclass,tsubclass):
        telid_to_type[tid] = (tcl,tsub)
        if type_to_telid.has_key( (tcl,tsub) ):
            type_to_telid[(tcl,tsub)].append(tid)
        else:
            type_to_telid[(tcl,tsub)] = [tid]

    return telid_to_type,type_to_telid


def loadLookupTableColumns( events, telarray ):

    try:
        tposx = telarray.data.field("POSX")
        tposy = telarray.data.field("POSY")
        telid = telarray.data.field("TELID")
        telMask   = events.data.field("TELMASK")
        telImpacts = events.data.field("HIL_TEL_IMPACT")
    except KeyError:
        print "One or more required columns was missing in the eventlist!"
        raise

    # sizes are already by telescope
    try:
        telSizes = events.data.field("HIL_TEL_SIZE") 
    except KeyError:
        telSizes = events.data.field("INT_TEL_SIZE") 

    nevents,ntels = telSizes.shape



    
    # apply some basic cuts to get rid of bad values
    valueMask =  np.isfinite(telImpacts) * np.isfinite(telSizes)
    valueMask *=  (telImpacts > 0.0) 
    valueMask *=  (telSizes > 0.0) 

    telMask *= valueMask  # mask off bad values

    return (telImpacts, np.log10(telSizes), telid,telMask )
