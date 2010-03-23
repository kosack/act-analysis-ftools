import pyfits
import numpy as np
import math
#from scipy import signal,ndimage

from kapteyn import wcs,maputils


class Histogram(object):
    """
    A simple N-dimensional histogram class with FITS IO 

    Uses numpy.histogram2d to generate histograms and read and write
    them to FITS files
    """
    
    
    def __init__(self, bins=None, range=None, name=None,
                 axisTypes=None, axisNames=None,initFromFITS=None):
        """ Initialize an unfilled histogram (need to call either
        fill() or loadFromFITS() to put data into it)

        `bins`: array listing binsize for each dimension
        `ranges`: array of ranges for each dimension
        `name`: name (use for FITS extension name) 
        """
        
        self.hist = np.zeros(bins)
        self._binLowerEdges= None  #TODO: should be a property, get only
        self._bins = bins
        self._range=range
        self.valueScale = None
        self.valueZero = None
        self.axisTypes=axisTypes
        self.axisNames=axisNames
        self.name = name
        self._proj = None

        if (initFromFITS):
            self.loadFromFITS(initFromFITS)

    @property
    def binLowerEdges(self):
        return self._binLowerEdges

    @property
    def bins(self):
        return self._bins
        
    @property
    def range(self):
        return self._range


    def fill(self, datapoints, **kwargs):
        """
        generate a histogram from data points
        
        Arguments:
        - `datapoints`: array of points (see numpy.histogramdd() documentation)
        """

        hist, binLowerEdges = np.histogramdd( datapoints, 
                                              bins=self._bins, 
                                              range=self._range, **kwargs)
        
        
        if (self._binLowerEdges == None):
            self._binLowerEdges = binLowerEdges

#        if ((binLowerEdges==self._binLowerEdges).all() == False):
#            raise exceptions.ArithmaticError("Bad Geometry!")

            
        self.hist += hist

        


    def asFITS(self):
        """ 
        return A FITS hdu, suitable for writing to disk
        
        to write it, just do 

        myhist.asFITS().writeto("outputfile.fits")

        """
        ohdu = pyfits.ImageHDU(data=self.hist.transpose())
        ohdu.name = self.name
        ndim = len(self._bins)

        # the lower-left edge of the first bin is (Xed[0],Yed[0]), which
        # is (0.5,0.5) in FITS pixel coordinates (the center of the bin is
        # at (1.0,1.0))


        # now note that this defines the first pixel in FITS coordinates
        # with the center (1.0,1.0). in integer python coordinates it is [0,0]

        # to transform a world value, you need to subtract 1.0 and round
        # down to get the bin number:
        #       ibin = round( Vpix-1.0 )
        # To get the value of ibin, you need to go the other way:
        #       Vworld[ibin] = proj.toworld( ibin+0.5 )
        

        for dim in xrange(ndim):
            width = self._range[dim][1] - self._range[dim][0]
            num = self._bins[dim]
            delta = width/float(num)
            bin0pix = 0.5 # lower-left corner of first bin
            bin0coord = self._range[dim][0]

            ohdu.header.update("CTYPE%d"%(dim+1),"LIN","Linear parameter") 
            ohdu.header.update("CDELT%d"%(dim+1), delta)
            ohdu.header.update("CRVAL%d"%(dim+1), bin0coord)
            ohdu.header.update("CRPIX%d"%(dim+1), bin0pix)
            
        if (self.valueScale):
            ohdu.header.update( "BSCALE", 1.0/float(self.valueScale) )

        if (self.valueZero):
            ohdu.header.update( "BZERO", float(self.valueZero) )

   
        return ohdu


    def loadFromFITS(self, inputFITS):
        """
        load a FITS histogram file
        
        Arguments:
        - `inputfits`: filename or fits HDU
        """

        if (type(inputFITS) == str):
            hdu = pyfits.open(inputFITS).pop(0)
        else:
            hdu = inputFITS
            
        self.hist = hdu.data.transpose()
        self._bins = self.hist.shape

        proj = wcs.Projection ( hdu.header )
        ndim = len(self._bins)
        
        edges = []
        self._range = []

        for dim in xrange(ndim):
            # note that histogramdd returns edges for 0-N+1 (including
            # the lower edge of the non-existant next bin), so we do
            # the same here to keep things the same
            a = np.zeros( (self._bins[dim]+1, ndim ) )
            a[:,dim] = np.arange( self._bins[dim]+1 )+0.5
            edges.append(proj.toworld( a )[:,dim])
            self._range.append( (edges[dim][0],edges[dim][-1]) )
            
        self._binLowerEdges = edges 
        self._range = np.array(self._range)
        
        if(hdu.header.get("BSCALE")):
            self.valueScale = hdu.header["BSCALE"]

        if(hdu.header.get("BSCALE")):
            self.valueZero = hdu.header["BZERO"]
    

    def getProjection(self ):
        """
        """
        
        if (self._proj == None):
            self._proj = wcs.Projection( self.asFITS().header )

        return self._proj

        
    def getValue(self, coords, outlierValue=None):
        """ Returns the values of the histogram at the given world
        coordinate(s) 
        
        Arguments:

        - `coords`: list of M coordinates of dimension N (where
          the N is the histogram dimension)

        - `outlierValue`: value for outliers, if None, coordinates
          outside the edges of the histogram will be given the edge
          value

        """

        if np.isnan(coords).any():
            raise ValueError("Bad coordinate value")

        world = np.array( coords, ndmin=2 )
        ndims = len(self._bins)
        maxbin = np.array(self.hist.shape)

        bins = np.array([np.digitize( world[:,ii], self._binLowerEdges[ii] )-1 
                         for ii in xrange(ndims)])

        print "DEBUG bins=",bins

        if (outlierValue==None):
            #extrapolate (simply for now, just takes edge value)
            bins[bins<0] = 0
            for ii in xrange(ndims):
                bins[ii][bins[ii]>=maxbin[ii]] = maxbin[ii]-1 
        else:
            if (bins>=maxbin).any() or (bins<0).any():
                return outlierValue
        
        print "DEBUG *bins=",bins

        return self.hist[tuple(bins)]
                
