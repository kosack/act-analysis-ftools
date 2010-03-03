import pyfits
import numpy as np
import math

from kapteyn import wcs,maputils


class FITSHistogram(object):
    """
    A simple N-dimensional histogram class with FITS IO 

    Uses numpy.histogram2d to generate histograms and read and write
    them to FITS files
    """
    
    
    def __init__(self, bins=(0), ranges=[[0,0]], name=None):
        """ Initialize an unfilled histogram (need to call either
        fill() or loadFromFITS() to put data into it)

        `bins`: array listing binsize for each dimension
        `ranges`: array of ranges for each dimension
        `name`: optional name/title of the histogram 
        """
        
        self.hist = None
        self.binLowerEdges= None
        self.bins = bins
        self.ranges=ranges
        self.valueScale = None
        self.axisTypes=None
        self.name = name


    def fill(self, datapoints, **kwargs):
        """
        generate a histogram from data points
        
        Arguments:
        - `datapoints`: array of points (see numpy.histogramdd() documentation)
        """

        hist, binLowerEdges = np.histogramdd( datapoints, 
                                              bins=self.bins, 
                                              range=self.ranges, *kwargs)
        
        if (self.hist != None and any(binLowerEdges==self.binLowerEdges) ==false):
            raise exceptions.ArithmaticError("Bad Geometry!")
        
        if (self.hist == None):
            self.hist = hist
        else:
            self.hist += hist

        self.binLowerEdges = binLowerEdges
        


    def asFITS(self):
        """ 
        return A FITS hdu, suitable for writing to disk
        
        to write it, just do 

        myhist.asFITS().writeto("outputfile.fits")

        """
        ohdu = pyfits.ImageHDU(data=self.hist.transpose())
        ndim = len(self.bins)

        for dim in xrange(ndim):
            width = self.ranges[dim][1] - self.ranges[dim][0]
            num = self.bins[dim]
            delta = width/float(num)
            bin0pix = 0.5 # lower-left corner of first bin
            bin0coord = self.ranges[dim][0]

            ohdu.header.update("CTYPE%d"%(dim+1),"LIN","Linear parameter") 
            ohdu.header.update("CDELT%d"%(dim+1), delta)
            ohdu.header.update("CRVAL%d"%(dim+1), bin0coord)
            ohdu.header.update("CRPIX%d"%(dim+1), bin0pix)
            
        if (self.valueScale):
            ohdu.header.update( "BSCALE", 1.0/float(self.valueScale) )
   
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
        self.bins = self.hist.shape

        proj = wcs.Projection ( hdu.header )
        ndim = len(self.bins)
        
        edges = []

        for dim in xrange(ndim):
            a = np.zeros( (self.bins[dim], ndim ) )
            a[:,dim] = np.arange( self.bins[dim] )+0.5
            edges.append(proj.toworld( a )[:,dim])
            
        self.binLowerEdges =  edges 

        


        
        
