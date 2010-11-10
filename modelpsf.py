# A very simple example code to fit a 2D gaussian model to 2D PSF
# datacube and plot the PSF. The real PSF probably needs a better
# model than this (double Gaussian at least, with a rotation angle)

import sys
import pyfits
from scipy import optimize
from pylab import *
from fitshistogram import Histogram

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def doublegaussian( h1,h2,cx1,cx2,cy1,cy2,wx1,wx2,wy1,wy2 ):
    return lambda x,y: (gaussian(h1,cx1,cy1,wx1,wy1)(x,y) 
                        + gaussian(h2,cx2,cy2,wx2,wy2)(x,y))

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    if (total==0):
        return 0,0,0,0,0
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fit_double_gaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    moms = moments(data)
    params = ravel(zip(moms,moms))
    errorfunction = lambda p: ravel(doublegaussian(*p)(*indices(data.shape)) -
                                    data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fit_gaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                    data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


if __name__ == '__main__':
    
    psfcube = sys.argv[1]    
    psf = pyfits.open( psfcube )[1]

    psfhist = Histogram( initFromFITS=psf )

    N = psf.data.shape[0]
    sig1 = list()
    sig2 = list()
    residuals = list()
    showplot = True
    allenergies = psfhist.binCenters(2)
    energies=list()
    
    for ii in range(N):

        if (psf.data[ii].sum() <= 1e-10): 
            print "skip"
            continue
        
        params = fit_gaussian(psf.data[ii])
        print ii,":",params
        fitfunc = gaussian( *tuple(params) )
        model = fitfunc(*indices(psf.data[ii].shape))    
        sig1.append( params[3] )
        sig2.append( params[4] )
        energies.append( allenergies[ii] )
        residuals.append( np.sum(psf.data[ii] - model) )

        if showplot:
            subplot(N,3,ii*3+1)
            pcolor( psf.data[ii] )
            contour(model, cmap=cm.copper)
            subplot(N,3,ii*3+2)
            pcolor( model )
            subplot(N,3,ii*3+3)
            pcolor( psf.data[ii]-model )
            colorbar()



    maxsig = np.maximum( sig1,sig2 )

    figure()
    title("PSF")
    scatter( energies, maxsig, color='b', label="sigmamax" )
    plot( energies, maxsig, color='b' )
    xlabel("Log10(E)")
    ylabel("Max sigma (bins)")


