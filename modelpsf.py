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

def fit_gaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, covar, info,msg,stat = optimize.leastsq(errorfunction, params, full_output=1)
    return p,sqrt(diag(covar))


if __name__ == '__main__':
    
    psfcube = sys.argv[1]    
    psf = pyfits.open( psfcube )[1]

    psfhist = Histogram( initFromFITS=psf )

    N = psf.data.shape[0]
    width = list()
    widthErr = list()
    residuals = list()
    showplot = False
    allenergies = psfhist.binCenters(2)
    dE = (allenergies[1]-allenergies[0])/2.0
    energies=list()

    for ii in range(N):

        if (psf.data[ii].sum() <= 20):  # require at least 20 counts
            print "skip"
            continue
        
        params,errors = fit_gaussian(psf.data[ii])
        print ii,":",params
        fitfunc = gaussian( *tuple(params) )
        model = fitfunc(*indices(psf.data[ii].shape))    
        width.append( maximum(params[3],params[4]) )
        widthErr.append( maximum( errors[3],errors[4] ))
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



    figure()
    title("PSF")
    errorbar( energies, width, yerr=widthErr, xerr=dE,fmt=None )
    ylim(0,5)
    xlabel("Log10(E)")
    ylabel("Max sigma (bins)")
    savefig("modelpsf.pdf")

    figure()
    scatter( arange(len(residuals)), residuals )
    xlabel("Log10(E)")

    show()
