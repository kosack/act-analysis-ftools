import pyfits
import numpy as np
from optparse import OptionParser
from kapteyn import wcs
import actutils
from fitshistogram import Histogram
import math
from pylab import *


# generates response functions from simulation data (e.g. effective
# area and Energy-resolution)

# currently this is very naive and assumes the thrown energy histogram
# is in exactly the same binning for all event lists, and that the
# external parameters (zenith angle, azimuth, offset) of the
# simulations are the same.

def normalizeToProb(x):
    """
    Makes x a probability distribution
    Arguments:
    - `x`:
    """
    x[x>0] /= np.sum(x)
    return x

def writeARF(EBinlow, EBinHigh, effectiveArea, outputFileName="spec_arf.fits"):
    """ 
    Write out effective area curve in ARF format
    
    - `EBinlow/High`: lower edges of bins (in linear energy units)
    - `effectiveArea`: effective area in cm^2 for each energy bin
    """

    colE_lo = pyfits.Column( name="ENERG_LO", format="E",
                             unit="TeV",array=EBinlow )
    colE_hi = pyfits.Column( name="ENERG_HI", format="E", 
                             unit="TeV", array=EBinHigh )
    colSpecResp = pyfits.Column( name="SPECRESP", format="E", 
                                 unit="cm2",array=effectiveArea )
    
    coldefs = pyfits.ColDefs( [colE_lo,colE_hi,colSpecResp])
    hdu = pyfits.new_table( coldefs )
    hdu.writeto( outputFileName, clobber=True )
    print "DEBUG:",hdu.header


if __name__ == '__main__':
    parser = OptionParser()
    (opts,args) = parser.parse_args()

    Ntrue_tot = None
    Nreco_tot = None
    NtrueAthrown_tot = None
    NrecoAthrown_tot = None
    Nthrown_tot = None
    energyResponseHist = Histogram( range=[[-1,2],[-1,2]], bins=[70,70],
                                    axisNames=["log10(Etrue)", "log10(Ereco)"])

    energyResolutionHist = Histogram( range=[[-1,2],[-1,1]], bins=[50,50],
                                      axisNames=["log10(Etrue)", 
                                                 "log10(Ereco/Etrue)"])
    count = 0

    for eventlist_filename in args:

        count += 1
        print eventlist_filename

        evfile = pyfits.open(eventlist_filename)
        events = evfile["EVENTS"]
        telarray = evfile["TELARRAY"]
        mcinfo = evfile["MCINFO"]
        mcenergy = evfile["MCENERGY"]
        tel2type,type2tel = actutils.getTelTypeMap( telarray )
        
        Emin = mcenergy.data.field("E_MIN")
        Emax = mcenergy.data.field("E_MAX")
        logEmin = np.log10(Emin)
        logEmax = np.log10(Emax)
        dlogE = logEmax - logEmin
        Nthrown = mcenergy.data.field("N")
        Nthrown_err = mcenergy.data.field("N_ERR")

        rmin = mcinfo.data.field("CORE_MIN")[0]
        rmax = mcinfo.data.field("CORE_MAX")[0]
        Athrown = math.pi * (rmax-rmin)**2

        nbins = len(Nthrown)
        histrange = (logEmin[0],logEmax[-1])

        print "Athrown=",Athrown
        print "Histogram:",nbins,"bins, range=",histrange
        
        # Effective area is Athrown*(Nreco/Nthrown)
        
        Etrue = events.data.field("MC_ENERGY")
        Ereco = events.data.field("ENERGY")
        Ntrue,bins = np.histogram( np.log10(Etrue), bins=nbins, 
                                   range=histrange )
        Nreco,bins = np.histogram( np.log10(Ereco), bins=nbins, 
                                   range=histrange )
        
        if (NrecoAthrown_tot == None):
            NrecoAthrown_tot = Nreco.copy()*Athrown
            NtrueAthrown_tot = Ntrue.copy()*Athrown
            Nthrown_tot = Nthrown.copy()
            Nreco_tot = Nreco.copy()
            Ntrue_tot = Ntrue.copy()
        else:
            NrecoAthrown_tot += Nreco*Athrown
            NtrueAthrown_tot += Ntrue*Athrown
            Nthrown_tot += Nthrown
            Nreco_tot += Nreco
            Ntrue_tot += Ntrue
        
        # fill 2D histograms

        energyResponseHist.fill( zip(np.log10(Etrue),np.log10(Ereco)) )
        energyResolutionHist.fill( zip(np.log10(Etrue), 
                                       np.log10(Ereco)-np.log10(Etrue)))
        evfile.close()


    print "Summed",count,"simulation runs"
    Aeff_reco = NrecoAthrown_tot/Nthrown_tot
    Aeff_true = NtrueAthrown_tot/Nthrown_tot

    writeARF( Emin, Emax, Aeff_reco*(100**2) )

    # Normalize the phonton distribution matrix (the integral along
    # the vertical axis should be 1.0, since it's a probability)
    np.apply_along_axis( normalizeToProb, arr=energyResolutionHist.hist, axis=1)

    # TODO: calculate statistical errors

    subplot(2,1,1)

    semilogy()
    plot( bins[0:-1], Aeff_reco, drawstyle="steps-post",
          label="Reco", color="red" )
    plot( bins[0:-1], Aeff_true, drawstyle="steps-post",
          label="True", color="blue" )
    legend(loc="lower right")
    title("Effective Area")
    xlabel("$Log_{10}(E)$")
    ylabel("$A_{\mathrm{eff}} (\mathrm{m}^2)$")
    grid()

    subplot(2,2,3)
    energyResponseHist.draw2D()
    l = energyResponseHist.binLowerEdges[0]
    plot( l,l, color="white") 
    colorbar()
    title("")
    

    # make this histogram normalized to have an integral of 1.0 along
    # the Y axis (so it is now a probability of reconstructing Ereco
    # for a given Etrue)
    subplot(2,2,4)
    energyResolutionHist.draw2D( vmax=0.25 )
    colorbar()
    title("Normalized Energy Resolution")
    l = energyResolutionHist.binLowerEdges[0]
    plot( l,zeros_like(l), color="black")



