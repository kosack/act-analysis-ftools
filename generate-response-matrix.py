import pyfits
import numpy as np
from optparse import OptionParser
import actutils
from fitshistogram import Histogram
import math
import pylab 


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
    
    # TODO: add all dependent parameters for effective areas (offset angle, azimuth, etc.)

    coldefs = pyfits.ColDefs( [colE_lo,colE_hi,colSpecResp])
    hdu = pyfits.new_table( coldefs )
    # put in required headers defined by OGIP:

    hdu.header.update("HDUCLASS", "OGIP", "Organization of definition" );
    hdu.header.update("HDUDOC", "CAL/GEN/92-019", "Document describing format" );
    hdu.header.update("HDUCLAS1", "RESPONSE" );
    hdu.header.update("HDUCLAS2", "EFFAREA" );
    hdu.header.update("HDUVERS", "1.0.0" );
    hdu.name = "EFF_AREA"
    hdu.writeto( outputFileName, clobber=True )
    print "DEBUG:",hdu.header


def writeRMF(responseHist, outputFileName="spec_rmf.fits"):

    # The redistribution matrix is just a set of 1D histogram (of
    # probability), one for each true-energy bin, that has N_CHAN
    # reco-energy bins.

    # `channel`: a reco-energy bin, one element of the redistribution
    # function, E_min to E_max

    # `energy bin`: bin in true energy, E_LO to E_HI



    Etrue_bins,Ereco_bins = responseHist.binLowerEdges
    Etrue_bins = 10**Etrue_bins
    Ereco_bins = 10**Ereco_bins
    nchan = len(Ereco_bins)-1
    nebin = len(Etrue_bins)-1

    chan = np.arange( 1, nchan )

    # first create the EBOUNDS extension (that defines the energy
    # ranges for each "channel"

    col_chan = pyfits.Column( name="CHANNEL", format='I',array=chan )
    col_emin = pyfits.Column( name="E_MIN", format='E', unit='TeV',
                              array=Ereco_bins[0:-1] )
    col_emax = pyfits.Column( name="E_MAX", format='E', unit='TeV',
                              array=Ereco_bins[1:] )

    ebounds = pyfits.new_table( [col_chan, col_emin, col_emax])

    
    # now create the RMF table. FOr now there is only a fixed
    # zenith/azimuth, etc, so we can just just a single group of
    # channels (for the real data, this will be more complicated)

    col_elo = pyfits.Column( name="E_LO", format='E', unit='TeV',
                              array=Etrue_bins[0:-1] )
    col_ehi = pyfits.Column( name="E_HI", format='E', unit='TeV',
                              array=Etrue_bins[1:] )
    col_ngrp = pyfits.Column( name="N_GRP", format='E', unit='TeV',
                              array=ones(nebins) )

    col_ngrp = pyfits.Column( name="N_GRP", format='I', unit='TeV',
                              array=ones(nebins) )

    col_fchan = pyfits.Column( name="F_CHAN", format='10I', unit='TeV',
                              array=ones(nebins) )

    col_nchan = pyfits.Column( name="N_CHAN", format='10I', unit='TeV',
                              array=ones(nebins) )

    # col_matrix = pyfits.Column( name="MATRIX", format='10I', unit='TeV',
    #                             array=ones(nebins) )



    ebounds.writeto( outputFileName, clobber=True )
    

    pass


if __name__ == '__main__':
    parser = OptionParser()
    (opts,args) = parser.parse_args()

    Ntrue_tot = None
    Nreco_tot = None
    NtrueAthrown_tot = None
    NrecoAthrown_tot = None
    Nthrown_tot = None
    energyResponseHist = Histogram( range=[[-1,2],[-1,2]], bins=[100,100],
                                    axisNames=["$\log_{10}(E_{true})$", 
                                               "$\log_{10}(E_{reco})$"])

    energyResolutionHist = Histogram( range=[[-1,2],[-1,1]], bins=[50,50],
                                      axisNames=["$\log_{10}(E_{true})$", 
                                                 "$\log_{10}(E_{reco}/E_{true})$"])
    count = 0

    for eventlist_filename in args:

        count += 1


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

        print eventlist_filename,": Athrown=",Athrown
        
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

    writeARF( Emin, Emax, Aeff_reco*(100**2), 
              outputFileName="spec_arf_reco.fits" )
    writeARF( Emin, Emax, Aeff_true*(100**2), 
              outputFileName="spec_arf_true.fits" )

    # Normalize the phonton distribution matrix (the integral along
    # the vertical axis should be 1.0, since it's a probability)
    np.apply_along_axis( normalizeToProb, arr=energyResolutionHist.hist, axis=1)

    # TODO: calculate statistical errors


    pylab.figure( figsize = (15.0,11.0), dpi=80 )
    pylab.subplot(2,1,1)

    pylab.semilogy()
    pylab.plot( bins[0:-1], Aeff_reco, drawstyle="steps-post",
          label="Reco", color="red" )
    pylab.plot( bins[0:-1], Aeff_true, drawstyle="steps-post",
          label="True", color="blue" )
    pylab.legend(loc="lower right")
    pylab.title("Effective Area")
    pylab.xlabel("$Log_{10}(E/\mathrm{TeV})$")
    pylab.ylabel("$A_{\mathrm{eff}} (\mathrm{m}^2)$")
    pylab.grid()


    # TODO: include background here, and 5 sigma requirement:
    # obstimeHrs = 50.0
    # obstime = obstimeHrs*60.0*60.0
    # nevents= 10.0
    # M2_CM2 = 100*2
    # sensitivity = nevents/(Aeff_reco*M2_CM2*obstime)
    # sensitivity[np.isinf(sensitivity)] = 0
    # pylab.subplot(2,2,2)
    # pylab.semilogy()
    # pylab.plot( bins[0:-1], sensitivity, drawstyle="steps-post" )
    # pylab.title("Sensitivity (%.1f hours, %d events)"%(obstimeHrs, nevents))
    # pylab.xlabel("$\log_{10}(E/\mathrm{TeV})$")
    # pylab.ylabel("$(dN/dE)_{min} (\mathrm{cm^{-2} s^{-1} TeV^{-1}})$")
    # pylab.grid()

    pylab.subplot(2,2,3)
    energyResponseHist.draw2D()
    l = energyResponseHist.binLowerEdges[0]

    pylab.plot( l,l, color="white") 
    pylab.colorbar()
    pylab.title("")
    

    # make this histogram normalized to have an integral of 1.0 along
    # the Y axis (so it is now a probability of reconstructing Ereco
    # for a given Etrue)
    pylab.subplot(2,2,4)
    energyResolutionHist.draw2D( vmax=0.25 )
    pylab.colorbar()
    pylab.title("Energy Response")
    l = energyResolutionHist.binLowerEdges[0]
    pylab.plot( l,np.zeros_like(l), color="black")
    pylab.grid(color='w')
    pylab.savefig("response.pdf", papertype="a4")


    writeRMF( energyResponseHist )

 
