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

if __name__ == '__main__':
    parser = OptionParser()
    (opts,args) = parser.parse_args()

    Ntrue_tot = None
    Nreco_tot = None
    NtrueAthrown_tot = None
    NrecoAthrown_tot = None
    Nthrown_tot = None
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
        Ntrue,bins = np.histogram( np.log10(Etrue), bins=nbins, range=histrange )
        Nreco,bins = np.histogram( np.log10(Ereco), bins=nbins, range=histrange )
        
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
        
        evfile.close()


    print "Summed",count,"simulation runs"
    Aeff_reco = NrecoAthrown_tot/Nthrown_tot
    Aeff_true = NtrueAthrown_tot/Nthrown_tot

    # todo: calculate statistical errors

    plot( bins[0:-1], Aeff_reco, drawstyle="steps-post",label="Reco", color="red" )
    plot( bins[0:-1], Aeff_true, drawstyle="steps-post",label="True", color="blue" )
    legend()
    title("Effective Area")
    xlabel("$Log_{10}(E)$")
    ylabel("$A_{\mathrm{eff}} (\mathrm{m}^2)$")
        
