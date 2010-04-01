import pyfits
import numpy as np
from optparse import OptionParser
from kapteyn import wcs
import actutils
from fitshistogram import Histogram

# generates response functions from simulation data (e.g. effective
# area and Energy-resolution)

if __name__ == '__main__':
    parser = OptionParser()
    (opts,args) = parser.parse_args()

    for eventlist_filename in args:

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
        nbins = len(Nthrown)
        histrange = (logEmin[0],logEmax[-1])

        print "Histogram:",nbins,"bins, range=",histrange
        
        # Effective area is Nreco/Nthrown
