import pyfits
import numpy as np
from optparse import OptionParser
import actutils
from fitshistogram import Histogram
import math
import pylab 
from scipy import optimize


# generates response functions from simulation data (e.g. effective
# area and Energy-resolution)

# currently this is very naive and assumes the thrown energy histogram
# is in exactly the same binning for all event lists, and that the
# external parameters (zenith angle, azimuth, offset) of the
# simulations are the same.

# TODO: use interpolation function for the histograms and then
# integrate the function over the bin width for the new bins (which
# can then be differently spaced)

def normalizeToProb(x):
    """
    Makes x a probability distribution
    Arguments:
    - `x`: array-like
    """
    x[x>0] /= np.sum(x)
    return x

def normalizeToMax(x):
    """
    Makes the maximum value equal to 1.0
    Arguments:
    - `x`: array-like
    """
    x[x>0] /= x.max()
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
    
    # TODO: add all dependent parameters for effective areas (offset
    # angle, azimuth, etc.)

    coldefs = pyfits.ColDefs( [colE_lo,colE_hi,colSpecResp])
    hdu = pyfits.new_table( coldefs )
    # put in required headers defined by OGIP:

    hdu.header.update("HDUCLASS", "OGIP", "Organization of definition" );
    hdu.header.update("HDUDOC", "CAL/GEN/92-019", "Document describing format" );
    hdu.header.update("HDUCLAS1", "RESPONSE" );
    hdu.header.update("HDUCLAS2", "SPECRESP" );
    hdu.header.update("HDUVERS", "1.0.0" );
    hdu.header.update("TELESCOP", "HESS", "Mission name" );
    hdu.header.update("INSTRUME", "HESS1", "Instrument name" );
    hdu.header.update("FILTER", "", "relic from X-ray analysis (required)" );
    hdu.name = "SPECRESP"
    hdu.writeto( outputFileName, clobber=True )



def writeRMF(responseHist, outputFileName="spec_rmf.fits"):
    """
    write out the photon distribution matrix in RMF format
    """
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
    ebounds.name="EBOUNDS"
    ebounds.header.update("HDUCLASS", "OGIP", "Organization of definition" );
    ebounds.header.update("HDUCLAS1", "RESPONSE", 
                         "dataset relates to spectral response" );
    ebounds.header.update("HDUCLAS2", "EBOUNDS", "dataset is a response matrix" );
    ebounds.header.update("HDUCLAS3", "REDIST", "photon redistribution matrix" );
    ebounds.header.update("HDUVERS", "1.3.0", 
                         "Version of format (OGIP memo CAL/GEN/92-002a)" );
    ebounds.header.update("HDUVERS1", "1.3.0", 
                          "Obsolete - included for backwards compatibility" );
    ebounds.header.update("CHANTYPE", "PI", "Required keyword, X-ray relic" );
    ebounds.header.update("TELESCOP", "HESS", "Mission name" );
    ebounds.header.update("INSTRUME", "HESS1", "Instrument name" );
    ebounds.header.update("DETCHANS", nchan )

    # now create the RMF matrix table. For now there is only a fixed
    # zenith/azimuth, etc, so we can just just a single group of
    # channels (for the real data, this will be more complicated)

    # see sec 3.2.1 of: http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html 

    col_elo = pyfits.Column( name="ENERG_LO", format='E', unit='TeV',
                              array=Etrue_bins[0:-1] )
    col_ehi = pyfits.Column( name="ENERG_HI", format='E', unit='TeV',
                              array=Etrue_bins[1:] )
    col_ngrp = pyfits.Column( name="N_GRP", format='I', unit='TeV',
                              array=np.ones(nebin) )

    # TODO: need to set TLMIN/TLMAX values for these.... how do we do
    # it with pyfits?
    col_fchan = pyfits.Column( name="F_CHAN", format='1I', unit='TeV',
                              array=np.ones(nebin) )
    col_nchan = pyfits.Column( name="N_CHAN", format='1I', unit='TeV',
                              array=np.ones(nebin) )

    
    col_matrix = pyfits.Column( name="MATRIX", format='{0:d}E'.format(nchan), 
                                array=responseHist.hist)


    matrix = pyfits.new_table( [col_elo, col_ehi,col_ngrp,
                                col_fchan, col_nchan]) #,col_matrix ] )

    # a hack to add TLMIN/MAX keywords to the F_CHAN and N_CHAN
    # columns (pyfits doesn't seem to support these)
    matrix.header.update("TLMIN4", 1,"Minimum channel number in group")
    matrix.header.update("TLMAX4", nchan, "Maximum channel number in group")
    matrix.header.update("TLMIN5", 1,"Minimum channel number in group")
    matrix.header.update("TLMAX5", nchan, "Maximum channel number in group")

    # set other header keywords
    matrix.name="MATRIX"
    matrix.header.update("HDUCLASS", "OGIP", "Organization of definition" );
    matrix.header.update("HDUCLAS1", "RESPONSE", 
                         "dataset relates to spectral response" );
    matrix.header.update("HDUCLAS2", "RSP_MATRIX", "dataset is a response matrix" );
    matrix.header.update("HDUCLAS3", "REDIST", "photon redistribution matrix" );
    matrix.header.update("HDUVERS", "1.3.0", 
                         "Version of format (OGIP memo CAL/GEN/92-002a)" );
    matrix.header.update("HDUVERS1", "1.3.0", 
                         "Obsolete - included for backwards compatibility" );
    matrix.header.update("CHANTYPE", "PI", "Required keyword, X-ray relic" );
    matrix.header.update("TELESCOP", "HESS", "Mission name" );
    matrix.header.update("INSTRUME", "HESSI", "Instrument name" );
    

    hdulist = pyfits.HDUList( hdus=[pyfits.PrimaryHDU(),ebounds,matrix] )
    hdulist.writeto( outputFileName, clobber=True )
    


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option( "-p","--plot", dest="plot", action="store_true",
                       default=False, help="generate plots")    
    (opts,args) = parser.parse_args()

    Ntrue_tot = None
    Nreco_tot = None
    NtrueAthrown_tot = None
    NrecoAthrown_tot = None
    Nthrown_tot = None
    psf = None

    energyResponseHist = Histogram( range=[[-1,2],[-1,2]], bins=[100,100],
                                    axisNames=["$\log_{10}(E_{true})$", 
                                               "$\log_{10}(E_{reco})$"])

    energyResponseHist.name="Photon Redistrubition Matrix"

    energyBiasHist = Histogram( range=[[-1,2],[-1,1]], bins=[50,50],
                                      axisNames=["$\log_{10}(E_{true})$", 
                                                 "$\log_{10}(E_{reco}/E_{true})$"])

    PSF_MAX_OFFS2 = 0.12
    PSF_ENERGY_RANGE = [-1,2] # in log10(TeV)
    maxpsfoffset = math.sqrt(PSF_MAX_OFFS2)
    
    psfCube = Histogram( range=[[-maxpsfoffset,maxpsfoffset],
                                [-maxpsfoffset,maxpsfoffset],
                                PSF_ENERGY_RANGE], 
                         bins=[80,80,10],
                         axisNames=["DETX","DETY","logE"])
    psfHist = Histogram( range=[[0,PSF_MAX_OFFS2],PSF_ENERGY_RANGE], 
                         bins=[50,10],
                         axisNames=["Offset^2","logE"])
    psfCube.name="PSF2D"
    psfHist.name="PSF1D"

    count = 0

    for ii,eventlist_filename in enumerate(args):

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

        detx = events.data.field("DETX")
        dety = events.data.field("DETY")

        rmin = mcinfo.data.field("CORE_MIN")[0]
        rmax = mcinfo.data.field("CORE_MAX")[0]
        Athrown = math.pi * (rmax-rmin)**2

        nbins = len(Nthrown)
        histrange = (logEmin[0],logEmax[-1])

        print "[{0:3d}/{1:3d}] ".format(ii,len(args)), \
            eventlist_filename, ": Athrown={0:.2f} km^2".format(Athrown/1e6)
        
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
        
        # PSF histogram:
            
        runhdr = events.header
        obspos = np.array([runhdr.get("RA_PNT",0.0), runhdr.get("DEC_PNT",0.0)])
        objpos = np.array([runhdr.get("RA_OBJ",0.0), runhdr.get("DEC_OBJ",0.0)])
        offX = actutils.angSepDeg( obspos[0], obspos[1], objpos[0], obspos[1] )
        offY = actutils.angSepDeg( obspos[0], obspos[1], obspos[0], objpos[1] )


        # fill 2D histograms

        energyResponseHist.fill( np.array(zip(np.log10(Etrue),
                                              np.log10(Ereco))) )
        energyBiasHist.fill( np.array(zip(np.log10(Etrue), 
                                       np.log10(Ereco)-np.log10(Etrue))))
        psfCube.fill( np.array(zip(detx+offX,dety+offY,np.log10(Ereco))) )
        theta2 = (detx+offX)**2+(dety+offY)**2
        psfHist.fill( np.array(zip( theta2,np.log10(Ereco))))
        evfile.close()
        
    


    print "Summed",count,"simulation runs"
    Aeff_reco = NrecoAthrown_tot/Nthrown_tot
    Aeff_true = NtrueAthrown_tot/Nthrown_tot

    writeARF( Emin, Emax, Aeff_reco*(100**2), 
              outputFileName="spec_arf_reco.fits" )
    writeARF( Emin, Emax, Aeff_true*(100**2), 
              outputFileName="spec_arf_true.fits" )

    # Normalize the PSFs to have the sum of 1.0 for each energy bin (axis=2)
#    psfCube.hist = np.apply_along_axis( normalizeToProb, axis=2, arr=psfCube.hist )
#    psfHist.hist = np.apply_along_axis( normalizeToProb, axis=1, arr=psfHist.hist )

    # write out the PSF datacube

    psfHDU = psfCube.asFITS()
    actutils.addInstrumentHeadersToHDU( psfHDU, telescope="HESSI", instrument="HESSI",
                                        filter="STD",
                                        hduclass=["OGIP","IMAGE","PSF","PREDICTED","NET"],
                                        hdudoc="CAL/GEN/92-020")
    psfHDU.header.update("ENERG_LO", 10**PSF_ENERGY_RANGE[0], "Minimum energy in cube (TeV)")
    psfHDU.header.update("ENERG_HI", 10**PSF_ENERGY_RANGE[1], "Maximum energy in cube (TeV)")
    #psfHDU.header.update("SUMRCTS", 1.0, "normalization of PSF")                            
    psfHDU.header.update("CBD10001", "THETA({0:.2f})deg".format(math.sqrt(offY**2+offX**2)),
                         "distance from optical axis")                            
    psfHDU.writeto("psf_cube.fits", clobber=True)


    psf1DHDU = psfHist.asFITS()
    psf1DHDU.writeto("psf_1d.fits", clobber=True)


    # Normalize the photon distribution matrix (the integral along
    # the vertical axis should be 1.0, since it's a probability)
    np.apply_along_axis( normalizeToProb, arr=energyBiasHist.hist, axis=1)
    np.apply_along_axis( normalizeToProb, arr=energyResponseHist.hist, axis=1)

    writeRMF( energyResponseHist )

    energyBiasHist.asFITS().writeto("energy_bias.fits",clobber=True)
    energyResponseHist.asFITS().writeto("energy_response.fits",clobber=True)

    # TODO: calculate statistical errors


    if (opts.plot):
        pylab.figure( figsize = (15.0,11.0), dpi=80 )
        pylab.subplot(2,2,1)

        pylab.semilogy()
        pylab.ylim( 1.0e-1, 1.0e8 )
        pylab.plot( bins[0:-1], Aeff_reco, drawstyle="steps-post",
              label="Reco", color="red" )
        pylab.plot( bins[0:-1], Aeff_true, drawstyle="steps-post",
              label="True", color="blue" )
        pylab.legend(loc="lower right")
        pylab.title("Effective Area")
        pylab.xlabel("$Log_{10}(E/\mathrm{TeV})$")
        pylab.ylabel("$A_{\mathrm{eff}} (\mathrm{m}^2)$")
        pylab.grid()

        pylab.subplot(2,2,2)
        pylab.semilogy()
        for ii in range(psfHist.hist.shape[1]):
            pylab.plot(psfHist.binCenters(0),psfHist.hist[:,ii], 
                       label="E={0:.2f}-{1:.2f}".format(10**psfHist.binLowerEdges[1][ii],
                                                        10**psfHist.binLowerEdges[1][ii+1]))
        pylab.title("PSF")
        pylab.xlabel("$\Theta^2 \mathrm{(deg^2)}$")
        pylab.ylabel("Counts")
        pylab.grid()
        pylab.legend(loc='upper right')
        




        pylab.subplot(2,2,3)
        energyResponseHist.draw2D()
        l = energyResponseHist.binLowerEdges[0]

        pylab.plot( l,l, color="white") 
        pylab.colorbar()
        pylab.title("Energy Response")


        # make this histogram normalized to have an integral of 1.0 along
        # the Y axis (so it is now a probability of reconstructing Ereco
        # for a given Etrue)
        pylab.subplot(2,2,4)
        energyBiasHist.draw2D( vmax=0.25 )
        pylab.colorbar()
        pylab.title("Energy Bias")
        l = energyBiasHist.binLowerEdges[0]
        pylab.plot( l,np.zeros_like(l), color="black")
        pylab.grid(color='w')
        pylab.savefig("response.pdf", papertype="a4")


            

            

        

 
