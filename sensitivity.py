import numpy as np
from pylab import *
import pyfits
from scipy import integrate,interpolate


#
# TODO: integrate in dLogE intead of dE to make it more well-behaved
#


def sourceSpectrum(E):
    return 3.0e-11* E**(-2.5)

def backgroundSpectrum(E):

    # cosmic ray spectrum at 1 TeV is approx 10-3 1/(m^2 sr s GeV)
    r_fov = 5.0
    A_fov = math.pi * r_fov**2 # in square degrees
    FOV_ster = A_fov * (math.pi/180.0)**2 # in steradians
    CM2_per_M2 = 100**2
    GeV_per_TeV = 1000

    N0 = 1.0e-6 *(1.0/CM2_per_M2) * GeV_per_TeV * FOV_ster

    return N0 * E**(-3.0)


if __name__ == '__main__':
    
    # TODO:  make the ARF file conform to OGIP standards (each row
    # has a vector of Aeffs)
    
    plot_rows = 3
    plot_cols = 2 
    curplot = 1       # current plot

    t_exp_hrs = 50.0


    arf_gamma = "gamma_arf.fits"   # true effective area
    arf_proton = "proton_arf.fits" # reco effective area
    
    print "EFFECTIVE AREAS:";
    print "    GAMMA: ", arf_gamma
    print "   PROTON: ", arf_proton
    
    # load the effective areas for protons and gamma-rays and set up
    # interpolating functions for them:
    
    table_gamma = pyfits.open(arf_gamma)["EFF_AREA"]
    table_proton = pyfits.open(arf_proton)["EFF_AREA"]

    Egmin = table_gamma.data.field("ENERG_LO")
    Egmax = table_gamma.data.field("ENERG_HI")
    Eg = (Egmax+Egmin)/2.0  
    Epmin = table_proton.data.field("ENERG_LO")
    Epmax = table_proton.data.field("ENERG_HI")
    Ep = (Epmax+Epmin)/2.0 
    Aeff_gamma = interpolate.interp1d(Eg,table_gamma.data.field("SPECRESP"),
                                      kind='linear',copy=False,
                                      bounds_error=False, fill_value=0)
    Aeff_backg = interpolate.interp1d(Ep,table_proton.data.field("SPECRESP"),
                                      kind='linear',  copy=False,
                                      bounds_error=False, fill_value=0)

    # set up plotting and calculating range
    E = logspace( -2,2,20)
    dE = E[1:]-E[0:-1]
    E = E[0:-1]

    figure( figsize=(10,10) )
    subplot(plot_rows,plot_cols,curplot)
    loglog( E, sourceSpectrum(E), color='b', label="Source" )
    loglog( E, backgroundSpectrum(E), color='r', label="Background" )
    legend()
    title("Intrinsic Spectra")
    ylabel("Flux")
    xlabel("E (TeV)")
    grid()
    curplot+=1


    subplot(plot_rows,plot_cols,curplot)
    loglog( Eg, Aeff_gamma(Eg), color="g", drawstyle="steps-mid" )
    loglog( Ep, Aeff_backg(Ep), color="g", drawstyle="steps-mid" )
    loglog( E, Aeff_gamma(E), color="b", label="Source" )
    loglog( E, Aeff_backg(E), color="r", label="Background" )
    title("Effective Area")
    ylabel("$A_{eff} (cm)$")
    xlabel("E (TeV)")
    grid()
    curplot+=1

    # differential rates: (dN/dt/dE)
    rate_gamma = lambda ee: sourceSpectrum(ee) * Aeff_gamma(ee) 
    rate_backg = lambda ee: backgroundSpectrum(ee) * Aeff_backg(ee) 
    
    subplot(plot_rows,plot_cols,curplot)
    loglog( E, rate_gamma(E),color='b', label="gamma" )
    loglog( E, rate_backg(E), color='r', label="background")
    title("Differential Rate")
    xlabel("E (TeV)")
    ylabel("dN/dt")
    grid()
    curplot+=1

    # integral number of gamma-ray events
    # N = \int( dN/(dt dE) dE ) * t_{exp}
    t_exp_sec = t_exp_hrs *60.0*60.0;
    N_gamma = zeros_like( E )
    N_backg = zeros_like( E )
    intflux = zeros_like( E )

    emax = np.max(E[rate_gamma(E) > 1e-10])
    emax=200
    print "E_max=",emax

    for ii in range( len(E) ):
        print "Integrating: E > {0:.2f}".format(E[ii])
        N_gamma[ii],err = integrate.quadrature( rate_gamma, E[ii], emax )
        N_backg[ii],err = integrate.quadrature( rate_backg, E[ii], emax ) 
        intflux[ii],err = integrate.quadrature( sourceSpectrum, E[ii], emax )

    N_gamma *= t_exp_sec 
    N_backg *= t_exp_sec
    N_exc = N_gamma - N_backg

    subplot(plot_rows,plot_cols,curplot)
    loglog( E, N_gamma,color='b', label="gamma" )
    loglog( E, N_backg, color='r', label="background")
    title("Integral Counts in {time} hrs".format(time=t_exp_hrs))
    xlabel("E (TeV)")
    ylabel("N")
    grid()
    curplot+=1

    # Significance:

    sig = N_gamma/np.sqrt(N_backg)
    sig[np.isfinite(sig)==False] = 0

    subplot(plot_rows,plot_cols,curplot)
    semilogx( E, sig, color='b' )
    xlabel("E (TeV)")
    ylabel("Significance (sigma)")
    curplot+=1

    # calculate minimum flux
    # scale to 5.0 sigma detection:

    sens = intflux.copy()
    minsigma = 5.0
    mask_sig = sig>minsigma
    sens[mask_sig] *= (minsigma/sig[mask_sig])

    # scale to at least 10 excess events:
#    mincounts = 10.0
#    mask_count = N_exc<mincounts
#    sens[mask_count] *= (mincounts/(N_exc[mask_count]+1e-20))
#    loglog( E, sens,color="g" )
#    ylim(1e-16, 1e-10)

    senscrab = sens/intflux

    subplot(plot_rows,plot_cols,curplot)
    loglog( E, senscrab, color="b" )
    ylabel("Minimum integral flux (C.U.)")
    xlabel("E (TeV)")
    title("$\mathrm{Sensitivity} (%.1f\sigma,%.1f \mathrm{hrs})$"%(minsigma,t_exp_hrs))
    grid()
    
