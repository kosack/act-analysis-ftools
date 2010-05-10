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
    
    pn = 3; pm=2 # plot grid MxN
    pi = 1       # current plot

    arf_gamma = "gamma_arf.fits"
    arf_proton = "proton_arf.fits"
    
    print "EFFECTIVE AREAS:";
    print "    GAMMA: ", arf_gamma
    print "   PROTON: ", arf_proton
    
    # load the effective areas for protons and gamma-rays and set up
    # interpolatin functions for them:
    
    table_gamma = pyfits.open(arf_gamma)["EFF_AREA"]
    table_proton = pyfits.open(arf_proton)["EFF_AREA"]

    Egmin = table_gamma.data.field("ENERG_LO")
    Egmax = table_gamma.data.field("ENERG_HI")
    Eg = (Egmax+Egmin)/2.0 # TODO: correctly weight center for log scale
    Epmin = table_proton.data.field("ENERG_LO")
    Epmax = table_proton.data.field("ENERG_HI")
    Ep = (Epmax+Epmin)/2.0 # TODO: correctly weight center for log scale
    Aeff_gamma = interpolate.interp1d(Eg,table_gamma.data.field("SPECRESP"),
                                      bounds_error=False, fill_value=0)
    Aeff_backg = interpolate.interp1d(Ep,table_proton.data.field("SPECRESP"),
                                      bounds_error=False, fill_value=0)

    # set up plotting and calculating range
    E = logspace( -2,2,40)
    dE = E[1:]-E[0:-1]
    E = E[0:-1]

    figure( figsize=(10,10) )
    subplot(pn,pm,pi)
    loglog( E, sourceSpectrum(E), color='b', label="Source" )
    loglog( E, backgroundSpectrum(E), color='r', label="Background" )
    legend()
    title("Intrinsic Spectra")
    ylabel("Flux")
    xlabel("E (TeV)")
    grid()
    pi+=1


    subplot(pn,pm,pi)
    loglog( Eg, Aeff_gamma(Eg), color="g", drawstyle="steps-mid" )
    loglog( Ep, Aeff_backg(Ep), color="g", drawstyle="steps-mid" )
    loglog( E, Aeff_gamma(E), color="b", label="Source" )
    loglog( E, Aeff_backg(E), color="r", label="Background" )
    title("Effective Area")
    ylabel("$A_{eff} (m)$")
    xlabel("E (TeV)")
    grid()
    pi+=1

    # differential rates: (dN/dt)
    rate_gamma = lambda Ex: sourceSpectrum(Ex) * Aeff_gamma(Ex) * Ex
    rate_backg = lambda Ex: backgroundSpectrum(Ex) * Aeff_backg(Ex) * Ex
    
    subplot(pn,pm,pi)
    loglog( E, rate_gamma(E),color='b', label="gamma" )
    loglog( E, rate_backg(E), color='r', label="background")
    title("Differential Rate")
    xlabel("E (TeV)")
    ylabel("dN/dt")
    grid()
    pi+=1

    # integral number of gamma-ray events
    t_exp_hrs = 50.0
    t_exp_sec = t_exp_hrs *60.0*60.0;
    N_gamma = zeros_like( E )
    N_backg = zeros_like( E )
    intflux = zeros_like( E )
    for ii in range( len(E) ):
        print "Integrating: E > {0:.2f}".format(E[ii])
        N_gamma[ii],err = integrate.quadrature( rate_gamma, E[ii], 200.0 )
        N_backg[ii],err = integrate.quadrature( rate_backg, E[ii], 200.0 ) 
        intflux[ii],err = integrate.quadrature( sourceSpectrum, E[ii], 200.0 )

    N_gamma *= t_exp_sec
    N_backg *= t_exp_sec
    N_exc = N_gamma - N_backg

    subplot(pn,pm,pi)
    loglog( E, N_gamma,color='b', label="gamma" )
    loglog( E, N_backg, color='r', label="background")
    title("Integral Counts in {time} hrs".format(time=t_exp_hrs))
    xlabel("E (TeV)")
    ylabel("N")
    grid()
    pi+=1

    # Significance:

    sig = N_gamma/np.sqrt(N_backg)
    sig[np.isfinite(sig)==False] = 0


    subplot(pn,pm,pi)
    semilogx( E, sig, color='b' )
    xlabel("E (TeV)")
    ylabel("Significance (sigma)")
    pi+=1

    # calculate minimum flux
    subplot(pn,pm,pi)
    loglog( E, intflux,color="r", label="Integral Flux" )
    sens = intflux.copy()

    # scale to at least 10 excess events:
    mask_count = N_exc>10.0
    sens[mask_count] *= (10.0/N_exc[mask_count])
    loglog( E, sens,color="g" )

    # scale to 5.0 sigma detection:
    mask_sig = sig>5.0
    mask = mask_sig
    sens[mask] *= (5.0/sig[mask])
    loglog( E, sens,color="b" )
    ylabel("Min Flux (dN/dE/dA/dt)")
    grid()


