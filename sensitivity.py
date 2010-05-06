import numpy as np
from pylab import *
import pyfits
from scipy import integrate


def sourceSpectrum(E):
    return 3.0e-11* E**(-2.5)

def backgroundSpectrum(E):

    # cosmic ray spectrum at 1 TeV is approx 10-3 1/(m^2 sr s GeV)
    r_fov = 5.0
    A_fov = math.pi * r_fov**2 # in square degrees
    FOV_ster = A_fov * (math.pi/180.0)**2 # in steradians
    CM2_per_M2 = 100**2
    GeV_per_TeV = 1000

    N0 = 1.0e-5 *(1.0/CM2_per_M2) * GeV_per_TeV * FOV_ster

    return N0 * E**(-2.7)


if __name__ == '__main__':
    
    # TODO: 
    #  - don't use Aeff histograms directly, interpolate (so you can use any set of E's)
    #  - make the ARF file conform to OGIP standards (each row has a vector of Aeffs)
    
    pn = 3
    pm = 2
    pi = 1

    arf_gamma = "gamma_arf.fits"
    arf_proton = "proton_arf.fits"
    
    print "EFFECTIVE AREAS:";
    print "    GAMMA: ", arf_gamma
    print "   PROTON: ", arf_proton
    
    # load the effective areas for protons and gamma-rays:
    
    table_reco_gamma = pyfits.open(arf_gamma)["EFF_AREA"]
    table_reco_proton = pyfits.open(arf_proton)["EFF_AREA"]

    Emin = table_reco_gamma.data.field("ENERG_LO")
    Emax = table_reco_gamma.data.field("ENERG_HI")
    E = (Emax+Emin)/2.0 # TODO: correctly weight center for log scale
    dE = Emax-Emin
    Aeff_gamma = table_reco_gamma.data.field("SPECRESP")
    Aeff_backg    = table_reco_proton.data.field("SPECRESP")

    figure( figsize=(10,10) )
    subplot(pn,pm,pi)
    loglog( E, sourceSpectrum(E), color='b', label="Source" )
    loglog( E, backgroundSpectrum(E), color='r', label="Background" )
    legend()
    grid()
    pi+=1


    subplot(pn,pm,pi)
    loglog( E, Aeff_gamma, color="b", label="Source" )
    loglog( E, Aeff_backg, color="r", label="Background" )
    title("Effective Area")
    grid()
    pi+=1

    # differential rates: (dN/dt)
    rate_gamma = sourceSpectrum(E) * Aeff_gamma * dE
    rate_backg = backgroundSpectrum(E) * Aeff_backg * dE
    
    subplot(pn,pm,pi)
    loglog( E, rate_gamma,color='b', label="gamma" )
    loglog( E, rate_backg, color='r', label="background")
    title("Differential Rate")
    xlabel("E (TeV)")
    ylabel("dN/dt")
    legend()
    grid()
    pi+=1

    # integral number of gamma-ray events
    # TODO: do with real interpolation/integration (integrate.quadriture)
    t_exp_hrs = 50.0
    t_exp_sec = t_exp_hrs *60.0*60.0;
    N_gamma = zeros_like( E )
    N_backg = zeros_like( E )
    intflux = zeros_like( E )
    for ii in range( len(E) ):
        N_gamma[ii] = np.sum(rate_gamma[ii:])
        N_backg[ii] = np.sum(rate_backg[ii:])
        intflux[ii],err = integrate.quad( sourceSpectrum, E[ii], 1000.0 )

    N_gamma *= t_exp_sec
    N_backg *= t_exp_sec
    N_exc = N_gamma - N_backg

    subplot(pn,pm,pi)
    loglog( E, N_gamma,color='b', label="gamma" )
    loglog( E, N_backg, color='r', label="background")
    title("Integral Counts in {time} hrs".format(time=t_exp_hrs))
    xlabel("E (TeV)")
    ylabel("N")
    legend()
    grid()
    pi+=1

    # Significance:
    subplot(pn,pm,pi)
    sig = N_gamma/np.sqrt(N_backg)
    sig[np.isfinite(sig)==False] = 0
    semilogx( E, sig, color='b' )
    xlabel("E (TeV)")
    ylabel("Significance")
    pi+=1


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
    grid()


