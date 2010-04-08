import pyfits
import numpy as np
from optparse import OptionParser
import sys
from numpy import random 
import math


# This is just a simple script to generate fake event lists for
# testing algorithms. Eventually it should be replaced with a fancier
# toy model that samples from physical parameter distributions.


def makePowerLaw(num, gamma=2.7,E0=1.0, Erange=[0.01,200.0]):
    """
    
    Arguments:
    - `num`:
    - `gamma`:
    - `E0`:
    """

    maxdist = ((Erange[0]/E0)**(-1.0/gamma))
    mindist = ((Erange[1]/E0)**(-1.0/gamma))
    unif = random.uniform(low=mindist,high=maxdist, size=num)
    E = unif**(-gamma)
    return E


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option( "-n","--num-events", dest="nevt", default=100, type="int",
                       help="number of events")
    parser.add_option( "-c","--center", dest="center", default="0.0,90.0",
                       help="center pos")
    parser.add_option( "-d","--offset", dest="offset", default="0.0,0.7",
                       help="center pos")
    parser.add_option( "-s","--sigma", dest="sigma", default="0.1",type="float",
                       help="gaussian sigma")
    parser.add_option( "-o","--output", dest="output", help="output file name")
    parser.add_option( "-b", "--num-background", dest="nbg", default=None, type="int",
                       help="number of background events")

    parser.add_option( "--index", dest="gamma_src", default=2.5, type="float",
                       help="spectral index of the source")
    parser.add_option( "--bgindex", dest="gamma_bg", default=3.0, type="float",
                       help="spectral index of the background")

    parser.add_option( "-f", "--fov", dest="fov", default=5.0, type="float",
                       help="field-of-view for background events" )
    (opts, args) = parser.parse_args()

    center = np.array(opts.center.split(",")).astype(float)
    if len(center) != 2:
        print "Center should be RA,Dec in deg"
        sys.exit(1)
    
    offset = np.array(opts.offset.split(",")).astype(float)
    if len(offset) != 2:
        print "offset should be RA,Dec in deg"
        sys.exit(1)

    point = center + offset

    nevt = int(opts.nevt)
    sig = float(opts.sigma)

    for opt,val in eval(opts.__str__()).items():
        print "%20s = %s" %(opt.upper(),val)

    random.seed()

    # signal events:
    X = random.randn(nevt)*sig
    Y = random.randn(nevt)*sig
    decs = Y + center[1]
    ras = X/math.cos(math.radians(center[1])) + center[0] # approx
                                                          # circular
                                                          # in sky
                                                          # degrees
    gam = float(opts.gamma_src)
    E = makePowerLaw( nevt, gamma=float(opts.gamma_src ))

    types = np.zeros(nevt)
    nbg = int(opts.nbg)

    # background events (if requested):
    if (opts.nbg):
        fov = float(opts.fov)/2.0
        Xb = random.randn(nbg)*fov
        Yb = random.randn(nbg)*fov
        decsb = Yb + point[1]
        rasb =  Xb/math.cos(math.radians(point[1])) + point[0]
        X = np.array(X.tolist() + Xb.tolist())
        Y = np.array(Y.tolist() + Yb.tolist())
        ras = np.array(ras.tolist() + rasb.tolist())
        decs = np.array(decs.tolist() + decsb.tolist())
        types = np.array(types.tolist() + np.ones(nbg).tolist()).astype(bool)

        Ebg = makePowerLaw( nbg, gamma=float(opts.gamma_bg ))
        E = np.append(E,Ebg)

    msw = np.zeros(nevt+nbg)
    msl = np.zeros(nevt+nbg)
    
    c1 = pyfits.Column( name='RA', format='D', array = ras )
    c2 = pyfits.Column( name='DEC', format='D', array = decs )
    c3 = pyfits.Column( name='DETX', format='D', array = X )
    c4 = pyfits.Column( name='DETY', format='D', array = Y )
    c5 = pyfits.Column( name='TYPE', format='B', array = types )
    c6 = pyfits.Column( name='HIL_MSW', format='D', array = msw )
    c7 = pyfits.Column( name='HIL_MSL', format='D', array = msl )
    c8 = pyfits.Column( name='ENERGY', format='D', array = E )

    evlist = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8])
    evlist.name = "EVENTS"
    evlist.header.update( "RA_PNT", point[0] )
    evlist.header.update( "DEC_PNT", point[1] )


    if (opts.output):
        evlist.writeto(opts.output, clobber=True)
