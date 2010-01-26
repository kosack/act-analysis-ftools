import pyfits
import numpy as np
from optparse import OptionParser
import sys
from numpy import random 
import math

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option( "-n","--num-events", dest="nevt", default=100,
                       help="number of events")
    parser.add_option( "-c","--center", dest="center", default="0.0,90.0",
                       help="center pos")
    parser.add_option( "-s","--sigma", dest="sigma", default="0.1",
                       help="gaussian sigma")
    parser.add_option( "-o","--output", dest="output", help="output file name")
    parser.add_option( "-b", "--num-background", dest="nbg", default=None, 
                       help="number of background events")
    parser.add_option( "-f", "--fov", dest="fov", default=5.0, 
                       help="field-of-view for background events" )
    (opts, args) = parser.parse_args()


    if (opts.center):
        center = np.array(opts.center.split(",")).astype(float)
        if len(center) != 2:
            print "Center should be RA,Dec"
            sys.exit(1)
    
    point = center + np.array((0,0.7))

    nevt = int(opts.nevt)
    sig = float(opts.sigma)

    print "CENTER: ",center
    print "   NUM: ",nevt

    random.seed()

    # signal events:
    X = random.randn(nevt)*sig
    Y = random.randn(nevt)*sig
    decs = Y + center[1]
    ras = X/math.cos(math.radians(center[1])) + center[0] # approx
                                                          # circular
                                                          # in sky
                                                          # degrees
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


    msw = np.zeros(nevt+nbg)
    msl = np.zeros(nevt+nbg)
    
    c1 = pyfits.Column( name='RA', format='D', array = ras )
    c2 = pyfits.Column( name='DEC', format='D', array = decs )
    c3 = pyfits.Column( name='DETX', format='D', array = X )
    c4 = pyfits.Column( name='DETY', format='D', array = Y )
    c5 = pyfits.Column( name='TYPE', format='B', array = types )
    c6 = pyfits.Column( name='HIL_MSW', format='D', array = msw )
    c7 = pyfits.Column( name='HIL_MSL', format='D', array = msl )

    evlist = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7])
    evlist.name = "EVENTS"
    # currently hard-coded offset of (0,0.7)
    evlist.header.update( "RA_PNT", center[0] )
    evlist.header.update( "DEC_PNT", center[1]+0.7 )


    if (opts.output):
        evlist.writeto(opts.output, clobber=True)
