#
# Convert HESS region file to DS9/fits region format
#

from __future__ import with_statement
import re
from optparse import OptionParser


if __name__ == "__main__":

    parser = OptionParser()
    (opts, args) = parser.parse_args()
    
    if ( len(args) < 1):
        print "hess2ds9.py <HESS excluion region file>"
        exit(1)

    filename=args.pop(0)

    with open(filename) as file:

        for line in file.readlines():
            line = line.strip("\n")
            line = re.sub( "#.*","" , line) # remove comments
            line = re.sub( "[\t ]+", "|", line )
            tokens = line.split("|")

            if (len(tokens)<6):
                continue 

            shape = tokens.pop(0)
            type = tokens.pop(0)
            sys = tokens.pop(0)
            lam = float(tokens.pop(0))
            bet = float(tokens.pop(0))
            name = tokens.pop(0)

            if (sys=="GAL"): 
                sys = "GALACTIC"
            elif (sys == "RADEC"): 
                sys="FK5"

            if (shape == "SEGMENT"):
                extra = "unknown"
                try:
                    r0 = float(tokens.pop(0))
                    r1 = float(tokens.pop(0))
                    phi0 = float(tokens.pop(0))
                    phi1 = float(tokens.pop(0))
                    extra = tokens.pop(0)
                except:
                    pass

                print "%s;circle(%f,%f,%f) # text={%s} move=0" % (sys, lam, bet, r1,name)


