#
# Convert HESS region file to DS9/fits region format
#
# TODO: for exclusion may need to make ellipse and divide by cos(dec)?
# TODO: use WCS transformation to get the right system (gal vs fk5)

from __future__ import with_statement
import re
from optparse import OptionParser
from kapteyn import wcs

if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option( "-t","--type", dest="type", default="ds9",
                       help="output type: fitsex (for fits exclusion file), "
                       +"ds9 for ds9 region file")

    parser.add_option( "-s","--scale", dest="scale", default=1.0,
                       help="scale factor for region radius (default 1.0)")
    
    (opts, args) = parser.parse_args()        

    if ( len(args) < 1):
        print "hess2ds9.py <HESS excluion region file>"
        exit(1)

    opts.scale = float(opts.scale)
    filename=args.pop(0)

    with open(filename) as file:

        gal2radec = wcs.Transformation( wcs.galactic, wcs.fk5 )

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
                # convert back to fk5
                (glam, gbet) = gal2radec( (lam,bet) )
                lam = glam
                bet = gbet
                sys="FK5"
            

            if (shape == "SEGMENT"):
                extra = "unknown"
                try:
                    r0 = float(tokens.pop(0))*opts.scale
                    r1 = float(tokens.pop(0))*opts.scale
                    phi0 = float(tokens.pop(0))
                    phi1 = float(tokens.pop(0))
                    extra = tokens.pop(0)
                except:
                    pass

                if opts.type == 'ds9':
                    print ("%s;circle(%f,%f,%f) # text={%s} move=0" 
                          % (sys, lam, bet, r1,name))
                elif opts.type=='fitsex':
                    print "-circle(%f,%f,%f)" % (lam, bet, r1)
                else:
                    raise Exception("Unknown type: "+opts.type)

            elif (shape == "BOX"):
                width = float(tokens.pop(0))
                height = float(tokens.pop(0))
                
                if opts.type == 'ds9':
                    print ("%s;Box(%f,%f,%f,%f,0.0) # text={%s} move=0" 
                           % (sys, lam, bet,width,height,name))
                elif opts.type=='fitsex':
                    print ("-box(%f,%f,%f,%f,0.0)" 
                           % (lam, bet, width,height))
                else:
                    raise Exception("Unknown type: "+opts.type)
