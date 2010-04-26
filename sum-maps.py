import os
import shutil
import numpy as np
from optparse import OptionParser
from math import floor

# sums a list of files in a fast manner. This program uses ftpixcalc
# to sum any number of images, up to 24 simultaneously, minimizing
# disk access.

parser = OptionParser()
parser.set_usage( "sum-maps.py [options] map1 map2 ...")
parser.add_option( "-o","--output", dest="output", 
                   default="sum.fits", 
                   help="output name base (rest of name will be appended)")

parser.add_option( "-v","--verbose", dest="verbose",action="store_true" )

(opts, args) = parser.parse_args()

alpha = np.array(list( 'ABCDEGHIJKLMNOPQRSUVWXYZ' )) # note no t or f
                                                     # since they are
                                                     # interpreted as
                                                     # true/false
bunchsize = len(alpha)
bunches = np.array_split( args, floor(len(args)/bunchsize)+1 )
tmpname = "temp-"+opts.output

if verbose: print args

tot = len(args) # first file is already summed
count = 0

for bunch in bunches:

    if opts.verbose:
        count += len(bunch)
        pct = count/float(tot)*100
        print "[%3d%%] SUMMING %d of %d total images ..." % (pct,len(bunch),tot) 

    expr = "+".join( alpha[0:len(bunch)] )
    vals = ""
    for var,fname in zip(alpha, bunch):
        vals += " %s=%s" % (var,fname)
        
    cmd = "ftpixcalc %s '%s' %s" % (tmpname,expr,vals)
    if  verbose: print "DEBUG: ",cmd
    os.system(cmd)
    os.rename( tmpname, opts.output )
        


