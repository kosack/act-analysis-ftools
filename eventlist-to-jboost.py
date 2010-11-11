from __future__ import with_statement
import pyfits
import numpy as np
import actutils
import sys
from optparse import OptionParser


"""
Converts an eventlist into an ASCII format that works as a JBoost
training dataset
"""


def generateSpecFile(hdu,cols,outputfile="training.spec"):
    """
    Generate a JBoost spec file from the columns in a fits data table
    
    Arguments:
    - `hdu`: a BINTABLE hdu
    - `cols`: list of columns to use for training
    """
    
    print "Generating Spec File"

    with open(outputfile, mode='w') as out:
        out.write("exampleTerminator=;\n")
        out.write("attributeTerminator=,\n")
        
        for col in cols:
            out.write("{0:20s} number\n".format(col) )


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.set_usage("eventlist-to-jboost.py [options] <eventlist> [<eventlist ..]")
    parser.add_option( "-o","--output", dest="output",  default="training", 
                       help="base output filename")    
    (opts,args) = parser.parse_args()
    if len(args)==0:
        parser.error("specify eventlist filename(s)")

    trainingColumns = ["HIL_MSW","HIL_MSL"]

    first = True

    with open(opts.output+".train",mode="w") as out:

        for ii,eventlist_filename in enumerate(args):

            print "[{0:3d}/{1:3d}] ".format(ii,len(args)), eventlist_filename

            evfile = pyfits.open(eventlist_filename)
            events = evfile["EVENTS"]
        
            if (first):
                generateSpecFile( hdu=events, cols=trainingColumns,
                                  outputfile=opts.output+".spec")
                first=False
                
            data = zip( *[events.data.field(trainingColumns[ii]) 
                          for ii in xrange(len(trainingColumns))] )
            
            for ii in xrange(len(data)):
                out.write(",".join(["{0}".format(val) for val in data[ii]]))
                out.write("\n")
                          

