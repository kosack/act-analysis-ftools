import pyfits
import sys

if __name__ == '__main__':
    
    filename = sys.argv.pop(1)
    outputname = sys.argv.pop(1)

    try:
        events = pyfits.open(filename)["EVENTS"]
    except:
        print "Couldn't EVENTS HDU in ",filename
        sys.exit( 1 )

    columns = ["RA","DEC","HIL_MSW"]
    headers = ["OBJECT","RUN_ID","RA_PNT","DEC_PNT","ALT_PNT", "AZ_PNT"]

    # check that events exist and necessary columns are there

    try:
        for col in columns:
            test = events.data.field(col)
    except:
        print "Couldn't find all the necessary columns in EVENTS for ",filename
        sys.exit(1)
        
    nevents = events.size()

    if nevents == 0:
        print "No events in file ",filename
        sys.exit(1)

    # check headers
    
    try:
        for key in headers:
            value = events.header[key]
    except:
        print "Couldn't find all necessary headers in ",filename
        sys.exit(1)

    # if we got this far, things are ok so write the output file:

    outf = open(outputname, 'w')
    outf.write("EVENTLIST: %s\n" % filename)
    outf.write("N_EVENTS: %d\n" % nevents)
    for key in headers:
        value = events.header[key]
        outf.write("%s: %s\n" % (key,value) )
        
    outf.close()
