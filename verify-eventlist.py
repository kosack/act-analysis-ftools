import pyfits
import sys

if __name__ == '__main__':
    
    filename = sys.argv.pop(1)

    try:
        events = pyfits.open(filename)["EVENTS"]
    except:
        print "Couldn't EVENTS HDU in ",filename
        sys.exit( 1 )

    columns = ["RA","DEC","HIL_MSW"]
    headers = ["RA_PNT","DEC_PNT","ALT_PNT"]

    # check that events exist and necessary columns are there

    try:
        for col in columns:
            test = events.data.field(col)
    except:
        print "Couldn't find all the necessary columns in EVENTS for ",filename
        sys.exit(1)
        
    nevents = len(test)

    if nevents == 0:
        print "No events in file ",filename
        sys.exit(1)

    print "EVENTLIST:",filename
    print "N_EVENTS: ",nevents

    # check headers
    
    try:
        for key in headers:
            value = events.header[key]
            print key,": ",value
    except:
        print "Couldn't find all necessary headers in ",filename
        sys.exit(1)


