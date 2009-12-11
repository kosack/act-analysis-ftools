#
#  Examples of manipulating CTA event lists with Python
#
#  Author: Karl Kosack <karl.kosack@cea.fr>
#

from pylab import *
import pyfits
import math
from matplotlib.patches import Circle

def display_telescope_pattern( evlist ):
    """ example of manipulating the telescope bitpattern """
    clf()
    title("test")
    # telmask is a boolean array:    
    telmask = evlist.data.field("TELMASK")

    # compare the sum of this array with the MULTIP column (they
    # should be the same):

    multip = evlist.data.field("MULTIP")
    ntels = array([ sum(mask) for mask in telmask ])

    error = sum(multip-ntels)
    print "ERROR: ", error
    
    if abs(error)>0:
        print "WARNING: muliplicity doesn't agree with telescope mask!"
    else:
        print "multiplicity and telescope mask agree for all events."

    # now select events with telescope 3 participating:

    totevts = len(telmask)
    print "TOTAL EVENTS:", 

    MAXTEL = evlist.header.get("N_TELS") 
    print "MAXTEL:",MAXTEL

    # ==================================================================
    # Calculate the single telescope participation fractions
    # ==================================================================

    particip_frac = zeros(MAXTEL) 
    for tel in range(MAXTEL): 
        n = sum(telmask[:,tel] == True) # count events with CT<tel>
        print "CT %d EVENTS:" % tel, n
        particip_frac[tel] = n

    particip_frac /= float(totevts)
    subplot(2,2,1)
    bar( arange(MAXTEL)-0.5, particip_frac, width=1 )
    title("1-Tel Participation frac")
    ylabel("Participation Fraction")
    xlabel("Telescope ID")

    # ==================================================================
    # Calculate the 2-telescope participation fractions
    # ==================================================================
    M = zeros( (MAXTEL,MAXTEL) )
    for a in range(MAXTEL):
        for b in range(MAXTEL):
            n = sum(telmask[:,a] & telmask[:,b] == True)
            M[a,b] = n

    M /= totevts

    subplot(2,2,2)
    spectral()
    pcolor(arange(MAXTEL+1)-0.5, arange(MAXTEL+1)-0.5, M)
    colorbar()
    xlabel("Telescope ID")
    ylabel("Telescope ID")
    title("2-Tel Participation frac")

    # ==================================================================
    # Plot the multiplicities
    # ==================================================================
    subplot(2,2,3)
    hist( multip, range=(0,MAXTEL+1), bins = MAXTEL, align='left' )
    xlabel("Multiplicity")
    ylabel("Frequency")

def get_cut_data(events, mswcut=(-2.0,0.7), mslcut=(-2.0,2.0)):
    """returns data after cuts are applied
    
    Arguments:
    - `events`: input fits HDU containing eventlist
    """

    print "Cutting data..."
    msw = events.data.field("HIL_MSW")
    msl = events.data.field("HIL_MSL")
    cutmask = ((msw > mswcut[0]) & (msw < mswcut[1]) 
               & (msl > mslcut[0]) & (msl < mslcut[1]))

    return events.data[cutmask]


def display_scaled_params( events ):

    msw = events.data.field("HIL_MSW")
    msl = events.data.field("HIL_MSW")
    
    cutdata = get_cut_data( events )
    
    print "Events [before cuts]: ",len(events.data)
    print "Events  [after cuts]: ",len(cutdata)

    fields = ["HIL_MSW","HIL_MSL","ENERGY","XMAX"]
    ranges = [(-2,2.4), (-4,4), (0,100), (0,600)]
    ii=1
    for field,range in zip(fields,ranges):
        subplot(2,2,ii)
        vals = events.data.field(field)
        cvals = cutdata.field(field) 
        n,bins,patches = hist( vals, range=range,
                               log=True, histtype='bar', 
                               bins=100,label='all', fc='y', lw=0,aa=False  )
        hist( cvals , log=True, histtype='bar', lw=0, aa=False,
              bins=100,label='cut', fc='g' )
        
        title(field)
        legend()
        ii+=1



def display_events(events):
    """ Make a simple histogram plot of the events in detector coordinates
    
    Arguments:
    - `eventlist`: input fits hdu
    """
    
    cutdata = get_cut_data(events)
    
    # note this is not quite right, since it's plotted at the bin
    # edges, not the centers, so there is a shift, but is is a simple
    # example.  Use SkyImage utility class to do proper plots with
    # arbitrary projections.

    subplot(2,2,1)
    try:
        X = events.data.field("DETX")
        Y = events.data.field("DETY")
    except:
        X = events.data.field("SKYX")
        Y = events.data.field("SKYY")

    (h2d,xe,ye) = histogram2d( X,Y, bins=(50,50), range=[(-3,3),(-3,3)] )
    pcolor( xe,-ye, h2d.transpose() )
    colorbar();
    title("No Cuts")

    plt = subplot(2,2,2)
    try:
        X = cutdata.field("DETX")
        Y = cutdata.field("DETY")
    except:
        X = cutdata.field("SKYX")
        Y = cutdata.field("SKYY")

    (h2d,xe,ye) = histogram2d( X,Y, bins=(50,50), range=[(-3,3),(-3,3)] )
    pcolor( xe,-ye, h2d.transpose() )
    colorbar();
    title("With Cuts")

    # get target position as an offset from the center:
    obs_ra  = events.header['RA_PNT']
    obs_dec = events.header['DEC_PNT']
    targ_ra  = events.header['RA_OBJ']
    targ_dec = events.header['DEC_OBJ']
    onpos = array( [(targ_ra - obs_ra) 
                    * math.cos(radians(obs_dec)),
                    targ_dec - obs_dec] )
    plt.add_artist( Circle( onpos, math.sqrt(0.05), fill=False, fc=None, 
                            ec='w', alpha=0.0, lw=2) )    


def display_rawsim(events):
    """displauy some raw simulation data (no scaled params)
    
    Arguments:
    - `events`:
    """
    Etrue = events.data.field('MC_ENERGY')
    hist( log10(Etrue), range=(-3,2.5), bins=50, log=True )
    xlabel("$\log(E_{true})$");
    

def hasSimulationData(events):
    """ Returns true if sim columns exist in the given data table
    
    Arguments:
    - `events`: HDU 
    """
    
    try:
        mc = events.data.field("MC_EVENTID")
    except:
        return False

    return True


def display_sim(events):
    """ display some simulation plots
    
    Arguments:
    - `events`: event-list HDU
    """
    Ereco = events.data.field('ENERGY')
    Etrue = events.data.field('MC_ENERGY')

    subplot(3,2,1)
    loglog( Etrue,Ereco, linestyle='None' )
    scatter( Etrue,Ereco )
    xlabel("$E_{true}$")
    ylabel("$E_{reco}$")
    
    subplot(3,2,2)
    (h2d,xe,ye) = histogram2d( log10(Etrue),log10(Ereco), 
                               bins=(100,100), range=[(-2,3),(-2,3)])
    extent = [xe[0], xe[-1], ye[0], ye[-1]]
    imshow( h2d[:,::-1].T, extent=extent, interpolation='nearest',
            aspect='auto')
    colorbar()
    
    subplot(3,2,3)
    Eres =  (Ereco-Etrue)/Etrue
    scatter( log(Etrue), Eres, alpha=0.3, edgecolors='none')
    ylim( (-2,2) )
    ylabel("($E_{reco}-E_{true})/E_{true}$")
    xlabel("$\log(E_{true})$")

    subplot(3,2,4)
    (bias,x,y) = histogram2d( log10(Etrue),Eres,
                              bins=(50,100), range=[(-4,5),(-2,2)])
    extent = [x[0], x[-1], y[0], y[-1]]
    imshow( bias[:,::-1].T, extent=extent, interpolation='nearest',
            aspect='auto')

    subplot(3,2,5)
    for ii in range(20,40,5):
        print ii
        plot( y[0:-1], (bias[ii,:]/max(bias[ii,:])) )
        xlim(-2,2)
        ylim(0,1)

def display_array(telarray):


    if (telarray.data == None):
        return

    X = telarray.data.field("POSX")
    Y = telarray.data.field("POSY")
    Z = telarray.data.field("POSZ")
    
    scatter( X,Y )
    title( "Telescope positions" )



def display(filename):
    """ perform all displays on the eventlist given """

    ff = pyfits.open(filename)
    hdu = ff['EVENTS']
    tel = ff['TELARRAY']

    figure()
    display_telescope_pattern( hdu )
    
    figure()
    display_scaled_params( hdu )

    figure()
    display_events(hdu)

    figure()
    display_array( tel )

    figure()
    if (hasSimulationData( hdu) ):
        display_sim( hdu )





def enableSim():
    def get_cut_data(hdu):
        return hdu

def main():

    print "Welcome. Run display(filename) in interactive mode to display "
    print "a fits eventlist"

if __name__ == "__main__":
    main()
