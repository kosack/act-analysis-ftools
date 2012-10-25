from pylab import * 
from sklearn import tree
from sklearn import cross_validation
import pyfits
from glob import glob
import os
from collections import defaultdict

DATADIR=os.path.expanduser("~kosack/Data/FITS/HESS/Simulations/Phase1b")
MAXTREEDEPTH = 6

def loadData( filepattern, cols=['HIL_MSW','HIL_MSL'], maxevents=50000 ):
    """ 
    returns array of MSW values for given data

    - filepattern: shell pattern for file(s) to open (e.g. proton-*.fits)
      . If it matches more than one file, all files will be opened and
      the result concatinated.

    - maxevents: keep loading data from files in the filepattern until
      this number of events is reached (discarding everything
      afterward)

      TODO: change to include columns=['HIL_MSW','HIL_MSL','ENERGY'],
      which returns not just the msw as now, but the appropriate
      column stack of the requested values

    """

    print filepattern
    data = defaultdict(list)
    count =0

    for thefile in glob( filepattern ) :
        print "\t reading",thefile
        fits = pyfits.open(thefile)['EVENTS']

        for column in cols:
            data[column].extend( fits.data.field(column) )
            
        del fits

        if len(data[cols[0]]) > maxevents:
            print "Stopping after",maxevents,"events"
            for column in cols:
                data[column] = np.array(data[column])[0:maxevents]
            break

    if (len(data[cols[0]])) == 0:
        raise ValueError("no data found")

    return data
    
#======================================================================
# Set up the data:

cols = ['HIL_MSW','HIL_MSL']

gammas = loadData(DATADIR+"/Gamma/0.7deg/run_00006*_eventlist.fits.gz",
                  maxevents=10000,cols=cols)
protons = loadData(DATADIR+"/Proton/run_000104*.fits.gz",
                   maxevents=30000,cols=cols)

testgammas = loadData(DATADIR+"/Gamma/0.7deg/run_00007*_eventlist.fits.gz",
                      maxevents=10000,cols=cols)
testprotons = loadData(DATADIR+"/Proton/run_00010[6789]*.fits.gz",
                       maxevents=30000,cols=cols)

gtype = ones_like(gammas[cols[0]])   # gammas are labeled 1
ptype = zeros_like(protons[cols[0]]) # protons are labeled 0

X_train_gammas = np.column_stack( gammas[col] for col in gammas )
X_train_protons = np.column_stack( protons[col] for col in protons )
X_train = np.concatenate( [X_train_gammas,X_train_protons] )

X_test_gammas = np.column_stack( testgammas[col] for col in testgammas )
X_test_protons = np.column_stack( testprotons[col] for col in testprotons )
X_test = np.concatenate( [X_test_gammas,X_test_protons] )

Y_train= np.concatenate([gtype,ptype]) # the array of gamma, proton labels

# Train the classifier:

print "Training using ", cols, "..."
classifier = tree.DecisionTreeClassifier(max_depth=MAXTREEDEPTH)
classifier.fit( X_train, Y_train )
print "done."

# Now plot the results of the classification

outfile = tree.export_graphviz(classifier,"classifier.dot")
outfile.close()

show()
