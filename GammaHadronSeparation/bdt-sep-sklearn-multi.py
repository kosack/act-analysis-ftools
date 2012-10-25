from pylab import * 
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn import cross_validation
import pyfits
from glob import glob
import os
from collections import defaultdict

DATADIR=os.path.expanduser("~kosack/Data/FITS/HESS/Simulations/Phase1b")
MAXTREEDEPTH = 10

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

def separation_plot(classfier, X_test, Y_test, bins=10):
    """
    
    Arguments:
    - `classfier`:
    - `X_signal`:
    - `Y_signal`:
    """
    
    probgamma = classifier.predict_proba(X_test[Y_test==0])[:,0] 
    probproton = classifier.predict_proba(X_test[Y_test==1])[:,0] 

    title("Separation Power: " + "+".join(cols))
    hist( probgamma, range=[0,1], bins=bins, color='g', 
          label="gammas", normed=True )
    hist( probproton, range=[0,1], bins=bins, color='r', alpha=0.5, 
          label="protons",normed=True)
    xlabel("Gamma probability")
    legend()
    figtext( 0.12,0.8, str(classifier), fontsize=8)

if __name__ == '__main__':
    
    
    #======================================================================
    # Set up the data:

    cols = ['HIL_MSW','HIL_MSL',"XMAX"]

    gammas = loadData(DATADIR+"/Gamma/0.7deg/run_00006*_eventlist.fits.gz",
                      maxevents=100000,cols=cols)
    protons = loadData(DATADIR+"/Proton/run_*.fits.gz",
                       maxevents=300000,cols=cols)


    gtype = ones_like(gammas[cols[0]])   # gammas are labeled 1
    ptype = zeros_like(protons[cols[0]]) # protons are labeled 0
    Y_all= np.concatenate([gtype,ptype]) # the array of gamma, proton labels

    X_gammas = np.column_stack( gammas[col] for col in gammas )
    X_protons = np.column_stack( protons[col] for col in protons )
    X_all = np.concatenate( [X_gammas,X_protons] )

    # separate into training and test sets

    X_train,X_test,\
        Y_train,Y_test = cross_validation.train_test_split( X_all,Y_all,
                                                            test_size=0.5)


    # Train the classifier:

    classifier = tree.DecisionTreeClassifier(max_depth=MAXTREEDEPTH)
    #classifier = RandomForestClassifier(n_estimators = 10)
    #classifier = GradientBoostingClassifier()

    print classifier
    print "Training using ", cols, "..."


    classifier.fit( X_train, Y_train )
    print "done."
    print "Score = ", classifier.score(X_test,Y_test)

    # Now plot the results of the classification

    if hasattr(classifier,"tree_") :
        outfile = tree.export_graphviz(classifier,"classifier.dot")
        outfile.close()


    figure()
    subplot(1,1,1)

    separation_plot( classifier, X_test, Y_test )

    show()
