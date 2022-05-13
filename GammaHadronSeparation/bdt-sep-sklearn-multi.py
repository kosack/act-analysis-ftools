from pylab import * 
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn import cross_validation
from sklearn.linear_model import Perceptron
import pyfits
from glob import glob
import os
from collections import defaultdict
from cPickle import dump,load

DATADIR=os.path.expanduser("~kosack/Data/FITS/HESS/Simulations/Phase1b")
MAXTREEDEPTH = 10

def loadData( filepattern, cols=['HIL_MSW','HIL_MSL'], maxevents=50000,
              verbose=False):
    """ 
    returns array of MSW values for given data

    - filepattern: shell pattern for file(s) to open
      (e.g. "proton-*.fits") . If it matches more than one file, files
      will be opened in sequence and the data will be accumulated from
      each.

    - maxevents: keep loading data from files in the filepattern until
      this number of events is reached (discarding everything
      afterward)
    """

    data = defaultdict(list)
    count =0

    print "Reading",filepattern,"..."

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
    Plot a the probability distributions for separation between gamma
    and hadron

    Arguments:
    - `classfier`: the classifier object used
    - `X_test`: the test data set 
    - `Y_test`: the test label set (array with actual event type)
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
    figtext( 0.2,0.8, str(classifier), fontsize=8, ha='left', va='top')
    suptitle("score={0}".format(classifier.score(X_test,Y_test)))




if __name__ == '__main__':
       
    #======================================================================
    # Set up the data:

    cols = ['HIL_MSW','HIL_MSL','XMAX','XMAX_ERR']

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


    # Construct several different types of classifier, loop over each
    # one, train it, and print out the results.  If the classifier
    # already has been trained and exists on disk (as a pickled
    # archive), load it instead of training a new one.

    classifiers = [tree.DecisionTreeClassifier(max_depth=MAXTREEDEPTH),
                   RandomForestClassifier(n_estimators = 10),
                   GradientBoostingClassifier()]
    
    for classifier in classifiers:

        name = str(classifier)
        name = name[0:name.find('(')]

        try: 
            # load up a saved classifier if it exists
            infile = open("{0}-{1}.pickle".format("-".join(cols),name),"r")
            print "Loading",infile,"..."
            classifier = load( infile ) 
            infile.close()
            
        except:
            # otherwise, train a new one
            print "Training",name,"using ", cols, "..."
            classifier.fit( X_train, Y_train )
            print "done."


        print "Score = ", classifier.score(X_test,Y_test)
        print classifier

        # Now plot the results of the classification

        if hasattr(classifier,"tree_") :
            outfile = tree.export_graphviz(classifier,"classifier.dot")
            outfile.close()


        figure()
        subplot(1,1,1)
        separation_plot( classifier, X_test, Y_test, bins=30 )
        savefig("classifier-{0}-{1}.pdf".format("-".join(cols),name))

        outfile =open("{0}-{1}.pickle".format("-".join(cols),name),"w")
        dump(classifier,outfile)
        outfile.close()

    show()
