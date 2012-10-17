from pylab import * 
from sklearn import tree
from sklearn import cross_validation
import pyfits
from glob import glob
import os

DATADIR=os.path.expanduser("~kosack/Data/FITS/HESS/Simulations/Phase1b")


def loadData( filepattern, maxevents=50000 ):
    """ returns array of MSW values for given data """

    print filepattern
    msw = np.array( [] )

    for thefile in glob( filepattern ) :
        print "\t reading",thefile
        fits = pyfits.open(thefile)['EVENTS']
        tmp_msw = fits.data.field("HIL_MSW")

        msw = np.concatenate( [msw,tmp_msw] )
        del fits

        if len(msw) > maxevents:
            msw = msw[0:maxevents]
            break

    if (len(msw)) == 0:
        raise ValueError("no data found")

    return msw
    


gmsw = loadData(DATADIR+"/Gamma/0.7deg/run_*_eventlist.fits.gz",maxevents=10000)
pmsw = loadData(DATADIR+"/Proton/run_*.fits.gz",maxevents=30000)


gtype = ones_like(gmsw) # gammas are labeled 1
ptype = zeros_like(pmsw) # protons are labeled 0

x= np.concatenate([gmsw,pmsw])[newaxis].T # the array of gamma, proton values
y= np.concatenate([gtype,ptype]) # the array of gamma, proton labels

print "Training..."
cl = tree.DecisionTreeClassifier(max_depth=4,
                                 min_samples_split=100,
                                 min_samples_leaf=100)
cl.fit(x,y)
print "done."

# Now plot the results of the classification

subplot(2,1,1)

title("Training: {0} gammas, {1} protons".format(len(gmsw),len(pmsw)))
ghist,xe,pat = hist(gmsw, range=[-3,10], bins=100,label="Gammas")
phist,xe,pat = hist(pmsw, range=[-3,10], bins=100, label="Protons")

# generate some test values to plot:

testx = linspace(-3,10,100)[newaxis].T
testy = cl.predict(testx) # predicts the class (1 or 0)

plot( testx,testy*ghist.max(), drawstyle='steps', lw=3,label="Class" )
xlabel("Mean-Reduced-Scaled-Width")
ylabel("Frequency")
legend()

subplot(2,1,2)
plot( testx,cl.predict_proba(testx)[:,1]*100,lw=3, color='g',label="gamma")
plot( testx,cl.predict_proba(testx)[:,0]*100,lw=1,color='r', label="proton")
xlabel("Mean-Reduced-Scaled-Width")
ylabel("Probability")
legend()



outfile = tree.export_graphviz(cl,"test.dot")
outfile.close()

