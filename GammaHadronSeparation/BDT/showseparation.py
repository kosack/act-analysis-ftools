from pylab import *

protons = loadtxt("classify_protons.dat")
gammas = loadtxt("classify_gammas.dat");

hrange = [min(protons.min(), gammas.min())-1.0,max(protons.max(), gammas.max())+1.0]
hbins = 20

pbins,edges = np.histogram(protons,bins=hbins, range=hrange, normed=True)
gbins,edges = np.histogram(gammas,bins=hbins, range=hrange, normed=True)


title("Gamma-Hadron sepration power")
plot( edges[0:-1], pbins, drawstyle='steps-pre', color='r', label="Protons",lw=2 )
plot( edges[0:-1], gbins, drawstyle='steps-pre', color='g', label="Gammas",lw=2 )
xlabel("$\zeta$")
ylabel("probability");
legend(loc="upper left")

savefig("separation.pdf")
