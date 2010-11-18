from pylab import *

protons = loadtxt("classify_protons.dat")
gammas = loadtxt("classify_gammas.dat");

hrange = [min(protons.min(), gammas.min())-1.0,max(protons.max(), gammas.max())+1.0]
hbins = 30

pbins,edges = np.histogram(protons,bins=hbins, range=hrange, normed=True)
gbins,edges = np.histogram(gammas,bins=hbins, range=hrange, normed=True)

pbins /= pbins.sum()
gbins /= gbins.sum()

subplot(2,1,1)
title("Gamma-Hadron sepration power")
plot( edges[0:-1], pbins, drawstyle='steps-pre', color='r', label="Protons",lw=2 )
plot( edges[0:-1], gbins, drawstyle='steps-pre', color='g', label="Gammas",lw=2 )
xlabel("$\zeta$")
ylabel("probability");
legend(loc="upper left")





eff_g = np.zeros_like(pbins)
eff_p = np.zeros_like(pbins)

ii = arange(len(pbins))

for bin in xrange(len(pbins)):
    eff_g[bin] = gbins[bin:].sum()
    eff_p[bin] = pbins[bin:].sum()

qfact = eff_g/sqrt(eff_p)


subplot(2,1,2)
plot( edges[0:-1], qfact, drawstyle="steps-pre" )
xlabel("$\zeta$")
ylabel("Q")
savefig("separation.pdf")
