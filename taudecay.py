import pythia8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
pythia = pythia8.Pythia()

pythia.readString("ProcessLevel:all = off")
pythia.readString("Random:setSeed = on")
pythia.readString("Random:seed = 0")
pythia.init()

UnstableParents = {15:"tau-",-15:"tau+"}

for parent in UnstableParents:
    # make parents unstable
    pythia.particleData.mayDecay(parent,True)

StableDaughters = {-12:"antinue", 12:"nue",
                   -14:"antinumu", 14:"numu",
                   -16:"antinutau", 16:"nutau",
                   -211:"pi-", 321:"k-"}

for daughter in StableDaughters:
    # make daughters stable
    pythia.particleData.mayDecay(daughter,False)

number_of_decays = 1000000000

# initial energy in GeV
momenta = [1.e6, 1.e7, 1.e8, 1.e9]

# dicts containing outgoing energies of nus
results = [{}, {}, {}, {}]

# histogram bins
bins = list(np.logspace(-3,0,101))

# auxiliary list
counts = [[],[],[],[]]

# weights for histogram
w = [[],[],[],[]]

for i in range(len(momenta)):
    for parent in UnstableParents:
        mass = pythia.particleData.m0(parent)
        energy = np.sqrt(mass*mass + momenta[i]*momenta[i])
        p4vec = pythia8.Vec4(momenta[i],0.,0.,energy)

        results[i][parent] = []

        for j in range(number_of_decays):
            # clean previous stuff
            pythia.event.reset()
            # pdgid, status, col, acol, p, m
            # status = 91 : make normal decay products
            # col, acol : color and anticolor indices
            pythia.event.append(parent,91,0,0,p4vec,mass)
            # go all the way
            pythia.forceHadronLevel()
            for k in range(pythia.event.size()):
                if (not pythia.event[k].isFinal()):
                    # if not done decaying, continue
                    continue
                if (pythia.event[k].id() == 16): # is desired neutrino
                    results[i][parent].append(pythia.event[k].e()/energy)

    # plot histograms
    counts[i], bins = np.histogram(results[i][15], bins=bins)
    for l in list(counts[i]): w[i].append(l/number_of_decays)
    plt.hist(bins[:-1], bins, weights=w[i], histtype='step', fill = False, label=r'$10^' + str(i+6) + '$ GeV')

plt.xscale('log')
plt.yscale('log')
bottom, top = plt.ylim()
plt.ylim(bottom=1e-5)
plt.ylim(top=1)

plt.xlabel(r"$E_{\nu_\tau}/E_\tau$ [GeV]", fontsize = 18)
plt.ylabel(r"$N_{\nu_\tau}/N_\tau$", fontsize = 13)

plt.legend(loc='best')
plt.show()
