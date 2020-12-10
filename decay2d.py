import traceback
import pythia8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import rc

try:
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

    number_of_decays = 10000

    # eini in GeV
    momenta = 1.e9

    # dicts containing outgoing energies of nus
    resultsx = {}
    resultsy = {}

    # histogram bins CHANGE
    bins = list(np.logspace(-3,0,60))

    # auxiliary list
    counts = [[],[],[],[]]

    # weights for histogram
    w = [[],[],[],[]]

    for parent in UnstableParents:
        mass = pythia.particleData.m0(parent)
        energy = np.sqrt(mass*mass + momenta*momenta)
        p4vec = pythia8.Vec4(momenta,0.,0.,energy)

        resultsx[parent] = []
        resultsy[parent] = []

        for i in range(number_of_decays):
            # clean previous stuff
            pythia.event.reset()
            # pdgid, status, col, acol, p, m
            # status = 91 : make normal decay products
            # col, acol : color and anticolor indices
            pythia.event.append(parent,91,0,0,p4vec,mass)
            # go all the way
            pythia.forceHadronLevel()
            for j in range(pythia.event.size()):
                if (not pythia.event[j].isFinal()):
                    # if not done decaying, continue
                    continue
                if (pythia.event[j].id() == 16): # is nutau
                    for sister in pythia.event[j].sisterList():
                        if(pythia.event[sister].id() == -14): # is desired sister
                            resultsx[parent].append(pythia.event[j].e()/energy)
                            resultsy[parent].append(pythia.event[sister].e()/energy)

    # plot histogramplt.figure(figsize=(6,5))
    H, xedges, yedges, fig = plt.hist2d(np.array(resultsx[15]),np.array(resultsy[15]),bins=[np.arange(0.,1.001,0.01),np.arange(0.,1.001,0.01)])

    # # Normalize H
    for i in range(len(H)):
        for j in range(len(H[i])):
            H[i][j] = H[i][j] / number_of_decays

    # cons of e function
    x = np.linspace(0, 1, 100)
    y = 1 - x

    plt.imshow(H, norm=matplotlib.colors.LogNorm())
    plt.plot(x, y, linewidth=2, label=r"$E_{\bar{\nu}_\mu} + E_{\nu_\tau} < E_\tau$", color="red")
    # plt.xscale('log')
    # plt.yscale('log')
    plt.legend(fontsize=14)
    plt.colorbar()
    plt.xlabel(r"$E_{\nu_\tau}/E_\tau$", fontsize = 14)
    plt.ylabel(r"$E_{\bar{\nu}_\mu}/E_\tau$", fontsize = 14)
    plt.show()

except:
    traceback.print_exc()