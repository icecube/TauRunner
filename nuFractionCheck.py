from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.colors
import numpy as np
import csv

def nuFractionCheck(eini, xvals_beacom, yvals_beacom, xvals_tr, yvals_tr):
    # Open and read csv files containing taurunner data (CHANGE FILE DIRECTORY)
    with open('checks/check2_confirm_' + eini + '.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            xvals_tr.append(float(row[0]))
            yvals_tr.append(float(row[1]))

    # TauRunner plot
    plt.plot(xvals_tr, yvals_tr, color='blue', label=r'$\nu_\tau + \bar{\nu}_\mu}$ TauRunner')
    
    # Open and read csv files containing Beacom data (CHANGE FILE DIRECTORY)
    with open('checks/check2_' + eini + '.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            xvals_beacom.append(float(row[0]))
            yvals_beacom.append(float(row[1]))

    # Empirical plot
    plt.plot(xvals_beacom, yvals_beacom, color='red', label=r'$\nu_\tau + \bar{\nu}_\mu}$ Beacom, Crotty, Kolb')
    bottom, top = plt.ylim()
    plt.ylim(bottom=0)
    plt.ylim(top=1.35)

    # NuTau flux is unity
    x = np.linspace(0, 90, 100)
    y = np.ones(100)
    plt.plot(x, y, linewidth=2, label=r"$\nu_\tau$", color="black")  

# Time to plot
nuFractionCheck('100000', [], [], [], [])
plt.legend()
nuFractionCheck('1000000', [], [], [], [])
nuFractionCheck('1000000000', [], [], [], [])
plt.xlabel('nadir angle (degrees)')
plt.ylabel('emergent flux / initial flux')
plt.show()
