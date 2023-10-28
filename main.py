"""
name,z,n,halflife(Seconds),d-halflife,alpha(keV),d-alpha
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize, stats

plt.rcParams['font.size'] = 12
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

from helper import *

data = list()
equal_z = dict()

## read data
with open('output.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        z, n, hl, hl_e  = float(row[0]), float(row[1]), float(row[2]), float(row[3])
        q, q_e, binding, binding_e = float(row[4]), float(row[5]), float(row[6]), float(row[7])

        nuclei = Nuclei(z, n, hl, hl_e, q, q_e, binding, binding_e)
        data.append(nuclei)

## separete data into equal z lists
for nuclei in data:
    try:
        equal_z[nuclei.z].append(nuclei)
    except KeyError:
        equal_z[nuclei.z] = [nuclei]

## sort equal z data based on z
equal_z = dict(sorted(equal_z.items()))

## sort equal z data based on hl
for z in equal_z:
    equal_z[z].sort()

## remove from equal z lists with only one element
equal_z_temp = equal_z.copy()
for z in equal_z_temp:
    if len(equal_z_temp[z]) == 1:
        equal_z.pop(z)

############ geiger law ############
if 1:
    cmap = plt.get_cmap('turbo', 12)
    for i, z in enumerate(equal_z):
        q_list = list()
        hl_list = list()

        for nuclei in equal_z[z]:
            q_list.append(nuclei.q)
            hl_list.append(nuclei.hl)
        
        q_list = np.asarray(q_list)**(-0.5)
        hl_list = np.asarray(hl_list)

        plt.semilogy(q_list, hl_list, label=f"Z={int(z)}",c=cmap(i), marker='.')


    plt.ylabel('Half-life (Seconds)')
    plt.xlabel(r'$\sqrt{\mathrm{Q (keVs)}}$')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid()
    plt.tight_layout()
    plt.savefig('geiger.png', dpi=300)
    plt.show()

############ root finding ############
if 1:
    x = []
    y = []
    for nuclei in data:
        hl = nuclei.hl
        Q = nuclei.q
        Z = nuclei.z - 2
        A = nuclei.z + nuclei.n - 4
    
        fitter = Fit(35000, hl, Q, Z, A)
        sol = fitter.find_root()
        
        x.append(A+4)
        y.append(sol.x[0])
    
    ## plot roots
    x = np.asarray(x)
    y = np.asarray(y)
    plt.plot(x, y*1e15, '+', markersize=5, c='b')

    ## fit func
    def fit_func(x, r0):
        return r0 * x**(1/3)

    popt, pcov = optimize.curve_fit(fit_func, x, y)
    print(f'a = {popt[0]} +- {np.sqrt(pcov[0,0])}')
    x_fit = np.linspace(0, 300, 300)
    y_fit = fit_func(x_fit, *popt)
    plt.plot(x_fit, y_fit*1e15, c='k')

    estimated = fit_func(x, *popt)
    estimated *= np.sum(y) / np.sum(estimated)
    chi2, p_value = stats.chisquare(y, estimated, 1)
    print(f'chi2 = {chi2}')
    print(f'p_value = {p_value}')

    plt.xlabel('A')
    plt.ylabel('R (fm)')
    plt.ylim(0, 10)
    plt.grid()
    plt.tight_layout()
    plt.savefig('radius.png', dpi=300)
    plt.show()