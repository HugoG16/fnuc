"""
name,z,n,halflife(Seconds),d-halflife,alpha(keV),d-alpha
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
fpath = Path(mpl.get_data_path(), "fonts/ttf/cmr10.ttf")
import matplotlib.font_manager as font_manager
font_path = fpath
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()
plt.rcParams['font.size'] = 12
plt.rcParams["axes.formatter.use_mathtext"] = True

from helper import *

data = list()
equal_z = dict()

## read data
with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        z, n, hl, hl_e, q, q_e = row[1], row[2], row[3], row[4], row[5], row[6]

        try:
            z, n, hl, hl_e, q, q_e = float(z), float(n), float(hl), float(hl_e), float(q), float(q_e)
        except:
            z, n, hl, hl_e, q, q_e = float(z), float(n), float(hl), 0, float(q), float(q_e)

        nuclei = Nuclei(z, n, hl, hl_e, q, q_e)
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
if 0:
    cmap = plt.get_cmap('turbo', 16)
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
    plt.xlabel('Q (keVs) ^ -1/2')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid()
    plt.tight_layout()
    # plt.savefig('geiger.png', dpi=300)
    plt.show()


############ fit func ############
## separete data into equal A lists
# equal_a = dict()
# for nuclei in data:
#     try:
#         equal_a[nuclei.a].append(nuclei)
#     except KeyError:
#         equal_a[nuclei.a] = [nuclei]

# print(equal_a)
# fitter = Fit(1, [1,1], [1,1], [1,1], [1,1], [0.1,0.1])
# fitter.fit()
