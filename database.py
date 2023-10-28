import csv
import numpy as np

z_list = []
n_list = []
hl_list = []
hl_e_list = []
q_list = []
q_e_list = []
binding_list = []
binding_e_list = []

with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        z, n, hl, hl_e, q, q_e = row[1], row[2], row[3], row[4], row[5], row[6]

        z_list.append(z)
        n_list.append(n)
        hl_list.append(hl)
        hl_e_list.append(hl_e)
        q_list.append(q)
        q_e_list.append(q_e)


binding_list = [None for i in z_list]
binding_e_list = [None for i in z_list]


nuclei_list = []
for i, z in enumerate(z_list):
    nuclei_list.append((z,n_list[i]))


with open('alldata.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        z, n, binding, binding_e = row[1], row[2], row[3], row[4]

        element = (str(int(z)-2), str(int(n)-2))

        if element in nuclei_list:
            index = [i for i, value in enumerate(nuclei_list) if value == element][0]
            binding_list[index] = binding
            binding_e_list[index] = binding_e


with open('output.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    for i, z in enumerate(z_list):
        row = [z_list[i], n_list[i], hl_list[i], hl_e_list[i], 
               q_list[i], q_e_list[i], binding_list[i], binding_e_list[i]]
        
        valid = all([(i is not None) and (i != '') for i in row])
        
        if valid:
            writer.writerow(row)
