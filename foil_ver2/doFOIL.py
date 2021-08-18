import FOIL
from time import sleep

#Phase_1
efficiency = 'efficiency.txt'
#Phase_2
NIST_file = 'NIST.txt'
JENDL = 'jendl-ad2017_300K.tar/jendl-ad2017_300K/'
density_file = 'density.csv'
Gammas_file = 'Gammas_molded.csv'
attenuation_folder = 'Attenuation'
Threshold_MIN = 0
Threshold_MAX = 20E+6

#計算準備##########################################################################################################################################
FOIL.Phase_1(efficiency)
FOIL.Phase_2(NIST_file,JENDL,density_file,Gammas_file,attenuation_folder,Threshold_MIN,Threshold_MAX)
result_list, cross_section_x, cross_section_y = FOIL.read_results(Threshold_MIN,Threshold_MAX)
###################################################################################################################################################

import numpy as np
from tqdm import tqdm

file = 'Fe_spectrum.txt'
f = open(file)
data = f.readlines()
energy = []
flax = []
DT = 10**9
for i in data:
    energy.append(float(i.split('\t')[0])*10E6)
    flax.append(float(i.split('\t')[1])*DT)
f.close()

Counts = []
Na = 6.022e23
ta = 8*3600
tb = 1*3600
tc = 15*3600
foil_thickness = 0.5 
foil_r = 1.5
count_threshold = 1
a=[]
for i in tqdm(range(len(result_list))):
    re_cross = FOIL.Calculate_interpolation(cross_section_x[i],cross_section_y[i],energy)
    N = foil_thickness*foil_r**2*3.14*Na*result_list[i][9]/result_list[i][8]
    R = 0
    for j in range(len(flax)):
        R += flax[j]*re_cross[j]*10e-24

    loss = (1-np.exp(-foil_thickness*result_list[i][9]*result_list[i][13]))/(foil_thickness*result_list[i][9]*result_list[i][13])*result_list[i][14]
    activate_num = N*R*result_list[i][9]/result_list[i][8]*result_list[i][7]*result_list[i][11]*result_list[i][12]/np.log(2)
    count = activate_num*loss*(1-np.exp(-np.log(2)*ta/result_list[i][12]))*np.exp(-np.log(2)*tb/result_list[i][12])*(1-np.exp(-np.log(2)*tc/result_list[i][12]))
    if count > count_threshold:
        Counts.append([result_list[i][0],result_list[i][1],result_list[i][10],(1-np.exp(-np.log(2)*ta/result_list[i][12]))*np.exp(-np.log(2)*tb/result_list[i][12])*(1-np.exp(-np.log(2)*tc/result_list[i][12])),count])
sleep(1)
for i in range(len(Counts)):
    print(Counts[i])