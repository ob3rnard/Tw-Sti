#!/usr/bin/python3

# Code in support of ePrint:2021/xxxx
# Copyright 2021, Tuong-Huy Nguyen
# GPL-3.0-only (see LICENSE file)


import os
from os import path
import sys
import matplotlib.pyplot as plt
from sage.misc.misc import cputime
import statistics
from sympy.ntheory.factor_ import totient
import numpy as np







## General data needed
data_dir = os.path.dirname(os.getcwd()) +'/data/'
list_m = []
for line in open(data_dir + 'list_ms_pcmp', 'r'):
    for m in line.split():
        list_m.append(m)


## orbits max (for all fields)
dmax = 10





print('Drawing...')


## For drawing, we choose the minimal tag giving the better GF. For urs/sat we choose d=1 whereas for su we choose d=dmax of the field

## Choice d=1 for URS and sat
chosen_orb = 1


t = cputime()

## Get a dictionnary of the form {m : {d : [list of chosen values urs_set, sat_set, su_set] } }
## In this script the list only contain 3 values (at most): min(urs),min(sat),min(su)
dict_m_orbmeans = {}
for m in list_m:
    dict_orb_means = {}
    for orb in range(1,dmax+1):
        file_path = data_dir + 'z' + str(m) + '/z' + str(m)+ '_d'  + str(orb) + '.gf'
        if path.exists(file_path):
            urs_iso_exp,urs_iso_tw,urs_noiso_exp,urs_noiso_tw,sat_iso_exp,sat_iso_tw,sat_noiso_exp,sat_noiso_tw,su_iso_exp,su_iso_tw,su_noiso_exp,su_noiso_tw = [],[],[],[],[],[],[],[],[],[],[],[]
            list_tags = [urs_iso_exp,urs_iso_tw,urs_noiso_exp,urs_noiso_tw,sat_iso_exp,sat_iso_tw,sat_noiso_exp,sat_noiso_tw,su_iso_exp,su_iso_tw,su_noiso_exp,su_noiso_tw]
            Yurs, Ysat, Ysu = [], [], []
            for line in open(file_path, 'r'):
                if not line.startswith('#'):
                    values = [float(s) for s in line.split()]
                    Yurs.append(min([values[i] for i in range(0,4)]))
                    Ysat.append(min([values[i] for i in range(4,8)]))

                    if len(values)>8:
                        Ysu.append(min([values[i] for i in range(8,12)]))
            

            if len(Ysu)>0:
                dict_orb_means[orb] = [statistics.mean(Yurs),statistics.mean(Ysat),statistics.mean(Ysu)]
            else:
                dict_orb_means[orb] = [statistics.mean(Yurs),statistics.mean(Ysat)]
    dict_m_orbmeans[m] = dict_orb_means


## Prepare coordinates for plot
##  X = list of x-coordinates, YYurs = list of y-coordinates for urs, YYsat = list pf y-coordinates for sat
X = [totient(m) for m,v in dict_m_orbmeans.items() if v != {}]
YYurs = [dict_m_orbmeans[m][chosen_orb][0] for m,v in dict_m_orbmeans.items() if v != {}]
YYsat = [dict_m_orbmeans[m][chosen_orb][1] for m,v in dict_m_orbmeans.items() if v != {}]


## Prepare average coordinates for plot 
minX, maxX = min(X), max(X)
means_X = []
means_urs_Y = []
means_sat_Y = []
for abs in range(minX,maxX+1):
    list_urs_y = [y for x,y in zip(X,YYurs) if x==abs]
    list_sat_y = [y for x,y in zip(X,YYsat) if x==abs]
    if len(list_urs_y) >0:
        means_X.append(abs)
        means_urs_Y.append(statistics.mean(list_urs_y))
        means_sat_Y.append(statistics.mean(list_sat_y))



## Dealing with the particular case of su coordinates, since they may not exist
## Scatter plot 
Xsu,YYsu = [],[]
for m in [m for m,v in dict_m_orbmeans.items() if v != {}]:
    if len(dict_m_orbmeans[m][chosen_orb]) > 2:
        Xsu.append(totient(m))
        chosen_orb_for_su = 1
        for orb in range(1,dmax+1):
            file_path = data_dir + 'z' + str(m) + '/z' + str(m)+ '_d'  + str(orb) + '.gf'
            if path.exists(file_path):
                chosen_orb_for_su = orb
        YYsu.append(dict_m_orbmeans[m][chosen_orb_for_su][2])



## Uncomment the following parts to have URS/2-sat/SU(when available) for all experiments
fig, ax = plt.subplots(figsize=(9.5, 4)) 


##Draw URS (and average if needed)
ax.scatter(X, YYurs, marker="s",s=10, color='b', label='URS (d=1)')
# ax.plot(means_X, means_urs_Y, color='b',linewidth=0.5)

## Draw sat (and average if needed)
ax.scatter(X, YYsat, marker='o',s=10, color='orange', label='2-saturated URS (d=1)')
# ax.plot(means_X, means_sat_Y, color='r',linewidth=0.5)

##Draw su (and average if needed)
ax.scatter(Xsu, YYsu, marker='^',s=10, facecolor='None', edgecolor='purple',label='S-units (d=dmax)') 
# ax.plot(means_su_X, means_su_Y, color='y',linewidth=0.5)



plt.ylabel('Approximation Factor (under GH)')
plt.xlabel('Cyclotomic Field degree')
leg = ax.legend(ncol=1, fontsize=14)
ax.set_axisbelow(True)
ax.xaxis.grid(color='lightgray', linestyle='dashed')
ax.set_axisbelow(True)
ax.yaxis.grid(color='lightgray', linestyle='dashed')
ax.set_ylim([0,80])

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles, labels, loc='upper left',fontsize=12)


plt.savefig(os.path.dirname(os.getcwd()) +'/figures/GF.png')
# plt.show()
plt.close(fig)
print("Done: minGF as fonction of field dimension in [t={:.2f}]".format(t))

## Draw the zoom (without URS)
fig, ax = plt.subplots(figsize=(9.5, 4)) 

## Draw sat and su
ax.scatter(X, YYsat, marker='o',s=20, color='orange', label='2-saturated URS (d=1)')
ax.scatter(Xsu, YYsu, marker='^',s=20, facecolor='None', edgecolor='purple',label='S-units (d=dmax)') 

plt.ylabel('Approximation Factor (under GH)')
plt.xlabel('Cyclotomic Field degree')
leg = ax.legend(ncol=1, fontsize=14)
ax.set_axisbelow(True)
ax.xaxis.grid(color='lightgray', linestyle='dashed')
ax.set_axisbelow(True)
ax.yaxis.grid(color='lightgray', linestyle='dashed')
ax.set_xlim([20,101])
ax.set_ylim([0,5])


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles, labels, loc='upper left',fontsize=12)

plt.savefig(os.path.dirname(os.getcwd()) +'/figures/zoomGF.png')
# plt.show()
plt.close(fig)
print("Done: zoom minGF as fonction of field dimension in [t={:.2f}]".format(t))
