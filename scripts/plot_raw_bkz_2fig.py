#!/usr/bin/python3

# Code in support of ePrint:2021/1384
# Copyright 2021, Tuong-Huy Nguyen
# GPL-3.0-only (see LICENSE file)

import os
from os import path
import sys
import matplotlib.pyplot as plt




if len(sys.argv) != 4:
    print("{} requires the number of orbits and two conductors for drawing".format(sys.argv[0][2:]) + "\n" + "Example: {} orb m1 m2"
          .format(sys.argv[0]));
    sys.exit(2);

## How many orbits ?
orb=int(sys.argv[1]);

## Input of 2 conductors for drawing next to one another
m1  = int(sys.argv[2])
m2 = int(sys.argv[3])
conductors = [m1,m2]

## Hardcoded parameter
dmax = 10


## General data needed
types = ["_su_iso_exp", "_su_noiso_exp", "_su_iso_tw", "_su_noiso_tw", "_urs_iso_exp", "_urs_noiso_exp", "_urs_iso_tw", "_urs_noiso_tw","_sat_iso_exp", "_sat_noiso_exp", "_sat_iso_tw", "_sat_noiso_tw"]
data_dir = os.path.dirname(os.getcwd()) +'/data/'


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4))

m = conductors[0]
for logemb in ['exp']:
    X, urs_raw, sat_raw, su_raw, urs_bkz, sat_bkz, su_bkz = [], [], [], [], [], [], []
    for set in ['urs', 'su', 'sat']:       
        file_path = data_dir + 'z' + str(m) + '/z' + str(m) + '_d' + str(orb) + '_' + set + '_noiso_' + logemb + '.gsn'
        if path.exists(file_path):
            if set == 'urs':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        X.append(values[0])
                        urs_raw.append(values[1]) 
                        urs_bkz.append(values[3]) 

            if set == 'su':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        su_raw.append(values[1]) 
                        su_bkz.append(values[3]) 

            if set == 'sat':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        sat_raw.append(values[1]) 
                        sat_bkz.append(values[3])
            

    if [len(urs_bkz),len(su_bkz),len(sat_bkz)] != [0,0,0]:
        ax1.set_ylim([1,3.5])
        if len(sat_bkz)>0:
            ax1.plot(X, sat_bkz, color='b', label='sat_bkz',linewidth=1) 
        if len(sat_raw)>0:
            ax1.plot(X, sat_raw, color='orange', label='sat_raw', linestyle='dashed',linewidth=0.8) 
        plt.xlim(xmin=0)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax1.text(0.05, 0.1, 'Q(z{}) with d={} and noiso/exp'.format(m,orb), transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props) 
        ax1.set_xlabel('index')
        ax1.xaxis.set_label_coords(0.5, -0.09)
        ax1.set_ylabel('ln(GSO)')
        ax1.yaxis.set_label_coords(-0.11, 0.5)
        ax1.set_axisbelow(True)
        ax1.xaxis.grid(color='lightgray', linestyle='dashed')
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(color='lightgray', linestyle='dashed')





m = conductors[1]
for logemb in ['exp']:
    X, urs_raw, sat_raw, su_raw, urs_bkz, sat_bkz, su_bkz = [], [], [], [], [], [], []
    for set in ['urs', 'su', 'sat']:       
        file_path = data_dir + 'z' + str(m) + '/z' + str(m) + '_d' + str(orb) + '_' + set + '_noiso_' + logemb + '.gsn'
        if path.exists(file_path):
            if set == 'urs':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        X.append(values[0])
                        urs_raw.append(values[1]) 
                        urs_bkz.append(values[3]) 

            if set == 'su':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        su_raw.append(values[1]) 
                        su_bkz.append(values[3]) 

            if set == 'sat':
                for line in open(file_path, 'r'):
                    if not line.startswith('#'):
                        values = [float(s) for s in line.split()]
                        sat_raw.append(values[1]) 
                        sat_bkz.append(values[3])
            

    if [len(urs_bkz),len(su_bkz),len(sat_bkz)] != [0,0,0]:

        ax2.autoscale()
        ax2.set_ylim([1,3.5])

        if len(sat_bkz)>0:
            ax2.plot(X, sat_bkz, color='b', label='sat_bkz',linewidth=1) 
        if len(sat_raw)>0:
            ax2.plot(X, sat_raw, color='orange', label='sat_raw', linestyle='dashed',linewidth=0.8) 
        plt.xlim(xmin=0)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax2.text(0.05, 0.1, 'Q(z{}) with d={} and noiso/exp'.format(m,orb), transform=ax2.transAxes, fontsize=10,
        verticalalignment='top', bbox=props) 
        ax2.set_xlabel('index')
        ax2.xaxis.set_label_coords(0.5, -0.09)
        leg = ax2.legend(ncol=1,fontsize=12)
        ax2.set_axisbelow(True)
        ax2.xaxis.grid(color='lightgray', linestyle='dashed')
        ax2.set_axisbelow(True)
        ax2.yaxis.grid(color='lightgray', linestyle='dashed')


    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper center',fontsize=12)

    plt.savefig(os.path.dirname(os.getcwd()) +'/figures/' + '/z{}-z{}_comparison_raw_bkz_d{}'.format(conductors[0], conductors[1],orb) + '.png')
    plt.show()
    
    plt.close(fig)
