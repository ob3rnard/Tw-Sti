#!/usr/bin/python3

# Code in support of ePrint:2021/1384
# Copyright 2021, Tuong-Huy Nguyen
# GPL-3.0-only (see LICENSE file)

import os
from os import path
import sys
import matplotlib.pyplot as plt
from sage.misc.misc import cputime


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
X, urs_raw_iso_tw, sat_raw_iso_tw, su_raw_iso_tw,urs_raw_noiso_tw, sat_raw_noiso_tw, su_raw_noiso_tw = [], [], [], [],[],[],[]
urs_raw_iso_exp, sat_raw_iso_exp, su_raw_iso_exp,urs_raw_noiso_exp, sat_raw_noiso_exp, su_raw_noiso_exp = [], [], [],[],[],[]


file_path = data_dir + 'z{}/z{}_d{}_sat_iso_tw.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            X.append(values[0])
            sat_raw_iso_tw.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_iso_exp.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_iso_exp.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_noiso_tw.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_noiso_tw.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_noiso_exp.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_noiso_exp.append(values[1]) 


if [len(urs_raw_iso_tw),len(su_raw_iso_tw),len(sat_raw_iso_tw)] != [0,0,0]:
    ax1.set_ylim([1,4])
    if len(sat_raw_iso_tw)>0:
        ax1.plot(X, sat_raw_iso_tw, color='b', label='sat_raw_iso_tw',linewidth=0.3) 

    if len(sat_raw_noiso_tw)>0:
        ax1.plot(X, sat_raw_noiso_tw, color='b', label='sat_raw_noiso_tw',linewidth=0.3,linestyle='dashed',marker='^',markersize=0.8) 

    if len(sat_raw_iso_exp)>0:
        ax1.plot(X, sat_raw_iso_exp, color='r', label='sat_raw_iso_exp',linewidth=0.3) 

    if len(sat_raw_noiso_exp)>0:
        ax1.plot(X, sat_raw_noiso_exp, color='r', label='sat_raw_noiso_exp',linewidth=0.3,linestyle='dashed',marker='s',markersize=0.8) 
    plt.xlim(xmin=0)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    ax1.text(0.5, 0.1, 'Q(z{}) with d={}'.format(m, orb), transform=ax1.transAxes, fontsize=10,
    verticalalignment='top', bbox=props) 
    ax1.set_xlabel('index')
    ax1.xaxis.set_label_coords(0.5, -0.07)
    ax1.set_ylabel('ln(GSO)')
    ax1.yaxis.set_label_coords(-0.11, 0.5)
    ax1.set_axisbelow(True)
    ax1.xaxis.grid(color='lightgray', linestyle='dashed')
    ax1.set_axisbelow(True)
    ax1.yaxis.grid(color='lightgray', linestyle='dashed')
    



m = conductors[1]
X, urs_raw_iso_tw, sat_raw_iso_tw, su_raw_iso_tw,urs_raw_noiso_tw, sat_raw_noiso_tw, su_raw_noiso_tw = [], [], [], [],[],[],[]
urs_raw_iso_exp, sat_raw_iso_exp, su_raw_iso_exp,urs_raw_noiso_exp, sat_raw_noiso_exp, su_raw_noiso_exp = [], [], [],[],[],[]


file_path = data_dir + 'z{}/z{}_d{}_sat_iso_tw.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            X.append(values[0])
            sat_raw_iso_tw.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_iso_exp.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_iso_exp.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_noiso_tw.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_noiso_tw.append(values[1]) 
file_path = data_dir + 'z{}/z{}_d{}_sat_noiso_exp.gsn'.format(str(m),str(m),str(orb),set)
if path.exists(file_path):
    for line in open(file_path, 'r'):
        if not line.startswith('#'):
            values = [float(s) for s in line.split()]
            sat_raw_noiso_exp.append(values[1]) 


if [len(urs_raw_iso_tw),len(su_raw_iso_tw),len(sat_raw_iso_tw)] != [0,0,0]:
    ax2.autoscale()
    ax2.set_ylim([1,4])

    if len(sat_raw_iso_tw)>0:
        ax2.plot(X, sat_raw_iso_tw, color='b',linewidth=0.3) 

    if len(sat_raw_noiso_tw)>0:
        ax2.plot(X, sat_raw_noiso_tw, color='b', linewidth=0.3,linestyle='dashed',marker='^',markersize=0.8) 

    if len(sat_raw_iso_exp)>0:
        ax2.plot(X, sat_raw_iso_exp, color='r',linewidth=0.3) 

    if len(sat_raw_noiso_exp)>0:
        ax2.plot(X, sat_raw_noiso_exp, color='r', linewidth=0.3,linestyle='dashed',marker='s',markersize=0.8) 
    plt.xlim(xmin=0)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  
    ax2.text(0.5, 0.1, 'Q(z{}) with d={}'.format(m, orb), transform=ax2.transAxes, fontsize=10,
    verticalalignment='top', bbox=props) 
    ax2.set_xlabel('index')
    ax2.xaxis.set_label_coords(0.5, -0.07)
    ax2.set_axisbelow(True)
    ax2.xaxis.grid(color='lightgray', linestyle='dashed')
    ax2.set_axisbelow(True)
    ax2.yaxis.grid(color='lightgray', linestyle='dashed')
        

    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper center',fontsize=12)

    
    plt.savefig(os.path.dirname(os.getcwd()) +'/figures/' + '/z{}-z{}_comparison_sets_d{}'.format(conductors[0], conductors[1],orb) + '.png')
    plt.show()
    plt.close(fig)
