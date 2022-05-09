#!/usr/bin/python3

# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard, Tuong-Huy Nguyen
# GPL-3.0-only (see LICENSE file)

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
# --------------------------------------------------------------------------------------

import matplotlib.pyplot as plt


# --------------------------------------------------------------------------------------
# Data folder 
data_dir   = os.getcwd() + "/data";
conductors_file = data_dir + "/list_ms_pcmp";
figure_dir = os.getcwd() + "/figures";


# --------------------------------------------------------------------------------------
# List of conductors
list_m = [];
for line in open(conductors_file, 'r'):
    for m in line.split():
        list_m.append(int(m));


# --------------------------------------------------------------------------------------
# Reading data
measure = "gf";
su_sets = ["low_dpw","cdw","urs","sat","su"]; # In the order of cols in the data files

print('Drawing...')

X = {"all":[], "su":[]};
Y = {_su:[] for _su in su_sets };

for m in list_m:
    data_file_path = data_dir + "/z{}/z{}.gf".format(m,m);
    if not (os.path.exists(data_file_path)):
        print("\x1b[33m[Warn] Skipping m={} (no gf file)\x1b[0m".format(m));
        continue;

    for _line in open(data_file_path, 'r'):
        if _line.startswith('#'):
            continue;
        
        values = { _su:s for _su,s in zip(["n"]+su_sets, _line.rstrip().split()) };
        assert(all(values[_k] != "NaN" for _k in values.keys() if _k != "su"));

        # Abscissae
        X["all"].append(int(values["n"]));
        if values["su"] != "NaN":
            X["su"].append(int(values["n"]));
        # Data points
        [Y[_su].append(float(values[_su])) for _su in su_sets if values[_su] != "NaN"];




# --------------------------------------------------------------------------------------
# Draw without CDW
fig, ax = plt.subplots(figsize=(9.5, 4)) 

# Draw URS
ax.scatter(X["all"], Y["urs"], marker="s",s=10, color='b',
           label='URS (d=1)');
# Draw sat 
ax.scatter(X["all"], Y["sat"], marker='o',s=10, color='orange',
           label='2-saturated URS (d=1)');
# Draw su 
ax.scatter(X["su"], Y["su"], marker='^',s=10, facecolor='None', edgecolor='purple',
           label='S-units (d=dmax)');

# Legends/...
plt.ylabel('Approximation Factor (under GH)');
plt.xlabel('Cyclotomic field degree');
leg = ax.legend(ncol=1, fontsize=14);
ax.set_axisbelow(True);
ax.xaxis.grid(color='lightgray', linestyle='dashed');
ax.set_axisbelow(True);
ax.yaxis.grid(color='lightgray', linestyle='dashed');
ax.set_ylim([0,40]);

handles, labels = ax.get_legend_handles_labels();
plt.legend(handles, labels, loc='upper left',fontsize=12);

# Save
plt.savefig(figure_dir + "/GF-mprandom2.png", bbox_inches='tight')
plt.close(fig)
print("Done: minGF as fonction of field dimension");



# --------------------------------------------------------------------------------------
# Draw the zoom (without URS)
fig, ax = plt.subplots(figsize=(9.5, 4)) 

# Draw sat and su
ax.scatter(X["all"], Y["sat"], marker='o', s=20, color='orange',
           label='2-saturated URS (d=1)');
ax.scatter(X["su"], Y["su"], marker='^', s=20, facecolor='None', edgecolor='purple',
           label='S-units (d=dmax)');

# Legends/...
plt.ylabel('Approximation Factor (under GH)');
plt.xlabel('Cyclotomic field degree');
leg = ax.legend(ncol=1, fontsize=14);
ax.set_axisbelow(True);
ax.xaxis.grid(color='lightgray', linestyle='dashed');
ax.set_axisbelow(True);
ax.yaxis.grid(color='lightgray', linestyle='dashed');
ax.set_xlim([20,101]);
ax.set_ylim([0,4]);

handles, labels = ax.get_legend_handles_labels();
plt.legend(handles, labels, loc='upper left', fontsize=12);

# Save
plt.savefig(figure_dir + "/zoomGF-mprandom2.png", bbox_inches='tight');
plt.close(fig);
print("Done: zoom minGF as fonction of field dimension"); 



# --------------------------------------------------------------------------------------
# Draw with CDW and CDW lower bound
fig, ax = plt.subplots(figsize=(9.5, 4)) 

# Draw URS 
ax.scatter(X["all"], Y["urs"], marker="s", s=10, color='b',
           label='URS (d=1)');
# Draw sat 
ax.scatter(X["all"], Y["sat"], marker='o', s=10, color='orange',
           label='2-saturated URS (d=1)');
# Draw su 
ax.scatter(X["su"], Y["su"], marker='^', s=10, facecolor='None', edgecolor='purple',
           label='S-units (d=dmax)') ;
# Draw CDW 
ax.scatter(X["all"], Y["cdw"], marker="p", s=3, color='red',
           label='CDW (d=1)');
# Draw asymptotic lower bound DPW
ax.scatter(X["all"], Y["low_dpw"], marker=".", s=3, color='green',
           label='CDW lower bound (d=1)');

# Legends/...
plt.ylabel('Approximation Factor (under GH)');
plt.xlabel('Cyclotomic field degree');
leg = ax.legend(ncol=1, fontsize=14);
ax.set_axisbelow(True);
ax.xaxis.grid(color='lightgray', linestyle='dashed');
ax.set_axisbelow(True);
ax.yaxis.grid(color='lightgray', linestyle='dashed');
ax.set_ylim([0,40]);

handles, labels = ax.get_legend_handles_labels();
plt.legend(handles, labels, loc='upper left', fontsize=12);

# Save
plt.savefig(figure_dir + "/GF_CDW-mprandom2.png", bbox_inches='tight');
plt.close(fig);
print("Done: minGF as fonction of field dimension (with CDW data)");


# --------------------------------------------------------------------------------------
# Draw with CDW and CDW lower bound, without URS
fig, ax = plt.subplots(figsize=(9.5, 4)) 

# # Draw URS 
# ax.scatter(X["all"], Y["urs"], marker="s", s=10, color='b',
#            label='URS (d=1)');
# Draw sat 
ax.scatter(X["all"], Y["sat"], marker='o', s=10, color='orange',
           label='2-saturated URS (d=1)');
# Draw su 
ax.scatter(X["su"], Y["su"], marker='^', s=10, facecolor='None', edgecolor='purple',
           label='S-units (d=dmax)') ;
# Draw CDW 
ax.scatter(X["all"], Y["cdw"], marker="p", s=3, color='red',
           label='CDW (d=1)');
# Draw asymptotic lower bound DPW
ax.scatter(X["all"], Y["low_dpw"], marker=".", s=3, color='green',
           label='CDW lower bound (d=1)');

# Legends/...
plt.ylabel('Approximation Factor (under GH)');
plt.xlabel('Cyclotomic field degree');
leg = ax.legend(ncol=1, fontsize=14);
ax.set_axisbelow(True);
ax.xaxis.grid(color='lightgray', linestyle='dashed');
ax.set_axisbelow(True);
ax.yaxis.grid(color='lightgray', linestyle='dashed');
ax.set_ylim([0,40]);

handles, labels = ax.get_legend_handles_labels();
plt.legend(handles, labels, loc='upper left', fontsize=12);

# Save
plt.savefig(figure_dir + "/GF_CDW_noURS-mprandom2.png", bbox_inches='tight');
plt.close(fig);
print("Done: minGF as fonction of field dimension (with CDW data)");


# --------------------------------------------------------------------------------------
# Draw the zoom (without URS, with CDW lower bound)
fig, ax = plt.subplots(figsize=(9.5, 4)) 

# Draw sat and su and cdw lower bound
ax.scatter(X["all"], Y["sat"], marker='o', s=20, color='orange',
           label='2-saturated URS (d=1)');
ax.scatter(X["su"], Y["su"], marker='^', s=20, facecolor='None', edgecolor='purple',
           label='S-units (d=dmax)');
ax.scatter(X["all"], Y["low_dpw"], marker=".", s=3, color='green',
           label='CDW lower bound (d=1)');

# Legends/...
plt.ylabel('Approximation Factor (under GH)');
plt.xlabel('Cyclotomic field degree');
leg = ax.legend(ncol=1, fontsize=14);
ax.set_axisbelow(True);
ax.xaxis.grid(color='lightgray', linestyle='dashed');
ax.set_axisbelow(True);
ax.yaxis.grid(color='lightgray', linestyle='dashed');
ax.set_xlim([20,101]);
ax.set_ylim([0,4]);

handles, labels = ax.get_legend_handles_labels();
plt.legend(handles, labels, loc='upper left', fontsize=12);

# Save
plt.savefig(figure_dir + "/zoomGF_with_cdw_lower-mprandom2.png", bbox_inches='tight');
plt.close(fig);
print("Done: zoom minGF (with CDW lower bound) as fonction of field dimension"); 


sys.exit(0);
