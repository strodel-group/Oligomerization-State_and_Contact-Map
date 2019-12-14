import matplotlib.pyplot as plt
from matplotlib import pyplot
from pylab import genfromtxt
import numpy as np
import copy as cp
import subprocess as sp
import os as os
import shutil as sh
import MDAnalysis as mdana
import sys
import networkx as nx
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use("TkAgg")

#input parameters

ref_structure=sys.argv[1]
traj=sys.argv[2]
inp = sys.argv[3]
sim_t =int(sys.argv[4])												#Simulation time in nano-seconds

#structure parameters

topology = md.load(ref_structure).topology
trajectory = md.load(traj, top=ref_structure)
frames=trajectory.n_frames											#Number of frames
chains=topology.n_chains											#Number of chains
atoms=int(topology.n_atoms/chains)									#Number of atoms in each monomer 
AminoAcids = int(topology.n_residues/chains)-2						#Number of residues per chain ('-2' avoid the N- and C- cap residues as individual residues)
adj=frames/sim_t													#Snaphots saved every "adj" ns


for i in range(len(inp)):
    data = genfromtxt(inp)
    plt.plot(data[:,0]/adj, data[:,1], linewidth=2, linestyle='-', color = '#ff0e15')

lty=np.arange(0.5,chains+2,1)
ltx=np.arange(0,(frames/adj)+(sim_t/10),(sim_t/4))

ticks_y = range(0,len(lty),1)
ticks_x = np.arange(0,(frames/adj)+(sim_t/10),(sim_t/4))
plt.yticks(list(ticks_y), labels=np.arange(chains+2))
plt.xticks(list(ticks_x), labels=np.arange(0,(frames/adj)+(sim_t/10),(sim_t/4)))

pyplot.ylabel('Oligomerization state')
pyplot.xlabel('Time(ns)')

plt.savefig('Oligomerisation-state.pdf', bbox_inches='tight', quality=200, dpi=400)
plt.show()
