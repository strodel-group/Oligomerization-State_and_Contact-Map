import matplotlib.pyplot as plt
from matplotlib import pyplot
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

#structure parameters

topology = md.load(ref_structure).topology
trajectory = md.load(traj, top=ref_structure)
frames=trajectory.n_frames								#Number of frames
chains=topology.n_chains								#Number of chains
atoms=int(topology.n_atoms/chains)						#Number of atoms in each monomer 
AminoAcids = int(topology.n_residues/chains)-2			#Number of residues per chain ('-2' avoid the N- and C- cap residues as individual residues)

isum=1
atoms_list=[]
residue_list=[]

for residue in topology.chain(0).residues:
    atoms_list.append(residue.n_atoms)
    residue_list.append(residue)
    ', '.join(map(lambda x: "'" + x + "'", str(residue_list)))
del residue_list[0]; del residue_list[-1] 


Norm_ContactMap = np.loadtxt(inp)
df=np.amax(Norm_ContactMap)
plt.imshow(Norm_ContactMap, cmap = 'Blues')

ticks = range(AminoAcids)
plt.xticks(list(ticks), labels=residue_list,rotation='vertical')
plt.yticks(list(ticks), labels=residue_list)
cbar=plt.colorbar(cmap='Blues', boundaries=np.arange(0.0, df+0.1, 0.1))
cbar=plt.cm.ScalarMappable(cmap='Blues')
plt.savefig('Contact-map.pdf')
plt.show()
