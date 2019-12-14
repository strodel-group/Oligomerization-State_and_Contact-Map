# Oligomerization-State_and_Contact-Map


Following dependencies need to be fulfilled to run the codes.\
  1)Python 2.7 or 3.6 (From Anaconda distribution)\
  2)MDanalysis        (conda install -c conda-forge mdanalysis)\
  3)MDtraj            (conda install -c omnia mdtraj)

# Calculation of intermolecular residue-residue contact map and oligomer state among the protein chains
The oligos-cmap.py code calculates the highest oligomer state formed by protein chains in the system within a specfic cutoff distance. It indicates the time (in the form of snaphots from simulation trajectory) at which dissociation or association of protein chains occurs and also identifies the protein chains involved in aggregated state.\
Copy/Download the codes to the analysis directory consisting of the reference PDB structure (.pdb file, other formats can be tested), the molecular dynamics trajectory file only with protein atoms (.xtc file, other formats can be tested) and the Minimum distance in nano-metres(nm) to consider association or aggregation of proteins.

Usage: python oligos-cmap.py <.pdb file> <.xtc file> <distance in nm> 

Output files:  oligomer-groups.dat        - protein chains involved in aggregated state\
               oligomer-states.dat        - Number of chains involved in aggregated state (quantitative value)\
             **oligo-highest-size.dat**   - Highest Oligomer size (maximum = number of protein chains)\
             **contact-map.dat**          - Average over the frequency of inter-residue contacts between protein chains\
               oligo-block-average.dat      - Moving average to smooth out fluctuations over simulation time dependent observables

# Plotting the Inter-residue Contact Map:

Usage: python plot-cmap.py <.pdb file> <.xtc file> **contact-map.dat** 

Output files: Contact-map.pdf             - Image file indicating average frequency of contacts between residues

# Plotting the Oligomerisation state:

Usage: python plot-oligostate.py <.pdb file> <.xtc file>  **oligo-highest-size.dat**  <simulation time in nanoseconds(ns)>

Output files: Oligomerisation-state.pdf   - Image file indicating evolution of aggregation kinetics over time
  
  
Note: **<filename.dat>**  files are output files.
