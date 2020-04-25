# Oligomerization-State_and_Contact-Map


Following dependencies need to be fulfilled to run the codes.
  - Python 2.7 or 3.6 (`From Anaconda distribution`)
  - MDanalysis        (`conda install -c conda-forge mdanalysis`)
  - MDtraj            (`conda install -c omnia mdtraj`)

# Calculation of intermolecular residue-residue contact map and oligomer state among the protein chains
The oligos-cmap.py code calculates the highest oligomer state formed by protein chains in the system within a specfic cutoff distance. It indicates the time (in the form of snaphots from simulation trajectory) at which dissociation or association of protein chains occurs and also identifies the protein chains involved in aggregated state.\
Copy/Download the codes to the analysis directory consisting of the reference PDB structure (.pdb file, other formats can be tested), the molecular dynamics trajectory file only with protein atoms (.xtc file, other formats can be tested) and the Minimum distance in nano-metres(nm) to consider association or aggregation of proteins.

Usage: `python2 oligos-cmap.py <.pdb file> <.xtc file> <distance in Angstrom>`\
Output files:\
             -   oligomer-groups.dat       - `groups the interacting protein chains.`\
             -   oligomer-states.dat       - `counts the number of chains in each oligomer group.`\
             - **oligo-highest-size.dat**  - `finds the maximum oligomer size formed per frame. (maximum size = total number of protein chains in simulation box)`\
             - **contact-map.dat**         - `saves the frequency of contacts between residues from different proteins.`\
             -   oligo-block-average.dat   - `creates a series of 25-frame moving average of  the maximum oligomer size.`

# Plotting the Inter-residue Contact Map:

Usage: `python3 plot-cmap.py <.pdb file> <.xtc file>` **contact-map.dat**

Output files: Contact-map.pdf             - Image file indicating average frequency of contacts between residues

# Plotting the Oligomerisation state:

Usage: `python3 plot-oligostate.py <.pdb file> <.xtc file> ` **oligo-highest-size.dat**  `<simulation time in nanoseconds(ns)>`

Output files: Oligomerisation-state.pdf   - Image file indicating evolution of aggregation kinetics over time
  
  
# Note: **<filename>.dat**  files are output files.
-----------------------------------------------------# Examples #-----------------------------------------------------\
Walkthrough the codes using the sample structure and trajectory files in example.zip\
Usage:
- python2 oligos-cmap.py protein_ref.pdb  protein_md.xtc 4
- python3 plot-cmap.py protein_ref.pdb  protein_md.xtc contact-map.dat
- python3 plot-oligostate.py protein_ref.pdb  protein_md.xtc oligo-highest-size.dat 1000
  
  
