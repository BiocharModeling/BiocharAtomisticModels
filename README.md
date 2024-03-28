# Biochar Atomistic Models
This new approach, implemented in Python and visualized using PyMOL, considers chemical constraints and characterization data. It incorporates polyaromatic structures and defects, randomly assigns functional groups and radicals, creates a cross-linked network between clusters, and captures porosity at a large-scale molecular level. 

The procedure starts with the ACSD obtained from BPCA data, followed by utilizing quantitative NMR data to integrate functional groups into the aromatic structures seamlessly. These functional groups include ether linkages, aliphatic chains (CH3, C2H5, and C3H7–the last for biochar at 400°C), carbonyl and carboxyl groups, and defective carbons. Additionally, ammonia, pyridine, hydroxyl, and pyrrolic groups are introduced to the structures to align with elemental analysis, FTIR, and XPS data. Finally, the H/C ratio establishes a cross-linked network between the modified clusters. Figure 1 illustrates this approach.

![plot](./Figures/1.png)
### Figure 1. Methodology to interpret characterization data and obtain atomistic representations of biochar. 
*Randomly select n-molecules to include the n-functional groups needed. Molecules can be selected multiple times, except for the inclusion of holes. 
Modifications are sequentially illustrated for the same molecule, but the process involves randomly selecting the molecule for each step. 

During the modification process, the clusters are strategically positioned in an "unoccupied" space to prevent molecule overlap and facilitate the modification procedure. As a result, the final coordinates of the modified clusters are converted from a PDB (Protein Data Bank) file to a LAMMPS data file with an atom-style charge, and the file conversion is carried out using OVITO; see Figure 2. 

![plot](./Figures/2.png)
### Figure 2. Approach to obtain the 3D representation of biochar. 

LAMMPS is used to construct a molecular dynamics-based biochar model. The simulation employs a canonical ensemble (NVT) at 298 K to mimic a compression process by reducing the simulation box size by < 1% in each time step until the target density is achieved. The temperature ensures a nonreactive environment by preventing bonds from breaking and forming. The C/H/O/N force field is applied with a 0.07 fs timestep.

###  How to use it? 
