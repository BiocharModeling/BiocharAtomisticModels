# BiocharAtomisticModels

This new approach, implemented in Python and visualized using PyMOL, considers chemical constraints and characterization data. It incorporates polyaromatic structures and defects, randomly assigns functional groups and radicals, creates a cross-linked network between clusters, and captures porosity at a large-scale molecular level. The process begins with the ACSD obtained from BPCA data (go to BPCA.py). Then, quantitative NMR data is employed to integrate functional groups into the aromatic structures. These functional groups encompass ether linkages, aliphatic chains (including CH3, C2H5, and C3H7–the last for biochar at 400°C), carbonyl functionalities, and carboxyl groups defective carbons. Additionally, ammonia, pyridine, hydroxyl, and pyrrolic groups are introduced to the structures to align with elemental analysis, FTIR, and XPS data. Defective carbons are represented as curvature structures rather than planar clusters, including non-hexagonal rings. Incorporating such defects is fundamental, not only due to their disruption of the overall structure of the biochar, which prevents the formation of planar regions but also because they amplify biochar's reactivity and adsorption behavior, which results from the introduction of heightened reactive defects and edge sites. Finally, the H/C ratio establishes a cross-linked network between the modified clusters. Figure 1 illustrates this approach.

![plot](./Figures/1.png)
### Figure 1. Methodology to interpret characterization data and obtain atomistic representations of biochar. 
*Randomly select n-molecules to include the n-functional groups needed. Molecules can be selected multiple times, except for the inclusion of holes. 
Modifications are sequentially illustrated for the same molecule, but the process involves randomly selecting the molecule for each step. 

During the modification process, the clusters are strategically positioned in an "unoccupied" space to prevent molecule overlap and facilitate the modification procedure. As a result, the final coordinates of the modified clusters are converted from a PDB (Protein Data Bank) file to a LAMMPS data file with an atom-style charge, and the file conversion is carried out using OVITO; see Figure 2. 

![plot](./Figures/2.png)
### Figure 2. Approach to obtain the solid representation of biochar. 
