# Biochar Atomistic Models
This new approach, implemented in Python and visualized using PyMOL, considers chemical constraints and characterization data. It incorporates polyaromatic structures and ring defects, randomly assigns functional groups and radicals, creates a cross-linked network between clusters, and captures porosity at a large-scale molecular level. Figure 1 illustrates this approach.

![plot](./Figures/1.png)
### Figure 1. Methodology to interpret characterization data and obtain atomistic representations of biochar. 
*Randomly select n-molecules to include the n-functional groups needed. Molecules can be selected multiple times, except for the inclusion of holes. 
Modifications are sequentially illustrated for the same molecule, but the process involves randomly selecting the molecule for each step. 

The structures are in a well-dispersed grid during construction stages to prevent molecule overlap and facilitate modification. After the cross-linking, the grid is reconstructed to place all the molecules closer, mimicking a gas phase. As a result, the final coordinates of the modified clusters are converted from a Protein Data Bank (PDB) file to a Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) data file with an atom-style charge, and the file conversion is carried out using OVITO; see Figure 6. LAMMPS is used to construct an MD-generated biochar representation

![plot](./Figures/2.png)
### Figure 2. Approach to obtain the 3D representation of biochar. 

LAMMPS is used to construct a molecular dynamics-based biochar model. We employed a canonical ensemble (NVT) at 298 K to mimic a compression process by reducing the simulation box size by < 1% in each time step until the target density is achieved. The temperature ensures a nonreactive environment by preventing bonds from breaking and forming. The C/H/O/N force field is applied with a 0.07 fs timestep.

You can use any other tool to obtain the 3D representation from the PBD file obtained after running Biochar.ipynb 

> [!Note]
> You will need at least the characterization data from the ultimate analysis, <sup>13</sup>C NMR (quantitative), and BPCA yields

> [!Tip]
> If you don't have the BPCA yields, you can use data from the literature or the data in the database file. However, this will increase the error by up to 30%

> [!Version used]
> Python 3.10.13 and PyMOL 2.5.5

### How to use it?

```
git clone https://github.com/BiocharModeling/BiocharAtomisticModels
cd BiocharAtomisticBuilder
```
- Install PyMol (https://github.com/schrodinger/pymol-open-source?tab=readme-ov-file)
- Connect PyMol to Jupyter Notebook (https://pymolwiki.org/index.php/Jupyter)

###  [How to run it? :computer:](https://www.youtube.com/watch?v=V4YJs-slrXk)

### Contributing :bulb:
We are open to contributing if you want to improve or suggest modifications. If you have issues running the code, don't hesitate to contact us; we will do our best to solve them. 

### Contact us :email:
**Email**: 
[Valentina Sierra](mailto:v.sierrajimenez@wsu.edu) <br>
**Address:** <br>
Biological Systems Engineering <br>
Washington State University <br> 
Pullman, WA, USA
