# 16S-central-folding

### Simulation code
The source code for the simulation is stored in `src/` with all the necessary parameters and coordinate files.
To compile, use a standard c++ compiler. There is no dependent library etc.
See `src/run_sample.sh` to run the program.  The code was entirely written by Dr. Natalia A. Denesyuk. 


### Supplementary Data
The supplemental dataset is available at https://doi.org/10.5281/zenodo.4304537

In a subdirectory `ContactMgConcentrations/`, numerical tables of <i>contact Mg concentrations</i> (See Eq. 1 in the paper for definition) are stored. 
Each filename indicates the bulk Mg concentration of the simulations. 
The first column of each table is nucleotide numbers, and the second column is the contact Mg concentration.

In `MgDensities` directory, data of the three-dimensional Mg densities are stored as OpenDx format.  
These can be visualized with standard software such as VMS and Chimera. A PDB file is also included for the reference structure.

