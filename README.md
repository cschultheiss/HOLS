# HOLS
R Software for performing HOLS checks:<br>
HOLS_procedure.R: implementation of the HOLS check <br>
Function HOLS.check() can be applied to data <br>
debias_z3.R procedures adapted from the debiased Lasso to orthogonalize z_j^3

To create the data for the figures in the paper:
- Figures 2, 3, 6: run the function stored in SEM_simulation.R with its default parameters
- Figure 5: run the function stored in mixgauss_simulation.R with its default parameters
- Figures 8 + 9 run the function stored in SEMHD_simulation.R with its default parameters
- Figures 10 - 12: run the function stored in block_simulation.R with its default parameters

To get Table 1 (in LaTeX encoding) for the analysis of the Sachs et al. dataset (https://www.science.org/doi/full/10.1126/science.1105809): run the function stored in cyto/cyto_table.R <br>
The files must be stored in a subfolder called Protein-signal <br>
Helper functions with explaining comments are stored in cyto/cyto_functions.R <br>

The file figures_execute.R contains a complete workflow to obtain each figure and the table <br>
The called plotting functions with explaining comments are stored in figures_fun.R
