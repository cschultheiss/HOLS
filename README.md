# HOLS
R Software for performing HOLS checks:
- HOLS_procedure.R : implementation of the HOLS check
- debias_z3.R : procedures adapted from the debiased Lasso to orthogonalize z_j^3

Create data for figures in paper (folder structure shall be adapated to user's needs):
- Figures 2, 3, 6: run SEM_simulation.R
- Figure 5: run mixgauss_simulation.R
- Figures 8 + 9: run SEMHD_simulation.R
- Figures 10 - 12: run block_simulation.R

The figures can be obtained with figures_execute.R. Note that the data are stored with time-dependent auto-created folder and file names. To get the figures, either the names or the R-script must be adapted to match eachother.


Create the table for the analys of the Sachs et al. dataset (https://www.science.org/doi/full/10.1126/science.1105809): store the files in a subfolder called Protein-signal, run cyto/cyto_table.R
