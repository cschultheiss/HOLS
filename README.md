# HOLS
R Software for performing HOLS checks:
- HOLS_procedure.R : implementation of the HOLS check
- debias_z3.R : procedures adapted from the debiased Lasso to orthogonalize z_j^3

Create data for figures in paper (folder structure shall be adapated to user's needs):
- Figures 2, 3, 6: run HA_simulation.R on this [commit](https://github.com/cschultheiss/HOLS/tree/268c18940cfbdb4ffc8aade09260fe606ebcb2d6) (R-version used for the paper: 4.1.0)
- Figure 5: run null_simulation.R on this [commit](https://github.com/cschultheiss/HOLS/tree/e38e8ad8973370ebcc64fbfa128fc1e41357ecb0) (R-version used for the paper: 4.1.0)
- Figures 8 + 9: run HA_simulation.R on this [commit](https://github.com/cschultheiss/HOLS/tree/429d1d05ccbdc48ef061951e4e165088bc5da88c) (R-version used for the paper: 4.1.1)
- Figures 10 - 12: run HA_simulation.R on this [commit](https://github.com/cschultheiss/HOLS/tree/40d2fd795b70980083b011fc2c53388019d66675) (R-version used for the paper: 4.1.0)

The figures can be obtained with figures_execute.R. Note that the data are stored with time-dependent auto-created folder and file names. To get the figures, either the names or the R-script must be adapted to match eachother.


Create the table for the analys of the Sachs et al. dataset (https://www.science.org/doi/full/10.1126/science.1105809): run cyto_analysis.R on this [commit](https://github.com/cschultheiss/HOLS/tree/bd0b055f14d09aab7bf601aba5014cce26b450c5)

