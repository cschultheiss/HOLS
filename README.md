# HOLS
R Software for performing HOLS checks:
- HOLS_procedure.R : implementation of the HOLS check
- debias_z3.R : procedures adapted from the debiased Lasso to orthogonalize z_j^3

Create data for figures (folder structure shall be adapated to user's needs):
- Figures 1 + 9: run null_simulation.R on commit   https://github.com/cschultheiss/HOLS/tree/e38e8ad8973370ebcc64fbfa128fc1e41357ecb0
- Figure 3 + 4: run HA_simulation.R on commit      https://github.com/cschultheiss/HOLS/tree/268c18940cfbdb4ffc8aade09260fe606ebcb2d6
- Figure 5: run nonlin_simulation.R on commit      https://github.com/cschultheiss/HOLS/tree/889dd281a1e88b53fd17cb364bf6aab08653c249
- Figure 6: run null_simulation.R on commit        https://github.com/cschultheiss/HOLS/tree/575ffc749c88acb52990c99ccb29afe2ca133a33
- Figures 7 + 8: run HA_simulation.R on commit     https://github.com/cschultheiss/HOLS/tree/40d2fd795b70980083b011fc2c53388019d66675

Create the tables for the analys of the Sachs et al. dataset (https://www.science.org/doi/full/10.1126/science.1105809): run cyto_analysis.R on commit https://github.com/cschultheiss/HOLS/tree/9a0feab83b2b6b6713156be0d97b634e52629e10
