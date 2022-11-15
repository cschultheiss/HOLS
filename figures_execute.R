source("figures_fun.R")

# folder/file names are manually set. Must be done accordingly when reproducing.

# creates Figures 2, 3, 6
HA_plot("results/SEM missing x3", exclude.chars = "+07", z.plot.ind = 1:6, z.plot.ind.label = (1:7)[-3], conf.ind = 2:3,
        beta0 = rep(0, 6), beta.OLS = sqrt(2.5) * c(0, -1/3, 2/3, 0, 0, 0), zlims.var = (0.1) * (1.1^(0:60)), which.line = 2:3,
        groups = list(2, 3, -c(2,3)), group.labels = c("2", "4", "U"))

# creates Figure 5
H0_plot("results/mix-Gauss null.RData")

# creates Figures 9, 8
HA_plot("results/SEM HD", z.plot = FALSE, groups = list(2, 3, -c(2,3)),
        group.labels = c("3", "4", "U"), hd = TRUE)

# creates Figures 10, 12, 11
HA_plot("results/block independent", z.plot.ind = c(1, 7, 9, 13, 14, 20, 22, 26), conf.ind = 1:13,
        beta0 = rep(0, 26), beta.OLS = c(0.3046434, -0.08445936, 0.01745333, 0.02742916, 0.1939396, 0.09398516,
                                         -0.5261166, -0.12003, -0.2272695, 0.03484701, 0.03762141, 0.05294573, 0.1150129,
                                         rep(0, 13)), zlims.var = (0.1) * (1.1^(0:60)), which.line = c(1, 7, 9, 13), colgroups = 2,
        groups = list(1, 9, 14:26), group.labels = c("1", "9", "U"))


# create the Table 1 in LaTeX format
source('cyto/cyto_table.R')
cyto_table()