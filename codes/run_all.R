## packages
remotes::install_github("ncn-foreigners/nonprobsvy@dev", force = TRUE)
library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(data.table)
library(xtable)
library(dplyr)

## seed
seed_for_sim <- 2024
cores <- 8

start_time <- Sys.time()

## run all
sims <- 500
source("codes/main_sim_dr_ipw_mi.R", echo=TRUE)
source("codes/dr.R", echo=TRUE)
source("codes/ipw.R", echo=TRUE)
source("codes/mi.R", echo=TRUE)
source("codes/var_sel.R", echo=TRUE)
source("codes/variance.R", echo=TRUE)

end_time <- Sys.time() - start_time
