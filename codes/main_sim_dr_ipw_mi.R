# TODO all works
##

## pmm paper 



# generate data -----------------------------------------------------------
## seed
set.seed(seed_for_sim)

## pop and sample sizes
N <- 100000
n_a <- c(500,1000)
n_b <- 500
n_a1 <- 0.7 * n_a
n_a2 <- 0.3 * n_a
## generate data
x1 <- rnorm(N, 2, 1)
x2 <- rnorm(N, 2, 1)
x3 <- rnorm(N, 2, 1)
e <- rnorm(N)
y1 <- 1 + 2*x1 + e
y2 <- -1 + x1 + x2 + x3 + e
y3 <- -10 + x1^2 + x2^2 + x3^2 +e
p1 <- plogis(1 + x1 + e)
p2 <- plogis(.5 + x1 + x2 + x3 + e)
p3 <- plogis(-2 + x1^2 + x2^2 + x3^2 + e)
strata <- x1 <= 2
pop <- data.frame(x1, x2, x3, y1, y2, y3, p1, p2, p3, strata)

y_formula_1 <- y1 ~ x1
y_formula_2 <- y2 ~ x1 + x2 + x3
y_formula_3_all <- y3 ~ I(x1^2) + I(x2^2) + I(x3^2)
y_formula_3_nn <- y_formula_3_all
y_formula_3_mis <- y3 ~ x1 + x2

pi_formula_1 <- ~ x1
pi_formula_2 <- ~ x1 + x2 + x3
pi_formula_3_all <- ~ I(x1^2) + I(x2^2) + I(x3^2)
pi_formula_3_mis <- ~ x1 + x2

y1_tar <- ~ y1
y2_tar <- ~ y2
y3_tar <- ~ y3

# main simulation ---------------------------------------------------------

## setup for parallel computation

cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())


results_simulation1 <- foreach(
  k=1:sims, .combine = rbind, .packages = c("survey", "nonprobsvy"), .options.snow = opts) %dopar% {

    ## nonprob sample
    pop1 <- subset(pop, strata == TRUE)
    pop2 <- subset(pop, strata == FALSE)
    sample_a_500 <- rbind(pop1[sample(1:nrow(pop1), n_a1[1]), ], pop2[sample(1:nrow(pop2), n_a2[1]), ])
    sample_a_1000 <- rbind(pop1[sample(1:nrow(pop1), n_a1[2]), ], pop2[sample(1:nrow(pop2), n_a2[2]), ])

    ## sample prob
    sample_b <- pop[sample(1:N, n_b), ]
    sample_b$w_b <- N / n_b
    svy_b <- svydesign(ids = ~ 1, weights = ~ w_b, data = sample_b)

    ## estimators
    ## true
    trues <- colMeans(pop[, c("y1", "y2", "y3")])
    ## naive
    naive_500 <- colMeans(sample_a_500[, c("y1", "y2", "y3")])
    naive_1000 <- colMeans(sample_a_1000[, c("y1", "y2", "y3")])

    ## mi glm
    mi_glm_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b)

    mi_glm_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b)

    # ipw mle

    ipw_mle_500_y1 <- nonprob(selection = pi_formula_1, target = y1_tar, data = sample_a_500, svydesign = svy_b)
    ipw_mle_500_y2 <- nonprob(selection = pi_formula_2, target = y2_tar, data = sample_a_500, svydesign = svy_b)
    ipw_mle_500_y3_all <- nonprob(selection = pi_formula_3_all, target = y3_tar, data = sample_a_500, svydesign = svy_b)
    ipw_mle_500_y3_mis <- nonprob(selection = pi_formula_3_mis, target = y3_tar, data = sample_a_500, svydesign = svy_b)

    ipw_mle_1000_y1 <- nonprob(selection = pi_formula_1, target = y1_tar, data = sample_a_1000, svydesign = svy_b)
    ipw_mle_1000_y2 <- nonprob(selection = pi_formula_2, target = y2_tar, data = sample_a_1000, svydesign = svy_b)
    ipw_mle_1000_y3_all <- nonprob(selection = pi_formula_3_all, target = y3_tar, data = sample_a_1000, svydesign = svy_b)
    ipw_mle_1000_y3_mis <- nonprob(selection = pi_formula_3_mis, target = y3_tar, data = sample_a_1000, svydesign = svy_b)


    # dr glm mle

    dr_glm_500_y1 <- nonprob(selection = pi_formula_1, outcome = y_formula_1, data = sample_a_500, svydesign = svy_b)
    dr_glm_500_y2 <- nonprob(selection = pi_formula_2, outcome = y_formula_2, data = sample_a_500, svydesign = svy_b)
    dr_glm_500_y3_all <- nonprob(selection = pi_formula_3_all, outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b) # to fix
    dr_glm_500_y3_mis <- nonprob(selection = pi_formula_3_mis, outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b)

    dr_glm_1000_y1 <- nonprob(selection = pi_formula_1, outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b)
    dr_glm_1000_y2 <- nonprob(selection = pi_formula_2, outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b)
    dr_glm_1000_y3_all <- nonprob(selection = pi_formula_3_all, outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b) # to fix
    dr_glm_1000_y3_mis <- nonprob(selection = pi_formula_3_mis, outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b)


    data.frame(
      k        = k,
      y        = c("y1", "y2", "y3_all", "y3_mis"),
      trues    = trues[c(1,2,3,3)],
      naive_500  = naive_500[c(1,2,3,3)],
      naive_1000 = naive_1000[c(1,2,3,3)],
      ## srs
      glm_500    = c(mi_glm_500_y1$output$mean, mi_glm_500_y2$output$mean, mi_glm_500_y3_all$output$mean, mi_glm_500_y3_mis$output$mean),
      glm_1000   = c(mi_glm_1000_y1$output$mean, mi_glm_1000_y2$output$mean, mi_glm_1000_y3_all$output$mean, mi_glm_1000_y3_mis$output$mean),
      ipw_500    = c(ipw_mle_500_y1$output$mean, ipw_mle_500_y2$output$mean, ipw_mle_500_y3_all$output$mean, ipw_mle_500_y3_mis$output$mean),
      ipw_1000   = c(ipw_mle_1000_y1$output$mean, ipw_mle_1000_y2$output$mean, ipw_mle_1000_y3_all$output$mean, ipw_mle_1000_y3_mis$output$mean),
      dr_500     = c(dr_glm_500_y1$output$mean, dr_glm_500_y2$output$mean, dr_glm_500_y3_all$output$mean, dr_glm_500_y3_mis$output$mean),
      dr_1000    = c(dr_glm_1000_y1$output$mean, dr_glm_1000_y2$output$mean, dr_glm_1000_y3_all$output$mean, dr_glm_1000_y3_mis$output$mean)
    )

  }

# processing results ------------------------------------------------------

## processing results
setDT(results_simulation1)

results_simulation1_process <- results_simulation1 |> melt(id.vars = 1:3)
#results_simulation1_process[, c("est", "sample", "ci"):=tstrsplit(variable, "_")]
#results_simulation1_process[sample %in% c("500b","1000b"), ci_boot := TRUE]
#results_simulation1_process[sample == "500b", sample := "500"]
#results_simulation1_process[sample == "1000b", sample := "1000"]

saveRDS(results_simulation1_process, file = "results/main_sim_dr_ipw_mi.RDS")
