### Doubly Robust Inference With Nonprobability
### Survey Samples

set.seed(seed_for_sim)

N <- 10000
n_b <- 1000
n_a <- 500

find_sigma_regsim <- function(sigma, rho) {
  y <- 2 + x1 + x2 + x3 + x4 + sigma * epsilon
  model_fit <- lm.fit(x = cbind(1, x1, x2, x3, x4), y = y)
  res <- cor(y, model_fit$fitted.values) - rho
  res
}

find_theta_logsim <- function(theta, n) {
  eta <- theta + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4
  rho <- exp(eta) / (1 + exp(eta))
  res <- sum(rho) - n
  res
}

z1 <- rbinom(N, 1, 0.5)
z2 <- runif(N, 0, 2)
z3 <- rexp(N)
z4 <- rchisq(N, 4)
epsilon <- rnorm(N)

x0 <- rep(1, N)
x1 <- z1
x2 <- z2 + 0.3 * x1
x3 <- z3 + 0.2 * (x1 + x2)
x4 <- z4 + 0.1 * (x1 + x2 + x3)
sigma30 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.3)$root
sigma60 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.6)$root
sigma80 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.8)$root
y30 <- 2 + x1 + x2 + x3 + x4 + sigma30 * epsilon
y60 <- 2 + x1 + x2 + x3 + x4 + sigma60 * epsilon
y80 <- 2 + x1 + x2 + x3 + x4 + sigma80 * epsilon
theta <- uniroot(find_theta_logsim, lower = -100, upper = 100, n = n_b)$root
rho <- exp(theta + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4) / (1 + exp(theta + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4))

p <- uniroot(f = function(x) max(x3 + x) - 50 * min(x3 + x), lower = -200, upper = 100)$root

sim_df <- data.frame(x0, x1, x2, x3, x4, y30, y60, y80, rho, srs = (x3 + p) / 10)


A <- Sys.time()

cl <- makeCluster(cores)
clusterExport(cl, c("N"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())


res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "sampling"),
               .errorhandling = "stop",
               .options.snow = opts) %dopar% {

                 flag_non <- rbinom(N, 1, prob = rho)
                 sample_nonprob <- sim_df[flag_non == 1, ]

                 sample_prob <- sim_df[sample(1:N, n_a), ]
                 sample_prob$w_b <- N / n_a
                 svy_prob <- svydesign(ids = ~ 1, weights = ~ w_b, data = sample_prob)

                 ################ ANALYTICAL #####################

                 # MI
                 # glm
                 # y_30
                 est_mi_glm_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_mi_nn_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn"
                 )

                 est_mi_pmm_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm"
                 )

                 # y_60
                 est_mi_glm_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_mi_nn_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn"
                 )

                 est_mi_pmm_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm"
                 )

                 # y_80
                 est_mi_glm_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_mi_nn_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn"
                 )

                 est_mi_pmm_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm"
                 )

                 # IPW
                 # y_30
                 est_ipw_mle_y30 <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_ipw_gee1_y30 <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_ipw_gee2_y30 <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )

                 # y_60
                 est_ipw_mle_y60 <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_ipw_gee1_y60 <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_ipw_gee2_y60 <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )


                 # y_90
                 est_ipw_mle_y80 <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_ipw_gee1_y80 <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_ipw_gee2_y80 <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )

                 # DR
                 # y_30
                 est_dr_mle_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_dr_gee1_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_dr_gee2_y30 <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )

                 # y_60
                 est_dr_mle_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_dr_gee1_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_dr_gee2_y60 <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )


                 # y_90
                 est_dr_mle_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob
                 )

                 est_dr_gee1_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1)
                 )

                 est_dr_gee2_y80 <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2)
                 )


                 ################ BOOTSTRAP #####################

                 # MI
                 # glm
                 # y_30
                 est_mi_glm_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_nn_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_pmm_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # y_60
                 est_mi_glm_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_nn_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_pmm_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # y_80
                 est_mi_glm_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_nn_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "nn",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_mi_pmm_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   method_outcome = "pmm",
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # IPW
                 # y_30
                 est_ipw_mle_y30_boot <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee1_y30_boot <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee2_y30_boot <- nonprob(
                   target = ~ y30,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # y_60
                 est_ipw_mle_y60_boot <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee1_y60_boot <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee2_y60_boot <- nonprob(
                   target = ~ y60,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )


                 # y_90
                 est_ipw_mle_y80_boot <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee1_y80_boot <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_ipw_gee2_y80_boot <- nonprob(
                   target = ~ y80,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # DR
                 # y_30
                 est_dr_mle_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee1_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee2_y30_boot <- nonprob(
                   outcome = y30 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 # y_60
                 est_dr_mle_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee1_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee2_y60_boot <- nonprob(
                   outcome = y60 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )


                 # y_90
                 est_dr_mle_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee1_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 1),
                   control_inference = controlInf(var_method = "bootstrap")
                 )

                 est_dr_gee2_y80_boot <- nonprob(
                   outcome = y80 ~ x1 + x2 + x3 + x4,
                   selection = ~ x1 + x2 + x3 + x4,
                   data = sample_nonprob,
                   svydesign = svy_prob,
                   control_selection = controlSel(est_method_sel = "gee", h = 2),
                   control_inference = controlInf(var_method = "bootstrap")
                 )


                 data.frame(
                   k = k,
                   y = c("y30", "y60", "y80"),
                   type = rep(c("analytical", "bootstrap"), each = 3),
                   dr_mle = c(
                     est_dr_mle_y30$output$SE, est_dr_mle_y30_boot$output$SE,
                     est_dr_mle_y60$output$SE, est_dr_mle_y60_boot$output$SE,
                     est_dr_mle_y80$output$SE, est_dr_mle_y80_boot$output$SE
                   ),
                   ipw_mle = c(
                     est_ipw_mle_y30$output$SE, est_ipw_mle_y30_boot$output$SE,
                     est_ipw_mle_y60$output$SE, est_ipw_mle_y60_boot$output$SE,
                     est_ipw_mle_y80$output$SE, est_ipw_mle_y80_boot$output$SE
                   ),
                   dr_gee1 = c(
                     est_dr_gee1_y30$output$SE, est_dr_gee1_y30_boot$output$SE,
                     est_dr_gee1_y60$output$SE, est_dr_gee1_y60_boot$output$SE,
                     est_dr_gee1_y80$output$SE, est_dr_gee1_y80_boot$output$SE
                   ),
                   ipw_gee1 = c(
                     est_ipw_gee1_y30$output$SE, est_ipw_gee1_y30_boot$output$SE,
                     est_ipw_gee1_y60$output$SE, est_ipw_gee1_y60_boot$output$SE,
                     est_ipw_gee1_y80$output$SE, est_ipw_gee1_y80_boot$output$SE
                   ),
                   dr_gee2 = c(
                     est_dr_gee2_y30$output$SE, est_dr_gee2_y30_boot$output$SE,
                     est_dr_gee2_y60$output$SE, est_dr_gee2_y60_boot$output$SE,
                     est_dr_gee2_y80$output$SE, est_dr_gee2_y80_boot$output$SE
                   ),
                   ipw_gee2 = c(
                     est_ipw_gee2_y30$output$SE, est_ipw_gee2_y30_boot$output$SE,
                     est_ipw_gee2_y60$output$SE, est_ipw_gee2_y60_boot$output$SE,
                     est_ipw_gee2_y80$output$SE, est_ipw_gee2_y80_boot$output$SE
                   ),
                   mi_glm = c(
                     est_mi_glm_y30$output$SE, est_mi_glm_y30_boot$output$SE,
                     est_mi_glm_y60$output$SE, est_mi_glm_y60_boot$output$SE,
                     est_mi_glm_y80$output$SE, est_mi_glm_y80_boot$output$SE
                   ),
                   mi_nn = c(
                     est_mi_nn_y30$output$SE, est_mi_nn_y30_boot$output$SE,
                     est_mi_nn_y60$output$SE, est_mi_nn_y60_boot$output$SE,
                     est_mi_nn_y80$output$SE, est_mi_nn_y80_boot$output$SE
                   ),
                   mi_pmm = c(
                     est_mi_pmm_y30$output$SE, est_mi_pmm_y30_boot$output$SE,
                     est_mi_pmm_y60$output$SE, est_mi_pmm_y60_boot$output$SE,
                     est_mi_pmm_y80$output$SE, est_mi_pmm_y80_boot$output$SE
                   ),
                   dr_mle_ci = c(
                     est_dr_mle_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_mle_y30$confidence_interval[, 2],
                     est_dr_mle_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_mle_y60$confidence_interval[, 2],
                     est_dr_mle_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_mle_y80$confidence_interval[, 2],
                     est_dr_mle_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_mle_y30_boot$confidence_interval[, 2],
                     est_dr_mle_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_mle_y60_boot$confidence_interval[, 2],
                     est_dr_mle_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_mle_y80_boot$confidence_interval[, 2]
                   ),
                   dr_gee1_ci = c(
                     est_dr_gee1_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_gee1_y30$confidence_interval[, 2],
                     est_dr_gee1_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_gee1_y60$confidence_interval[, 2],
                     est_dr_gee1_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_gee1_y80$confidence_interval[, 2],
                     est_dr_gee1_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_gee1_y30_boot$confidence_interval[, 2],
                     est_dr_gee1_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_gee1_y60_boot$confidence_interval[, 2],
                     est_dr_gee1_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_gee1_y80_boot$confidence_interval[, 2]
                   ),
                   dr_gee2_ci = c(
                     est_dr_gee2_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_gee2_y30$confidence_interval[, 2],
                     est_dr_gee2_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_gee2_y60$confidence_interval[, 2],
                     est_dr_gee2_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_gee2_y80$confidence_interval[, 2],
                     est_dr_gee2_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_dr_gee2_y30_boot$confidence_interval[, 2],
                     est_dr_gee2_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_dr_gee2_y60_boot$confidence_interval[, 2],
                     est_dr_gee2_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_dr_gee2_y80_boot$confidence_interval[, 2]
                   ),
                   ipw_mle_ci = c(
                     est_ipw_mle_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_mle_y30$confidence_interval[, 2],
                     est_ipw_mle_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_mle_y60$confidence_interval[, 2],
                     est_ipw_mle_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_mle_y80$confidence_interval[, 2],
                     est_ipw_mle_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_mle_y30_boot$confidence_interval[, 2],
                     est_ipw_mle_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_mle_y60_boot$confidence_interval[, 2],
                     est_ipw_mle_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_mle_y80_boot$confidence_interval[, 2]
                   ),
                   ipw_gee1_ci = c(
                     est_ipw_gee1_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_gee1_y30$confidence_interval[, 2],
                     est_ipw_gee1_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_gee1_y60$confidence_interval[, 2],
                     est_ipw_gee1_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_gee1_y80$confidence_interval[, 2],
                     est_ipw_gee1_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_gee1_y30_boot$confidence_interval[, 2],
                     est_ipw_gee1_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_gee1_y60_boot$confidence_interval[, 2],
                     est_ipw_gee1_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_gee1_y80_boot$confidence_interval[, 2]
                   ),
                   ipw_gee2_ci = c(
                     est_ipw_gee2_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_gee2_y30$confidence_interval[, 2],
                     est_ipw_gee2_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_gee2_y60$confidence_interval[, 2],
                     est_ipw_gee2_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_gee2_y80$confidence_interval[, 2],
                     est_ipw_gee2_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_ipw_gee2_y30_boot$confidence_interval[, 2],
                     est_ipw_gee2_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_ipw_gee2_y60_boot$confidence_interval[, 2],
                     est_ipw_gee2_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_ipw_gee2_y80_boot$confidence_interval[, 2]
                   ),
                   mi_glm_ci = c(
                     est_mi_glm_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_glm_y30$confidence_interval[, 2],
                     est_mi_glm_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_glm_y60$confidence_interval[, 2],
                     est_mi_glm_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_glm_y80$confidence_interval[, 2],
                     est_mi_glm_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_glm_y30_boot$confidence_interval[, 2],
                     est_mi_glm_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_glm_y60_boot$confidence_interval[, 2],
                     est_mi_glm_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_glm_y80_boot$confidence_interval[, 2]
                   ),
                   mi_nn_ci = c(
                     est_mi_nn_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_nn_y30$confidence_interval[, 2],
                     est_mi_nn_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_nn_y60$confidence_interval[, 2],
                     est_mi_nn_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_nn_y80$confidence_interval[, 2],
                     est_mi_nn_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_nn_y30_boot$confidence_interval[, 2],
                     est_mi_nn_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_nn_y60_boot$confidence_interval[, 2],
                     est_mi_nn_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_nn_y80_boot$confidence_interval[, 2]
                   ),
                   mi_pmm_ci = c(
                     est_mi_pmm_y30$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_pmm_y30$confidence_interval[, 2],
                     est_mi_pmm_y60$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_pmm_y60$confidence_interval[, 2],
                     est_mi_pmm_y80$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_pmm_y80$confidence_interval[, 2],
                     est_mi_pmm_y30_boot$confidence_interval[, 1] < mean(sim_df$y30) & mean(sim_df$y30) < est_mi_pmm_y30_boot$confidence_interval[, 2],
                     est_mi_pmm_y60_boot$confidence_interval[, 1] < mean(sim_df$y60) & mean(sim_df$y60) < est_mi_pmm_y60_boot$confidence_interval[, 2],
                     est_mi_pmm_y80_boot$confidence_interval[, 1] < mean(sim_df$y80) & mean(sim_df$y80) < est_mi_pmm_y80_boot$confidence_interval[, 2]
                   )
                 )
               }

stopCluster(cl)

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:3)
# results_simulation1_process[, c("est", "var", "ci"):=tstrsplit(variable, "_")]

saveRDS(results_simulation1_process, file = "results/variance.rds")

## coverage
# tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100,
#                                                  se = sd(value)*100,
#                                                  rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100),
#                                     keyby=.(sample, y, est, var)] |>
#   melt(id.vars = c(1, 4,2,3)) |>
#   transform(sample=paste(sample, variable, sep = "_")) |>
#   transform(variable=NULL) |>
#   dcast(... ~ sample, value.var = "value") |>
#   {\(x) x[order(y, est, var)]}()
#
# tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100),
#                                     keyby=.(sample, y, est, var)]  |>
#   dcast(... ~ sample, value.var = "ci")
#
# tab1[tab2, on = c("y", "est", "var")] |>
#   setcolorder(c("y", "est", "var", "b1_bias", "b1_se", "b1_rmse", "b1", "b2_bias", "b2_se", "b2_rmse", "b2")) |>
#   {\(x) x[,y:=NULL][]}()
