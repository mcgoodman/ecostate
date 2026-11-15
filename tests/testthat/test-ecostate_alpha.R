# EBS test inputs (subset of taxa for speed) ---------------

data("eastern_bering_sea", package = "ecostate")
ebs <- eastern_bering_sea

taxa_all <- names(ebs$P_over_B)
years <- sort(unique(ebs$Survey$Year))
type <- sapply(taxa_all, switch, "Detritus" = "detritus", "Chloro" = "auto", "hetero")

EE <- B <- U <- setNames(rep(NA, length(taxa_all)), taxa_all)
B[] <- aggregate(ebs$Survey$Mass, list(Taxon = factor(ebs$Survey$Taxon, levels = taxa_all)), mean, drop = FALSE)[,2]
EE[] <- 1
U[] <- 0.2

X <- array(2, dim = c(length(taxa_all), length(taxa_all)), dimnames = list(taxa_all, taxa_all))
X[,'Chloro'] <- 91

taxa <- c("Pollock", "Krill", "Arrowtooth", "Cod")

fit_alpha_model <- function(control, ...) {
  ecostate(
    taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey,
    PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions,
    EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
    control = control, ...
  )
}

# Tests ---------------------------------------------------

test_that("alpha process errors are penalized through derived epsilon", {

  fit <- fit_alpha_model(
    control = ecostate_control(
      process_error = "alpha", random = NULL, profile = NULL,
      nlminb_loops = 0, getsd = FALSE, start_tau = 0.1
    ),
    fit_eps = "Pollock"
  )
  obj <- fit$obj

  # alpha_ti[1,] is mapped off because initial biomass is governed by delta_i
  i_alpha <- which(names(obj$par) == "alpha_ti")
  expect_length(i_alpha, length(years) - 1)

  set.seed(101)
  par <- obj$par
  par[i_alpha] <- rnorm(length(i_alpha), mean = 0, sd = 0.05)
  rep0 <- obj$report(obj$par)
  rep1 <- obj$report(par)

  # Penalty is evaluated on the derived deviation log(B / Bhat) for t >= 2
  i <- match("Pollock", taxa)
  tau <- exp(fit$internal$parhat$logtau_i[["Pollock"]])
  eps_t <- log(rep1$B_ti[-1, i] / rep1$Bhat_ti[-1, i])
  expect_equal(rep1$loglik2_ti[-1, i], dnorm(eps_t, 0, tau, log = TRUE))
  expect_equal(rep1$loglik2_ti[1, ], rep(0, length(taxa)))

  # ... so it responds to alpha, rather than being a constant dnorm(0, 0, tau)
  expect_false(isTRUE(all.equal(sum(rep0$loglik2_ti), sum(rep1$loglik2_ti))))
  expect_false(all(rep1$loglik2_ti[-1, i] == dnorm(0, 0, tau, log = TRUE)))
})

test_that("dev_penalty uses derived epsilon under alpha", {

  fit <- fit_alpha_model(
    control = ecostate_control(
      process_error = "alpha", random = NULL, profile = NULL,
      nlminb_loops = 0, getsd = FALSE, start_tau = 0.1, dev_penalty = TRUE
    ),
    fit_eps = "Pollock"
  )
  obj <- fit$obj

  set.seed(101)
  par <- obj$par
  i_alpha <- which(names(obj$par) == "alpha_ti")
  par[i_alpha] <- rnorm(length(i_alpha), mean = 0, sd = 0.05)
  rep <- obj$report(par)

  i <- match("Pollock", taxa)
  eps_t <- log(rep$B_ti[-1, i] / rep$Bhat_ti[-1, i])
  expect_equal(rep$dev_penalty, sum(eps_t)^2)
})

test_that("alpha simulator draws new process errors", {

  fit <- fit_alpha_model(
    control = ecostate_control(
      process_error = "alpha", random = NULL, profile = NULL,
      nlminb_loops = 0, getsd = FALSE, start_tau = 0.1
    ),
    fit_eps = "Pollock"
  )
  i <- match("Pollock", taxa)

  set.seed(1)
  sim1 <- fit$simulator(simulate_random = TRUE)
  set.seed(2)
  sim2 <- fit$simulator(simulate_random = TRUE)
  sim_fixed <- fit$simulator(simulate_random = FALSE)

  # Fixed random effects reproduce the fitted biomass
  expect_equal(sim_fixed$B_ti, fit$rep$B_ti)
  # Simulated random effects change process errors and biomass for Pollock
  expect_false(isTRUE(all.equal(sim1$epsilon_ti[-1, i], sim2$epsilon_ti[-1, i])))
  expect_false(isTRUE(all.equal(sim1$B_ti[, i], fit$rep$B_ti[, i])))
  # Taxa without process error keep zero deviations
  expect_true(all(sim1$epsilon_ti[, -i] == 0))
})

test_that("SEM forces epsilon process errors without modifying caller's control", {

  control <- ecostate_control(process_error = "alpha", nlminb_loops = 0, getsd = FALSE)
  sem <- "eps_Pollock <-> eps_Pollock, 0, tau_plk, 0.5"

  expect_warning(
    fit <- fit_alpha_model(control = control, sem = sem),
    "SEM forces innovation"
  )

  expect_equal(fit$internal$control$process_error, "epsilon")
  expect_equal(control$process_error, "alpha")
  expect_true(is.finite(fit$obj$fn(fit$obj$par)))
})
