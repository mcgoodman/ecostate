
library("dsem")

# Load data -----------------------------------------------

data(gulf_of_alaska)
goa <- gulf_of_alaska
#years

goa$taxa <- gsub(" ", ".", goa$taxa)
names(goa$type) <- names(goa$B) <- names(goa$P_over_B) <- names(goa$Q_over_B) <- names(goa$EE) <- goa$taxa
colnames(goa$Diet_proportions) <- rownames(goa$Diet_proportions) <- goa$taxa
goa$catch_data$Taxon <- gsub(" ", ".", goa$catch_data$Taxon)
goa$biomass_data$Taxon <- gsub(" ", ".", goa$biomass_data$Taxon)

names(goa$stanza_groups) <- names(goa$Amax) <- gsub(" ", ".", names(goa$stanza_groups))
goa$stanza_groups <- gsub(" ", ".", goa$stanza_groups)
names(goa$K) <- names(goa$d) <- names(goa$Amat) <- names(goa$Wmat) <- names(goa$Wmatslope) <- names(goa$agecomp_data) <- gsub(" ", ".", names(goa$K))
names(goa$Leading) <- gsub(" ", ".", names(goa$Leading))

attach(goa)
#years = 1959:2023

# Inputs, Parameters to estimate --------------------------

data( goa_mle )
attach(goa_mle)
SpawnX = c( "Walleye.pollock" = 4 / (5-1/0.999), "Sablefish" = 4 / (5-1/0.999) )

settings_dgmrf <- stanza_settings(
  taxa = taxa, stanza_groups = stanza_groups,
  K = K, Wmat = Wmat, Amat = Amat, d = d, Amax = Amax, SpawnX = SpawnX,
  STEPS_PER_YEAR = STEPS_PER_YEAR, comp_weight = comp_weight, Leading = Leading,
  Wmatslope = Wmatslope
)

control <- ecostate_control(
  nlminb_loops = 0,
  n_steps = n_step,
  profile = NULL,
  random = NULL,
  derived_quantities = c(),
  getsd = FALSE,
)

test_that(
  "basic stanzas work", {

  # Fix SD of biomass and recruitment errors (penalized likelihood)
  sem = "
    eps_Euphausiids <-> eps_Euphausiids, 0, NA, 1,
    eps_Large.copepods <-> eps_Large.copepods, 0, NA, 1,
    phi_Walleye.pollock <-> phi_Walleye.pollock, 0, NA, 1,
    phi_Sablefish <-> phi_Sablefish, 0, NA, 1
  "

  fit <- ecostate(
    taxa = taxa, years = years, type = type,
    catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
    PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
    fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem,
    log_prior = log_prior, settings = settings_dgmrf, control = control
  )

  obj = fit$obj
  start_val = obj$fn( goa_mle$mle_opt$par )

  ## match objective function
  expect_equal(goa_mle$mle_opt$obj, start_val, tolerance = 1e-4)

})


test_that("epsilon on multi-stanza taxa behave correctly", {

  stanza_taxa <- c("Walleye.pollock.juv", "Walleye.pollock.adult")

  # Biomass process errors for pollock stanzas
  sem = "
    phi_Walleye.pollock <-> phi_Walleye.pollock, 0, NA, 1,
    phi_Sablefish <-> phi_Sablefish, 0, NA, 1
    eps_Walleye.pollock.juv <-> eps_Walleye.pollock.juv, 0, NA, 0.1
    eps_Walleye.pollock.adult <-> eps_Walleye.pollock.adult, 0, NA, 0.1
  "

  fit <- ecostate(
    taxa = taxa, years = years, type = type,
    catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
    PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
    fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem,
    log_prior = log_prior, settings = settings_dgmrf, control = control
  )

  obj <- fit$obj
  i_eps <- which(names(obj$par) == "epsilon_ti")
  eps_index <- matrix(NA, nrow = length(fit$internal$years), ncol = length(taxa), dimnames = list(NULL, taxa))
  eps_index[, stanza_taxa] <- i_eps
  rep0 <- obj$report(obj$par)

  t <- 20
  for (tx in stanza_taxa) {

    par <- obj$par
    par[eps_index[t, tx]] <- 0.2
    rep1 <- obj$report(par)
    i <- match(tx, taxa)
    s2 <- match(i, rep1$stanza_data$stanzainfo_s2z[, "s"])

    # Positive biomass process error raises the taxon's biomass
    expect_gt(log(rep1$B_ti[t, i] / rep0$B_ti[t, i]), 0.05)
    # by lowering its total mortality within the stanza update
    expect_lt(rep1$TotalZ_ts2[t, s2], rep0$TotalZ_ts2[t, s2])

  }

})