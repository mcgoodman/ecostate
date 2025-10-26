
#' @title Conditional simulation from a GMRF
#'
#' @description
#' Generates samples from a Gaussian Markov random field (GMRF) conditional upon
#' fixed values for some elements. Duplicated from package \code{tinyVAST}.
#'
#' @param Q precision for a zero-centered GMRF.
#' @param observed_idx integer vector listing rows of \code{Q} corresponding to
#'        fixed measurements
#' @param x_obs numeric vector with fixed values for indices \code{observed_idx}
#' @param n_sims integer listing number of simulated values
#' @param what Whether to simulate from the conditional GMRF, or predict the mean
#'        and precision
#'
#' @return
#' A matrix with \code{n_sims} columns and a row for every row of \code{Q} not in
#' \code{observed_idx}, with simulations for those rows
#'
#' @export
conditional_gmrf <-
  function( Q,
            observed_idx,
            x_obs,
            n_sims = 1,
            what = c("simulate","predict") ){
    
    # Required libraries
    #library(Matrix)
    what = match.arg(what)
    
    # Error checks
    if( !all(observed_idx %in% seq_len(nrow(Q))) ){
      stop("Check `observed_idx` in `conditional_gmrf")
    }
    if( length(observed_idx) != length(x_obs) ){
      stop("Check length of `observed_idx` and `x_obs`")
    }
    if( any(is.na(x_obs)) ){
      stop("`x_obs` cannot include NA values")
    }
    
    # Calculate conditional mean and variance
    predict_conditional_gmrf <- function(Q, observed_idx, x_obs) {
      all_idx <- seq_len(nrow(Q))
      unobserved_idx <- setdiff(all_idx, observed_idx)
      
      # Partition Q
      #Q_oo <- Q[observed_idx, observed_idx, drop = FALSE]
      Q_ou <- Q[observed_idx, unobserved_idx, drop = FALSE]
      Q_uo <- Matrix::t(Q_ou)
      Q_uu <- Q[unobserved_idx, unobserved_idx, drop = FALSE]
      
      # Compute conditional mean and covariance
      #mu_cond <- -Q_uu_inv %*% Q_uo %*% x_obs
      mu_cond <- -1 * Matrix::solve(Q_uu, Q_uo %*% x_obs)
      
      out = list( mean = as.vector(mu_cond),
                  Q_uu = Q_uu,
                  unobserved_idx = unobserved_idx )
      return(out)
    }
    
    res <- predict_conditional_gmrf(Q, observed_idx, x_obs )
    if( what == "predict" ){
      return(res$mean)
    }else{
      y = c(res$mean + rgmrf0(n = n_sims, Q = res$Q_uu))
      return(y)
    }
  }

# Unexported
# Update first n values in atomics in list x with those in list y
substitute_pars <- function(x, y) {
  
  ind <- intersect(names(x), names(y))
  
  for (i in ind) {
    
    yd <- dim(y[[i]])
    
    if (is.null(yd)) {
      # Vectors
      x[[i]][seq_len(length(y[[i]]))] <- y[[i]]
    } else {
      # Matrices / arrays
      yind <- lapply(yd, seq_len)
      args <- c(list(x[[i]]), yind, list(value = y[[i]]))
      x[[i]] <- do.call("[<-", args)
    }
  }
  return(x)
}

#' @title Project ecostate model to future times (EXPERIMENTAL)
#'
#' @description
#' Projects a fitted model forward in time, optionally 
#' conditioning on covariates and/or catches.
#'
#' @inheritParams predict.tinyVAST
#' @param model fitted model from \code{ecostate}
#' @param nsim Number simulations
#' @param extra_years a vector of extra years, matching values in \code{catch} and 
#'        \code{covariates}
#' @param catch data.frame formatted as for \code{catch} argument in \code{ecostate},
#'        with column names "Year", "Mass", and "Taxon".
#' @param covariates matrix with rownames corresponding to \code{extra_years}
#'        and columns corresponding to SEM covariates. Covariates corresponding to 
#'        missing covariate columns and/or NA values within included covariate columns
#'        will be estimateed conditional on included covariates.
#' @param future_var logical indicating whether to simulate future process errors
#'        from GMRFs, or just compute the predictive mean
#' @param past_var logical indicating whether to re-simulate past process errors
#'        from predictive distribution of random effects, thus changing the boundary
#'        condition of the forecast
#' @param parm_var logical indicating whether to re-sample fixed effects from their
#'        predictive distribution, thus changing the GMRF for future process errors
#' @param obj_new Optional; a re-built model object with parameter blocks covering
#'        the full range of years in the simulation. Intended for internal use.
#'
#' @return
#' A list (or list of lists, if nsim > 1) of estimated parameters and 
#' derived quantities with structure resembling \code{model$internal$parhat}
project <- function(
    model, 
    nsim = 1,
    extra_years,
    catch,
    covariates, 
    future_var = TRUE, 
    past_var = FALSE, 
    parm_var = FALSE, 
    obj_new
) {
  
  years_all <- c(model$internal$years, extra_years)
  
  if (!identical(min(years_all):max(years_all), years_all)) {
    stop("extra_years must be an integer sequence starting from the year following the last year in the model.")
  }
  
  if (!missing(catch)) {
    catch_full <- rbind(model$internal$catch, catch)
  } else {
    catch_full <- model$internal$catch
  }
  
  if (isFALSE(parm_var) & isFALSE(past_var)) {
    parlist <- model$internal$parhat
  } else {
    message(paste(
      "resmapling of fixed effects and and past random effects", 
      "not yet implemented"
    ))
    parlist <- model$internal$parhat
  }
  
  # Populate / format covariates matrix
  if (!missing(covariates)) {
    
    cov_hind <- parlist$covariates
    
    if (!all(colnames(covariates) %in% colnames(cov_hind))) {
      message(paste(
        "ignoring covariates not included in the model:",
        paste(names(covariates)[!(names(covariates) %in% colnames(cov_hind))], collapse = ", ")
      ))
      covariates <- covariates[,colnames(covariates) %in% colnames(cov_hind), drop = FALSE]
    }
    
    if (!all(rownames(covariates) %in% as.character(extra_years))) {
      message("ignoring years in covariates matrix not in extra_years")
      covariates <- covariates[rownames(covariates) %in% as.character(extra_years),,drop = FALSE]
    }
    
    cov_new <- matrix(
      NA, nrow = length(extra_years), ncol = ncol(cov_hind), 
      dimnames = list(extra_years, colnames(cov_hind))
    )
    
    for (i in seq_along(colnames(covariates))) {
      cov_new[rownames(covariates),colnames(covariates)[i]] <- covariates[,colnames(covariates)[i]]
    }
    
    cov_full <- rbind(cov_hind, cov_new)
    
  } else if (length(parlist$covariates)) {
    
    cov_full <- parlist$covariates
    
  }
  
  if (missing(obj_new)) {
    # Control settings 
    ctrl_new <- model$internal$control
    ctrl_new$nlminb_loops <- ctrl_new$newton_loops <- 0
    ctrl_new$getsd <- FALSE
    
    # Some objects needed to update the model are not returned by
    # ecostate(), or may not always evaluate correctly (e.g. values from call()). 
    # Use update() to recycle most arguments by calling in the parent env., 
    # but pull specific values to use from the function environment
    local <- environment()
    years_value <- get("years_all", envir = local)
    ctrl_value <- get("ctrl_new", envir = local)
    catch_value <- get("catch_full", envir = local)
    cov_value <- get("cov_full", envir = local)
    
    # Call now contains the values for these objects
    new_call <- bquote(update(
      model, 
      years = .(years_value), 
      catch = .(catch_value), 
      covariates = .(cov_value), 
      control = .(ctrl_value)
    ))
    
    # Build
    obj_new <- eval(new_call, envir = parent.frame())
  }
  
  # Replace values in years prior to projection with fitted estimates
  par_new <- obj_new$internal$parhat 
  par_new <- substitute_pars(par_new, parlist)
  
  # Predict random effects in future years conditional on fitted estimates
  if (!is.null(model$sem)) {
    
    sem_path <- eval(model$call$sem, envir = parent.frame())
    
    # Store back to model call, to avoid problems with recursion when nsim > 1
    model$call$sem <- bquote(.(sem_path))
    
    sem <- parse_ecostate_sem(
      sem_path, covariates = cov_full, taxa = model$internal$taxa, 
      years = years_all, fit_eps = c(), fit_nu = c(), 
      settings = model$internal$settings,
      control = model$internal$control
    )$model
    
    # SEM precision matrix
    variables <- unique(c(sem$first, sem$second))
    sem_mat <- make_matrices(par_new$beta, sem, years_all, variables)
    Q <- t(sem_mat$IminusP_kk) %*% sem_mat$invV_kk %*% sem_mat$IminusP_kk
    
    # Matrix to bundle process errors and covariates
    Xit <- matrix(NA, nrow = length(years_all), ncol = length(variables))
    colnames(Xit) <- variables
    
    # Subtract estimated covariate means
    if (!is.null(dim(cov_full))) {
      Xit[,colnames(cov_full)] <- cov_full
      Xit[,colnames(cov_full)] <- sweep(Xit[,colnames(cov_full), drop = FALSE], 2, par_new$mu) 
    }
    
    taxa <- model$internal$taxa
    yr_ind <- seq_along(model$internal$years)
    
    # Pull out process error values from matrices
    for (i in seq_len(ncol(Xit))) {
      
      if (gsub("eps_", "", colnames(Xit)[i]) %in% taxa) {
        Xit[yr_ind,i] <- parlist$epsilon_ti[,which(taxa %in% gsub("eps_", "", colnames(Xit)[i]))]
      } else if (grepl("nu_", colnames(Xit)[i])) {
        if (gsub("nu_", "", colnames(Xit)[i]) %in% taxa) {
          Xit[yr_ind,i] <- parlist$nu_ti[,which(taxa %in% gsub("nu_", "", colnames(Xit)[i]))]
        } else if (all(strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]] %in% taxa)) {
          pred_prey <- strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]]
          Xit[yr_ind,i] <- parlist$nu_tij[, pred_prey[1], pred_prey[2]]
        } 
      } else if (gsub("phi_", "", colnames(Xit)[i]) %in% model$internal$settings$unique_stanza_groups) {
        Xit[yr_ind,i] <- parlist$phi_tg2[,which(model$internal$settings$unique_stanza_groups %in% gsub("phi_", "", colnames(Xit)[i]))]
      }
    }
    
    # Stack columnwise
    Xvec <- c(Xit)
    
    # Treat non-NA indices as fixed and condition the GMRF on them
    X_fixed <- which(!is.na(Xvec))
    
    # Conditional simulation from GMRF
    Xit[is.na(Xit)] <- conditional_gmrf(
      Q, 
      observed_idx = X_fixed, 
      x_obs = Xvec[X_fixed], n_sims = 1, 
      what = ifelse(isTRUE(future_var), "simulate", "predict")
    )
    
    # Add back in estimated covariate means
    if (!is.null(dim(cov_full))) {
      Xit[,colnames(cov_full)] <- sweep(Xit[,colnames(cov_full), drop = FALSE], 2, par_new$mu, FUN = "+")
    }
    
    # Unpack samples from GMRF
    for (i in seq_len(ncol(Xit))) {
      if (gsub("eps_", "", colnames(Xit)[i]) %in% taxa) {
        par_new$epsilon_ti[, which(taxa %in% gsub("eps_", "", colnames(Xit)[i]))] <- Xit[,i]
      } else if (grepl("nu_", colnames(Xit)[i])) {
        if (gsub("nu_", "", colnames(Xit)[i]) %in% taxa) {
          par_new$nu_ti[, which(taxa %in% gsub("nu_", "", colnames(Xit)[i]))] <- Xit[,i]
        } else if (all(strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]] %in% taxa)) {
          pred_prey <- strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]]
          par_new$nu_tij[, pred_prey[1], pred_prey[2]] <- Xit[,i]
        } 
      } else if (gsub("phi_", "", colnames(Xit)[i]) %in% model$internal$settings$unique_stanza_groups) {
        par_new$phi_tg2[, which(model$internal$settings$unique_stanza_groups %in% gsub("phi_", "", colnames(Xit)[i]))] <- Xit[,i]
      } else if (colnames(Xit)[i] %in% colnames(par_new$covariates)) {
        par_new$covariates[,which(colnames(par_new$covariates) == colnames(Xit)[i])] <- Xit[,i]
      }
    }
    
  } else if (isTRUE(future_var)) {
    
    # If model does not use SEM,
    # implement i.i.d normal simulation for process errors
    error("Projection not yet implemented for non-SEM models")
    
  }
  
  # Simulate
  # Random is false because random simulation handled above
  if (nsim == 1) {
    sim <- obj_new$simulator(par_new, simulate_random = FALSE)
  } else {
    # Call project() recursively but provide obj_new to avoid rebuilding
    sim <- vector("list", length = nsim)
    for (i in seq_len(nsim)) {
      sim[[i]] <- project(
        model = model, obj_new = obj_new, nsim = 1,
        extra_years = extra_years, catch = catch, covariates = covariates, 
        future_var = future_var, past_var = past_var, parm_var = parm_var
      )
    }
  }
  
  return(sim)
  
}
