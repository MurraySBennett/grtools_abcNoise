
#' Fit a hierarchy of traditional GRT models to identification data
#' 
#' Fits a hierarchy of traditional GRT models to data from a 2x2 identification 
#' experiment, using the BFGS optimization method (See Ashby & Soto, 2015). It
#' then selects the best-fitting model using the AIC. 
#' 
#' @param cmat A 4x4 confusion matrix (see Details).
#' @param rand_pert Maximum value of a random perturbation added to the starting
#'   parameters. Defaults to 0.3. With a value of zero, the optimization is started exactly at the 
#'   default starting parameters (see Details). As the value of \code{rand_pert} is 
#'   increased, the starting parameters become closer to be "truly random."
#' @param n_reps Number of times the optimization algorithm should be run, each time
#' with a different value for the starting parameters. The function will return the
#' model with a highest log-likelihood from all the runs. The value of \code{n_reps}
#' defaults to ten.
#' @param control A list of optional control parameters for the \code{optim} function. See 
#'   \code{\link[stats]{optim}}. Note that the parameter \code{ndeps} entered 
#'   here should be a single number instead of the vector that is usually passed 
#'   to \code{optim}. This single value is repeated inside \code{grt_hm_fit} to 
#'   create the appropriate vectors.
#' @param noise_models A list of noise model specifications. Each element of the list should be a character vector. 
#'   The first element specifies the noise model type ("none", "uniform", or "differential"). 
#'   If "differential" is chosen, the second element should specify the noise ratio. Defaults to \code{list(c("none"), c("uniform"), c("differential", 2))}.
#' @param alpha A numeric value between 0 and 1 representing the influence of the noise model on the overall likelihood. 
#'   0 indicates no influence of the noise model (GRT model only), and 1 indicates noise model only. Defaults to 0.
#'
#' @return An object of class "\code{grt_hm_fit}."
#'   
#'   The function \code{summary} is used to obtain a summary of results from the
#'   model fit and selection process, including the best-fitting model and
#'   conclusions about perceptual separability and perceptual independence
#'   (decisional separability is assumed by all models)
#'   
#'   The function \code{\link[=plot.grt_hm_fit]{plot}} is used to print a
#'   graphical representation of the best-fitting model.
#'   
#' @details A 2x2 identification experiment involves two dimensions, A and B,
#' each with two levels, 1 and 2. Stimuli are represented by their level in each
#' dimension (A1B1, A1B2, A2B1, and A2B2) and so are their corresponding correct
#' identification responses (a1b1, a1b2, a2b1, and a2b2).
#' 
#' The data from a single participant in the experiment should be ordered in a 
#' 4x4 confusion matrix with rows representing stimuli and columns representing 
#' responses. Each cell has the frequency of responses for the stimulus/response
#' pair. Rows and columns should be ordered in the following way:
#' 
#' \itemize{ \item{Row 1: Stimulus A1B1} \item{Row 2: Stimulus A2B1} 
#' \item{Row 3: Stimulus A1B2} \item{Row 4: Stimulus A2B2} \item{Column
#' 1: Response a1b1} \item{Column 2: Response a2b1} \item{Column 3: Response a1b2} 
#' \item{Column 4: Response a2b2} }
#' 
#' The default starting parameters for the optimization algorithm are the
#' following: \itemize{ \item{Means:}{ A1B1=(0,0), A2B1=(1,0), A1B2=(1,0),
#' A2B1=(1,1)} \item{Variances:}{ All set to one} \item{Correlations:}{ All set
#' to zero} }
#' 
#' Decisional separability is assumed for all models (i.e., decision bounds are
#' fixed and orthogonal to the dimension they divide)
#' 
#' Note that a random value will be added to the default starting parameters if 
#' \code{rand_pert} is given a value higher than zero.
#' 
#' @references Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal
#'   detection theory. In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A.
#'   Eidels (Eds.), \emph{Oxford handbook of computational and mathematical
#'   psychology} (pp. 13-34). Oxford University Press: New York, NY.
#'   
#' @examples 
#' # Create a confusion matrix
#' # Inside the c(...) below, we enter the data from row 1 in the 
#' # matrix, then from row 2, etc.
#' cmat <- matrix(c(140, 36, 34, 40,
#'                  89, 91, 4, 66,
#'                  85, 5, 90, 70,
#'                  20, 59, 8, 163),
#'                  nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Perform model fit and selection
#' hm_fit_results <- grt_hm_fit(cmat)
#' 
#' # See a summary of the fitting and selection results
#' summary(hm_fit_results)
#' 
#' # plot a graphical representation of the best-fitting model
#' plot(hm_fit_results)
#' 
#' @export
#' 
grt_hm_fit <- function(
    cmat, rand_pert=0.3, n_reps=10, control=list(),
    noise_models = list(
      c("none"), c("uniform"), c("differential", 2),
      alpha = 0
    )
  ) {
  results_list <- list()
  for (noise_model in noise_models) {
    fitted_models <- fit_grt_models(cmat, rand_pert = rand_pert, n_reps = n_reps, control = control, noise_model = noise_model, alpha = alpha) 
    if (noise_model[1] == "none"){
      results_list[[paste0("GRT-only")]] <- list(
        table = order_aic(fitted_models),
        best_model = extract_best_model(fitted_models, order_aic(fitted_models))
      )
    } else {
      results_list[[paste0("GRT-", paste(noise_model, collapse = "-"))]] <- list(
        table = order_aic(fitted_models),
        best_model = extract_best_model(fitted_models, order_aic(fitted_models))
      )
    }
  }
  
  class(results_list) <- "grt_hm_fit"
  return(results_list)
}

########################################################
# Function that actually performs maximum-likelihood estimation
# for all models:

fit_grt_models <- function(cmat, rand_pert=0, n_reps=1, control=control, noise_model = "none", alpha=0){
  
  # if ndeps was not selected by the user, assign a default value
  if (is.null(control[["ndeps"]])) {
    control$ndeps <- 1e-1
  }
  
  create_objective <- function(fn) {
    if (noise_model[1] == "none") {
      return(function(par) fn(par, cmat, alpha = alpha, noise_model = noise_model))
    } else {
      return(function(par) fn(par, cmat, alpha = alpha, noise_model = noise_model))    
    }
  }
  
  # fit model 1
  start_params <- c(1, 1, -.5, -.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf, -Inf, -Inf, -Inf)
  high_params <- c(Inf, Inf, Inf, Inf)
  min_nll <- Inf
  for(i in 1:n_reps) {
    init_par <- rand_start(start_params, low_params, high_params, rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod1), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model1 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 2
  start_params <- c(1,0,1,1,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod2), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model2 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 3
  start_params <- c(1,0,1,1,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod3), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model3 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 4
  start_params <- c(1,1,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod4), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model4 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 5
  start_params <- c(1,0,1,1,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod5), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model5 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 6
  start_params <- c(1,0,0,1,1,1,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod6), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model6 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 7
  start_params <- c(1,0,1,1,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod7), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model7 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 8
  start_params <- c(1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod8), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model8 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 9
  start_params <- c(1,0,1,1,0,0,0,0,-0.5,-0.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod9), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model9 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 10
  start_params <- c(1,0,0,1,1,1,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod10), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model10 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 11
  start_params <- c(1,0,1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf)
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod11), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model11 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  # fit model 12
  start_params <- c(1,0,0,1,1,1,0,0,0,0,-.5,-.5)
  ctrl <- control
  ctrl$ndeps <- rep(ctrl$ndeps, times=length(start_params))
  low_params <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-1,-1,-1,-1,-Inf,-Inf)
  high_params <- c(Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,Inf,Inf) 
  min_nll <- Inf
  for(i in 1:n_reps){
    init_par <- rand_start(start_params,low_params,high_params,rand_pert)
    candidate <- optim(par=init_par, fn=create_objective(negloglik_mod12), data=cmat, 
                       method="L-BFGS-B", lower=low_params, upper=high_params, 
                       control=ctrl)
    if(candidate$value < min_nll) {
      mle_model12 <- candidate
      min_nll <- candidate$value
    }
  } 
  
  fitted_models <- list(mle_model1, mle_model2, mle_model3, mle_model4, mle_model5, 
                     mle_model6, mle_model7, mle_model8, mle_model9, mle_model10,
                     mle_model11, mle_model12)
  return(fitted_models)
}


extract_best_model <- function(fitted_models, ordered_aic) {
  best_model_index <- which(sapply(fitted_models, function(model) identical(model$value, -ordered_aic$`log-likelihood`[1])))
  return(fitted_models[[best_model_index]])
}


##############################################
# Function that computes AIC and ranks models according to it:

order_aic <-function(fitted_models) {
  aic_list <- rep(0,12)
  L <- rep(0,12)
  model <- c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
  for(i in 1:12){
    L[i] <- -fitted_models[[i]]$value
    m <- length(fitted_models[[i]]$par)
    aic_list[i] <- -2*L[i]+2*m+(2*m^2+2*m)/(16-m-1)
  }
  
  aic_exp <- rep(0,12)    
  aic_weight <- rep(0,12)
  aic_exp <- exp(-(aic_list-min(aic_list))/2)
  aic_weight <- aic_exp/sum(aic_exp)
  
  ordered_aic <- data.frame(model,L,aic_list,aic_weight)
  colnames(ordered_aic) <- c("model","log-likelihood", "AIC", "AIC weight")
  ordered_aic <- ordered_aic[order(-aic_weight),]
  ordered_aic[4] <- prettyNum(round(ordered_aic[4], digits=3)[[1]], nsmall=2)
  return(ordered_aic)
}
