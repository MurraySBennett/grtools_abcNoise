
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
#'   If "differential" is chosen, the second element should specify the noise ratio.
#'   Defaults to \code{list(c("none"))} if no argument is provided.
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
# Function to extract the best model based on AIC
extract_best_model <- function(fitted_models, ordered_table) {
  best_model_index <- which(names(fitted_models) == ordered_table$model[1])
  if (length(best_model_index) > 0) {
    return(fitted_models[[best_model_index[1]]])
  } else {
    return(NULL)
  }
}

# Function to order models by AIC
order_aic <- function(fitted_models) {
  aic_values <- sapply(fitted_models, function(model) {
    if (!is.null(model)) {
      return(2 * length(model$par) + 2 * model$value)
    } else {
      return(Inf)
    }
  })
  model_names <- names(fitted_models)
  ordered_table <- data.frame(model = model_names, AIC = aic_values)
  ordered_table <- ordered_table[order(ordered_table$AIC), ]
  return(ordered_table)
}

# Function to fit GRT models
fit_grt_models <- function(cmat, rand_pert = 0.3, n_reps = 10, control = list(), noise_model = NULL, alpha = 0) {
  fitted_models <- list()
  model_list <- c(
    "{PI, PS, DS}",
    "{PI, PS(A), DS}",
    "{PI, PS(B), DS}",
    "{1_RHO, PS, DS}",
    "{1_RHO, PS(A), DS}",
    "{PI, DS}",
    "{1_RHO, PS(B), DS}",
    "{PS, DS}",
    "{PS(A), DS}",
    "{1_RHO, DS}",
    "{PS(B), DS}",
    "{DS}"
  )
  for (i in 1:length(model_list)) {
    model_name <- model_list[i]
    negloglik_func <- switch(model_name,
      "{PI, PS, DS}" = negloglik_mod1,
      "{PI, PS(A), DS}" = negloglik_mod2,
      "{PI, PS(B), DS}" = negloglik_mod3,
      "{1_RHO, PS, DS}" = negloglik_mod4,
      "{1_RHO, PS(A), DS}" = negloglik_mod5,
      "{PI, DS}" = negloglik_mod6,
      "{1_RHO, PS(B), DS}" = negloglik_mod7,
      "{PS, DS}" = negloglik_mod8,
      "{PS(A), DS}" = negloglik_mod9,
      "{1_RHO, DS}" = negloglik_mod10,
      "{PS(B), DS}" = negloglik_mod11,
      "{DS}" = negloglik_mod12
    )

    # Define starting parameters, lower bounds, and upper bounds based on the model
    if (model_name == "{PI, PS, DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0)
      high_params <- c(1, 1, 1, 1)
    } else if (model_name == "{PI, PS(A), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, 0, 0)
      high_params <- c(1, 1, 1, 1, 1, 1)
    } else if (model_name == "{PI, PS(B), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, 0, 0)
      high_params <- c(1, 1, 1, 1, 1, 1)
    } else if (model_name == "{1_RHO, PS, DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, -1, 0)
      high_params <- c(1, 1, 1, 1)
    }else if (model_name == "{1_RHO, PS(A), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, -1)
      high_params <- c(1, 1, 1, 1, 1)
    }else if (model_name == "{PI, DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, 0, 0)
      high_params <- c(1, 1, 1, 1, 1, 1)
    }else if (model_name == "{1_RHO, PS(B), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, -1)
      high_params <- c(1, 1, 1, 1, 1)
    }else if (model_name == "{PS, DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, -1, -1, -1, -1)
      high_params <- c(1, 1, 1, 1, 1, 1)
    }else if (model_name == "{PS(A), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, -1, -1, -1, -1)
      high_params <- c(1, 1, 1, 1, 1, 1, 1, 1)
    }else if (model_name == "{1_RHO, DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, 0, 0, -1)
      high_params <- c(1, 1, 1, 1, 1, 1, 1)
    }else if (model_name == "{PS(B), DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, -1, -1)
      high_params <- c(1, 1, 1, 1, 1, 1)
    }else if (model_name == "{DS}") {
      start_params <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5)
      low_params <- c(0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 0, 0)
      high_params <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    }
    
    best_fit <- NULL
    best_value <- Inf
    
    for (j in 1:n_reps) {
      random_params <- rand_start(start_params, low_params, high_params, rand_pert)
      fit <- optim(par = random_params, fn = negloglik_func, data = pmatrix(cmat), control = control, alpha = alpha, noise_model = noise_model)
      
      if (fit$value < best_value) {
        best_value <- fit$value
        best_fit <- fit
      }
    }
    fitted_models[[model_name]] <- best_fit
  }
  return(fitted_models)
}

grt_hm_fit <- function(
    cmat, rand_pert = 0.3, n_reps = 10, control = list(),
    noise_models = NULL, alpha = 0
) {
  if (is.null(noise_models)) {
    noise_models <- list(c("none"))
  }
  
  results_list <- list()
  
  for (noise_model in noise_models) {
    print(paste0("Fitting GRT with ", noise_model, " noise."))
    fitted_models <- fit_grt_models(cmat, rand_pert = rand_pert, n_reps = n_reps, control = control, noise_model = noise_model, alpha = alpha)
    ordered_table <- order_aic(fitted_models)
    best_model <- extract_best_model(fitted_models, ordered_table)
    
    # Construct the best_model list with means, covmat, etc.
    if (!is.null(best_model)) {
      o_match <- ordered_table[1, ]$model
      model_list <- c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", 
                      "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", "{PI, DS}", 
                      "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", 
                      "{1_RHO, DS}", "{PS(B), DS}", "{DS}")
      model_num <- pmatch(o_match, model_list)
      w <- best_model$par
      
      best_model_list <- list()
      best_model_list$means <- matrix(0, 4, 2, byrow = TRUE)
      best_model_list$covmat <- list()
      best_model_list$a1 <- 0
      best_model_list$a2 <- 0
      best_model_list$model <- paste("GRT-", o_match, sep = "")
      
      switch(model_num,
             mod1 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- diag(2)
               best_model_list$covmat[[2]] <- diag(2)
               best_model_list$covmat[[3]] <- diag(2)
               best_model_list$covmat[[4]] <- diag(2)
               best_model_list$a1 <- w[3]
               best_model_list$a2 <- w[4]
             },
             mod2 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 2] <- w[3]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[4]
               best_model_list$covmat[[1]] <- diag(2)
               best_model_list$covmat[[2]] <- diag(2)
               best_model_list$covmat[[3]] <- diag(2)
               best_model_list$covmat[[4]] <- diag(2)
               best_model_list$a1 <- w[5]
               best_model_list$a2 <- w[6]
             },
             mod3 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[4, 1] <- w[4]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- diag(2)
               best_model_list$covmat[[2]] <- diag(2)
               best_model_list$covmat[[3]] <- diag(2)
               best_model_list$covmat[[4]] <- diag(2)
               best_model_list$a1 <- w[5]
               best_model_list$a2 <- w[6]
             },
             mod4 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[4]
               best_model_list$a2 <- w[5]
             },
             mod5 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 2] <- w[3]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[4]
               best_model_list$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[6]
               best_model_list$a2 <- w[7]
             },
             mod6 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[3, 2] <- w[4]
               best_model_list$means[4, 1] <- w[5]
               best_model_list$means[4, 2] <- w[6]
               best_model_list$covmat[[1]] <- diag(2)
               best_model_list$covmat[[2]] <- diag(2)
               best_model_list$covmat[[3]] <- diag(2)
               best_model_list$covmat[[4]] <- diag(2)
               best_model_list$a1 <- w[7]
               best_model_list$a2 <- w[8]
             },
             mod7 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[4, 1] <- w[4]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[9]
               best_model_list$a2 <- w[10]
             },
             mod8 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- matrix(c(1, w[3], w[3], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[4], w[4], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[7]
               best_model_list$a2 <- w[8]
             },
             mod9 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 2] <- w[3]
               best_model_list$means[4, 1] <- w[1]
               best_model_list$means[4, 2] <- w[4]
               best_model_list$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[9]
               best_model_list$a2 <- w[10]
             },
             mod10 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[3, 2] <- w[4]
               best_model_list$means[4, 1] <- w[5]
               best_model_list$means[4, 2] <- w[6]
               best_model_list$covmat[[1]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[8]
               best_model_list$a2 <- w[9]
             },
             mod11 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[3, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[4, 1] <- w[4]
               best_model_list$means[4, 2] <- w[2]
               best_model_list$covmat[[1]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[5], w[5], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[6], w[6], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[7]
               best_model_list$a2 <- w[8]
             },
             mod12 = {
               best_model_list$means[2, 1] <- w[1]
               best_model_list$means[2, 2] <- w[2]
               best_model_list$means[3, 1] <- w[3]
               best_model_list$means[3, 2] <- w[4]
               best_model_list$means[4, 1] <- w[5]
               best_model_list$means[4, 2] <- w[6]
               best_model_list$covmat[[1]] <- matrix(c(1, w[7], w[7], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[2]] <- matrix(c(1, w[8], w[8], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[3]] <- matrix(c(1, w[9], w[9], 1), 2, 2, byrow = TRUE)
               best_model_list$covmat[[4]] <- matrix(c(1, w[10], w[10], 1), 2, 2, byrow = TRUE)
               best_model_list$a1 <- w[11]
               best_model_list$a2 <- w[12]
             }
               )
               best_model_list$convergence <- best_model$convergence
               best_model_list$message <- best_model$message
               best_model_list$predicted <- as.vector(matrix_predict(
                 best_model_list$means, best_model_list$covmat, diag(2), 
                 matrix(c(best_model_list$a1, best_model_list$a2), 2, 1)
               ))
               best_model_list$observed <- as.vector(pmatrix(cmat))
               
               best_model <- best_model_list
             }
             
             if (noise_model[1] == "none") {
               results_list[["GRT-only"]] <- list(
                 table = ordered_table,
                 best_model = best_model
               )
             } else {
               results_list[[paste0("GRT-", paste(noise_model, collapse = "-"))]] <- list(
                 table = ordered_table,
                 best_model = best_model
               )
             }
    }
    
    class(results_list) <- "grt_hm_fit"
    return(results_list)
  }