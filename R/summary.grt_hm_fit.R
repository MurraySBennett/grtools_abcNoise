#' @export
summary.grt_hm_fit <- function(hm_list) {
  
  # Find the best model across all noise models
  best_model_info <- NULL
  best_aic <- Inf
  best_model_name <- NULL
  
  for (noise_model_name in names(hm_list)) {
    current_aic <- hm_list[[noise_model_name]]$table[1, "AIC"] # Assuming AIC is in the "AIC" column
    if (current_aic < best_aic) {
      best_aic <- current_aic
      best_model_info <- hm_list[[noise_model_name]]$best_model
      best_model_name <- noise_model_name
    }
  }
  
  if (!is.null(best_model_info) && !is.null(best_model_info$convergence)) {
    if (best_model_info$convergence == 0) {
      cat("The optimization algorithm was successful\n\n")
    } else {
      cat("The optimization algorithm may have failed\n")
      cat("The following message was produced by the optim() function:\n")
      cat(paste("\t", best_model_info$message, "\n\n"))
    }
  } else {
    cat("Could not determine optimization convergence.\n\n")
  }
  
  # get check table to print conclusions
  check_table <- data.frame( model=c("{PI, PS, DS}", "{PI, PS(A), DS}", "{PI, PS(B), DS}", "{1_RHO, PS, DS}", "{1_RHO, PS(A), DS}", 
                                     "{PI, DS}", "{1_RHO, PS(B), DS}", "{PS, DS}", "{PS(A), DS}", "{1_RHO, DS}", "{PS(B), DS}", "{DS}") )
  check_table$PS_A <- c("yes", "yes", "no", "yes", "yes", "no", "no", "yes", "yes", "no", "no", "no")
  check_table$PS_B <- c("yes", "no", "yes", "yes", "no", "no", "yes", "yes", "no", "no", "yes", "no")
  check_table$PI <- c("yes", "yes", "yes", "no", "no", "yes", "no", "no", "no", "no", "no", "no" )
  
  cat("Summary of measures of fit for all models:\n")
  cat("(Models are ranked according to AIC)\n\n")
  
  for (noise_model_name in names(hm_list)) {
    cat(paste("Noise Model: ", noise_model_name, "\n"))
    names(hm_list[[noise_model_name]]$table) <- c("Model", "  Log-likelihood", "      AIC", "AIC weight")
    print(hm_list[[noise_model_name]]$table, row.names=F)
    cat("\n")
  }
  
  if (!is.null(best_model_info)){
    best_model_str <- best_model_info$par #get the best model parameters.
    for (noise_model_name in names(hm_list)){
      table = hm_list[[noise_model_name]]$table
      for (i in 1:nrow(table)){
        if(identical(best_model_info$par, hm_list[[noise_model_name]]$best_model$par)){
          best_model_str = table[i,1]
          break
        }
      }
    }
    cat(paste("Best fitting model:", best_model_str, "\n"))
    cat(paste("Perceptual Separability of A?:", check_table[check_table$model==best_model_str, 2], "\n"))
    cat(paste("Perceptual Separability of B?:", check_table[check_table$model==best_model_str, 3], "\n"))
    cat(paste("Perceptual Independence?:", check_table[check_table$model==best_model_str, 4], "\n"))
  }
}