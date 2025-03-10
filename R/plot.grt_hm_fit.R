#' Plot a \code{grt_hm_fit} object
#' 
#' Plot the object returned by \code{\link{grt_hm_fit}}
#' 
#' @param model A \code{grt_hm_fit} object
#' @param labels Optional names for the labels of dimensions A and B
#' @param ellipse_width Parameter controlling the width of the drawn ellipses
#' @export

plot.grt_hm_fit <- function(
    model, labels = c("dim A", "dim B"),
    marginals = TRUE, scatter = TRUE,
    show_assumptions = FALSE, show_labels = TRUE, subtitle = NA,
    ellipse_width = 0.8
) {
    if (inherits(model, "grt_hm_fit")) {
        n_models <- length(model)
        par(mfrow = c(1, n_models)) # Arrange plots in a row

        for (noise_model_name in names(model)) {
            plot_single_model(model[[noise_model_name]]$best_model, labels, marginals, scatter, show_assumptions, show_labels, noise_model_name, ellipse_width, data = model[[noise_model_name]]$best_model$data, noise_model_name = noise_model_name) # Pass noise_model_name
        }
    } else if (is.list(model)) {
        # Single model: plot directly
        noise_model_name <- names(model)
        plot_single_model(model, labels, marginals, scatter, show_assumptions, show_labels, subtitle, ellipse_width, data = model$data, noise_model_name = noise_model_name) # Pass noise_model_name
    } else {
        stop("Input 'model' must be a grt_hm_fit object or a single model list.")
    }
}

# Helper function to plot a single model
plot_single_model <- function(model, labels, marginals, scatter, show_assumptions, show_labels, subtitle, ellipse_width, data, noise_model_name = NULL) {
    plot.new()

    if (show_assumptions) {
        held_assumptions <- c(
            grepl("PS\\(A\\))", model$model),
            grepl("PS\\(B\\))", model$model),
            grepl("PI", model$model)
        )
        fig_title <- ""
        if (any(held_assumptions)) {
            assumptions <- c(
                paste0("PS(", labels[1], ")"),
                paste0("PS(", labels[2], ")"),
                "PI"
            )
            counter <- 1
            for (held in held_assumptions) {
                if (held) {
                    held_assumptions[counter] <- FALSE
                    fig_title <- paste0(
                        fig_title,
                        assumptions[counter]
                    )
                    if (any(held_assumptions)) {
                        fig_title <- paste0(fig_title, " | ")
                    }
                }
                counter <- counter + 1
            }
        }
    }

    # Extract noise model from noise_model_name
    # if (!is.null(noise_model_name) && noise_model_name != "GRT-only") {
    #     noise_model_str <- sub("^GRT-", "", noise_model_name)
    #     if (is.na(subtitle) || subtitle == "") {
    #         subtitle <- paste("Noise Model:", noise_model_str)
    #     } else {
    #         subtitle <- paste(subtitle, "| Noise Model:", noise_model_str)
    #     }
    # }
  
  # first plot the main panel
  if (!marginals && !scatter) {
    par(
      mar = c(5, 5, 5, 5) + 0.1,
      fig = c(0, 1, 0, 1),
      lty = 1, new = TRUE
    )
  } else {
    par(mar = c(4, 4, 2, 2) + 0.1, fig = c(0.2, 1, 0.2, 1), lty = 1, new = TRUE)
  }
  
# get range of values
  buffer <- 1.5
  ranx <- c(min(model$means[, 1] - buffer), max(model$means[, 1] + buffer))
  rany <- c(min(model$means[, 2] - buffer), max(model$means[, 2] + buffer))
  
  range_x <- ranx[2] - ranx[1]
  range_y <- rany[2] - rany[1]
  if (range_x > range_y) {
    diff <- range_x - range_y
    rany[1] <- rany[1] - diff / 2
    rany[2] <- rany[2] + diff / 2
  } else if (range_y > range_x) {
    diff <- range_y - range_x
    ranx[1] <- ranx[1] - diff / 2
    ranx[2] <- ranx[2] + diff / 2
  }
  
  # draw model$means of distributions
  plot(model$means[, 1], model$means[, 2],
       pch = 3, xlim = ranx, ylim = rany,
       xlab = ifelse(show_labels, labels[1], ""),
       ylab = ifelse(show_labels, labels[2], ""),
       cex.lab = ifelse(marginals, 1.5, 2)
  )
  
  if (show_assumptions) {
    title(
      main = fig_title,
      sub = ifelse(!is.na(subtitle), subtitle, ""),
      cex.main = 2
    )
  }
  
  # draw contours of distributions
  ellipse <- function(s, t) {
    u <- c(s, t) - center
    u %*% sigma.inv %*% u / 2
  }
  n <- 200
  x <- 1:200 / 10 - 10
  y <- 1:200 / 10 - 10
  for (i in 1:4) {
    center <- model$means[i, ]
    sigma.inv <- solve(model$covmat[[i]])
    z <- mapply(ellipse, as.vector(rep(x, n)), as.vector(outer(rep(0, n), y, `+`)))
    contour(x, y, matrix(z, n, n), levels = ellipse_width, drawlabels = FALSE, add = TRUE)
  }
  
  if (marginals) {
    # add marginal distributions at the bottom
    par(mar = c(1, 3.7, 1, 1.7) + 0.1, fig = c(0.2, 1, 0, 0.2), new = TRUE)
    for (i in 1:4) {
      x <- 1:100 * (ranx[2] - ranx[1]) / 100 + ranx[1]
      y <- dnorm(x, mean = model$means[i, 1], sd = sqrt(model$covmat[[i]][1, 1]))
      if (i > 2) {
        par(new = TRUE, lty = 2)
      } else {
        par(new = TRUE, lty = 1)
      }
      plot(x, y, type = "l", axes = FALSE, ylab = "", xlab = "", xlim = ranx)
    }
    par(new = TRUE)
    Axis(side = 1)
    
    # add marginal distributions to the left
    par(mar = c(3.7, 1, 1.7, 1) + 0.1, fig = c(0, 0.2, 0.2, 1), new = TRUE)
    for (i in 1:4) {
      x <- 1:100 * (rany[2] - rany[1]) / 100 + rany[1]
      y <- dnorm(x, mean = model$means[i, 2], sd = sqrt(model$covmat[[i]][2, 2]))
      if (i == 2 || i == 4) {
        par(new = TRUE, lty = 2)
      } else {
        par(new = TRUE, lty = 1)
      }
      plot(y, x, type = "l", axes = FALSE, ylab = "", xlab = "", ylim = rany)
    }
    par(new = TRUE)
    Axis(side = 2)
  }
  
  if (scatter) {
    # add scatterplot if there are predicted and observed values
    if (any(names(model) == "predicted") && any(names(model) == "observed")) {
      par(mar = c(1.5, 1.5, 1, 1), fig = c(0, 0.33, 0, 0.33), new = TRUE)
      plot(model$predicted, model$observed, pch = 21, cex = 0.3, col = "gray40", bg = "gray40", bty = "n", axes = FALSE)
      abline(a = 0, b = 1, lty = 1)
      axis(side = 1, at = c(0, 1), mgp = c(3, 0.5, 0))
      axis(side = 2, at = c(0, 1), mgp = c(3, 0.5, 0))
    }
  }
}