

noise_uniform <- function(data) {
  return(sum(data) * log(1/nrow(data)))
}

noise_differential <- function(data, noise_ratio = 2) {
  # noise_ratio: controls the relative probability of noise in hl/lh vs. hh/ll
  total_responses <- sum(data)
  hh_ll_responses <- data[1, 1] + data[4, 4] # hh and ll
  hl_lh_responses <- data[1, 3] + data[3, 1] + data[2, 4] + data[4, 2] + data[2, 2] + data[3, 3] # hl and lh
  
  # Probabilities based on noise_ratio
  hl_lh_prob <- noise_ratio / (2 * noise_ratio + 2)
  hh_ll_prob <- 1 / (2 * noise_ratio + 2)
  
  log_likelihood <- hh_ll_responses * log(hh_ll_prob) + hl_lh_responses * log(hl_lh_prob)
  
  return(log_likelihood)
}
