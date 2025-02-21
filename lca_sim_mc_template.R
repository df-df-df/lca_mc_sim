################################################################################
# LCA Template Script (Generic)
#
# This script demonstrates how to:
#   1) Prepare data for LCA (latent class analysis),
#   2) Create composite variables,
#   3) Fit multiple LCA models and store fit indices,
#   4) Run a Monte Carlo simulation to assess model stability,
#   5) Optionally perform parallel processing,
#   6) Log results and visualize outcomes.
#
# Author: df-df-df-df
# Date:   2025-02-21
################################################################################

##############################
### (0) Load / Install Required Packages
##############################
# Uncomment and run if needed:
# install.packages("poLCA")
# install.packages("mclust")
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("ggplot2")

library(poLCA)       # Latent Class Analysis
library(mclust)      # For adjustedRandIndex (if needed)
library(foreach)     # For optional parallel processing
library(doParallel)  # For optional parallel processing
library(ggplot2)     # For visualization

##############################
### (1) Configuration & Logging
##############################

# A configuration object to store settings for LCA
lca_config <- list(
  max_iterations         = 1000,
  convergence_criterion  = 1e-10,
  n_random_starts        = 20,
  parallel               = FALSE,
  verbose                = TRUE,
  save_intermediate      = FALSE
)

# Simple logging function to append run info to a log file
log_results <- function(results, filename = "lca_analysis_log.txt") {
  sink(filename, append = TRUE)
  cat("\n=== Analysis Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
  cat("Number of replications:", results$parameters$n_rep, "\n")
  cat("Computation time (mins):", round(as.numeric(results$computation_time), 2), "\n")
  cat("Error count:", results$error_count, "\n")
  cat("Error rate (%):", round(results$error_count/results$parameters$n_rep * 100, 2), "\n")
  sink()
}

##############################
### (2) Helper Functions
##############################

# (2.1) Prepare LCA Data: Subset columns and recode them to positive integers.
prepare_lca_data <- function(df, indicators) {
  # Basic input validation:
  if (!is.data.frame(df)) stop("Input 'df' must be a data frame")
  if (!all(indicators %in% names(df))) {
    missing <- setdiff(indicators, names(df))
    stop("The following indicators are not found in the data: ", paste(missing, collapse = ", "))
  }
  
  sub_df <- df[, indicators, drop = FALSE]
  
  cat("\nData preparation summary (after recoding):\n")
  for (colnm in indicators) {
    # Convert to numeric and recode to ensure lowest value is 1
    sub_df[[colnm]] <- as.numeric(sub_df[[colnm]])
    min_val <- min(sub_df[[colnm]], na.rm = TRUE)
    if (min_val < 1) {
      sub_df[[colnm]] <- sub_df[[colnm]] - min_val + 1
    }
    sub_df[[colnm]] <- round(sub_df[[colnm]])
    rng <- range(sub_df[[colnm]], na.rm = TRUE)
    cat(colnm, ": range [", paste(rng, collapse = ", "), "]\n")
  }
  return(sub_df)
}

# (2.2) Run LCA Models: Fit models for k = 1..max_classes using poLCA.
run_lca_models <- function(data, indicators, max_classes = 4, config = lca_config) {
  gc()  # Garbage collection for memory management
  model_list <- vector("list", max_classes)
  fit_list <- data.frame(
    nclass       = integer(),
    logLik       = numeric(),
    AIC          = numeric(),
    BIC          = numeric(),
    entropy      = numeric(),
    n_parameters = numeric(),
    class_sizes  = character(),
    stringsAsFactors = FALSE
  )
  
  for (k in 1:max_classes) {
    cat("\nFitting", k, "-class model...\n")
    f <- as.formula(paste0("cbind(", paste(indicators, collapse = ", "), ") ~ 1"))
    
    model_k <- tryCatch({
      poLCA(f, data = data, nclass = k, nrep = config$n_random_starts,
            verbose = FALSE, na.rm = FALSE)
    }, error = function(e) {
      cat("Error in fitting", k, "-class model:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(model_k)) {
      model_list[[k]] <- NA
      fit_list <- rbind(fit_list, data.frame(
        nclass = k, logLik = NA, AIC = NA, BIC = NA,
        entropy = NA, n_parameters = NA, class_sizes = NA,
        stringsAsFactors = FALSE
      ))
    } else {
      logLik_val  <- model_k$llik
      AIC_val     <- model_k$aic
      BIC_val     <- model_k$bic
      class_probs <- model_k$P
      entropy_val <- NA  # Not computed here
      n_params    <- model_k$npar
      class_size_str <- paste(round(class_probs * 100, 1), collapse = "/")
      
      model_list[[k]] <- model_k
      fit_list <- rbind(fit_list, data.frame(
        nclass       = k,
        logLik       = logLik_val,
        AIC          = AIC_val,
        BIC          = BIC_val,
        entropy      = entropy_val,
        n_parameters = n_params,
        class_sizes  = class_size_str,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(list(models = model_list, fit_indices = fit_list))
}

# (2.3) Simulate Dataset from an LCA Model (Parametric Bootstrap).
simulate_dataset <- function(N, pi_est, probs_list) {
  K <- length(pi_est)
  latent_class <- sample(seq_len(K), size = N, replace = TRUE, prob = pi_est)
  
  sim_data <- data.frame(matrix(ncol = length(probs_list), nrow = N))
  colnames(sim_data) <- names(probs_list)
  
  for (j in seq_along(probs_list)) {
    prob_mat <- probs_list[[j]]
    n_cat <- ncol(prob_mat)
    responses <- numeric(N)
    for (i in seq_len(N)) {
      this_class <- latent_class[i]
      pvec <- prob_mat[this_class, ]
      responses[i] <- sample(seq_len(n_cat), size = 1, prob = pvec)
    }
    sim_data[[j]] <- responses
  }
  return(list(data = sim_data, latent = latent_class))
}

# (2.4) Run Monte Carlo Simulation
run_monte_carlo <- function(best_model, indicator_names, n_rep = 10, max_classes = 4, seed = 42042) {
  set.seed(seed)
  start_time <- Sys.time()
  
  pi_est     <- as.numeric(best_model$P)
  probs_list <- best_model$probs
  N          <- best_model$Nobs
  
  class_counts    <- integer(max_classes)
  error_count     <- 0
  error_messages  <- character()
  
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3)
  
  for (rep_i in seq_len(n_rep)) {
    setTxtProgressBar(pb, rep_i)
    
    sim_out          <- simulate_dataset(N, pi_est, probs_list)
    sim_data         <- sim_out$data
    sim_data_prepped <- prepare_lca_data(sim_data, indicator_names)
    
    sim_lca <- tryCatch({
      run_lca_models(sim_data_prepped, indicator_names, max_classes = max_classes)
    }, error = function(e) {
      error_count <<- error_count + 1
      error_messages <<- c(error_messages, paste("Replication", rep_i, ":", conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(sim_lca)) next
    
    valid_fit <- sim_lca$fit_indices[!is.na(sim_lca$fit_indices$BIC), ]
    if (nrow(valid_fit) == 0) {
      error_count <- error_count + 1
      next
    }
    
    # Choose the number of classes by minimum BIC
    bic_vals  <- valid_fit$BIC
    chosen_k  <- valid_fit$nclass[which.min(bic_vals)]
    class_counts[chosen_k] <- class_counts[chosen_k] + 1
  }
  
  close(pb)
  computation_time <- difftime(Sys.time(), start_time, units = "mins")
  
  cat("\nMonte Carlo Simulation (n_rep =", n_rep, ") results:\n")
  for (k in 1:max_classes) {
    cat(k, "-class chosen in", class_counts[k], "out of", n_rep,
        "replications (", round(100 * class_counts[k] / n_rep, 1), "% )\n")
  }
  cat("Replications with errors or no valid model:", error_count, "\n")
  
  results <- list(
    class_counts     = class_counts,
    error_count      = error_count,
    error_messages   = error_messages,
    computation_time = computation_time,
    parameters = list(n_rep = n_rep, max_classes = max_classes, seed = seed)
  )
  return(results)
}

# (2.5) Parallel Processing Option
run_monte_carlo_parallel <- function(best_model, indicator_names, n_rep = 10,
                                     max_classes = 4, seed = 123, n_cores = NULL) {
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results <- foreach(i = 1:n_rep, .combine = rbind, .packages = c("poLCA")) %dopar% {
    set.seed(seed + i)
    sim_out <- simulate_dataset(best_model$Nobs, as.numeric(best_model$P), best_model$probs)
    sim_data <- prepare_lca_data(sim_out$data, indicator_names)
    
    sim_lca <- tryCatch({
      run_lca_models(sim_data, indicator_names, max_classes = max_classes)
    }, error = function(e) { NULL })
    
    if (is.null(sim_lca)) {
      return(c(rep(NA, max_classes)))
    }
    valid_fit <- sim_lca$fit_indices[!is.na(sim_lca$fit_indices$BIC), ]
    if (nrow(valid_fit) == 0) {
      return(c(rep(NA, max_classes)))
    }
    bic_vals <- valid_fit$BIC
    chosen_k <- valid_fit$nclass[which.min(bic_vals)]
    
    counts <- rep(0, max_classes)
    counts[chosen_k] <- 1
    return(counts)
  }
  
  stopCluster(cl)
  class_counts <- colSums(results, na.rm = TRUE)
  return(class_counts)
}

# (2.6) Visualization
plot_monte_carlo_results <- function(results) {
  df <- data.frame(
    Classes = factor(1:length(results$class_counts)),
    Count   = results$class_counts
  )
  ggplot(df, aes(x = Classes, y = Count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = "Monte Carlo Simulation Results",
         x     = "Number of Classes",
         y     = "Frequency")
}

##############################
### (3) Data Preparation & Composite Variables (Example)
##############################
# Example of how you might create composite variables in your own data.
# In practice, rename the placeholders (var1, var2, etc.) to the columns in your dataset.

# Suppose your dataset is named 'my_data' and you've already loaded it.
# Example:
# my_data <- read.csv("my_data.csv")

cat("\n# Creating composite variables (example) #\n")
my_data$Composite_1 <- rowMeans(
  my_data[, c("var1", "var2")],
  na.rm = TRUE
)
my_data$Composite_2 <- rowMeans(
  my_data[, c("var3", "var4")],
  na.rm = TRUE
)
my_data$Composite_3 <- rowMeans(
  my_data[, c("var5", "var6")],
  na.rm = TRUE
)
my_data$Composite_4 <- rowMeans(
  my_data[, c("var7", "var8")],
  na.rm = TRUE
)

##############################
### (4) Define LCA Indicators (Example)
##############################
cat("\n# Defining new LCA indicators (composites + selected items) #\n")
new_lca_indicators <- c(
  "Composite_1",
  "Composite_2",
  "Composite_3",
  "Composite_4",
  # ... possibly more columns ...
  "var9",
  "var10"
)

##############################
### (5) Prepare Data & Fit LCA Models
##############################
cat("\n# Preparing data for LCA #\n")
lca_data_new <- prepare_lca_data(my_data, new_lca_indicators)

cat("\n# Running LCA models (1..4 classes) #\n")
lca_results_new <- run_lca_models(lca_data_new, new_lca_indicators, max_classes = 4, config = lca_config)

cat("\nFit Indices for LCA Models:\n")
print(lca_results_new$fit_indices)

# Pick the "best" model by minimum BIC:
best_model_index_new <- which.min(lca_results_new$fit_indices$BIC)
best_model_new <- lca_results_new$models[[best_model_index_new]]

# You can store predicted class membership as a new column in your dataset:
my_data$lca_class <- best_model_new$predclass

##############################
### (6) Print Correlation Matrix & Conditional Probabilities
##############################
cat("\nPearson Correlation Matrix for LCA Indicators:\n")
simple_cor_new <- cor(lca_data_new[, new_lca_indicators], use = "pairwise.complete.obs")
print(simple_cor_new)

cat("\nConditional probabilities for the best LCA model:\n")
for (i in seq_along(best_model_new$probs)) {
  indicator_name <- names(best_model_new$probs)[i]
  cat("\nIndicator:", indicator_name, "\n")
  print(round(best_model_new$probs[[i]], 3))
}

##############################
### (7) Monte Carlo Simulation (Example: 10 replications)
##############################
cat("\n# Running Parametric Monte Carlo Simulation (for demonstration) #\n")
mc_results <- run_monte_carlo(
  best_model     = best_model_new,
  indicator_names = new_lca_indicators,
  n_rep          = 500,    # Increase later as needed
  max_classes    = 4,
  seed           = 42042
)

print(mc_results)

# Optionally, log the results to a file:
log_results(mc_results)

# Optionally, visualize the results:
plot_monte_carlo_results(mc_results)
