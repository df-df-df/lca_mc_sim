LCA Simulation and Monte Carlo Template
=========================================

This repository provides a Latent Class Analysis (LCA) pipeline in R, featuring:
  - Data preparation (with recoding to positive integers)
  - (Optional) creation of composite variables
  - Fitting multiple LCA models (e.g., 1 to 4 classes)
  - Selecting the best model by BIC (or another criterion)
  - Running a Monte Carlo simulation to assess model stability via parametric bootstrapping
  - (Optional) parallel processing for faster simulations
  - Logging of results
  - Visualization (e.g., bar plot of Monte Carlo class‐selection counts)

---------------------------
File Overview
---------------------------
1. lca_mc_sim_example_script.R
   - A self‐contained demonstration of the entire workflow using randomly generated data.
   - Great for testing that everything works on your machine.

2. lca_mc_sim_template.R
   - A template script designed to be adapted to your own data and variables.
   - Contains helper functions for:
         * Data preparation
         * Model fitting (using the poLCA package)
         * Monte Carlo simulation
         * Parallelization
         * Logging and visualization

---------------------------
Required R Packages
---------------------------
The scripts rely on the following R packages:
  - poLCA (https://cran.r-project.org/package=poLCA)
  - mclust (https://cran.r-project.org/package=mclust) [optional for adjusted Rand Index]
  - foreach (https://cran.r-project.org/package=foreach) [for parallel loops]
  - doParallel (https://cran.r-project.org/package=doParallel) [for parallel loops]
  - ggplot2 (https://cran.r-project.org/package=ggplot2) [for plotting]

Install them with:
   install.packages(c("poLCA", "mclust", "foreach", "doParallel", "ggplot2"))

---------------------------
How to Use
---------------------------
1. Clone or Download the Repository
   - If you’re on GitHub, either clone the repository via Git or click “Download ZIP” and extract it locally.

2. Check Your Environment
   - In R or RStudio, ensure you can load and run the scripts:
         * lca_mc_sim_example_script.R
         * lca_mc_sim_template.R
   - If any required packages are missing, install them as shown above.

3. Run the Example Script (to get the gist of the script)
   - lca_mc_sim_example_script.R generates random dummy data, fits a 1–4 class LCA, runs a Monte Carlo simulation, and prints/logs/plots the results.
   - To run it, open R or RStudio and type:
         source("lca_mc_sim_example_script.R")
   - Observe the output:
         * Fit indices for 1–4 class models
         * The chosen best model by BIC
         * Monte Carlo simulation results (class counts and error rate)
         * A plot (if plotting is enabled)
         * A log file (lca_analysis_log.txt) containing basic run statistics

4. Adapt the Template Script for Real Data
   - lca_mc_sim_template.R is a general script with placeholders:
         * Data loading: Replace the dummy data section with your own data (e.g., my_data <- read.csv("my_data_file.csv")).
         * Composite creation: Adjust or remove the example composite variables as needed.
         * Indicators list: Update the new_lca_indicators vector to match the columns you want to include in the LCA.
         * Max classes: Increase or decrease max_classes (e.g., 2–6) based on your requirements.
         * Monte Carlo reps: Increase n_rep for a more robust simulation (e.g., 100, 500, or even 1000, depending on your available compute resources).
   - Then run:
         source("lca_mc_sim_template.R")
   - to execute your real analysis.

5. Optional: Parallel Monte Carlo
   - If you have multiple cores, you can speed up the Monte Carlo phase by using the parallel function:
         results_parallel <- run_monte_carlo_parallel(
                              best_model      = your_best_model,
                              indicator_names = your_indicator_names,
                              n_rep           = 100,
                              max_classes     = 4,
                              seed            = 123,
                              n_cores         = 4
                            )
   - This returns a vector of class counts across replications instead of a full results list, but you can adapt it for your needs.

6. Logging and Visualization
   - Logging: The log_results() function writes a brief summary of each simulation run to lca_analysis_log.txt. You can change the filename or disable logging as needed.
   - Plotting: The plot_monte_carlo_results() function (using ggplot2) produces a bar plot of how many times each class solution was chosen across the Monte Carlo replications.

---------------------------
Troubleshooting
---------------------------
1. Missing Packages:
   Ensure you have installed all the listed R packages. Verify with library(packageName) that they load properly.

2. Errors in Data Preparation:
   The function prepare_lca_data() expects columns to exist in your dataset. If it can’t find a specified column, it will stop. Check that your column names match exactly.

3. Model Fitting Issues:
   If a particular class solution fails (e.g., due to singularities or empty categories), the script catches the error but continues running. You might see messages about errors or empty solutions.

4. Large n_rep or Large max_classes:
   Increasing these can significantly increase runtime. Consider using parallel processing if you have enough CPU cores.

5. Plot Doesn’t Display:
   Some R environments don’t automatically show plots unless you explicitly print them. Try:
         print(plot_monte_carlo_results(mc_results))
   or use RStudio.
