##--------------------------------------------##
##    Katrina Sky Munsterman                  ##
##    kmunster@umich.edu                      ##
##    www.github.com/kmunsterman              ##
##    question 2 script                       ##
##--------------------------------------------##

### Load packages
library(remotes)
remotes::install_github(repo = "Allgeier-Lab/arrR", ref = "fixed-population")

library(arrR)
library(tidyverse)
library(rslurm)

#### Load values ####

parameters <- arrR::default_parameters
starting_values <- arrR::default_starting

#### Setup environment ####

# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

#### Setup experiment ####

## time ##

# one iterations equals 120 minutes
min_per_i <- 120

# run the model for 10 years
years <- 10
max_i <- (60 * 24 * 365 * years) / min_per_i
save_each <- max_i / 120

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# fish population values ##

starting_values$pop_sd_size <- 0

pop_n_values <- c(5, 32, 10, 64, 20, 128)
pop_mean_size_values <- c(18, 10, 18, 10, 18, 10)
pop_ldie_values <- c(18.5, 10.5, 18.5, 10.5, 18.5, 10.5)

input_df_size <- data.frame(pop_n = pop_n_values, pop_mean_size = pop_mean_size_values,
                            pop_ldie = pop_ldie_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 9))

## movement values ##

parameters$pop_reserves_max <- 0.1

sd <- (rep(c(2, 2, 2, 0.2, 0.01, 0.001, 0.01, 0.001, 0.0001), times = 1, each = 1))
pop_reserves_consump_values <- c(0.99, 0.05, 0.01)
distribution_values <- c("rbeta", "rnorm", "rbeta")

input_df_move <- expand.grid(pop_reserves_consump = pop_reserves_consump_values,
                             distribution = distribution_values) %>% 
  dplyr::mutate(distribution = as.character(distribution)) %>%
  cbind(sd) %>%
  dplyr::slice(rep(1:dplyr::n(), 1))

# create a dataframe with each treatment combination to run 5 simulations each 
input_df <- cbind(input_df_size, input_df_move) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 50))

#### Run simulations ####

foo <- function(pop_n, pop_mean_size, pop_ldie, pop_reserves_consump, distribution, sd, class_width = 5) {
  
  starting_values$pop_n <- pop_n
  starting_values$pop_mean_size <- pop_mean_size
  parameters$pop_ldie <- pop_ldie
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1,
                                         reef = reef_matrix, starting_values = starting_values, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor,
                                       starting_values = starting_values,
                                       parameters = parameters,
                                       use_log = TRUE, verbose = FALSE)
  
  foo_helper <- get(distribution, mode = "function")
  
  values <- foo_helper(nrow(input_fishpop), pop_reserves_consump, sd)
  
  for (i in 1:length(values)) {
    
    while(values[i] > 1 | values[i] < 0) {
      
      values[i] <- foo_helper(nrow(input_fishpop), pop_reserves_consump, sd)
      
    }
    
  }
  
  threshold_mat <- cbind(input_fishpop$id, exp(-9.984 * values), values)
  
  result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                 parameters = parameters, threshold_mat = threshold_mat, 
                                 movement = "behav", max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = save_each, 
                                 verbose = FALSE) %>%
    arrR::filter_mdlrn(filter = c(max_i/2, max_i), reset = TRUE)
  
  foo_mutate <- function(df) { 
    df %>%       
      dplyr::mutate(bio = dplyr::case_when(
        pop_n == 5 ~ "low",
        pop_n == 32 ~ "low", 
        pop_n == 10 ~ "med", 
        pop_n == 64 ~ "med",
        pop_n == 20 ~ "high",
        pop_n == 128 ~ "high"),
        size_struc = dplyr::case_when(
          pop_mean_size == 18 ~ "large",
          pop_mean_size == 10 ~ "small"),
        mean_mvmnt = dplyr::case_when(
          pop_reserves_consump == 0.99 ~ "low",
          pop_reserves_consump == 0.05 ~ "med",
          pop_reserves_consump == 0.01 ~ "high"),
        intra_mvmnt = dplyr::case_when(
          distribution == "rbeta" & sd == 2 ~ "bold",
          distribution == "rnorm" ~ "normal",
          distribution == "rbeta" & sd %in% c(0.01, 0.001, 0.0001) ~ "shy"))
  }
  
  intra_behav <- 
    dplyr::select(result$fishpop, id, behavior) %>%
    dplyr::mutate(behavior = dplyr::case_when(behavior == 1 ~ "shelter",
                                              behavior %in% c(2, 3) ~ "forage")) %>%
    dplyr::group_by(id, behavior) %>% 
    dplyr::count(behavior) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(freq = n/ sum(n)) %>%
    dplyr::ungroup() %>%
    foo_mutate()
  
  intra_fish <- 
    dplyr::group_by(result$fishpop, timestep) %>% 
    dplyr::summarise(biomass = sum(weight)) %>%
    foo_mutate()
  
  intra_sg_agbg <- 
    dplyr::filter(result$seafloor, timestep == max_i, reef == 0) %>% 
    dplyr::summarise(tot_ag_prod = sum(ag_production), tot_bg_prod = sum(bg_production)) %>%
    dplyr::mutate(tot_prod = tot_ag_prod + tot_bg_prod) %>%
    dplyr::mutate(ag_scaled = tot_ag_prod/2475/1825, bg_scaled =  tot_bg_prod/2475/1825) %>%
    dplyr::mutate(tot_scaled = tot_prod/2475/1825) %>%
    foo_mutate()
  
  intra_sg_agbg_dist <- 
    dplyr::filter(result$seafloor, timestep == max_i, reef == 0) %>%
    dplyr::mutate(dist = sqrt(x ^ 2 + y ^ 2), 
                  dist_class = cut(dist, 
                                   breaks = seq(from = 0, to = max(dist) + class_width,
                                                by = class_width), ordered_result = TRUE)) %>%
    dplyr::select(dist_class, ag_production, bg_production) %>%
    tidyr::pivot_longer(cols = c(ag_production, bg_production), 
                        names_to = "part", values_to = "production") %>%
    dplyr::mutate(distance = dplyr::case_when(dist_class == "(0,5]" ~ "near",
                                              dist_class != "(0,5]"  ~ "far")) %>%
    dplyr::group_by(part, distance) %>%
    dplyr::summarise(tot_value_scaled = (sum(production))/1825, .groups = "drop") %>%
    dplyr::group_by(distance) %>%
    dplyr::mutate(tot_prod_scaled = sum(tot_value_scaled)) %>%
    foo_mutate()
  
  return(list(intra_behav = intra_behav, intra_fish = intra_fish, intra_sg_agbg = intra_sg_agbg, intra_sg_agbg_dist = intra_sg_agbg_dist))
  
}

#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

# create .sh script
sbatch_intra <- rslurm::slurm_apply(f = foo, params = input_df,
                                    global_objects = globals, jobname = "behav_q2",
                                    nodes = nrow(input_df), cpus_per_node = 1,
                                    slurm_options = list("account" = "kmunster0",
                                                         "partition" = "standard",
                                                         "time" = "00:30:00", ## hh:mm::ss
                                                         "mem-per-cpu" = "10G"),
                                    pkgs = c("arrR", "dplyr"),
                                    rscript_path = rscript_path,
                                    submit = FALSE)
