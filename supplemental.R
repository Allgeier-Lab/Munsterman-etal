##--------------------------------------------##
##    Katrina Sky Munsterman                  ##
##    kmunster@umich.edu                      ##
##    www.github.com/kmunsterman              ##
##    supplemental tests                      ##
##--------------------------------------------##

#### Load packages ####

library(remotes)
remotes::install_github(repo = "Allgeier-Lab/arrR", ref = "fixed-population")

library(arrR)
library(tidyverse)
library(rslurm)
library(purrr)

##### PARAMETER A, B SENSITIVITY #####

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
max_i <- (60 * 24 * 365 * years) / min_per_i #43800
save_each <- max_i / 120 # saving every 30 days

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

## fish population values ##

starting_values$pop_sd_size <- 0
pop_n_values <- c(10)

## size structure values ##

pop_mean_size_values <- c(18)
pop_ldie_values <- c(18.5)

input_df_size <- data.frame(pop_n = pop_n_values, pop_mean_size = pop_mean_size_values,
                            pop_ldie = pop_ldie_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 5))

## movement values ##

parameters$pop_reserves_max <- 0.1

pop_reserves_consump_values <- c(0.05, 0.01, 0.99, 0.05, 0.05)
pop_reserves_thres_mean_values <- c(0.5, 0.5, 0.5, 0.99, 0.0001)

input_df_move <- data.frame(pop_reserves_consump = pop_reserves_consump_values,
                            pop_reserves_thres_mean = pop_reserves_thres_mean_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 1))

# create a dataframe with each treatment combination to run 25 simulations each 
input_df <- cbind(input_df_size, input_df_move) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 25))

#### Run simulations ####

foo <- function(pop_n, pop_mean_size, pop_ldie, 
                pop_reserves_consump, pop_reserves_thres_mean, class_width = 5) {
  
  starting_values$pop_n <- pop_n
  starting_values$pop_mean_size <- pop_mean_size
  parameters$pop_ldie <- pop_ldie
  parameters$pop_reserves_consump <- pop_reserves_consump
  parameters$pop_reserves_thres_mean <- pop_reserves_thres_mean
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1,
                                         reef = reef_matrix, 
                                         starting_values = starting_values, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor,
                                       starting_values = starting_values,
                                       parameters = parameters,
                                       use_log = TRUE, verbose = FALSE)
  
  result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                 parameters = parameters, movement = "behav",
                                 max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = save_each, 
                                 verbose = FALSE) %>%
    arrR::filter_mdlrn(filter = c(max_i/2, max_i), reset = TRUE)
  
  foo_mutate <- function(df) { 
    df %>%  
      dplyr::mutate(
        pop_reserves_consump = pop_reserves_consump,
        pop_reserves_thres_mean =  pop_reserves_thres_mean)
  }
  
  sens_sg_agbg <- 
    dplyr::filter(result$seafloor, timestep == max_i, reef == 0) %>% 
    dplyr::summarise(tot_ag_prod = sum(ag_production), tot_bg_prod = sum(bg_production)) %>%
    dplyr::mutate(tot_prod = tot_ag_prod + tot_bg_prod) %>%
    dplyr::mutate(ag_scaled = tot_ag_prod/2475/1825, bg_scaled = tot_bg_prod/2475/1825) %>%
    dplyr::mutate(tot_scaled = tot_prod/2475/1825) %>%
    foo_mutate()
  
  return(list(sens_sg_agbg = sens_sg_agbg))
  
}


#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

# create .sh script
sbatch_sens <- rslurm::slurm_apply(f = foo, params = input_df,
                                   global_objects = globals, jobname = "behav_sens",
                                   nodes = nrow(input_df), cpus_per_node = 1,
                                   slurm_options = list("account" = "kmunster0",
                                                        "partition" = "standard",
                                                        "time" = "00:15:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "10G"),
                                   pkgs = c("arrR", "dplyr"),
                                   rscript_path = rscript_path,
                                   submit = FALSE)

#### Results ####

behav_sens <- rslurm::get_slurm_out(sbatch_sens, outtype = "raw")

#### Summarize data ####

sens <- purrr::map(behav_sens[1:125], "sens_sg_agbg")
sensdf <- sens %>%
  map_df(as_tibble)

sens1 <- sensdf %>%
  mutate(
    treatment = case_when(
      pop_reserves_consump == 0.05 & pop_reserves_thres_mean == 0.5 ~ "Med",
      pop_reserves_consump == 0.01 & pop_reserves_thres_mean == 0.5 ~ "A High",
      pop_reserves_consump == 0.99 & pop_reserves_thres_mean == 0.5 ~ "A Low",
      pop_reserves_consump == 0.05 & pop_reserves_thres_mean == 0.99 ~ "B High",
      pop_reserves_consump == 0.05 & pop_reserves_thres_mean == 0.0001 ~ "B Low"
    )
  )

sens2 <- sens1 %>%
  group_by(treatment) %>%
  mutate(ag_mean = mean(ag_scaled), bg_mean = mean(bg_scaled), tot_mean = mean(tot_scaled)) %>%
  select(treatment, ag_mean, bg_mean, tot_mean) %>%
  unique()

med_row <- sens2[sens2$treatment == "Med", ]

percent_diff <- sens2
percent_diff$ag_mean_diff <- 100 * (percent_diff$ag_mean - med_row$ag_mean) / med_row$ag_mean
percent_diff$bg_mean_diff <- 100 * (percent_diff$bg_mean - med_row$bg_mean) / med_row$bg_mean
percent_diff$tot_mean_diff <- 100 * (percent_diff$tot_mean - med_row$tot_mean) / med_row$tot_mean

diff <- percent_diff %>%
  filter(treatment != "Med") %>%
  select(treatment, ag_mean_diff, bg_mean_diff, tot_mean_diff)

diff_long <- diff %>%
  pivot_longer(
    cols = c(ag_mean_diff, bg_mean_diff, tot_mean_diff),
    names_to = "part", 
    values_to = "percent_diff"
  )

#### Plot data ####

diff_long <- diff_long %>%
  mutate(part = factor(part, levels=c("tot_mean_diff", "bg_mean_diff", "ag_mean_diff")))

MyPalette <- c("#BD582C", "#865640", "#9B8357")

ggplot(diff_long, aes(x = treatment, y = percent_diff, fill = part)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = MyPalette) +
  labs(title = "",
       x = "",
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")


##### TIME WITHIN EACH BEHAVIOUR STATE #####

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
max_i <- (60 * 24 * 365 * years) / min_per_i #43800

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

## fish population values ##

starting_values$pop_sd_size <- 0
pop_n_values <- c(10)

## size structure values ##

pop_mean_size_values <- c(18)
pop_ldie_values <- c(18.5)

input_df_size <- data.frame(pop_n = pop_n_values, pop_mean_size = pop_mean_size_values,
                            pop_ldie = pop_ldie_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 3))

## movement values ##

parameters$pop_reserves_max <- 0.1

pop_reserves_consump_values <- c(0.99, 0.05, 0.01)
pop_reserves_thres_mean_values <- c(0.0001, 0.5, 0.99)

input_df_move <- data.frame(pop_reserves_consump = pop_reserves_consump_values,
                            pop_reserves_thres_mean = pop_reserves_thres_mean_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 1))

# final df

input_df <- cbind(input_df_size, input_df_move)

#### Run simulations ####

foo <- function(pop_n, pop_mean_size, pop_ldie, 
                pop_reserves_consump, pop_reserves_thres_mean, class_width = 5) {
  
  starting_values$pop_n <- pop_n
  starting_values$pop_mean_size <- pop_mean_size
  parameters$pop_ldie <- pop_ldie
  parameters$pop_reserves_consump <- pop_reserves_consump
  parameters$pop_reserves_thres_mean <- pop_reserves_thres_mean
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1,
                                         reef = reef_matrix, 
                                         starting_values = starting_values, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor,
                                       starting_values = starting_values,
                                       parameters = parameters,
                                       use_log = TRUE, verbose = FALSE)
  
  result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                 parameters = parameters, movement = "behav",
                                 max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = 1, 
                                 verbose = FALSE) %>%
    arrR::filter_mdlrn(filter = c(max_i/2, max_i), reset = TRUE)
  
  foo_mutate <- function(df) { 
    df %>%       
      dplyr::mutate(bio = dplyr::case_when(
        pop_reserves_consump_values == 0.99 ~ "low",
        pop_reserves_consump_values == 0.05 ~ "med",
        pop_reserves_consump_values == 0.01 ~ "high"))
    
  }
  
  behavior <- 
    dplyr::filter(result$fishpop, timestep >=43789) %>% # 12 timesteps (final day in 120 minute iterations)
    dplyr::select(id, behavior, timestep) %>%
    dplyr::mutate(behavior = dplyr::case_when(behavior == 1 ~ "shelter",
                                              behavior %in% c(2, 3) ~ "forage")) %>%
    foo_mutate()
  
  return(behavior)
  
}

#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

# create .sh script
sbatch_time <- rslurm::slurm_apply(f = foo, params = input_df,
                                   global_objects = globals, jobname = "behav_time",
                                   nodes = nrow(input_df), cpus_per_node = 1,
                                   slurm_options = list("account" = "kmunster0",
                                                        "partition" = "standard",
                                                        "time" = "00:15:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "50G"),
                                   pkgs = c("arrR", "dplyr"),
                                   rscript_path = rscript_path,
                                   submit = FALSE)

#### Results ####

behav_time <- rslurm::get_slurm_out(sbatch_time, outtype = "table")

#######################################################################

##### WATER COLUMN NUTRIENTS AND DETRITUS OVER TIME #####

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
max_i <- (60 * 24 * 365 * years) / min_per_i #43800
save_each <- max_i / 120 # saving every 30 days

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

## fish population values ##

starting_values$pop_sd_size <- 0
pop_n_values <- c(10, 64)

## size structure values ##

pop_mean_size_values <- c(18, 10)
pop_ldie_values <- c(18.5, 10.5)

input_df_size <- data.frame(pop_n = pop_n_values, pop_mean_size = pop_mean_size_values,
                            pop_ldie = pop_ldie_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 1))

## movement values ##

parameters$pop_reserves_max <- 0.1

pop_reserves_consump_values <- c(0.05)
pop_reserves_thres_mean_values <- c(0.5)

input_df_move <- data.frame(pop_reserves_consump = pop_reserves_consump_values,
                            pop_reserves_thres_mean = pop_reserves_thres_mean_values) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 2))

# create a dataframe with each treatment combination to run 50 simulations each 
input_df <- cbind(input_df_size, input_df_move) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 50))

#### Run simulations ####

foo <- function(pop_n, pop_mean_size, pop_ldie, 
                pop_reserves_consump, pop_reserves_thres_mean, class_width = 5) {
  
  starting_values$pop_n <- pop_n
  starting_values$pop_mean_size <- pop_mean_size
  parameters$pop_ldie <- pop_ldie
  parameters$pop_reserves_consump <- pop_reserves_consump
  parameters$pop_reserves_thres_mean <- pop_reserves_thres_mean
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1,
                                         reef = reef_matrix, 
                                         starting_values = starting_values, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor,
                                       starting_values = starting_values,
                                       parameters = parameters,
                                       use_log = TRUE, verbose = FALSE)
  
  result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                 parameters = parameters, movement = "behav",
                                 max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = save_each, 
                                 verbose = FALSE)
  
  foo_mutate <- function(df) { 
    df %>%       
      dplyr::mutate(size_struc = dplyr::case_when(
        pop_mean_size == 18 ~ "large",
        pop_mean_size == 10 ~ "small"))
  }
  
  
  nuts_det <- 
    dplyr::filter(result$seafloor) %>%
    dplyr::select(timestep, nutrients_pool, detritus_pool) %>%
    dplyr::group_by(timestep) %>%
    dplyr::mutate(wc_scaled = (sum(nutrients_pool))/2500, det_scaled = (sum(detritus_pool))/2500) %>%
    select(wc_scaled, det_scaled) %>%
    unique() %>%
    foo_mutate()
  
  return(nuts_det)
  
}

#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

# create .sh script
sbatch_sf <- rslurm::slurm_apply(f = foo, params = input_df,
                                 global_objects = globals, jobname = "seafloor",
                                 nodes = nrow(input_df), cpus_per_node = 1,
                                 slurm_options = list("account" = "kmunster0",
                                                      "partition" = "standard",
                                                      "time" = "00:15:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "25G"),
                                 pkgs = c("arrR", "dplyr"),
                                 rscript_path = rscript_path,
                                 submit = FALSE)

#### Results ####

seafloor <- rslurm::get_slurm_out(sbatch_sf, outtype = "table")

#### Summarize data ####

sfdf <- seafloor %>%
  group_by(timestep, size_struc) %>%
  summarise(water_column = mean(wc_scaled), detritus = mean(det_scaled)) %>%
  mutate(day = timestep / 12) %>%
  unique() %>%
  pivot_longer(cols = c(water_column, detritus), 
               names_to = "seafloor", 
               values_to = "amount")

#### Plot data ####

MyPalette <- c( "#5b3e31", "#3b719f")

sf_plot <- ggplot(sfdf, aes(x = day, y = amount, color = seafloor)) +
  geom_point(size = 4) +
  scale_colour_manual(values = MyPalette) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")

seafloor_plot <- sf_plot + facet_wrap(~ size_struc) + theme(strip.text = element_text(size = 16)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
