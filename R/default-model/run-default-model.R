source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)
fs = .Platform$file.sep


# Setup  ------------------------------------------------------------------
# target="targets_test_config"
# target = "targets_with_or"
target = "targets_paper_resubmission"

# set seed once randomly
# seed_wppl = as.numeric(Sys.time())
seed_wppl = "1644922777.24488"
save_seed_subdir = TRUE

#run_sweep = TRUE
run_sweep = FALSE

# Run model ---------------------------------------------------------------
for(i in seq(1,4)){
  if(i == 1){
    params <- configure("config.yml", c("pl", target))
  } else if(i==2){
    params <- configure("config.yml", c("speaker", target))
  } else if(i==3){
    params <- configure("config.yml", c("speaker_literal", target))
  } else if (i==4){
    params <- configure("config.yml", c("priorN", target))
  }
  params$seed_webppl = seed_wppl

  
  # Setup -------------------------------------------------------------------
  if(save_seed_subdir){
    base_dir <- file.path(params$base_dir, seed_wppl, fsep=fs)
    params$seed_dir <- base_dir
  } else {
    base_dir <- params$base_dir
    params$seed_dir <- base_dir
  }
  alphas <- params$alpha
  thetas <- params$theta
  if(run_sweep) {
    pars <- configure("config.yml", c("sweep"), load_default=FALSE)
    alphas <- pars$alpha
    thetas <- pars$theta
  }
  for(alpha in alphas) {
    for(theta in thetas) {
      params$alpha <- alpha
      params$theta <- theta
      subdir <- paste(alpha, theta, sep="--")
      target_dir <- file.path(base_dir, subdir, fsep=fs)
      if(!dir.exists(target_dir)) dir.create(target_dir, recursive=T)
      
      params$target_dir <- target_dir
      params$target <- file.path(
        target_dir, paste("results-", params$fn_suffix, ".rds", sep=""), fsep=fs
      )
      params$target_params <- str_replace(params$target, "results", "params")
      params$utts_path <- file.path(params$base_dir, params$utts_fn, fsep=fs)
      
      params$plot_dir = paste(target_dir, "figs", sep=fs)
      if(!dir.exists(params$plot_dir)) dir.create(params$plot_dir)
  
      ##---- Generate/Retrieve utterances ----##
      if(!file.exists(params$utts_path)){
        utterances <- generate_utts(params)
      } else {
        utterances <- readRDS(params$utts_path)
        print(paste("utterances read from:", params$utts_path))
      }
      params$utterances <- utterances
   
      # Run Model ---------------------------------------------------------------
      posterior <- run_webppl(params$model_path, params)
      
      # structure + save data
      if(params$level_max == "speaker") {
        speaker <- posterior %>% structure_speaker_data(params)
        if(params$save) {
          save_data(posterior$all_ids %>% rename(bn_id=value),
                  str_replace(params$target, "results", "sample-ids"))
        }
      } else if(params$level_max == "priorN"){
          data <- structure_bns(posterior, params)
      } else {
        data <- posterior %>% structure_listener_data(params)
      }
    }
  }
}
