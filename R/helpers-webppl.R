# formats webppl distributions of P(state|utterance), i.e. for listeners + prior
webppl_distrs_to_tibbles <- function(posterior){
  posterior_tibbles <- map2(posterior, names(posterior), function(x, y){
    x <- x %>% rowid_to_column() 
    bn_probs <- x %>% dplyr::select("probs", "rowid")
    data_tibble <- x$support %>% rowid_to_column() %>% 
                    unnest(cols = c(table.probs, table.support)) %>%
                    as_tibble() %>% 
                    left_join(bn_probs, by = "rowid") %>% 
                    mutate("rowid" = as.character(rowid)) %>%
                    add_column(level=y) %>% 
                    rename(prob=probs, val=table.probs, cell=table.support)
    return(data_tibble)             
  })
  df <- bind_rows(posterior_tibbles)
  return(df)
}

structure_bns <- function(posterior, params){
  data.long <- posterior$bns %>% rowid_to_column(var = "rowid") %>%
    unnest(c(table.probs, table.support)) %>%
    rename(val=table.probs, cell=table.support) %>%
    add_column(level=params$level_max)

  if(params$add_accept_conditions){
    df_wide <- data.long %>% spread(key=cell, val=val)
    df <- acceptability_conditions(df_wide)
    data.long <- df %>% group_by(rowid, r, level) %>%
      pivot_longer(cols = c(`AC`, `A-C`, `-AC`, `-A-C`),
                   names_to = "cell", values_to = "val")
  }
  if(params$save){
    data.long %>% save_data(params$target)
    params %>% save_data(params$target_params)
  }
  return(data.long)
}

run_webppl <- function(path_wppl_file, params){
  if(params$verbose){
    print(paste('model file read from:', path_wppl_file))
    print(paste('packages loaded from:' ,params$packages))
  }
  data <-   webppl(program_file = path_wppl_file,
                   data = params,
                   data_var = "data",
                   random_seed = params$seed_webppl,
                   packages=params$packages
                  )
  # data is a list of lists
  data <- data %>% map(function(x){as_tibble(x)})
  return(data)
}

structure_listener_data <- function(posterior, params){
  df_long <- posterior %>% webppl_distrs_to_tibbles()
  if(params$add_accept_conditions){
    df_wide <- df_long %>% spread(key=cell, val=val)
    df <- acceptability_conditions(df_wide)
    df_long <- df %>% group_by(rowid, r, level) %>%
      pivot_longer(cols=c(AC, `A-C`, `-AC`, `-A-C`),
                   names_to="cell", values_to="val")
  }
  if(params$save){
    df_long %>% save_data(params$target)
    params %>% save_data(params$target_params)
  }
  return(df_long)
}


# summarise webppl distributions ------------------------------------------
webppl_speaker_distrs_to_tibbles <- function(posterior){
  speaker <- posterior[names(posterior) != "bns"] 
  posterior_tibbles <- map2(speaker, names(speaker), function(x, y){
    data_tibble <- x %>% rowid_to_column() %>% unnest(cols = c(probs, support)) %>% 
      rename(utterance=support) %>% 
      add_column(level=y)
    return(data_tibble)             
  })
  speaker <- bind_rows(posterior_tibbles) 
  bns_unique <- posterior$bns %>% rowid_to_column() %>%
    unnest(cols = c(table.probs, table.support)) %>% 
    rename(cell=table.support, val=table.probs) %>%
    spread(key=cell, val=val) %>%
    nest(data = c(r, `-A-C`, `-AC`, `A-C`, `AC`))

  # add AC, A-C, -AC, -A-C cell entries to speaker data
  bns <- bns_unique[speaker$rowid,]$data
  speaker_wide <- speaker %>% add_column(bn=bns) %>% unnest(cols = c(bn)) %>%
    spread(key=utterance, val=probs, fill=0)
  
  return(speaker_wide)
}

structure_speaker_data <- function(posterior, params, tbls.map=NA){
  posterior_dist = posterior$distributions
  speaker_wide <- webppl_speaker_distrs_to_tibbles(posterior_dist)
  bns = posterior_dist$bns %>% rowid_to_column() %>% group_by(rowid) %>%
    unnest(c(table.probs, table.support)) %>%
    pivot_wider(names_from=table.support, values_from=table.probs) %>%
    select(-r);
  df.wide = left_join(speaker_wide %>% select(-`AC`, -`A-C`, -`-AC`, -`-A-C`),
                      bns, by = "rowid")
  df <- acceptability_conditions(df.wide)
  df <- df %>% group_by(rowid) %>%
    pivot_longer(cols=c(-rowid, -bn_id, -starts_with("cell."), 
                        -p_delta, -p_rooij, -p_diff,
                        -level, -r, -`AC`, -`A-C`, -`-AC`, -`-A-C`),
                  names_to = "utterance", values_to = "probs")
  if(!is.na(tbls.map)){
    sp = left_join(df, tbls.map, by="bn_id")
  } else { sp = df}

  if(params$save){
    sp %>% save_data(params$target)
    params %>% save_data(params$target_params)
  }
  return(sp)
}
