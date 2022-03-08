library(viridisLite)
library(yaml)

LABELS_R = c(
  "A || C" = expression(paste("A,C independent")),
  "A implies C" = expression(paste(A%->%C)),
  "A implies -C" = expression(paste(A%->%C)),
  "-A implies C" = expression(paste(A%->%C)),
  "-A implies -C" = expression(paste(A%->%C)),
  "C implies A" = expression(paste(C%->%A)),
  "C implies -A" = expression(paste(C%->%A)),
  "-C implies A" = expression(paste(C%->%A)),
  "-C implies -A" = expression(paste(C%->%A))
)
LABELS_r = c(
  "A implies C" = parse(text=TeX('$A^+C^+$')),
  "C implies A" = parse(text=TeX('$C^+A^+$')),
  "-A implies -C" = parse(text=TeX('$A^-C^-$')),
  "-C implies -A" = parse(text=TeX('$C^-A^-$')),
  "A implies -C" = parse(text=TeX('$A^+C^-$')),
  "-A implies C" = parse(text=TeX('$A^-C^+$')),
  "C implies -A" = parse(text=TeX('$C^+A^-$')),
  "-C implies A" = parse(text=TeX('$C^-A^+$')),
  "A || C" = expression(paste("A,C independent"))
)
LABELS_levels <- 
  c(`prior` = expression(paste("Prior")),
    `literal-speaker` = expression(paste("Literal speaker")),
    `pragmatic-speaker` = expression(paste("Pragmatic speaker")))

viridis1 <- viridis(20)
viridis2 <- plasma(20)
viridis3 <- cividis(20)
COLS_r <- c(
  "C implies A" = viridis2[13],
  "A implies C" = viridis2[10],
  "-A implies -C" = viridis2[4],
  "-C implies -A" = viridis2[7],
  "A || C" = viridis3[18],
  "-C implies A" = viridis1[16],
  "-A implies C" = viridis1[14],
  "C implies -A" = viridis1[12],
  "A implies -C" = viridis1[18]
)
  

save_data <- function(data, target_path){
  data %>% write_rds(target_path)
  print(paste("saved to:", target_path))
}

filter_vars <- function(df_long, vars){
  df <- df_long %>% mutate(keep=TRUE)
  for(var in vars){
    if(str_detect(var, "^-")){
      # negative variable, starts with -
      df <- df %>% mutate(keep=case_when(!keep ~ keep, TRUE ~ str_detect(cell, var)))
    }
    else {
      token <- paste("-", var, sep="")
      df <- df %>% mutate(keep=case_when(!keep ~ keep, TRUE ~ !str_detect(cell, token)))
    }
  }
  return(df)
}

# Utterances --------------------------------------------------------------
generate_utts <- function(params){
  utterances <- run_webppl(
    here("model", "default-model", "utterances.wppl"),
    params
  )
  utterances <- utterances %>% map(function(x){x %>% pull(value)}) %>% unlist()
  if(params$save) utterances %>% save_data(params$utts_path)
  return(utterances)
}

# instead of all different utterances, chunk them into categories (for plotting)
chunk_utterances <- function(data, utts_kept=c()){
  levels = c("likely + literal", "conditional", "literal", "conjunction");
  s = paste(utts_kept, collapse="");
  if(str_detect(s, ">") || str_detect(s, "if")){
    levels = c("likely + literal", "other conditional", "literal", "conjunction",
               utts_kept);
  }
  data = data %>% mutate(
    utterance = case_when(
      utterance %in% utts_kept ~ utterance,
      startsWith(utterance, "likely") ~ "likely + literal",
      str_detect(utterance, ">") ~ levels[[2]],
      str_detect(utterance, "and") ~ "conjunction",
      TRUE ~ "literal"
    ),
    utterance = str_replace_all(utterance, "-", "¬"),
    utterance = str_replace(utterance, ">", "->"),
    utterance = factor(utterance, levels=
                         c(map(utts_kept, function(s){
                           s <- str_replace_all(s, "-", "¬")
                           return(str_replace(s, ">", "->"))
                         }),
                         levels)
    )
  );
  return(data)
}

add_graph <- function(data) {
  data = data %>%
    mutate(r = as.character(r),
           r_graph = case_when(str_detect(r, "A implies") ~ "A implies C",
                                str_detect(r, "C implies") ~ "C implies A",
                                T ~ r),
           relation = case_when(
             r == "A || C" ~ "independent",
             startsWith(r, "A") | startsWith(r, "-A") ~ "A%->%C",
             startsWith(r, "C") | startsWith(r, "-C") ~ "C%->%A"),
           relation = as.factor(relation),
           r = factor(r, levels = names(LABELS_R)), 
           r_graph=factor(r_graph,
                           levels = c("A || C", "A implies C", "C implies A")))
  return(data)
}

#@arg tables: with columns AC, A-C, -AC, -A-C
table_to_utts = function(tables, theta){
  tbls = tables %>% add_probs() %>%
    mutate(p_a = p_a >= theta, 
           p_c = p_c >= theta, 
           p_na = p_na >= theta,
           p_nc = p_nc >= theta,
           p_c_given_a = p_c_given_a >= theta,
           p_c_given_na = p_c_given_na >= theta,
           p_nc_given_a = p_nc_given_a >= theta,
           p_nc_given_na = p_nc_given_na >= theta,
           p_a_given_c = p_a_given_c >= theta,
           p_a_given_nc = p_a_given_nc >= theta,
           p_na_given_c = p_na_given_c >= theta,
           p_na_given_nc = p_na_given_nc >= theta,
           p_likely_a = p_likely_a >= 0.5, 
           p_likely_na = p_likely_na >= 0.5, 
           p_likely_c = p_likely_c >= 0.5, 
           p_likely_nc = p_likely_nc >= 0.5, 
           p_ac = AC >= theta,
           p_anc = `A-C` >= theta,
           p_nac = `-AC` >= theta,
           p_nanc = `-A-C` >= theta) %>% 
    mutate_if(is.logical, as.numeric)
  return(tbls)
}

# Probabilities -----------------------------------------------------------
# takes the expected value of column 'p' with probability in column 'prob'
# @args:
#   df_wide; tibble with one bn per row, at least columns: p, prob, level
#   value_str: str describing value, e.g. 'P(A)' for expected val of P(A)
expected_val <- function(df_wide, value_str){
  evs <- df_wide %>% mutate(ev_prod=p * prob)
  evs <- evs %>% group_by(level)
  evs <- evs %>% summarise(ev=sum(ev_prod), .groups="drop") %>% add_column(p=value_str) %>% ungroup()
  
  # fill non-existent levels for plotting
  levels <- evs$level 
  if(is.na(match("prior", levels))){
    evs <- evs %>% add_row(level="prior", ev=0, p=value_str)}
  if(is.na(match("LL", levels))){
    evs <- evs %>% add_row(level="LL", ev=0, p=value_str)}
  if(is.na(match("PL", levels))){
    evs <- evs %>% add_row(level="PL", ev=0, p=value_str)}
  return(evs)
}

compute_cond_prob <- function(distr_wide, prob){
  if(prob=="P(C|A)"){
    distr <- distr_wide %>% mutate(p=`AC`/(`AC`+`A-C`))
  } else if(prob=="P(A|C)"){
    distr <- distr_wide %>% mutate(p=`AC`/(`-AC`+`AC`))
  } else if(prob=="P(C|-A)"){
    distr <- distr_wide %>% mutate(p=`-AC`/(`-AC`+`-A-C`))
  } else if(prob=="P(A|-C)"){
    distr <- distr_wide %>% mutate(p=`A-C`/(`A-C`+`-A-C`))
  
  } else if(prob=="P(-C|A)"){
    distr <- distr_wide %>% mutate(p=`A-C`/(`AC`+`A-C`))
  } else if(prob=="P(-A|C)"){
    distr <- distr_wide %>% mutate(p=`-AC`/(`-AC`+`AC`))
  } else if(prob=="P(-C|-A)"){
    distr <- distr_wide %>% mutate(p=`-A-C`/(`-AC`+`-A-C`))
  } else if(prob=="P(-A|-C)"){
    distr <- distr_wide %>% mutate(p=`-A-C`/(`A-C`+`-A-C`))
  
  }  else{
    stop("not implemented.")
  }
  return(distr)
}
# @arg df: data frame containing columns `AC`, `A-C`, `-AC`
add_probs <- function(df){
  df <- df %>% mutate(p_a = AC + `A-C`, p_c = AC + `-AC`,
                      p_na = `-AC` + `-A-C`, p_nc = `A-C` + `-A-C`) %>%
    mutate(p_c_given_a = AC / p_a,
           p_c_given_na = `-AC` / p_na,
           p_a_given_c = AC / p_c, 
           p_a_given_nc = `A-C` / p_nc, 
           p_nc_given_a = `A-C`/p_a,
           p_nc_given_na = `-A-C`/p_na,
           p_na_given_c = `-AC`/p_c,
           p_na_given_nc = `-A-C`/p_nc,
           p_likely_a = p_a,
           p_likely_na=p_na,
           p_likely_c = p_c,
           p_likely_nc=p_nc
    )
  return(df)
}

# other functions ---------------------------------------------------------
# if identical fields in different sections,the one that occurs first in list
# 'sections' is taken
configure <- function(filename, sections, load_default=TRUE) {
  config = yaml.load_file(filename, eval.expr=T)
  if(load_default) sections = c(sections, "default")
  data <- list()
  for(sec in sections) {
    data <- c(data, config[[sec]])
  }
  return(data)
}

# plotting functions ------------------------------------------------------
plot_speaker_conditions <- function(data) {
  df <- data %>%
    mutate(p=round(as.numeric(p), 2),
           utterance=as.character(utterance), 
           r_graph = factor(r_graph, levels = c("A || C", "A implies C", 
                                                  "C implies A")))
  p <- df %>%
    ggplot(aes(y=utterance, x=p, fill = r_graph)) +
    geom_bar(stat="identity", 
             position = position_dodge(preserve="single")) +
    labs(x="proportion", y="best utterance") + theme_minimal() +
    facet_wrap(~speaker_condition, labeller=label_parsed) +
    theme(axis.text.y=element_text(), legend.position="top",
          legend.key.size = unit(0.75,"line"), 
          ) +
    scale_fill_viridis(discrete=T, name="causal relation", labels = LABELS_R) +
    guides(fill=guide_legend(reverse = TRUE))
  return(p)
}

data_cp_plots <- function(params, data=NA){
  if(is.na(data)) data <- read_rds(params$target)
  data <- data %>% ungroup() %>% group_by(bn_id, level) %>%
    select(-p_delta, -p_rooij, -p_diff)
  data.wide <- data %>%
    pivot_wider(names_from = "cell", values_from = "val") %>% 
    compute_cond_prob("P(-C|-A)") %>%  rename(`-C_-A` = p) %>% 
    compute_cond_prob("P(A|C)") %>% rename(`A_C` = p)
  
  # Expected values for P(-C|-A) and P(A|C) and for causal nets
  ev_nc_na = data.wide %>% rename(p=`-C_-A`) %>% expected_val("P(-C|-A)")
  ev_a_c = data.wide %>% rename(p=`A_C`) %>% expected_val("P(A|C)") 
  ev_probs <- bind_rows(ev_a_c, ev_nc_na) %>%
    mutate(level=factor(level, levels=c("PL", "LL", "prior")),
           p=str_replace_all(p, "-", "¬"), 
           p=str_replace_all(p, "A", "a"),
           p=str_replace_all(p, "C", "c")) %>%
    rename(val_type=p) %>% add_column(val="p")
  
  ev_rels = data.wide %>% ungroup() %>% group_by(level, r) %>% 
    summarise(ev=sum(prob), .groups="drop_last") %>%
    mutate(level=factor(level, levels=c("PL", "LL", "prior"))) %>% 
    add_graph() %>% rename(val_type=r) %>% add_column(val="relations")
  
  data <- bind_rows(ev_probs, ev_rels) %>% 
    mutate(level=factor(level, levels=c("PL", "LL", "prior"), 
                        labels=c("Pragmatic listener", "Literal listener", "Prior")))
  return(data)  
}

plot_cp_relations <- function(cp.relations){
  p.relations <- cp.relations %>% 
    ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar(position=position_stack(), stat="identity") +
    scale_fill_manual(name="causal relation", labels=LABELS_r,
                      values = COLS_r) +
    labs(x="Degree of belief", y="") +
    theme_minimal() +
    theme(legend.position="top", legend.key.size = unit(0.75,"line")) +
    guides(fill=guide_legend(reverse = TRUE))
  return(p.relations)
}

plot_cp_probs <- function(data.cp){
  p.probs <- data.cp %>% filter(val=="p") %>% 
    ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    scale_fill_viridis(discrete=TRUE, name="Probability") +
    labs(x="Degree of belief", y="") +
    theme_minimal() +
    theme(legend.position="top", legend.key.size = unit(0.75,"line")) +
    guides(fill=guide_legend(reverse = TRUE))
  return(p.probs)
}

# Acceptability/Assertability conditions ----------------------------------
# p_rooij: (P(e|i) - P(e|¬i)) / (1-P(e|¬i))
# p_delta: P(e|i) - P(e|¬i)
acceptability_conditions <- function(data_wide){
  df <- data_wide %>% compute_cond_prob("P(C|A)") %>% rename(p_c_given_a=p) %>% 
    compute_cond_prob("P(C|-A)") %>% rename(p_c_given_na=p) %>%
    mutate(p_delta=p_c_given_a - p_c_given_na,
           p_nc_given_na=1-p_c_given_na,
           p_rooij=p_delta/(1-p_c_given_na),
           pc=`AC` + `-AC`,
           p_diff=round(p_c_given_a - pc, 5)) %>%
    select(-p_nc_given_na, -p_c_given_a, -p_c_given_na, -pc)
  return(df)
}