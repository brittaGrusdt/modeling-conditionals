library(latex2exp)
library(RColorBrewer)
library(yaml)

fs = .Platform$file.sep

LABELS_R_TEX = c(
  "A || C" = "$A\\indep C$",
  "A implies C" = "$A\\rightsquigarrow C$",
  "-A implies C" = "$A\\rightsquigarrow C$",
  "C implies A" = "$C\\rightsquigarrow A$",
  "-C implies A" = "$C\\rightsquigarrow A$"
)
LABELS_cp_probs = c(
  "P(A|C)" = "$\\mathbb{E}[P^{(s)}(a\\mid c)]$", 
  "P(-C|-A)" = "$\\mathbb{E}[P^{(s)}(\\neg c \\mid \\neg a)]$",
  "P(C|A)" = "$\\mathbb{E}[P^{(s)}(c \\mid a)]$"
)
LABELS_r = c(
  "A implies C" = parse(text=TeX('$A^+C^+$')),
  "C implies A" = parse(text=TeX('$C^+A^+$')),
  "-A implies C" = parse(text=TeX('$A^-C^+$')),
  "-C implies A" = parse(text=TeX('$C^-A^+$')),
  "A || C" = expression(paste("A,C independent"))
)
LABELS_sp_levels <-
  c(`prior` = expression(paste("Prior")),
    `literal-speaker` = expression(paste("Literal speaker")),
    `pragmatic-speaker` = expression(paste("Pragmatic speaker")))

brewer_pairs = brewer.pal(5, name = 'Set2')
COLS_r <- c(
  "C implies A" = brewer_pairs[1],
  "A implies C" = brewer_pairs[2],
  "A || C" = brewer_pairs[5],
  "-C implies A" = brewer_pairs[3],
  "-A implies C" = brewer_pairs[4])

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
  levels = c("likely + literal", "conditional", "literal", "conjunction", 
             "disjunction");
  s = paste(utts_kept, collapse="");
  if(str_detect(s, ">") || str_detect(s, "if")){
    levels = c("likely + literal", "other conditional", "literal", "conjunction",
               "disjunction", utts_kept);
  }
  data = data %>% mutate(
    utterance = case_when(
      utterance %in% utts_kept ~ utterance,
      startsWith(utterance, "likely") ~ "likely + literal",
      str_detect(utterance, ">") ~ levels[[2]],
      str_detect(utterance, "and") ~ "conjunction",
      str_detect(utterance, "or") ~ "disjunction",
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
           r = factor(r, levels = names(LABELS_R_TEX)), 
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
  evs <- evs %>% summarise(ev=sum(ev_prod), .groups="drop") %>% 
    add_column(p=value_str) %>% ungroup()
  
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

compute_cond_prob <- function(distr_wide, prob, name_p=NULL){
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
  if(!is.null(name_p)){
    distr <- distr %>% rename(!!name_p := p)
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
           r_graph = factor(r_graph, 
                            levels = c("A || C", "A implies C", "C implies A")))
  p <- df %>%
    ggplot(aes(y=utterance, x=p, fill = r_graph, color = r_graph)) +
    geom_bar_pattern(
      aes(pattern = r_graph, pattern_linetype = r_graph),
      pattern_density = 0.1, pattern_colour = 'white', pattern_fill = 'white',
      position = position_dodge2(preserve = "single"), stat = "identity",
      pattern_key_scale_factor = 0.25
    ) +
    facet_wrap(~speaker_condition) +
    scale_fill_brewer(name = "causal relation $R$", palette = 'Set2', labels = LABELS_R_TEX) +
    scale_color_brewer(name = "causal relation $R$", palette = 'Set2', labels = LABELS_R_TEX) +
    scale_pattern_linetype_discrete(name = "causal relation $R$", labels = LABELS_R_TEX) +
    scale_pattern_discrete(name = "causal relation $R$", labels = LABELS_R_TEX) +
    theme(axis.text.y=element_text(), legend.spacing.x = unit(1.25, "line"), 
          panel.spacing = unit(1.5, "lines")) +
    labs(x = "proportion", y="") +
    guides(fill=guide_legend(reverse = TRUE, title.position = "left", spacing.x = 1), 
           color=guide_legend(reverse = TRUE, title.position = "left", spacing.x = 1), 
           pattern=guide_legend(reverse = TRUE, title.position = "left", spacing.x = 1),
           pattern_linetype=guide_legend(reverse = TRUE, title.position = "left", spacing.x = 1))
  return(p)
}

data_cp_plots <- function(params, data=NA){
  if(is.na(data)) data <- read_rds(params$target)
  data <- data %>% ungroup() %>% group_by(bn_id, level) %>%
    select(-p_delta, -p_rooij, -p_diff)
  data.wide <- data %>%
    pivot_wider(names_from = "cell", values_from = "val") %>% 
    compute_cond_prob("P(-C|-A)", name_p = "-C_-A") %>% 
    compute_cond_prob("P(A|C)", name_p = "A_C") %>% 
    compute_cond_prob("P(C|A)", name_p = "C_A")
  
  # Expected values for P(-C|-A) and P(A|C) and for causal nets
  ev_nc_na = data.wide %>% rename(p=`-C_-A`) %>% expected_val("P(-C|-A)")
  ev_a_c = data.wide %>% rename(p=`A_C`) %>% expected_val("P(A|C)") 
  ev_c_a = data.wide %>% rename(p=`C_A`) %>% expected_val("P(C|A)") 
  ev_probs <- bind_rows(ev_a_c, ev_nc_na, ev_c_a) %>%
    mutate(level=factor(level, levels=c("PL", "LL", "prior"))) %>% 
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

plot_cp_relations <- function(cp.relations, fn, w, h){
  p.relations <- cp.relations %>% 
    mutate(val_type = factor(val_type, 
                             levels=c("A || C", "C implies A", "-C implies A", 
                                      "-A implies C", "A implies C")),
           r_graph = factor(r_graph, levels = c("A || C", "C implies A", "A implies C"))
           ) %>% 
    ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar_pattern(aes(pattern = r_graph, pattern_linetype = r_graph),
                     pattern_density = 0.1, pattern_colour =  'white',
                     pattern_fill = 'white',
                     pattern_key_scale_factor = 0.25,
                     position=position_stack(), stat="identity") +
    scale_fill_manual(name="instance $r$", labels=LABELS_r, values = COLS_r,
                      guide = guide_legend(reverse = TRUE, 
                                           override.aes = 
                                             list(pattern = "none", 
                                                  pattern_linetype = "none"))) +
    scale_pattern_linetype_discrete(name="relation $R$", labels=LABELS_R_TEX,
                                    guide = guide_legend(reverse = TRUE)) + 
    scale_pattern_discrete(name = "relation $R$", labels = LABELS_R_TEX,
                           guide = guide_legend(reverse = TRUE)) +
    labs(x="Degree of belief", y="") +
    theme(legend.key.size = unit(0.75,"line"), legend.box = 'vertical',
          legend.spacing.x = unit(1.25, "line"))
  return(p.relations)
}

plot_cp_probs <- function(dat, fn, w, h){
  p.probs <- dat %>% filter(val=="p") %>% 
    ggplot(aes(y=level, x=ev, fill = val_type, color = val_type)) + 
    geom_bar_pattern(aes(pattern = val_type, pattern_linetype = val_type),
                     pattern_density = 0.1, pattern_colour =  'white', 
                     pattern_fill = 'white', pattern_key_scale_factor = 0.25, 
                     stat="identity", position=position_dodge()) +
    scale_pattern_linetype_discrete(name = "value ", labels = LABELS_cp_probs) +
    scale_pattern_discrete(name = "value ", labels = LABELS_cp_probs) +
    scale_fill_brewer(palette = 'Set2', name="value ", labels = LABELS_cp_probs) +
    scale_color_brewer(palette = 'Set2', name="value ", labels = LABELS_cp_probs) +
    labs(x="Degree of belief", y="") +
    theme(legend.key.size = unit(0.75,"line"), legend.spacing.x = unit(1.25, "line")) +
    guides(fill=guide_legend(reverse = TRUE), color=guide_legend(reverse = TRUE), 
           pattern_linetype = guide_legend(reverse = TRUE),
           pattern = guide_legend(reverse = TRUE))
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
