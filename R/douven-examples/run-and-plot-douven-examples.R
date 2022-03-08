library(ggplot2)
library(here)
library(rwebppl)
library(tidyverse)
source("R/helper-functions.R")
source("R/helpers-webppl.R")

# Run Model ---------------------------------------------------------------
get_douven_data <- function(params){
  if(!dir.exists(params$target_dir)) {
    dir.create(params$target_dir, recursive=TRUE)
  }
  data.all <- run_webppl(params$model_path, params) 
  vals = names(data.all)

  data.with_cells = data.all[vals[vals %in% c("prior", "LL", "PL")]] %>%
    structure_listener_data(params)
  
  data = data.with_cells %>%
    distinct_at(vars(c(rowid, bn_id, r_graph, prob, level)))

  bns = data.all$bns %>% unnest(c(table.probs, table.support)) %>% 
    rename(cell=table.support, val=table.probs)
  return(list(data=data, bns=bns, all=data.all)) #, data_with_cells=data.with_cells))
}

# condition each state (a joint prob.distribution) on var_conditioned
condition_results <- function(data, var_conditioned){
  if(!is.na(var_conditioned)){
    df <- data %>% 
      filter_vars(var_conditioned) %>%
      group_by(bn_id, keep) %>% 
      mutate(val=val/sum(val), val=case_when(keep ~ val, TRUE ~ 0)) %>%
      ungroup() %>% select(-keep) %>% 
      mutate(bn_id =  paste(bn_id, var_conditioned, sep="-"))
  }
  return(df)
}

get_marginal = function(bns, vars){
  df = bns %>% filter_vars(vars) %>% filter(keep) %>% 
    select(-keep) %>% group_by(bn_id, r_graph) %>%
    summarize(p=sum(val), .groups = "drop_last") %>% 
    add_column(marginal=str_to_lower(paste(vars, collapse="")))
  return(df)
}

add_conditioned_bn <- function(results, conditioned_on){
  df <- bind_rows(results %>% mutate(conditioned = FALSE), 
                  results %>% mutate(bn_id=paste(bn_id, conditioned_on, sep="-"),
                                     level=paste(level, "observed", sep="-"),
                                     conditioned = TRUE))
  return(df)
}

get_ev_marginal = function(df){
  df.evs = df %>% 
    group_by(conditioned, level, bn_id, marginal, example, r_graph) %>% 
    summarize(ev = sum(prob * p), .groups = "drop_last")
  return(df.evs)
}

# 1. Skiing ---------------------------------------------------------------
params.skiing <- configure("config.yml", c("skiing"))
skiing <- get_douven_data(params.skiing)
skiing.results = skiing$data
skiing.bns = skiing$bns

# condition each state on listener's observation which results in "new" states
skiing.bns_observed <- condition_results(skiing.bns, params.skiing$condition_on)
skiing.bns_all = bind_rows(skiing.bns, skiing.bns_observed)

# compute P(e) for each state
skiing.bns_marginal = skiing.bns_all %>% get_marginal(c("E")) %>% 
  add_column(example="skiing")

# add conditioned Bayes nets
skiing.data = add_conditioned_bn(skiing.results, params.skiing$condition_on)

df.skiing = left_join(skiing.data, skiing.bns_marginal, by=c("bn_id", "r_graph"))
df.skiing.evs = df.skiing %>% get_ev_marginal()

# ev for causal nets does not change based on listener's observation,
# only the expected value of P(e)

# 2. Garden Party ---------------------------------------------------------
params.gp <- configure("config.yml", c("garden_party"))
garden_party <- get_douven_data(params.gp)
gp.results <- garden_party$data
gp.bns <- garden_party$bns

gp.bns_observed <- condition_results(gp.bns, params.gp$condition_on)
gp.bns_all <- bind_rows(gp.bns, gp.bns_observed)

# compute P(e) for each state
gp.bns_marginal = gp.bns_all %>% get_marginal(c("D")) %>% 
  add_column(example="gardenParty")

# add conditioned Bayes nets
garden_party.data = add_conditioned_bn(gp.results, params.gp$condition_on)
df.garden_party = left_join(garden_party.data, gp.bns_marginal, 
                            by=c("bn_id", "r_graph"))
df.garden_party.evs = df.garden_party %>% get_ev_marginal()
df.garden_party.evs

# 3. Sundowners -----------------------------------------------------------
# condition on has NA-value (no observation made here)
params.sundowners <- configure("config.yml", c("sundowners"))
sundowners <- get_douven_data(params.sundowners)
sundowners.results <- sundowners$data
sundowners.bns <- sundowners$bns

sundowners.bns_marginal_r = sundowners.bns %>% get_marginal(c("R")) %>% 
  add_column(example="sundowners")
sundowners.bns_marginal_rs = sundowners.bns %>% get_marginal(c("R", "S")) %>% 
  add_column(example="sundowners")
sundowners.bns_marginals <- bind_rows(sundowners.bns_marginal_r,
                                      sundowners.bns_marginal_rs)

df.sundowners <- left_join(sundowners.results, sundowners.bns_marginals, 
                           by=c("bn_id", "r_graph")) %>% 
  add_column(conditioned=FALSE)
df.sundowners.evs = df.sundowners %>% get_ev_marginal()

# Joint sundowners+skiing+garden party data
douven.data <- bind_rows(df.skiing.evs, df.garden_party.evs, df.sundowners.evs)

plot_douven_cases <- function(data){
  # for labels
  df = tribble(~limits, ~labels, ~example,
               "R||S", "R,S indep.", "sundowners", 
               "R > W > S", expression(R %->% W %->% S), "sundowners",
               "E || S>C", expression(paste("E, S indep.,", S %->% C)), "skiing",
               "E>S>C", expression(E %->% S %->% C), "skiing",
               "D || G>-S", expression(paste("D, G indep.,",G %->%"","¬S")), "gardenParty", 
               "D>G>-S", expression(paste(D %->% G %->%"", "¬S")), "gardenParty")
  ex.in = data$example %>% unique()
  max.x = round(data$ev %>% max(), 1)
  df = df %>% filter(example %in% ex.in)
  
  y_limits <- c("PL-observed", "PL", "LL", "prior")
  y_labels <- c(
    paste(strwrap("Pragmatic interpretation with listener's observation",
                  width=28), collapse="\n"),
    paste(strwrap("Pragmatic interpretation w/o listener's observation",
                  width=28), collapse="\n"),
    paste(strwrap("Literal interpretation", width=28), collapse="\n"),
    "Prior Belief"
  )
  
  p <- data %>% mutate(level=as.factor(level)) %>%
    ggplot() +
    geom_bar(mapping = aes(y=level, x=ev, fill=r_graph), stat="identity",
             position="stack") +
    # facet_wrap(~example) + 
    facet_wrap(~marginal, labeller=
                 as_labeller(c(`r`="Sundowners: P(r)", `rs`="Sundowners: P(r,s)",
                               `e`="Skiing: P(e)", `d`="Garden Party: P(d)"))) +
    scale_y_discrete(
      name = "",
      limits = y_limits[y_limits %in% (data$level %>% unique())],
      labels = y_labels[y_limits %in% (data$level %>% unique())]
    ) +
    scale_fill_viridis(name="causal relation", limits = df$limits, 
                       labels = df$labels, discrete = TRUE) +
    scale_x_continuous(breaks = seq(0, max.x, by=0.1)) +
    labs(x="Expected value", title="") +
    theme_minimal() +
    theme(legend.position="top",
          panel.spacing = unit(2, "lines"), 
          axis.text = element_text(size=12),
          legend.key.size = unit(0.75,"line"),
          axis.text.y = element_text(hjust = 0, size=12),
          legend.text = element_text(size=12))

  if(data$example %>% unique() != "sundowners") {
    p <- p + guides(fill=guide_legend(reverse = TRUE))
  }
  return(p)
}

for(ex in c("gardenParty", "sundowners","skiing")) {
  p <- douven.data %>% filter(example == !!(ex)) %>% plot_douven_cases()
  ggsave(here("data", "douven-examples", paste(ex, ".png", sep="")), p,
         width=7, height=4, dpi="print")
}
message(paste('saved plots to', here("data", "douven-examples")))
