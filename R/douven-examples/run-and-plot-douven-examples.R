library(ggplot2)
library(here)
library(rwebppl)
library(tidyverse)
library(tikzDevice)
library(viridis)
source("R/helper-functions.R")
source("R/helpers-webppl.R")

theme_set(theme_minimal(base_size=12) + theme(legend.position = "top"))
read_data = TRUE

# Run Model ---------------------------------------------------------------
get_douven_data <- function(params, read = TRUE){
  if(!dir.exists(params$target_dir)) {
    dir.create(params$target_dir, recursive=TRUE)
  }
  fn_data = paste(params$target_dir, paste(params$context, "-all.rds", sep=""), sep=fs)
  if(read && file.exists(fn_data)) {
    data.all <- readRDS(fn_data)
    params$save <- FALSE
  } else {
    data.all <- run_webppl(params$model_path, params) 
    if(params$save) save_data(data.all, fn_data)
  }
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
skiing <- get_douven_data(params.skiing, read = read_data)
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
garden_party <- get_douven_data(params.gp, read = read_data)
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
sundowners <- get_douven_data(params.sundowners, read = read_data)
sundowners.results <- sundowners$data
sundowners.bns <- sundowners$bns

sundowners.speaker = sundowners$all$speaker %>% 
  unnest(c(speaker.probs, speaker.support)) %>% 
  rename(prob = speaker.probs, u = speaker.support) %>% 
  add_column(level = "speaker")

sundowners.surprise_u = left_join(sundowners.speaker, 
          sundowners.results %>% filter(level == "prior") %>% 
            rename(prior_bn_id = prob) %>% dplyr::select(bn_id, prior_bn_id))

sundowners.surprise_u %>% group_by(u) %>% 
  summarize(surprise_u = sum(prob * prior_bn_id))


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


# Join data from all three examples ---------------------------------------
douven.data <- bind_rows(df.skiing.evs, df.garden_party.evs, df.sundowners.evs)

plot_dir = here("data", "douven-examples")
dat <- douven.data %>% 
  mutate(state = factor(r_graph, levels = c("R||S", "R > S", "E || S>C", 
                                            "E>S>C", "D || G>-S", "D>G>-S"), 
                        labels = c("$s_{ind}$", "$s_{dep}$", "$s_{ind}$", 
                                   "$s_{dep}$", "$s_{ind}$", "$s_{dep}$")),
         r_graph = factor(r_graph, levels = c("R||S", "R > S", "E || S>C", 
                                              "E>S>C", "D || G>-S", "D>G>-S"), 
                          labels = c("$R\\indep S$", 
                                     "$R\\rightarrow S$",
                                     "$E\\indep S$", "$E\\rightarrow S$",
                                     "$D \\indep G$", "$D \\rightarrow G$")))
max_x = dat %>%  group_by(example, level) %>% 
  summarize(max = sum(round(ev, 1) + 0.1), .groups = "drop_last")
df.pragmatic = dat %>% filter(level %in% c("PL", "LL", "prior")) %>%
  mutate(level = factor(level, levels=c("PL", "LL", "prior"), 
                        labels=c("Pragmatic interpretation", 
                                 "Literal interpretation", "Prior belief")))
df.obs = dat %>% filter(level == "PL-observed" | level == "PL") %>% 
  mutate(level=factor(level, levels=c("PL-observed", "PL"), 
                      labels=c("Pragmatic interpretation\nwith listener's observation", 
                               "Pragmatic interpretation\n w/o listener's observation")))

# PLOTS -------------------------------------------------------------------
# 1.Skiing Example --------------------------------------------------------
ex <- "skiing"
ex.max_x <- max_x %>% filter(example==ex & level %in% c("prior", "LL", "PL")) %>%
  pull(max) %>% max()
df = df.pragmatic %>% filter(example == ex)
tikz(paste(plot_dir, paste(ex, "pragmatic.tex", sep="-"), sep=fs),
     width = 5, height = 2.5, standAlone = FALSE,
     packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                  "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
df %>% ggplot() +
  geom_bar(mapping = aes(y=level, x=ev, fill=state), stat="identity", position="stack") +
  facet_wrap(~marginal, labeller=as_labeller(c(`e`="$\\mathbb{E}[P^{(s)}(e)]$"))) +
  scale_fill_viridis(name="state", discrete = TRUE) +
  scale_x_continuous(breaks = seq(0, ex.max_x, by=0.1)) +
  labs(x="", y="", title="") +
  theme(panel.spacing = unit(2, "lines"), legend.key.size = unit(0.75,"line"), 
        legend.spacing.x = unit(1.25, "line"), legend.position = "top") +
  guides(fill=guide_legend(reverse = TRUE))
dev.off()

# interpretation with listener's observation
ex.max_x <- max_x %>% filter(example==ex & level %in% c("PL-observed", "PL")) %>%
  pull(max) %>% max()
df = df.obs %>% filter(example == ex)
tikz(paste(plot_dir, paste(ex, "with-observation.tex", sep="-"), sep=fs), 
     width = 5, height = 2.5, standAlone = FALSE, 
     packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                  "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
df %>% ggplot() + 
  geom_bar(mapping = aes(y=level, x=ev, fill = state), stat="identity", 
           position = "stack") + 
  facet_wrap(~marginal, labeller=as_labeller(c(`e`="$\\mathbb{E}[P^{(s)}(e\\mid c)]$"))) + 
  scale_fill_viridis(name="state", discrete = TRUE) +
  scale_x_continuous(breaks = seq(0, ex.max_x, by=0.1)) +
  labs(x="", y="", title="") +
  theme(panel.spacing = unit(2, "lines"), legend.key.size = unit(0.75,"line"),
        legend.spacing.x = unit(1.25, "line"), legend.position="top") + 
  guides(fill=guide_legend(reverse = TRUE))
dev.off()

# 2.Garden Party Example --------------------------------------------------
ex <- "gardenParty"
ex.max_x <- max_x %>% filter(example==ex & level %in% c("prior", "LL", "PL")) %>%
  pull(max) %>% max()
df = df.pragmatic %>% filter(example == ex)
tikz(paste(plot_dir, paste(ex, "pragmatic.tex", sep="-"), sep=fs),
     width = 5, height = 2.5, standAlone = FALSE,
     packages = c("\\usepackage{tikz, amssymb, amsmath}",
                  "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
df %>% ggplot() +
  geom_bar(mapping = aes(y=level, x=ev, fill=state), stat="identity", position="stack") +
  facet_wrap(~marginal, labeller=as_labeller(c(`d`="$\\mathbb{E}[P^{(s)}(d)]$"))) +
  scale_fill_viridis(name="state", discrete = TRUE) +
  scale_x_continuous(breaks = seq(0, ex.max_x, by=0.1)) +
  labs(x="", y="", title="") +
  theme(panel.spacing = unit(2, "lines"), legend.key.size = unit(0.75,"line"),
        legend.spacing.x = unit(1.25, "line")) +
  guides(fill=guide_legend(reverse = TRUE))
dev.off()
# interpretation with listener's observation
ex.max_x <- max_x %>% filter(example==ex & level %in% c("PL-observed", "PL")) %>%
  pull(max) %>% max()
df = df.obs %>% filter(example == ex)
tikz(paste(plot_dir, paste(ex, "with-observation.tex", sep="-"), sep=fs), 
     width = 5, height = 2.5, standAlone = FALSE, 
     packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                  "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
df %>% ggplot() + 
  geom_bar(mapping = aes(y=level, x=ev, fill = state), stat="identity", 
           position = "stack") + 
  facet_wrap(~marginal, labeller=as_labeller(c(`d`="$\\mathbb{E}[P^{(s)}(d\\mid s)]$"))) + 
  scale_fill_viridis(name="state", discrete = TRUE) +
  scale_x_continuous(breaks = seq(0, ex.max_x, by=0.1)) +
  labs(x="", y="", title="") +
  theme(panel.spacing = unit(2, "lines"), legend.key.size = unit(0.75,"line"),
        legend.spacing.x = unit(1.25, "line"), legend.position="top") + 
  guides(fill=guide_legend(reverse = TRUE))
dev.off()

# 3.Sundowners Example ----------------------------------------------------
ex <- "sundowners"
ex.max_x <- max_x %>% filter(example==ex & level %in% c("prior", "LL", "PL")) %>%
  pull(max) %>% max()

df = df.pragmatic %>% filter(example == ex)
# add inferred relations
df = bind_rows(df, 
               sundowners.results %>% group_by(level, r_graph) %>% rename(ev=prob) %>% 
                 mutate(r_graph = factor(r_graph, levels = c("R||S", "R > S", "E || S>C", 
                                                             "E>S>C", "D || G>-S", "D>G>-S"), 
                                         labels = c("$R\\indep S$", 
                                                    "$R\\rightarrow S$",
                                                    "$E\\indep S$", "$E\\rightarrow S$",
                                                    "$D \\indep G$", "$D \\rightarrow G$")),
                        marginal = "causal relation",
                        level=factor(level, levels=c("PL", "LL", "prior"), 
                                     labels=c("Pragmatic interpretation", 
                                              "Literal interpretation", "Prior belief")), 
                        conditioned = FALSE, example = ex) %>% 
                 dplyr::select(-rowid)
)

tikz(paste(plot_dir, paste(ex, "pragmatic.tex", sep="-"), sep=fs),
     width = 6, height = 2.5, standAlone = FALSE,
     packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                  "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
df %>% ggplot() +
  geom_bar(mapping = aes(y=level, x=ev, fill=r_graph), stat="identity", position="stack") +
  facet_wrap(~marginal, labeller=as_labeller(
    c(`r`="$\\mathbb{E}[P^{(s)}(r)]$", `rs`="$\\mathbb{E}[P^{(s)}(r,s)]$", 
      `causal relation` = "$\\mathbb{E}[P^{(s)}(\\mathcal{R})]$")
  )) +
  scale_fill_viridis(name="causal relation $\\mathcal{R}$", discrete = TRUE, 
                     labels = LABELS_R_TEX) +
  scale_x_continuous(labels = number_format(accuracy = 0.1)) + 
  labs(x="", y="", title="") +
  theme(panel.spacing = unit(2, "lines"), legend.key.size = unit(0.75,"line"),
        legend.spacing.x = unit(1.25, "line")) +
  guides(fill=guide_legend(reverse = TRUE))

dev.off()

message(paste('saved plots to', plot_dir))
