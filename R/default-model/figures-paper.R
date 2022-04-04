library(config)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(grid)
library(latex2exp)
library(scales)
library(tidyverse)
library(tikzDevice)
library(unikn)
library(viridis)
library(viridisLite)
library(xtable)
source("R/helper-functions.R")
source("R/default-model/helpers-tables.R")

# target <- "test-config"
target <- "paper-resubmission-config"
# target <- "with-or"

# seed <- ""
seed <- "1644922777.24488"

# figures for specific param config or plot for all combinations
theta <- 0.9
alpha <- 3
pars <- paste(alpha, theta, sep="--")
# pars <- "all"

plot_tables = TRUE

theme_set(theme_minimal(base_size=12) + theme(legend.position = "top"))
# funcitons ---------------------------------------------------------------
plot_accept_conditions <- function(dat, cond_tex_str){
  data <- dat %>% 
    mutate(condition = case_when(r == "A || C" & level=="prior" ~ 
                                   as.double(scientific(condition, digits=1)),
                                 T ~ condition),
           r = factor(r, levels = names(LABELS_r)), 
           level_f = factor(level, levels = names(LABELS_sp_levels)),
           r_f = r) 
  # important to use function levels
  levels(data$r_f) = LABELS_r
  levels(data$level_f) = LABELS_sp_levels
  
  plots = list()
  levels = data$level %>% unique()
  for(i in seq(1, length(levels))){
    df <- data %>% filter(level == levels[[i]])
    p <- df %>% 
      ggplot(aes(y=condition, x=r)) +
      geom_boxplot(outlier.size = 0.1) +
      facet_wrap(~r_f, scales="free", ncol=4, labeller = label_parsed) +
      labs(y=TeX(cond_tex_str), x="", title=LABELS_sp_levels[levels[[i]]]) +
      theme_minimal() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    if(levels[[i]] == "pragmatic-speaker") {
      p <- p + scale_y_continuous(limits = c(df$condition %>% min(), 1))
    }
    if(levels[[i]] == "prior"){
      p <- p + theme(legend.position="bottom",
                     legend.key.size = unit(0.75,"line"))
    } else {
      p <- p + theme(legend.position="none")
    }
    plots[[i]] = p
  }
  return(plots)
}

# setup -------------------------------------------------------------------
fs = .Platform$file.sep
data_dir = here("data", "default-model", target)
if(seed != "") data_dir <- paste(data_dir, fs, seed, sep="")
if(pars == "all") {
  subdirs <- list.files(data_dir, recursive = F, full.names = F, 
                        pattern = ".*--.*")
} else {
  subdirs <- c(pars)
}


# Run plots ---------------------------------------------------------------
csv_data <- tibble()
for(par_config in subdirs) {
  params = read_rds(file.path(data_dir, par_config, "params-prior-ll-pl.rds",
                              fsep=fs))
  theta = params$theta
  alpha = params$alpha
  pars <- str_split(par_config, "--")[[1]]
  if(pars[1] != alpha | pars[2] != theta) stop("error with subdir names + params")
  data <- read_rds(params$target) %>% add_graph()

  # get params for other conditions than for default prior-ll-pl
  params.speaker <- read_rds(paste(params$target_dir, "params-speaker.rds", sep=fs))
  params.prior <- read_rds(paste(params$target_dir, "params-priorN.rds", sep=fs))
  params.sp_literal <- read_rds(paste(params$target_dir, "params-speaker-literal.rds", sep=fs))

  plot_dir = params$plot_dir
  
  # Figure 1 ----------------------------------------------------------------
  tables = data %>% filter(level == "prior") %>% 
    ungroup() %>% dplyr::select(r, relation, bn_id, cell, val) %>% 
    group_by(bn_id) %>% pivot_wider(names_from="cell", values_from="val") %>% 
    mutate(vs=list(c("AC", "A-C", "-AC", "-A-C")),
           ps=list(c(`AC`, `A-C`, `-AC`, `-A-C`))) 
  tables_df = tables %>%
    select(-`AC`, -`A-C`, -`-AC`, -`-A-C`) %>% 
    mutate(bn_id.tmp=bn_id) %>% 
    separate(bn_id, into=c("r", "t0", "t1", "t2", "t3"), sep="_") %>% 
    unite("table_id", "t0", "t1", "t2", "t3", sep="_") %>% 
    rename(bn_id=bn_id.tmp)
  
  # only plot tables once (identical for the same seed, indep. of theta/alpha)
  if(plot_tables && which(subdirs == par_config) == 1) {
    table_dir <- paste(params$seed_dir, "figs-tables", sep = fs)
    if(!dir.exists(table_dir)) dir.create(table_dir)
    plot_tables_relations("", table_dir, w=3, h=3.5, tables_df)
  }
  
  # Informativeness of utterances (Table Appendix)
  tables_utts = table_to_utts(tables, theta) %>% ungroup() %>% 
    dplyr::select(-vs, -ps, -bn_id, -AC, -`A-C`, -`-AC`, -`-A-C`) %>%
    rowid_to_column() %>% group_by(rowid) %>%
    arrange(relation, r)
  
  df.tbls_utts = tables_utts %>% 
    pivot_longer(cols=c(-rowid, -r, -relation), 
                 names_to="utterance", values_to="assertable") %>% 
    group_by(utterance, relation) %>% 
    summarize(n_assertable = sum(assertable), .groups = "drop_last") %>%
    pivot_wider(names_from="relation", values_from="n_assertable") %>% 
    mutate(order = case_when(str_detect(utterance, "likely") ~ 1,
                             str_detect(utterance, "given") ~ 2,
                             utterance %in% c("p_a", "p_c", "p_na", "p_nc") ~ 3,
                             T ~ 4)) %>% 
    arrange(order) %>% dplyr::select(-order)
  sums = rowSums(df.tbls_utts %>% ungroup() %>% dplyr::select(-utterance))
  df.tbls_utts <- df.tbls_utts %>% add_column(total = sums)
  tbls.utts.tex = xtable(df.tbls_utts, digits=c(0,0,0,0,0,0))
  if(par_config == "3--0.9") print(tbls.utts.tex, include.rownames = FALSE)
  
  # Fig.Appendix 
  # tables where combined conditional 'A->C and ¬A->C' is assertable
  # vs. tables where 'C' is assertable
  literal_combined_ifs = tables %>%
    mutate(a=AC+`A-C`, ifac=AC/a>=theta, na=`-AC`+`-A-C`, ifnac=`-AC`/na>=theta,
           c=AC+`-AC` >= theta) %>% 
    group_by(r, ifac, ifnac, c)
  # for all states where both conditionals are true, C is also true?
  if(nrow(literal_combined_ifs %>% filter(ifac & ifnac & !c))){
    stop("both conditionals true, but not C.")
  }
  literal_combined_ifs <- literal_combined_ifs %>% 
    dplyr::count() %>% filter(ifac &  ifnac) %>% ungroup() %>% 
    mutate(ratio=n/sum(n)) %>% dplyr::select(-ifac, -ifnac, -c) %>% 
    add_column(utterance = "combined ifs")
  
  literal_c = tables %>%
    filter(AC+`-AC`>=theta) %>% group_by(r) %>% 
    dplyr::count() %>% ungroup() %>% 
    mutate(ratio=n/sum(n)) %>% 
    add_column(utterance = "C")
  
  df.literal <- bind_rows(literal_combined_ifs, literal_c) %>% 
    add_graph() %>% group_by(relation, utterance) %>% 
    summarize(n=sum(n), ratio=sum(ratio), .groups = "drop_last") %>% 
    mutate(utterance = factor(utterance, levels = c("C", "combined ifs")))
  
  p.combined_ifs = df.literal %>% ggplot(aes(x=relation, y=n)) + 
    geom_bar(aes(fill=utterance), stat="identity", position=position_dodge()) + 
    theme_minimal() +
    theme(legend.position="top") +
    scale_x_discrete(labels = function(label) parse(text=label)) + 
    labs(x="causal relation", 
         y = paste(strwrap("Number of states", 
                           width=28), collapse="\n")) +
    scale_fill_viridis(discrete=TRUE,
      name = "assertable utterance", option="D",
      labels = c(C="C", 
                 `combined ifs`= 
                   expression((A%->%C)~symbol("\331")~(symbol("\330")~A%->%C))))
  
  # p.combined_ifs
  ggsave(paste(plot_dir, "combined-ifs-c-true.png", sep=fs), p.combined_ifs,
         width=7, height=2.5)
  
  # Figure 2 ----------------------------------------------------------------
  # speaker results unique sampled states from n samples 
  # (rowid is unique, bn_id not necessarily)
  data.speaker.all <- read_rds(params.speaker$target) %>%
    select(-level, -p_delta, -p_diff) %>% ungroup()
  
  # Fig.Appendix: expected value to choose each utterance
  speaker.evs = data.speaker.all %>% group_by(r, utterance) %>% 
    summarize(probs=mean(probs), .groups="keep") %>% arrange(desc(probs)) %>%
    chunk_utterances() %>% summarize(probs=sum(probs), .groups="drop_last")
  
  p.speaker_evs = speaker.evs %>% add_graph() %>% 
    group_by(relation, utterance) %>%
    summarize(probs=mean(probs), .groups = "drop_last") %>% 
    ggplot(aes(x=relation, y=probs, fill=utterance)) + 
    geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
    theme_minimal() + theme(legend.position="top") +
    scale_fill_viridis(name="utterance type", discrete = TRUE) +
    labs(x="causal relation", y=str_wrap("expected utterance choice probability", width=20)) +
    scale_x_discrete(labels = function(label) parse(text=label))
  
  # p.speaker_evs
  ggsave(paste(plot_dir, "speaker-evs.png", sep=fs), p.speaker_evs,
         width=7, height=2.5)

  df.csv = tibble(val = speaker.evs %>% filter(utterance == "conditional" & 
                                                 r=="A || C") %>% pull(probs),
                  key = "sp_ev_any_conditional_indep")
  
  data.speaker = data.speaker.all %>% 
    distinct_at(vars(c(bn_id, utterance)), .keep_all = T)
  # bn_ids of sampled states with nb of occurrence in overall set of samples
  bn_ids = read_rds(str_replace(params.speaker$target, "results", "sample-ids")) %>% 
    group_by_all() %>% summarize(n_sampled=n(), .groups="drop_last")
  
  # note:for some states there might be several best utterance choices for the speaker!
  # (depending on informativeness of utterances)
  data.speaker.best <- data.speaker %>% group_by(rowid) %>%
    mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
    unnest(c(u_best)) %>% select(-p_best) %>% 
    filter(utterance == u_best)
  
  data.speaker.best = left_join(data.speaker.best, bn_ids, by=c("bn_id"))
  
  sp.best.conditions <- data.speaker.best %>%
    mutate(pa=`AC` + `A-C`, pc=`AC` + `-AC`,
           certainA = pa >= theta | pa <= 1-theta,
           certainC = pc >= theta | pc <= 1-theta,
           certain=certainA & certainC,
           uncertain=!certainA & !certainC, 
           speaker_condition = case_when(
             certain ~ "certain",
             uncertain ~ "uncertain",
             T ~ "xor"
           ), sp_condition=speaker_condition) %>%
    select(-certainA, -certainC, -certain, -uncertain) %>%
    chunk_utterances()
  
  sp.best.conditions.chunked <- sp.best.conditions %>%
    add_graph() %>% 
    group_by(speaker_condition, r_graph, utterance) %>%
    dplyr::select(r, r_graph, utterance, speaker_condition, n_sampled) %>%
    summarize(n_sampled=sum(n_sampled), .groups="drop_last") %>%
    mutate(N=sum(n_sampled), p=n_sampled/N) %>%
    mutate(speaker_condition=factor(speaker_condition,
                                    levels=c("certain", "uncertain", "xor")))
  levels(sp.best.conditions.chunked$speaker_condition) <-
    c("certain" = "(i) A certain, C certain",
      "uncertain" = "(ii) A uncertain, C uncertain",
      "xor" = "(iii) A xor C certain")
  
  tikz(paste(plot_dir, "speaker_freq_best_un_certain_other.tex", sep=fs), 
       width = 7.5, height = 2.5, standAlone = FALSE,
       packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                    "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
  p = plot_speaker_conditions(sp.best.conditions.chunked)
  p
  dev.off()
  
  # 1.check A,C uncertain + independent + best utterance is conditional
  # almost true states?
  df.unc_ind_best_if <- sp.best.conditions %>% 
    filter(sp_condition=="uncertain" & utterance == "conditional" & r=="A || C")
  key = "both_uncertain_nb_ind_best_conditional"
  if(nrow(df.unc_ind_best_if) == 0){
    df.csv <- bind_rows(df.csv, tibble(key=key, val = 0))  
  } else {
    check = df.unc_ind_best_if %>%
      mutate(pa=`AC`+`A-C`, pc=`AC`+`-AC`, pna=`-AC`+`-A-C`, pnc=`A-C`+`-A-C`) %>%
      select(bn_id, pa, pc, pna, pnc) %>% 
      pivot_longer(cols=c(-bn_id), names_to="p", values_to="val") %>%
      arrange(desc(val)) %>% 
      distinct_at(vars(c("bn_id")), .keep_all = TRUE)
    
    ratio.almost_true = nrow(check %>% filter(val>=0.85)) / nrow(check)
    print(paste('ratio pragmatic speaker r=ind. + literal almost true (p>=0.85):',
                ratio.almost_true))
    
    df.csv <- bind_rows(df.csv, tibble(key=key, val = nrow(df.unc_ind_best_if)))
    df.csv <- bind_rows(df.csv, tibble(key = "both_uncertain_literal_almost_true",
                                       val = ratio.almost_true))
    check %>% filter(val<0.85)
  } 
  
  # 2. check A certain, C certain and independent, best utterance is not a 
  # conjunction, are these states where no conjunction is applicable in the first
  # place?
  check_certain_ind_best_not_conj = sp.best.conditions %>% add_graph() %>% 
    filter(sp_condition == "certain" & r == "A || C" & 
             !str_detect(u_best, "and")) %>% 
    dplyr::select(rowid, bn_id, AC, `A-C`, `-AC`, `-A-C`) %>% 
    pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), names_to="cell", 
                 values_to="val") %>% 
    filter(val >= theta)
  # should be TRUE (1)
  df.csv = bind_rows(df.csv, tibble(key="both_certain_ind_best_not_conj_since_not_assertable",
                                    val = nrow(check_certain_ind_best_not_conj) == 0))
  
  # Figure 3+4 ----------------------------------------------------------------
  data.cp = data_cp_plots(params) 
  cp.relations = data.cp %>% filter(val=="relations")
  p.probs <- data.cp %>% filter(val=="p") %>% 
    plot_cp_probs(paste(plot_dir, "cp-evs-probs.tex", sep=fs), w=7, h=2.5)
  p.relations <- cp.relations %>% plot_cp_relations()
  # p.probs
  # p.relations
  ggsave(paste(plot_dir, "evs-relations.png", sep=fs), p.relations, width=7, height=2.5)
  
  # Figure 5 ----------------------------------------------------------------
  cond_all = c("p_delta", "p_rooij", "p_diff")
  cond <- "p_rooij"
  
  if(cond == "p_rooij") cond_tex_str <- "$\\Delta^{*}  P$"
  if(cond == "p_delta") cond_tex_str <- "$\\Delta  P$"
  if(cond == "p_diff") cond_tex_str <- "$P(C\\mid A)-P(C)$"
  
  prior <-  read_rds(params.prior$target) %>%
    pivot_wider(names_from = "cell", values_from = "val") %>% 
    mutate(level = "prior") %>%
    pivot_longer(cols=c(p_delta, p_rooij, p_diff),names_to="condition",
                 values_to = "val")
  
  ## a)speaker.literal by filtering samples from prior s.t. those remain where
  ## A > C is assertable (just to check differences)
  # ifac_not_assertable = read_rds(params.speaker$target) %>%
  #   filter(utterance == "A > C" & probs == 0) %>% 
  #   pull(rowid) %>% unique()
  # 
  # speaker.lit = anti_join(read_rds(params.speaker$target),
  #                         tibble(rowid=ifac_not_assertable)) 
  
  # b) speaker literal: samples from prior conditioned on A > C being assertable. 
  # here only total of 18 utterances used, 2 are not applicable due to condition
  # We use this version!
  speaker.lit <- read_rds(params.sp_literal$target)
  speaker.literal <- speaker.lit %>% 
    pivot_longer(cols=c("p_delta", "p_rooij", "p_diff"), names_to = "condition_name",
                 values_to="condition") %>% 
    filter(!condition_name %in% cond_all[cond_all != cond]) %>% 
    dplyr::select(-condition_name) %>% 
    mutate(utterance = paste("utt", utterance, sep="_"),
           level="literal-speaker") %>% 
    group_by(rowid) %>% 
    mutate(p_best=max(probs), u_best = probs == p_best)
  
  # pragmatic speaker condition: literal speaker's best utterance is A > C
  speaker.literal.best = speaker.literal %>%
    filter(u_best) %>% dplyr::select(-p_best, -u_best)
  speaker.pragmatic <- speaker.literal.best %>%
    filter(utterance == "utt_A > C") %>% 
    mutate(level = "pragmatic-speaker")
  speakers = bind_rows(speaker.literal, speaker.pragmatic) %>% 
    group_by(rowid,level)
  
  prior_vals = prior %>% filter(condition == "p_delta") %>% 
    dplyr::select(-condition) %>%  rename(condition=val)
  
  df_accept <- bind_rows(prior_vals, speakers) %>%
    mutate(level=factor(level, levels = c("prior", "literal-speaker",
                                          "pragmatic-speaker"))) %>% 
    group_by(rowid, level)
  plots.all = plot_accept_conditions(df_accept, cond_tex_str)
  p <- plot_grid(plots.all[[3]], plots.all[[2]], plots.all[[1]], ncol=1)
  ggsave(paste(plot_dir, "accept-conditions.png", sep=fs), p, width=8, height=10)
  
  
  # ---%%%%%% additional checks on tables of pragmatic speaker condition %%%%%%---
  prag.ind = speaker.pragmatic %>% filter(r == "A || C") 
  sp.prag.ratio.ind = nrow(prag.ind) / nrow(speaker.pragmatic)
  df.csv <- bind_rows(df.csv, tibble(key = "ratio_pragmatic_speaker_cond_r_ind", 
                                     val = sp.prag.ratio.ind))
  
  # if applicable check states where pragmatic speaker's best utterance is A->C
  # and A,C independent
  if(nrow(prag.ind) > 0){
    prag.ind.probs = prag.ind %>% dplyr::select(-utterance) %>% group_by(bn_id) %>% 
      mutate(pa=`AC`+`A-C`, pc=`AC`+`-AC`, pna=`-AC`+`-A-C`, pnc=`A-C`+`-A-C`) %>%
      select(bn_id, pa, pc, pna, pnc) %>% 
      pivot_longer(cols=c(-bn_id), names_to="p", values_to="val") %>%
      arrange(desc(val)) %>% 
      distinct_at(vars(c("bn_id")), .keep_all = TRUE)
    df.almost_true = prag.ind.probs %>% filter(val >= 0.85) 
    ratio.almost_true = nrow(df.almost_true) / prag.ind.probs %>% nrow()
    df.csv <- bind_rows(df.csv, tibble(key="ratio_pragmatic_speaker_cond_r_ind_lit_almost_true", 
                                       val = ratio.almost_true))
    prior.unc = prior %>%
      mutate(conj=AC>theta | `A-C`>theta | `-AC`>theta | `-A-C`>theta,
             lit=(AC+`A-C` > theta) | (AC+`-AC`>theta) |
                 (AC+`A-C` < 1-theta) | (AC+`-AC`<1-theta)) %>%
      filter(!conj &!lit) %>% ungroup() %>%
      filter(AC/(AC+`A-C`) > theta)
    ids = prior.unc %>% filter(r=="A || C") %>% pull(bn_id)
    ratio = nrow(prag.ind %>% filter(bn_id %in% ids)) / nrow(prag.ind)
    df.csv <- bind_rows(df.csv, tibble(key="ratio_ind_tables_ifac_true_but_no_lit_no_conj", 
                                       val = ratio))
  }
  
  # some more checks
  ratio.sp_prag = speaker.pragmatic %>% filter(condition > theta) %>% nrow() / (speaker.pragmatic %>% nrow())
  ratio.sp_lit  = speaker.literal %>% filter(condition > theta) %>% nrow() / (speaker.literal %>% nrow())
  ratio.prior = prior_vals %>% filter(condition > theta) %>% ungroup() %>% nrow() / (prior_vals  %>% nrow())
  df.csv <- bind_rows(df.csv, tibble(key="ratio_pragmatic_speaker_cond_accept_cond_larger_theta", val = ratio.sp_prag))
  df.csv <- bind_rows(df.csv, tibble(key="ratio_literal_speaker_cond_accept_cond_larger_theta", val = ratio.sp_lit))
  df.csv <- bind_rows(df.csv, tibble(key="ratio_prior_accept_cond_larger_theta", val = ratio.prior))
  
  ratio.prior_neg = prior_vals %>% filter(condition < 0) %>% ungroup() %>% nrow() / (prior_vals  %>% nrow())
  df.csv <- bind_rows(df.csv, tibble(key="ratio_prior_accept_cond_neg", val = ratio.prior_neg))
  
  # Figures Appendix -----------------------------------------------------------
  # speaker results for states from literal speaker condition where the 
  # speaker's best utterance is NOT A->C
  df.sp_lit_best_not_ifac = speaker.literal.best %>%
    filter(utterance != "utt_A > C") %>% 
    mutate(utterance = str_replace(utterance, "utt_", ""), utt=utterance) %>%
    chunk_utterances() %>% add_graph()
  
  tikz(paste(plot_dir, "literal-speaker_deltap_best_not_ac.tex", sep=fs), 
       width = 7, height = 2.5, standAlone = FALSE,
       packages = c("\\usepackage{tikz, amssymb, amsmath}", 
                    "\\newcommand{\\indep}{\\rotatebox[origin=c]{90}{$\\models$}}"))
  
  p <- df.sp_lit_best_not_ifac %>% 
    mutate(utterance = as.character(utterance), 
           utterance = case_when(utterance == "conditional" ~ "(other) conditional",
                                 T ~ utterance),
           utterance = as.factor(utterance)) %>%
    ggplot(aes(x=utterance, y=condition, color=r_graph)) + 
    geom_boxplot() +
    labs(x="", y="$\\Delta^{*}  P^{(s)}$") +
    theme(axis.text.y=element_text(), legend.key.size = unit(0.75,"line"),
          legend.spacing.x = unit(1.25, "line")) +
    scale_color_viridis(name="causal relation $\\mathcal{R}$", discrete=TRUE,
                        labels=LABELS_R_TEX) +
    coord_cartesian(ylim=c(0,1)) # zoom in for positive values
  p
  dev.off()
  
  # get ranges
  df.sp_lit_best_not_ifac %>% add_graph() %>% group_by(utterance, r_graph) %>% 
    summarize(min=min(condition), max=max(condition))
  
  # literal-speaker condition, but best utterance is another conditional than A->C
  df.lit_sp = df.sp_lit_best_not_ifac %>% filter(utterance=="conditional") %>% 
    mutate(conj=AC>theta | `A-C`>theta | `-AC`>theta | `-A-C`>theta,
           lit=(AC+`A-C` > theta) | (AC+`-AC`>theta) |
             (AC+`A-C` < 1-theta) | (AC+`-AC`<1-theta)) 
  
  nb_lit_or_conj_assertable = nrow(df.lit_sp %>% filter(conj | lit))
  # should be 0
  df.csv <- bind_rows(df.csv, tibble(key="nb_tables_literal_speaker_cond_other_conditional_best_utt_but_conj_or_lit_assertable",
                                     val = nb_lit_or_conj_assertable))
  
  df.csv = df.csv %>% add_column(alpha=alpha, theta=theta)
  csv_data = bind_rows(csv_data, df.csv)
}

if(pars == "all") {
  write_csv(csv_data, paste(params$seed_dir, "results-across-params.csv", sep=fs))
} 
csv_data$key %>% unique
csv_data %>% filter(key == "sp_ev_any_conditional_indep") %>% arrange(alpha)
csv_data %>% filter(key == "both_uncertain_nb_ind_best_conditional" & val!=0)
csv_data %>% filter(key == "both_certain_ind_best_not_conj_since_not_assertable" &
                      val != 1)
csv_data %>% filter(key == "ratio_pragmatic_speaker_cond_r_ind" & val != 0)
csv_data %>% 
  filter(key %in% c("ratio_literal_speaker_cond_accept_cond_larger_theta",
                    "ratio_pragmatic_speaker_cond_accept_cond_larger_theta")) %>% 
  arrange(key, theta) %>% dplyr::select(val, key, theta) %>% distinct()

csv_data %>% filter(key == "ratio_prior_accept_cond_larger_theta") %>% 
  arrange(key, theta) %>% dplyr::select(val, key, theta) %>% distinct()
csv_data %>% filter(key == "ratio_prior_accept_cond_neg") %>% 
  dplyr::select(val, key) %>% distinct()

csv_data %>% 
  filter(key == "nb_tables_literal_speaker_cond_other_conditional_best_utt_but_conj_or_lit_assertable") %>% 
  filter(val != 0)



# conditional vs. or ------------------------------------------------------
df.or = tables %>% compute_cond_prob("P(C|A)", "p_if") %>% 
  mutate(p_or = AC + `-AC` + `-A-C`) %>% 
  mutate(or_true = p_or >= theta, if_true = p_if >= theta) %>% 
  dplyr::select(bn_id, p_if, p_or, or_true, if_true) %>% ungroup() 

df.or %>% summarize(n_if=sum(if_true), n_or = sum(or_true))
# when p_if is true, p_or is true as well
df.or %>% filter(if_true) %>% summarize(s=sum(or_true)/n())
# when p_or is true, p_if may be false
df.or %>% filter(or_true) %>% summarize(s=sum(if_true)/n())
# or is 
df.or %>% filter(or_true & !if_true)

