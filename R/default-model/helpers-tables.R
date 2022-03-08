library(dplyr)
library(here)
library(latex2exp)
source(here("R", "helper-functions.R"))

unnest_tables <- function(tables){
  tables <- tables %>% rowid_to_column()
  tables_long <- tables %>% unnest(cols=c(vs, ps)) %>% rename(cell=vs, val=ps)
  return(tables_long)
}

# plot densities of generated tables for different causal nets
# @arg tables: long format
plot_tables_relations <- function(tables_path, plot_dir, w, h, tables=NA){
  if(is.na(tables)){
    tables <- readRDS(tables_path)
  }
  tables.wide <- tables %>% unnest_tables() %>% ungroup() %>%
    select(-table_id) %>% group_by(rowid) %>% 
    pivot_wider(names_from = cell, values_from = val)
  tables.evs = tables.wide %>% group_by(r) %>%
    summarize(AC = round(mean(AC), 2), `A-C` = round(mean(`A-C`), 2), 
              `-AC` = round(mean(`-AC`), 2) , `-A-C` = round(mean(`-A-C`), 2))
  
  tables.long <- tables.wide %>% 
    pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), names_to = "cell",
                 values_to = "val") %>% 
    group_by(bn_id, r) %>% 
    mutate(cell=factor(cell, levels=c("AC", "A-C", "-AC", "-A-C")),
           r=as.factor(r))
  labels_cells <- c(`AC` = TeX("$P(w_{AC})$"), `A-C` = TeX("$P(w_{A})$"),
                    `-AC` = TeX("$P(w_{C})$"),`-A-C` = "P("~w[symbol("\306")]~")") 
  all_plots = list()
  relations <- tables.long$r %>% unique
  df_meta = list()
  for(i in seq(1, length(relations))) {
    r <- relations[i]
    df_meta[[i]] = 
      data.frame(label = tables.evs %>% filter(r == !!r) %>%
                   dplyr::select(-r) %>% as.numeric(),
                 cell=c(1,2,3,4),
                 x=rep(0.8, 4), y=rep(200,4)
                 #lab = c("P(a,c)", "P(a,¬c)", "P(¬a,c)", "P(¬a,¬c)")
                ) %>% 
      mutate(cell = factor(cell))
    if(r == "A || C"){
      df_meta[[i]][["y"]] <- rep(500, 4)
    }
    levels(df_meta[[i]]$cell) <- labels_cells
  }
  names(df_meta) <- relations

  for(i in seq(1, length(relations))) {
    r <- relations[i] %>% as.character()
    tbls <- tables.long %>% filter(r == relations[[i]])
  
    if(nrow(tbls) > 0) {
      levels(tbls$cell) <- labels_cells
      p <- tbls  %>% group_by(bn_id) %>% 
        ggplot(aes(x=val)) +
        geom_histogram(fill="gray") + 
        facet_wrap(~cell, ncol = 2, #scales = "free",
                   # labeller = labeller(cell = c(`AC`="P(a,c)", `A-C`="P(a,¬c)",
                   #                             `-AC`="P(¬a,c)", `-A-C`="P(¬a,¬c)"))
                   labeller = label_parsed) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        labs(x="probability", y="density") +
        theme_minimal() +
        theme(legend.position = "none") +
        ggtitle(LABELS_r[[`r`]]) +
        geom_text(data = df_meta[[`r`]], mapping = aes(x=x, y=y, label=label))
      all_plots[[i]] = p
    
      save_to = paste(plot_dir, paste("tables_", str_replace_all(r, " ", ""),
                                      ".png", sep=""),
                      sep=.Platform$file.sep)
      ggsave(save_to, p, width=w, height=h)
      print(paste('saved to', save_to))
    }
  }
  return(all_plots)
}
