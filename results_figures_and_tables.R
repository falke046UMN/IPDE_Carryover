library(kableExtra)
library(ggplot2)
library(stringr)
library(arrow)
library(dplyr)
library(tidyr)

result_list <- c("28days", "56days", "0.61beta_sd", "indif_0.025", "indif_0.1", "fixed_n", "prop_prob", "prior_0.3", "prior_0.6", "prior_0.9", "main")

## Functions #####################

results_to_latex <- function(input, midrule_space = 4, col_1, file){
  input[,1] <- ifelse(1:nrow(input) %% midrule_space == 1, 
                      paste0("\\multirow{", midrule_space,"}{2em}{",input[[1]],"}"),
                      " ")
  body <- kable(input, format = "latex", booktabs = TRUE, escape = FALSE,linesep = c(rep("",midrule_space-1), "\\midrule"))
  body <- str_replace(body, pattern = "^\n.*\n.*\n.*\n.*", replacement = "")
  body <- str_replace(body, pattern = "\n.*$", replacement = "\n\\\\end{tabular*}")
  
  header <- paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}clllllllllll@{\\extracolsep\\fill}}
\\toprule
&  & \\multicolumn{3}{c}{Final Selection (\\%)$^1$} & \\multicolumn{3}{c}{Allocations (\\%)} & &  &  &  \\\\
\\multirow{2}{2em}{",col_1,"} & \\multirow{2}{2em}{Method} & \\multirow{2}{2em}{MTD} & \\multirow{2}{2em}{Above MTD} & \\multirow{2}{2em}{Below MTD} & \\multirow{2}{2em}{MTD} & \\multirow{2}{2em}{Above MTD} & \\multirow{2}{2em}{Below MTD} & \\multirow{2}{2em}{Toxic~\\%$^2$} & \\multirow{2}{2em}{DLTs} & \\multirow{2}{2em}{Trial Size} & \\multirow{2}{2em}{Days}  \\\\ \\\\ 
\\midrule")
  out <- paste0(header, body)
  
  write(out, file = file, append = FALSE)
}

ipde_summarize <- .%>%
  summarize(pcs = round(mean(mtd_estimate == mtd)*100),
            pcs_o = round(mean(mtd_estimate > mtd)*100),
            pcs_u = round(mean(mtd_estimate < mtd)*100),
            pca = round(mean(assignments_to_mtd/i)*100),
            pca_o = round(mean(assignments_over_mtd/i)*100),
            pca_u = round(mean(assignments_below_mtd/i)*100),
            tox_pct = round(mean(toxic_assignments/i)*100),
            dlts = round(mean(dlts), 1),
            size = round(mean(n), 1), 
            days = round(mean(days)))

#####################

for(result in result_list){
  
  if(result %in% c("prior_0.3", "prior_0.6", "prior_0.9")){
    file_list_summmary <- grep(list.files("output/Data/output_main_summary//", recursive = TRUE, full.names = TRUE), invert = TRUE, pattern = "decision_rule=alphacrm", value = TRUE)
    file_list_summmary <- c(file_list_summmary, grep(list.files(paste0("output/Data/output_", result, "_summary//"), recursive = TRUE, full.names = TRUE), pattern = "decision_rule=alphacrm", value = TRUE))
    data <- open_dataset(file_list_summmary, format = "parquet")
    
    file_list_full <- grep(list.files("output/Data/output_main_full//", recursive = TRUE, full.names = TRUE), invert = TRUE, pattern = "decision_rule=alphacrm", value = TRUE)
    file_list_full <- c(file_list_full, list.files(paste0("output/Data/output_", result, "_full//"), recursive = TRUE, full.names = TRUE))
    full_data <- open_dataset(file_list_full, format = "parquet")
  } else {
    data <- open_dataset(paste0("output/Data/output_", result, "_summary//"))
    full_data <- open_dataset(paste0("output/Data/output_", result, "_full//"))
  }
  
  Results_parameter <- data%>%
    group_by(parameter,decision_rule)%>%
    ipde_summarize()%>%
    collect()%>%
    mutate(decision_rule = factor(x = decision_rule, levels = c("crm0", "crm", "alphacrm","pocrm"), labels = c("CRM", "AIDE-CRM","$\\alpha$-CRM","IPDE-POCRM")),
           parameter = substr(parameter, start = 7, stop = 12))%>%
    arrange(parameter,decision_rule)
  results_to_latex(Results_parameter, col_1 = "Alpha", file = paste0("output/Tables/alpha_table_", result,".tex"))
  
  Results_scenario <- data%>%
    group_by(scenario,decision_rule)%>%
    ipde_summarize()%>%
    collect()%>%
    mutate(decision_rule = factor(x = decision_rule,  levels = c("crm0", "crm", "alphacrm","pocrm"), labels = c("CRM", "AIDE-CRM","$\\alpha$-CRM","IPDE-POCRM")),
           scenario = substr(scenario, start = 5, stop = 10))%>%
    arrange(scenario,decision_rule)
  results_to_latex(Results_scenario, col_1 = "Scenario", file = paste0("output/Tables/scenario_table_", result,".tex"))
  
  #### Plot of results
  
  plot_results <- data%>%
    group_by(scenario,parameter,decision_rule)%>%
    ipde_summarize()%>%
    select(scenario,parameter,decision_rule, days, pcs, dlts)%>%
    collect()%>%
    tidyr::pivot_longer(cols = c("days", "pcs", "dlts"), names_to = "outcome")%>%
    mutate(outcome = factor(outcome, 
                             levels = c("pcs", "dlts", "days"),
                             labels = c("Correct Selection Percent", "DLTs per Trial", "Days per Trial")),
           decision_rule = factor(decision_rule, 
                                   levels = c("crm0", "crm", "alphacrm", "pocrm"), 
                                   labels = c("CRM", "AIDE-CRM", "\u03b1-CRM", "IPDE-POCRM")),
           scenario = paste("Scenario", substr(scenario, start = 5, stop = 10)),
           parameter = substr(parameter, start = 7, stop = 12))
  
  plot_results_graph <- ggplot(plot_results, aes(x = parameter, y = value, color = decision_rule, group = decision_rule))+
    geom_line()+
    facet_grid(rows = vars(outcome), cols = vars(scenario), scales = "free_y")+
    labs(x = "Cumulative Toxicity Strength (\u03b1)", y = "", color = "Decison Method")+
    theme_bw()+
    theme(strip.text.y = element_text(size = 7))+
    ylim(c(0,NA))
  
  ggsave(plot_results_graph, filename = paste0("results_figure_", result,".tiff"), device = "tiff", dpi = 600, width = 8, height = 4.8, scale = 1, units = "in", path = "output/Figures/")
  ggsave(plot_results_graph, filename = paste0("results_figure_", result,".png"), device = "png", dpi = 600, width = 8, height = 4.8, scale = 1, units = "in", path = "output/Figures/")
}

file.copy(from = "output/Figures/results_figure_main.tiff", to = paste0("output/Figures/falke_4_t.tiff"), overwrite = TRUE)
file.copy(from = "output/Figures/results_figure_main.png", to = paste0("output/Figures/falke_4_t.png"), overwrite = TRUE)

### Text values
data <- open_dataset(paste0("output/Data/output_main_summary//"))
  
Results_detailed <- data%>%
  group_by(scenario,parameter,decision_rule)%>%
  ipde_summarize()%>%
  collect()

# Efficiency

data%>%
  filter(decision_rule == "crm0", scenario == "mtd_5")%>%
  summarise(mean(mtd_estimate == 0)*100)%>%
  collect()

Results_detailed%>%
  pivot_wider(id_cols = c("scenario","parameter"), 
              names_from = "decision_rule", 
              values_from = "days")%>%
  arrange(scenario,parameter)%>%
  transmute(aide_dif = 100-crm/crm0*100,
            alpha_dif = 100-alphacrm/crm0*100,
            pocrm_dif = 100-pocrm/crm0*100)%>%
  print(n = 24)


data%>%
  group_by(decision_rule)%>%
  summarise(mean_days = mean(days),
            mean_n = mean(n))%>%collect()

data%>%
  filter(decision_rule  != "crm0", scenario %in% c("mtd_1","mtd_5"))%>%
  group_by(decision_rule,scenario)%>%
  summarise(mean_days = mean(days),
            mean_n = mean(n))%>%collect()

# Accuracy
Results_detailed%>%
  pivot_wider(id_cols = c("scenario","parameter"), 
              names_from = "decision_rule", 
              values_from = "pcs")%>%
  arrange(scenario,parameter)%>%
  transmute(aide_dif = crm - crm0,
            alpha_dif = alphacrm - crm0,
            pocrm_dif = pocrm - crm0)%>%
  print(n = 24)

Results_detailed%>%
  group_by(decision_rule)%>%
  summarise(overall_pcs = mean(pcs))

# Safety
Results_detailed%>%
  group_by(parameter,decision_rule)%>%
  summarise_at(c("dlts","pca"),mean)

Results_detailed%>%
  filter(parameter == "alpha_0.9")%>%
  group_by(decision_rule)%>%
  summarise_at(.vars = "dlts", .funs = mean)%>%
  pivot_wider(names_from = "decision_rule", 
              values_from = "dlts")%>%
  transmute(aide_excess_dlts = crm - crm0,
            alpha_excess_dlts = alphacrm - crm0,
            pocrm_excess_dlts = pocrm - crm0)

# Sensitivity Analyses

## 28 Days
out <- open_dataset(paste0("output/Data/output_28days_summary//"))%>%
  group_by(decision_rule)%>%
  summarise(mean_days = mean(days),
            mean_n = mean(n))%>%
  collect()
(1-mean(out$mean_days[1:3])/out$mean_days[4])*100
(1-mean(out$mean_n[1:3])/out$mean_n[4])*100

## 56 Days
out <- open_dataset(paste0("output/Data/output_56days_summary//"))%>%
  group_by(decision_rule)%>%
  summarise(mean_days = mean(days),
            mean_n = mean(n))%>%collect()
(1-mean(out$mean_days[1:3])/out$mean_days[4])*100
(1-mean(out$mean_n[1:3])/out$mean_n[4])*100

## Indif interval 0.1
out <- cbind(open_dataset(paste0("output/Data/output_indif_0.1_summary////"))%>%
  group_by(decision_rule)%>%
  summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect(),
pcs_main = open_dataset(paste0("output/Data/output_main_summary/////"))%>%
  group_by(decision_rule)%>%
  summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%pull(pcs))
out$dif <- out$pcs_main - out$pcs

## Indif interval 0.025
out <- cbind(open_dataset(paste0("output/Data/output_indif_0.025_summary//"))%>%
               group_by(decision_rule)%>%
               summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect(),
             pcs_main = open_dataset(paste0("output/Data/output_main_summary/////"))%>%
               group_by(decision_rule)%>%
               summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%pull(pcs))
out$dif <- out$pcs_main - out$pcs

## Prop prob
out <- cbind(open_dataset(paste0("output/Data/output_prop_prob_summary//"))%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect(),
             open_dataset(paste0("output/Data/output_main_summary/////"))%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_main = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%ungroup()%>%select(pcs_main))
out$dif <- out$pcs_main - out$pcs

# Fixed n
out <- open_dataset(paste0("output/Data/output_fixed_n_summary/"))%>%
  group_by(decision_rule)%>%
  summarise(mean_dosings = mean(i))%>%collect()
mean(out$mean_dosings[1:3] - out$mean_dosings[4])

out <- cbind(open_dataset(paste0("output/Data/output_fixed_n_summary/"))%>%
               filter(decision_rule  != "crm0")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs = round(mean(mtd_estimate == mtd)*100),
                         dlts = mean(dlts))%>%collect(),
             open_dataset(paste0("output/Data/output_main_summary/////"))%>%
               filter(decision_rule  != "crm0")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_main = round(mean(mtd_estimate == mtd)*100),
                         dlts_main = mean(dlts))%>%collect()%>%ungroup()%>%select(pcs_main, dlts_main))
mean(out$pcs_main - out$pcs)
mean(out$dlts_main - out$dlts)

# Priors
out <- cbind(open_dataset(paste0("output/Data/output_prior_0.3_summary//"))%>%
               filter(decision_rule  == "alphacrm")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_0.3 = round(mean(mtd_estimate == mtd)*100))%>%collect(),
             open_dataset(paste0("output/Data/output_prior_0.6_summary//"))%>%
               filter(decision_rule  == "alphacrm")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_0.6 = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%ungroup()%>%select(pcs_0.6),
             open_dataset(paste0("output/Data/output_prior_0.9_summary//"))%>%
               filter(decision_rule  == "alphacrm")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_0.9 = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%ungroup()%>%select(pcs_0.9),
             open_dataset(paste0("output/Data/output_main_summary//"))%>%
               filter(decision_rule  == "alphacrm")%>%
               group_by(decision_rule, parameter)%>%
               summarise(pcs_main = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%ungroup()%>%select(pcs_main))
mean(out$pcs_main - out$pcs)

# Patient heterogeneity
out <- cbind(open_dataset(paste0("output/Data/output_0.61beta_sd_summary//"))%>%
               group_by(decision_rule)%>%
               summarise(pcs = round(mean(mtd_estimate == mtd)*100))%>%collect(),
             open_dataset(paste0("output/Data/output_main_summary//"))%>%
               group_by(decision_rule)%>%
               summarise(pcs_main = round(mean(mtd_estimate == mtd)*100))%>%collect()%>%ungroup()%>%select(pcs_main))
out$dif <- out$pcs_main - out$pcs
