library(kableExtra)
library(ggplot2)
library(stringr)
library(arrow)
library(dplyr)
library(tidyr)

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
            dlts = round(mean(dlts),1),
            size = round(mean(n),1), 
            days = round(mean(days)))

#####################

result_list <- c("_28days", "_56days", "_0.61beta_sd", "_main")

for(result in result_list){
  data <- open_dataset(paste0("output_summary", result, "//"))
  full_data <- open_dataset(paste0("output_full", result, "//"))
  
  if(result == "_main"){result <- ""}
  
  Results_parameter <- data%>%
    group_by(parameter,decision_rule)%>%
    ipde_summarize()%>%
    collect()%>%
    mutate(decision_rule = factor(x = decision_rule, levels = c("crm0", "crm", "alphacrm","pocrm"), labels = c("CRM", "AIDE-CRM","$\\alpha$-CRM","IPDE-POCRM")),
           parameter = substr(parameter, start = 7, stop = 12))%>%
    arrange(parameter,decision_rule)
  results_to_latex(Results_parameter, col_1 = "Alpha", file = paste0("Figures/fig_alpha_table", result,".tex"))
  file.copy(from = paste0("Figures/fig_alpha_table", result,".tex"), to = paste0("../Manuscript 1/Figures/fig_alpha_table", result,".tex"), overwrite = TRUE)
  
  
  Results_scenario <- data%>%
    group_by(scenario,decision_rule)%>%
    ipde_summarize()%>%
    collect()%>%
    mutate(decision_rule = factor(x = decision_rule,  levels = c("crm0", "crm", "alphacrm","pocrm"), labels = c("CRM", "AIDE-CRM","$\\alpha$-CRM","IPDE-POCRM")),
           scenario = substr(scenario, start = 5, stop = 10))%>%
    arrange(scenario,decision_rule)
  results_to_latex(Results_scenario, col_1 = "Scenario", file = paste0("Figures/fig_scenario_table", result,".tex"))
  
  Results_detailed <- data%>%
    group_by(scenario,parameter,decision_rule)%>%
    ipde_summarize()%>%
    collect()
  
  if(FALSE){
    for(scenar in 0:5){
      Results_detailed%>%
        filter(scenario == paste0("mtd_",scenar))%>%
        group_by(parameter,decision_rule)%>%
        summarise_at(.vars = 2:11,
                     .funs = mean)%>%
        mutate(decision_rule = factor(x = decision_rule, levels = c("crm0", "crm", "alphacrm","pocrm"), labels = c("CRM", "AIDE-CRM","$\\alpha$-CRM","IPDE-POCRM")),
               parameter = substr(parameter, start = 7, stop = 12))%>%
        arrange(parameter,decision_rule)%>%
        results_to_latex(col_1 = "Alpha", file = paste0("Figures/fig_scenario_",scenar,"_table", result,".tex"))
      file.copy(from = paste0("Figures/fig_scenario_",scenar,"_table", result,".tex"), to = paste0("../Manuscript 1/Figures/fig_scenario_",scenar,"_table", result,".tex"), overwrite = TRUE)
    }
  }
  
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
  
  ggsave(plot_results_graph, filename = paste0("falke_4_t", result,".tiff"), device = "tiff", dpi = 600, width = 8, height = 4.8, scale = 1, units = "in", path = "Figures/")
  ggsave(plot_results_graph, filename = paste0("falke_4_t", result,".png"), device = "png", dpi = 600, width = 8, height = 4.8, scale = 1, units = "in", path = "Figures/")
}

### Text values

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
