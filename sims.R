task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job_name <- as.character(Sys.getenv("SLURM_JOB_NAME"))

source(file = "IPDEsim.R")

options(mc.cores = 24)

scenarios <- list(
  mtd_0 = logistic(1:5, x0 = 0, k = 1/2),
  mtd_1 = logistic(1:5, x0 = 1, k = 1/2),
  mtd_2 = logistic(1:5, x0 = 2, k = 1/2),
  mtd_3 = logistic(1:5, x0 = 3, k = 1/2),
  mtd_4 = logistic(1:5, x0 = 4, k = 1/2),
  mtd_5 = logistic(1:5, x0 = 5, k = 1/2))

if(job_name == "IPDEsim_main"){
simulate(folder_name = "output_main",
         iters = 4000,
         max_n = 30,
         task_id = task_id,
         decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
         true_dlts = scenarios,
         alphas = c(0, 0.3, 0.6, 0.9))
}

if(job_name == "IPDEsim_28days"){
  simulate(folder_name = "output_28days",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           accrual_mean = 28)
}

if(job_name == "IPDEsim_56days"){
  simulate(folder_name = "output_56days",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           accrual_mean = 56)
}

if(job_name == "IPDEsim_0.61beta_sd"){
  simulate(folder_name = "output_0.61beta_sd",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           response_sim_args = list(beta_sd = 0.61))
}

if(job_name == "IPDEsim_indif_0.1"){
  simulate(folder_name = "output_indif_0.1",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           decision_rule_args = list(indif_interval = 0.1))
}

if(job_name == "IPDEsim_indif_0.025"){
  simulate(folder_name = "output_indif_0.025",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           decision_rule_args = list(indif_interval = 0.025))
}

if(job_name == "IPDEsim_fixed_n"){
  simulate(folder_name = "output_fixed_n",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           count_ipde_dosings = FALSE)
}

if(job_name == "IPDEsim_prop_prob"){
  simulate(folder_name = "output_prop_prob",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           response_sim = prop_prob_sim)
}

if(job_name == "IPDEsim_prior_0.3"){
  simulate(folder_name = "output_prior_0.3",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("alphacrm"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9), 
           decision_rule_args = list(alpha_a = 0.6, alpha_b = 1.4))
}

if(job_name == "IPDEsim_prior_0.6"){
  simulate(folder_name = "output_prior_0.6",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("alphacrm"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9), 
           decision_rule_args = list(alpha_a = 1.2, alpha_b = 0.8))
}

if(job_name == "IPDEsim_prior_0.9"){
  simulate(folder_name = "output_prior_0.9",
           iters = 4000,
           max_n = 30,
           task_id = task_id,
           decision_rules = c("alphacrm"),
           true_dlts = scenarios,
           alphas = c(0, 0.3, 0.6, 0.9),
           decision_rule_args = list(alpha_a = 1.8, alpha_b = 0.2))
}
