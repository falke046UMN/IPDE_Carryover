task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source(file = "IPDEsim.R")

options(mc.cores = 24)

logistic <- function(x, x0, k, phi = 0.3){signif(1/(1+exp(-(log(phi/(1-phi)) + k*(x-x0)))))}

scenarios <- list(
  mtd_0 = logistic(1:5, x0 = 0, k = 1/2),
  mtd_1 = logistic(1:5, x0 = 1, k = 1/2),
  mtd_2 = logistic(1:5, x0 = 2, k = 1/2),
  mtd_3 = logistic(1:5, x0 = 3, k = 1/2),
  mtd_4 = logistic(1:5, x0 = 4, k = 1/2),
  mtd_5 = logistic(1:5, x0 = 5, k = 1/2))

simulate(iters = 1000,
         max_n = 30,
         task_id = task_id,
         decision_rules = c("pocrm", "alphacrm", "crm", "crm0"),
         true_dlts = scenarios,
         alphas = c(0, 0.3, 0.6, 0.9),
         return_value = FALSE)
