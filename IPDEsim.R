
suppressPackageStartupMessages({
suppressWarnings({
library(magrittr)
library(ggplot2)
library(parallel)
library(rstan)
library(dplyr)  
library(arrow)  
library(data.table)
})})

##### Functions #####

#### Main functions

new_ipdeTrial <- function(max_n = 30, max_d = 5, t_max = 28, accrual_dist = poisson_accrual, accrual_mean = 14, tox_dist = unif_tox, ipde = TRUE){
  t <- accrual_dist(accrual_mean = accrual_mean, max_n = max_n)
  if(ipde){t <- c(t, t + t_max)}
  order_t <- order(t)
  N <- length(t)
  out <- list(pid = rep(sprintf(fmt = "P%03d", 1:max_n), ifelse(ipde, 2, 1))[order_t],
              ipde = rep(c(0,1), each = max_n)[order_t],
              t = t[order_t],
              t_eval = rep(NA, N),
              y = rep(NA, N),
              d1 = rep(NA, N),
              dlt_rate = rep(NA, N),
              beta_shift = rep(rnorm(max_n), 2)[order_t],
              min_d = rep(1, N),
              max_d = rep(max_d,N),
              n = 0,
              i = 0,
              mtd_estimate = NA,
              parameter_estimates = list()
              )
  class(out) <- "ipdeTrial"
  return(out)
}

length.ipdeTrial <- function(trial){length(trial$d1)}

as.data.frame.ipdeTrial <- function(trial){
  data.frame(pid = trial$pid,
             ipde = trial$ipde,
             d1 = trial$d1,
             y = trial$y,
             t = trial$t,
             t_eval = trial$t_eval,
             dlt_rate = trial$dlt_rate,
             beta_shift = trial$beta_shift,
             min_d = trial$min_d,
             max_d = trial$max_d)
}

print.ipdeTrial <- function(trial, remove_na = FALSE){
#  print(noquote(paste0("n = ",trial$n)))
#  print(noquote(paste0("i = ",trial$i)))
  temp <- as.data.frame(trial)
  if(remove_na){temp <- temp[!is.na(temp$d1),]}
  temp[is.na(temp)] <- "-"
  print(temp[,!(colnames(temp) %in% c("beta_shift", "min_d", "max_d"))])
  print(noquote(paste0("MTD estimate = ",trial$mtd_estimate)))
  print(noquote(paste0("True MTD = ",trial$mtd)))
}

  
evalble <- function(trial, i = NULL, done = FALSE){
  out <- as.data.frame(trial)
  if(done){return(out[!is.na(trial$t_eval),])}
  out <- out[!is.na(trial$t_eval) & (out$t_eval <= trial$t[i]),]
  return(out)
}

last_assign <- function(trial){
  max(which(!is.na(trial$d1)))
}

rewind <- function(trial, i, assigned = FALSE){
  N <- length(trial$y)
  ia <- i
  if(assigned){ia <- i+1}
  
  trial$n <- length(na.omit(unique(trial$pid[1:(ia-1)])))
  trial$d1[ia:N] <- NA
  trial$y[i:N] <- NA
  trial$t_eval[i:N] <- NA
  trial$dlt_rate[i:N] <- NA
  trial$min_d[i:N] <- trial$min_d[1]
  trial$max_d[i:N] <- trial$max_d[1]
  return(trial)
}

sim_trials <- function(max_n = 24,
                       t_max = 28,
                       phi = 0.3,
                       d1_mg = c(15, 20, 30, 35, 45),
                       ipde = TRUE,
                       count_ipde_dosings = TRUE,
                       tox_dist = unif_tox,
                       accrual_dist = poisson_accrual,
                       accrual_mean = 14,
                       response_sim = new_sim,
                       response_sim_args = list(),
                       decision_rule,
                       decision_rule_args = list(),
                       stopadjust_rule = no_stopadjust,
                       final_decision_rule,
					             crm_decision_model,
					             alphacrm_decision_model){
  # iterations = number of trials to simulate
  # max_n = Maximum possible number of participants
  # t_max = Last day for DLT evaluation
  # G = Within-cycle toxicity function. should sample number in (0,1) and will be scaled to (0,t_max) 
  # d1_mg = Dose levels in mg
  # phi = target DLT rate
  # patient_max_cohorts = max number of cohorts a patient may be assigned to i.e. number of times each patient may be escalated
  # response_sim = name of function simulating patient responses
  # response_sim_args = list of arguments for response_sim
  # decision_rule = name of function giving next dose level(s) given trial data
  # decision_rule_args = list of arguments for decision rule
  # stopadjust_rule = name of function containing rules for stopping for safety or eliminating doses for safety
  
  max_d <-  length(d1_mg)
  results <- list()
  trial <- new_ipdeTrial(max_n, max_d, t_max, accrual_dist, accrual_mean, tox_dist, ipde = ipde)
  for(i in 1:length(trial)){
    if(count_ipde_dosings & (sum(!is.na(trial$t_eval)) >= max_n)){break}
    
    # Patients with a DLT are not reassigned
    if(1 %in% trial$y[trial$pid == trial$pid[i]]){next}
    
    # Pass the trial data and parameters to decision rule and get d1, the next dose level of agent 1
    if(nrow(evalble(trial, i)) < 3){
      d1 <- 1
    } else if(nrow(evalble(trial, i)) == nrow(evalble(trial, last_assign(trial)))){
      d1 <- trial$d1[last_assign(trial)]
      trial$parameter_estimates[[i]] <- trial$parameter_estimates[[last_assign(trial)]]
    } else{
      decision <- do.call(what = decision_rule, args = c(decision_rule_args, 
                                                    list(trial = trial, 
                                                         i = i, 
                                                         phi = phi, 
                                                         d1_mg = d1_mg, 
                                                         crm_decision_model = crm_decision_model, 
                                                         alphacrm_decision_model = alphacrm_decision_model)))
      trial$parameter_estimates[[i]] <- decision$posteriors
      d1 <- as.integer(truncat(trial$d1[last_assign(trial)] + sign(decision$d1 - trial$d1[last_assign(trial)]), trial$min_d[i], trial$max_d[i]))
      # Must evaluate three participants at dose level before escalating
      prev_d1 <- trial$d1[last_assign(trial)]
      if(sum(evalble(trial, i)$d1 %in% prev_d1) < 3){d1 <- pmin(d1, prev_d1)}
    }
    if(trial$ipde[i] & (d1 <= trial$d1[trial$pid == trial$pid[i]][1])){next}
    trial$d1[i] <- d1
    trial$n <- length(unique(trial$pid[!is.na(trial$y)])) + 1
    
    # Simulate the response of patient
    trial$dlt_rate[i] <- do.call(what = response_sim, args = c(response_sim_args, 
                                                          list(trial = trial, 
                                                               i = i, 
                                                               d1_mg = d1_mg)))
    trial$y[i] <- rbinom(n = 1, size = 1, prob = trial$dlt_rate[i])
    trial$t_eval[i] <- ifelse(trial$y[i] == 1, 
                              trial$t[i] + tox_dist(t_max = t_max, max_n = max_n, phi = phi), 
                              trial$t[i] + t_max)
    
    # Determine if trial is terminated for safety
    if(nrow(evalble(trial, i)) > 3){
      if(trial$parameter_estimates[[i]]$tox_probs[1] > 0.95){
        trial$mtd_estimate <- 0
        trial$mtd <- max(c(0,which(response_sim_args$true_dlt <= phi+0.0001)))
        return(trial)
      }
    }
  }
  decision <- do.call(what = final_decision_rule, args = c(decision_rule_args, 
                                                      list(trial = trial, 
                                                           i = i, 
                                                           phi = phi, 
                                                           d1_mg = d1_mg, 
                                                           crm_decision_model = crm_decision_model, 
                                                           alphacrm_decision_model = alphacrm_decision_model)))
  trial$mtd_estimate <- decision$d1
  trial$mtd <- max(c(0,which(response_sim_args$true_dlt <= phi+0.0001)))
  return(trial)
}


#### Response simulation functions

new_sim <- function(trial, i, true_dlt =  c(0.10, 0.19, 0.30, 0.42, 0.54), d1_mg, cohort = NULL, patient_pool = NULL, alpha = 0, beta_sd = 0){
  d1 <- trial$d1[i]
  d1_mg_true <- d1_mg[d1]
  if(trial$ipde[i] == 1){
    prev_d1 <- trial$d1[which(trial$pid == trial$pid[i])[1]]
    d1_mg_true <- d1_mg_true + alpha*d1_mg[prev_d1]
  }
  
  dlt_rate <- approx(x = c(d1_mg, tail(d1_mg, n = 1)+10*diff(tail(d1_mg, n = 2))), 
                     y = c(true_dlt, tail(true_dlt, n = 1)+10*diff(tail(true_dlt, n = 2))), 
                     rule = 2, 
                     xout = d1_mg_true)$y
  dlt_rate <- dlt_rate^exp(trial$beta_shift[i]*beta_sd)
  dlt_rate <- truncat(dlt_rate,0,1)
  
  return(dlt_rate)
}

prop_prob_sim <- function(trial, i, true_dlt = c(0.10, 0.19, 0.30, 0.42, 0.54), d1_mg, cohort = NULL, patient_pool = NULL, alpha = 0, beta_sd = 0){
  d1 <- trial$d1[i]
  dlt_rate <- true_dlt[d1]
  if(trial$ipde[i] == 1){
    prev_d1 <- trial$d1[which(trial$pid == trial$pid[i])[1]]
    dlt_rate <- dlt_rate + alpha*true_dlt[prev_d1]
  }
  
  dlt_rate <- dlt_rate^exp(trial$beta_shift[i]*beta_sd)
  dlt_rate <- truncat(dlt_rate,0,1)
  
  return(dlt_rate)
}

#### Patient accrual functions

fixed_accrual <- function(accrual_mean, max_n){
  return((0:(max_n-1))*accrual_mean)
}

poisson_accrual <- function(accrual_mean, max_n){
  out <- rexp(n = max_n-1, rate = 1/accrual_mean)
  out <- floor(cumsum(c(0,out)))
  return(out)
}

#### Time to toxicity functions

unif_tox <- function(t_max, max_n, ...){
  out <- runif(n = max_n, min = 1, max = (t_max+1))
  out <- floor(out)
  return(out)
}

weibull_tox <- function(t_max, max_n, phi){
  scale <- t_max/qweibull(phi, shape = 4, scale = 1)
  probs <- diff(c(0, pweibull(q = 2:(t_max+1), shape = 4, scale = scale)))
  out <- sample(1:t_max, size = max_n, prob = probs, replace = TRUE)
  return(out)
}

#### Helper functions

logistic <- function(x, x0, k, phi = 0.3){signif(1/(1+exp(-(log(phi/(1-phi)) + k*(x-x0)))))}

truncat <- function(x, a, b){
  # truncate x to the interval [a,b]
  pmax(pmin(x, b),a)
}

lad <- function(y, x = 1:length(y), target){
  ads <- abs(y - target)
  lads <- sapply(ads,all.equal,min(ads))=="TRUE"
  if(sum(lads)==1){
    out <- x[lads]
  } else if(any(y[lads] <= target)){ 
    out <- max(x[lads & (y <= target)]) 
  } else {
    out <- min(x[lads])
  }
  
  return(as.integer(out))
}


lee_cheung_skeleton <- function(indif_interval, phi, D){
  prior_mtd <- ceiling((D+1)/2)
  D <- as.integer(D)
  out <- numeric(D)
  
  out[prior_mtd] <- phi
  
  if(prior_mtd < D){
    for(i in (prior_mtd+1):D){
      out[i] <- exp(log(phi + indif_interval)/log(phi - indif_interval)*log(out[i-1]))
    }
  }
  
  if(prior_mtd > 1L){
    for(i in (prior_mtd-1):1){
      out[i] <- exp(log(phi - indif_interval)/log(phi + indif_interval)*log(out[i+1]))
    }
  }
  return(out)
}


#### Decision Rules

crm_decision <- function(trial, i, phi = 0.3, indif_interval = 0.05, final_decision = FALSE, mcmc_iters = 5000, crm_decision_model, ...){
  D <- trial$max_d[1]
  current_d1 <- trial$d1[last_assign(trial)]
  trial <- evalble(trial, i, done = final_decision)

  d <- lee_cheung_skeleton(indif_interval, phi, D)
  
  trial$d1 <- factor(trial$d1, levels = 1:D)
  n <- tapply(trial$y, trial$d1, length, default = 0)
  y <- tapply(trial$y, trial$d1, sum, default = 0)
  
  stan_data <- list(y = y, n = n, d = d, D = D)
  stan_out <- sampling(object = crm_decision_model,
                       data = stan_data, 
                       chains = 1, 
                       iter = mcmc_iters,
                       warmup = 1500,
                       verbose = FALSE, 
                       cores = 1,
                       refresh = 0,
                       init = "random",
                       pars = "theta")
  tox_probs <- sapply(d, \(x) mean(x^exp(extract(stan_out)$theta) > phi))
  p_est <- sapply(d, \(x) mean(x^exp(extract(stan_out)$theta)))
  theta_posterior <- mean(extract(stan_out)$theta)
  rm(stan_out)
  
  mtd_est <- which.min((p_est - phi)^2)
  
  return(list(d1 = mtd_est, posteriors = list(ps = p_est, theta = theta_posterior, tox_probs = tox_probs)))
}


pocrm_decision <- function(trial, i, phi = 0.3, indif_interval = 0.05, d1_mg = NULL, final_decision = FALSE, mcmc_iters = 5000, crm_decision_model, ...){
  D <- trial$max_d[1]
  current_d1 <- trial$d1[last_assign(trial)]
  trial <- evalble(trial, i, done = final_decision)
  N <- nrow(trial)
  
  d1_ipde <- factor(paste(trial$d1, trial$ipde), levels = c(paste(1:D,0), paste(2:D,1)))
  n <- tapply(trial$y, d1_ipde, length, default = 0)
  y <- tapply(trial$y, d1_ipde, sum, default = 0)
  
  d <- lee_cheung_skeleton(indif_interval/2, phi, 2*D-1)
  d_noipde <- lee_cheung_skeleton(indif_interval, phi, D)
  
  diffs <- apply(expand.grid(log(d1_mg[-1]),log(d1_mg[-1])), 1, diff)
  diffs <- sort(diffs[!duplicated(round(diffs, 10)) & round(diffs, 10) >= 0 & diffs < log(2)])
  M <- length(diffs)+1
  
  d_mat <- matrix(nrow = M, ncol = 2*D-1)
  d_mat[1,] <- c(d_noipde, d_noipde[-1])
  for(m in 1:(M-1)){
    m_rank <- rank(c(log(d1_mg), log(d1_mg)[-1] + diffs[m] + 0.0001))
    d_mat[m+1,] <-  d[m_rank]
  }
  
  tox_probs <- matrix(nrow = D, ncol = M)
  p_est <- matrix(nrow = D, ncol = M)
  theta_posteriors <- c()
  stan_data <- list(y = y, n = n, D = 2*D-1)
  for(m in 1:M){
    stan_data$d <- d_mat[m,]
    stan_out <- sampling(object = crm_decision_model, 
                         data = stan_data, 
                         chains = 1, 
                         iter = mcmc_iters, 
                         verbose = FALSE, 
                         cores = 1,
                         refresh = 0,
                         init = "random",
                         pars = "theta")
    tox_probs[,m] <- sapply(d_mat[m,1:D], \(x) mean(x^exp(extract(stan_out)$theta) > phi))
    p_est[,m] <- sapply(d_mat[m,1:D], \(x) mean(x^exp(extract(stan_out)$theta)))
    theta_posteriors[m] <- mean(extract(stan_out)$theta)
    rm(stan_out)
  }
  model_probs <- numeric(M)
  theta_priors <- rnorm(2000, mean = 0, sd = 2)
  for(m in 1:M){
    model_probs[m] <- mean(sapply(theta_priors, \(x) exp(sum(trial$y*exp(x)*log(d_mat[m,as.numeric(d1_ipde)])+(1-trial$y)*log(1-d_mat[m,as.numeric(d1_ipde)]^exp(x))))))
  }
  model_probs <- model_probs/sum(model_probs)
  tox_probs <- c(tox_probs %*% model_probs)
  p_est <- c(p_est %*% model_probs)
  
  mtd_est <- which.min((p_est - phi)^2)
  
  return(list(d1 = mtd_est, posteriors = list(thetas = theta_posteriors, model_probs = model_probs, tox_probs = tox_probs)))
}

alphacrm_decision <- function(trial, i, phi = 0.3, indif_interval = 0.05, prior_mtd = NULL, d1_mg = NULL, final_decision = FALSE, final_mtd = FALSE, patient_max_cohorts = NULL, mcmc_iters = 5000, alphacrm_decision_model, alpha_a = 1, alpha_b = 1, ...){
  D <- trial$max_d[1]
  current_d1 <- trial$d1[last_assign(trial)]
  trial <- evalble(trial, i, done = final_decision)
  
  d <- lee_cheung_skeleton(indif_interval, phi, D)
  
  trial$d1 <- factor(trial$d1, levels = 1:D)
  n <- tapply(trial$y[trial$ipde == 0], trial$d1[trial$ipde == 0], length, default = 0)
  y <- tapply(trial$y[trial$ipde == 0], trial$d1[trial$ipde == 0], sum, default = 0)
  n_ipde <- tapply(trial$y, trial$d1, length, default = 0) - n
  y_ipde <- tapply(trial$y, trial$d1, sum, default = 0) - y

  stan_data <- list(D = D, y = y, n = n,  y_ipde = y_ipde, n_ipde = n_ipde, d = d, d1_mg = d1_mg, alpha_a = alpha_a, alpha_b = alpha_b)
  stan_out <- sampling(object = alphacrm_decision_model, 
                       data = stan_data, 
                       chains = 1, 
                       iter = mcmc_iters, 
                       verbose = FALSE, 
                       cores = 1,
                       refresh = 0,
                       init = "random",
                       pars = c("theta", "alpha"))
  tox_probs <- sapply(d, \(x) mean(x^exp(extract(stan_out)$theta) > phi))
  p_est <- sapply(d, \(x) mean(x^exp(extract(stan_out)$theta)))
  theta_posterior <- mean(extract(stan_out)$theta)
  alpha_summary <- summary(extract(stan_out)$alpha)
  rm(stan_out)
  
  mtd_est <- which.min(abs(p_est - phi))
  
  return(list(d1 = mtd_est, posteriors = list(ps = p_est, theta = theta_posterior,  alpha_summary = alpha_summary, tox_probs = tox_probs)))
}


##### MTD selection algorithms

logistic_mtd <- function(trial, phi, ...){
  trial <- trim(trial)
  if(var(trial$d1) == 0){ return(trial$d1[1])}
  suppressWarnings(reg <- glm(y ~ d1, family = binomial, data = trial))
  p_hats <- data.frame(d1 = sort(unique(na.omit(trial$d1))))
  p_hats$y <- predict(reg, newdata = p_hats, type = "response")
  mtd <- lad(x = p_hats$d1, y = p_hats$y, target = phi)
  return(mtd)
}

alphacrm_mtd <- function(...){
  args <- list(...)
  args$final_decision <- FALSE
  mtd <- do.call(alphacrm_decision, args)
  return(mtd)
}

crm_mtd <- function(...){
  args <- list(...)
  args$final_decision <- FALSE
  mtd <- do.call(crm_decision, args)
  return(mtd)
}

pocrm_mtd <- function(...){
  args <- list(...)
  args$final_decision <- FALSE
  mtd <- do.call(pocrm_decision, args)
  return(mtd)
}



simulate <- function(iters = 1000, task_id = NULL, folder_name = "output", phi = 0.3, decision_rules, true_dlts, alphas = 0, ...){
  if(!exists("crm_decision_model")){
  crm_decision_model <- stan_model(file = "stan_files/crm_decision.stan",
                                   auto_write = FALSE,
                                   save_dso = FALSE, 
                                   allow_optimizations = TRUE)}
  if(!exists("alphacrm_decision_model")){
  alphacrm_decision_model <- stan_model(file = "stan_files/alphacrm_decision.stan",
                                        auto_write = FALSE,
                                        save_dso = FALSE, 
                                        allow_optimizations = TRUE)}
  
  iter_groups <- aggregate(1:iters ~ cut_interval(1:iters, length = 501), FUN = c)[[2]]
  if(iters <= 501){iter_groups <- list(1:iters)}
  if(is.null(names(true_dlts))){names(true_dlts) <- paste0("scenario",1:length(true_dlts))}
  par <- "alpha"
  par_vals <- alphas
  
  sim_schedule <- expand.grid(decision_rule = 1:length(decision_rules), par = 1:length(par_vals), true_dlt = 1:length(true_dlts))
  sched_rows <- (task_id-1)*length(decision_rules) + 1:length(decision_rules)
  
  for(i in sched_rows){
    print(i)
    sched <- sim_schedule[i,]
    
    args <- list(...)
    args$phi <- phi
    args$ipde <- ifelse(decision_rules[sched$decision_rule] == "crm0", FALSE, TRUE)
    args$crm_decision_model <- crm_decision_model
    args$alphacrm_decision_model <- alphacrm_decision_model
    args$response_sim_args$true_dlt <- true_dlts[[sched$true_dlt]]
    args$response_sim_args[[par]] <- par_vals[sched$par]
    args$decision_rule <- switch(decision_rules[sched$decision_rule], crm = crm_decision, crm0 = crm_decision, alphacrm = alphacrm_decision, pocrm = pocrm_decision)
    args$final_decision_rule <- switch(decision_rules[sched$decision_rule], crm = crm_mtd, crm0 = crm_mtd, alphacrm = alphacrm_mtd, pocrm = pocrm_mtd)
    
    for(iter_group in iter_groups){
      trials <- mclapply(iter_group, 
                         FUN = do.call,
                         what = sim_trials,
                         args = args,
                         mc.preschedule = FALSE)
    
      summary_stats <- data.frame(
        scenario = names(true_dlts)[sched$true_dlt],
        mtd = max(c(0, which(args$response_sim_args$true_dlt <= (phi + 0.0001)))),
        parameter = paste0(par,"_",par_vals[sched$par]), 
        decision_rule = decision_rules[sched$decision_rule],
        task_id = task_id,
        sched_rows = i,
        iter = iter_group,
        key = paste0(i, "-",iter_group))
      summary_stats$n <- sapply(trials, "[[", "n")
      summary_stats$i <- sapply(trials, function(trial){sum(!is.na(trial$y), na.rm = TRUE)})
      summary_stats$mtd_estimate <- sapply(trials, "[[", "mtd_estimate")
      summary_stats$mtd <- sapply(trials, "[[", "mtd")
      summary_stats$days <- sapply(trials, function(trial){max(trial$t_eval[!is.na(trial$t_eval)], na.rm = TRUE)})
      summary_stats$ipdes <- sapply(trials, function(trial){sum(trial$ipde[!is.na(trial$t_eval)], na.rm = TRUE)})
      summary_stats$dlts <- sapply(trials, function(trial){sum(trial$y, na.rm = TRUE)})
      summary_stats$assignments_to_mtd <- sapply(trials, function(trial){sum(trial$d1[!is.na(trial$t_eval)] == trial$mtd, na.rm = TRUE)})
      summary_stats$assignments_over_mtd <- sapply(trials, function(trial){sum(trial$d1[!is.na(trial$t_eval)] > trial$mtd, na.rm = TRUE)})
      summary_stats$assignments_below_mtd <- sapply(trials, function(trial){sum(trial$d1[!is.na(trial$t_eval)] < trial$mtd, na.rm = TRUE)})
      summary_stats$toxic_assignments <- sapply(trials, function(trial){sum(trial$dlt_rate > (phi + 0.000001), na.rm = TRUE)})
      all_n <- sapply(trials,length)
      trials <- rbindlist(lapply(trials, "[", c("pid", "ipde", "t", "t_eval", "y", "d1", "dlt_rate", "beta_shift", "min_d", "max_d")))
      trials[,iter := c(unlist(mapply(FUN = rep, iter_group, each = all_n, SIMPLIFY = TRUE)))]
      trials[,key := paste0(i,"-",iter)]
        
      folders <- paste0("sched_row=",i,"/",
                         "decision_rule=", summary_stats$decision_rule[1],"/",
                         "scenario=", summary_stats$scenario[1],"/",
                         "parameter=", summary_stats$parameter[1],"/")
      dir.create(paste0("./output/Data/", folder_name, "_summary/",folders), recursive = TRUE)
      dir.create(paste0("./output/Data/", folder_name, "_full/",folders), recursive = TRUE)
      
      filename <- paste0("iterations_",min(iter_group),"_",max(iter_group),".parquet")
      
      write_parquet(summary_stats, sink = paste0("./output/Data/", folder_name, "_summary/", folders, filename))
      write_parquet(trials, sink = paste0("./output/Data/", folder_name, "_full/", folders, filename))
      
      rm(crm_decision_model)
      rm(alphacrm_decision_model)
      gc(verbose = FALSE)
      crm_decision_model <- stan_model(file = "stan_files/crm_decision.stan",
                                       auto_write = FALSE,
                                       save_dso = FALSE, 
                                       allow_optimizations = TRUE)
        
      alphacrm_decision_model <- stan_model(file = "stan_files/alphacrm_decision.stan",
                                            auto_write = FALSE,
                                            save_dso = FALSE, 
                                            allow_optimizations = TRUE)
    }
  }
}


print.ipdeResults <- function(results){
  missing <- sum(sapply(results$trials, class)!="ipdeTrial")
  all_n <- length(results$scenario)
  iterations <- max(results$iter)
  scenarios <- sort(unique(results$scenario))
  parameters <- sort(unique(results$parameter))
  decision_rules <- sort(unique(results$decision_rule))
  
  cat(paste0("\nTrials per scenario: ", iterations, "\n",
             "Total trials: ",all_n,"\n",
             "Missing: ", missing,"\n\n",
             "Dose response curves: ", paste0(scenarios, collapse = ", "),"\n",
             "Decision rules: ", paste0(decision_rules, collapse = ", "),"\n",
             "Parameters: ", paste0(parameters, collapse = ", "),"\n "))
}


plot.ipdeTrial <- function(trial){
  D <- trial$max_d[1]
  mtd <- trial$mtd
  trial <- as.data.frame(trial)
  trial$d1 <- factor(trial$d1, levels = 1:D)
  trial$assignment <- 1:nrow(trial)
  ggplot(trial[!is.na(trial$d1),], 
         aes(x = assignment, 
             y = d1, 
             color = factor(ipde),
             shape = factor(y)
         ))+
    geom_hline(yintercept = factor(mtd), color = "gold")+
    geom_point(size = 3)+
    scale_color_manual(values = c('#94b8b8', 'red'))+
    labs(color = "IPDE", x = "Assignment", y = "d1", shape = "DLT")
  
}


