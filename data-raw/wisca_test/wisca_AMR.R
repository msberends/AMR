
# Copyright (c) [2022] [Larisse Bolton (author), Aislinn Cook (contributor)]
# Adapted from: Bielicki JA, Sharland M, Heath PT, et al. Evaluation of the coverage of 3 antibiotic regimens for
# neonatal sepsis in the hospital setting across Asian countries. JAMA Netw Open.
# 2020:3(2):e1921124. doi:10.1001.jamanetworkopen.2019.21124


wisca_funct <- function(wisca_in,analysis){

#function used to gerenrate coverage estimate parameters from input data
priors <- function(wisca_prior_in){
  #INCIDENCE
  #prior Dirichlet (Gamma) parameters
  A <- rep(1, times = length(unique(wisca_prior_in$mo)))
  Multiplier <- rep(1, times = length(unique(wisca_prior_in$mo)))
  multinomial_obs <- round((wisca_prior_in$n.x)*Multiplier,0)
  
  #posterior Dirichlet (Gamma) parameters
  gamma_A <- A + multinomial_obs
  B <- rep(1, times = length(gamma_A))
  
  ##SUSCEPTIBILITY 
  #prior Beta parameters
  beta_A <- rep(1, times = length(unique(wisca_prior_in$mo)))
  beta_B <- rep(1, times = length(unique(wisca_prior_in$mo)))
  
  #Binomial distribution parameters
  r <- wisca_prior_in$S_n #number of pathogens that tested sensitive
  n_obs <- wisca_prior_in$tested_n #number of pathogens tested
  # n_tested <- round(n_obs*(wisca_calc_input$prop.sens))
  diff_nr <- n_obs - r #difference between number of pathogens that were tested and those that tested sensitive
  
  #posterior beta parameters
  post_beta_1 <- beta_A + r
  post_beta_2 <- beta_B + diff_nr
  
  priors_df <- bind_cols(gamma_A,B,post_beta_1,post_beta_2)
  names(priors_df) <- c("gamma_A","B","post_beta_1","post_beta_2")
  return(priors_df)
}

#simulations per bug
#function to generate simulations
coverage <- function(wisca_sim_in){
  sim_coverage <- wisca_sim_in %>%
    mutate(random_incidence =  runif(n = nrow(wisca_sim_in), min = 0, max = 1),
           random_susceptibility = runif(n = nrow(wisca_sim_in), min = 0, max = 1)) %>%
    mutate(simulation_bug = qgamma(random_incidence, shape = gamma_A, scale = B),
           simulation_suscep = qbeta(p = random_susceptibility, 
                                     shape1 = post_beta_1, shape2 = post_beta_2)) %>%
    mutate(sim_inc_total = sum(simulation_bug),
           simulated_incidence = simulation_bug/sim_inc_total) %>%
    mutate(simulated_bugcover = simulated_incidence *simulation_suscep,
           coverage = sum(simulated_bugcover))
  sim_coverage_out <- slice_head(sim_coverage, n = 1)$coverage  
  return(sim_coverage_out)
}

##will run this for every regimen
set.seed(1243)
cover <- data.frame()

for (gg in 1:length(unique(wisca_in$keyantimicrobials))){
   simulation_nr <- 1000
  #generate 1000 coverage estimates and determine mean and CI
  wisca <- subset(wisca_in, keyantimicrobials == unique(wisca_in$keyantimicrobials)[gg])
    
  if (analysis %in% names(wisca)){
    cover_level <- data.frame()
    for (ii in 1:length(levels(wisca[[analysis]]))){
      wisca_prep <- subset(wisca, wisca[[analysis]] == levels(wisca[[analysis]])[ii])
      params_priors <- priors(wisca_prior_in = wisca_prep)
      
    coverage_simulation_total <- replicate(n = simulation_nr,coverage(wisca_sim_in = params_priors))
    
    av_coverage <- mean(coverage_simulation_total)
    ci_coverage_lower <- quantile(coverage_simulation_total, probs = 0.025)
    ci_coverage_upper <- quantile(coverage_simulation_total, probs = 0.975)
    
    combine_out <- cbind.data.frame(unique(wisca_in$keyantimicrobials)[gg],
                                    levels(wisca[[analysis]])[ii],
                                    av_coverage,
                                    ci_coverage_lower,
                                    ci_coverage_upper,
                                    row.names = NULL)
    names(combine_out) <- c("Regimen",str_to_sentence(analysis),"Coverage","Lower_CI","Upper_CI")
    cover_level <- rbind.data.frame(cover_level,combine_out, row.names = NULL)
    }
    cover <- rbind.data.frame(cover,cover_level, row.names = NULL)
    
  
    } else {
    params_priors <- priors(wisca_prior_in = wisca)
    
    coverage_simulation_total <- replicate(n = simulation_nr,coverage(wisca_sim_in = params_priors))
    
    av_coverage <- mean(coverage_simulation_total)
    ci_coverage_lower <- quantile(coverage_simulation_total, probs = 0.025)
    ci_coverage_upper <- quantile(coverage_simulation_total, probs = 0.975)
    
    combine_out <- cbind.data.frame(unique(wisca_in$keyantimicrobials)[gg],
                                    av_coverage,
                                    ci_coverage_lower,
                                    ci_coverage_upper,
                                    row.names = NULL)
    names(combine_out) <- c("Regimen","Coverage","Lower_CI","Upper_CI")
    cover <- rbind.data.frame(cover,combine_out, row.names = NULL)
    
    }
  
}
cover <- cover %>%
  mutate(across(.cols = Coverage:Upper_CI, ~ round(.x*100,1))) %>%
  rowwise() %>%
  mutate(Regimen_full = ifelse(str_detect(Regimen, fixed("+")),
                               str_c(sapply(unlist(unique(str_split(Regimen, fixed("+")))),ab_name,USE.NAMES = FALSE), collapse = " + "), 
                               ab_name(Regimen))) %>%
  relocate(Regimen_full,.before = Regimen)

return(cover)

}
