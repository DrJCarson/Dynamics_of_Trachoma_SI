load(file = "Chain1_Batch1.rdata")

# Load libraries
library(MASS)
library(Rcpp)
library(MCMCpack)

# Source c++ code
Rcpp::sourceCpp('Trachoma_cpp_rif.cpp')

# Source MCMC R code
source("Trachoma_MCMC_rif.R")

acc_hist_ber <- 0
acc_hist_jal <- 0

acc_tr_ber <- 0
acc_tr_jal <- 0

acc_i <- 0

acc_id_1 <- 0
acc_id_2 <- 0

acc_d_1 <- 0
acc_d_2 <- 0

acc_rid <- 0

acc_reinf <- 0

acc_re_hist_ber <- 0
acc_re_hist_jal <- 0

acc_iffbs_mean <- numeric(mcmc_length)
acc_iffbs_ind <- numeric(individuals_count_total)

acc_missing_ages <- numeric(length(missing_ages))

runtime <- system.time({

  for (i in 1:mcmc_length) {
    
    if (!(i%%10)) {
      
      print(i)
      
    }
    
    
    ## Update missing ages
    for (j in missing_ages) {
      
      if (village_num[j] == 1) {
        
        age_prop <- sample(as.integer(age_raw[which(village_num == 1 & !is.na(age_raw))]), size = 1)
        
        if (as.integer(states[j, 1] / disease_comps) == 0) {
          
          mh_age <- dpois(0, (age_prop + 0.5) * berending_history) / dpois(0, (age[j] + 0.5) * berending_history)
          
        } else {
          
          mh_age <- ppois(0, (age_prop + 0.5) * berending_history, lower.tail = F) / ppois(0, (age[j] + 0.5) * berending_history, lower.tail = F)
          
        }
        
      } else {
        
        age_prop <- sample(as.integer(age_raw[which(village_num == 2 & !is.na(age_raw))]), size = 1)
        
        if (as.integer(states[j, 1] / disease_comps) == 0) {
          
          mh_age <- dpois(0, (age_prop + 0.5) * jali_history) / dpois(0, (age[j] + 0.5) * jali_history)
          
        } else {
          
          mh_age <- ppois(0, (age_prop + 0.5) * jali_history, lower.tail = F) / ppois(0, (age[j] + 0.5) * jali_history, lower.tail = F)
          
        }
        
      }
      
      if (runif(1) < mh_age) {
        
        age[j] <- as.integer(age_prop)
        acc_missing_ages[which(missing_ages == j)] <- acc_missing_ages[which(missing_ages == j)] + 1
        
      }
      
      missing_ages_chain[i, which(missing_ages == j)] <- age[j]
      
    }

    ## Update initial state probabilities
    berending_initial <- init_gibbs(states[, 1], 1, village_num, villages_count, disease_groups, disease_comps, prior_initial)
    jali_initial <- init_gibbs(states[, 1], 2, village_num, villages_count, disease_groups, disease_comps, prior_initial)
    
    berending_initial_chain[i,,] <- berending_initial
    jali_initial_chain[i,,] <- jali_initial
    
    
    ## Update observation parameters
    obs_parms <- obs_gibbs(obs_pcr, obs_att, obs_cex, states[, obs_days + 1], disease_groups, disease_comps, prior_pcr_specificity_shape_1, prior_pcr_specificity_shape_2, prior_pcr_sensitivity_shape_1, prior_pcr_sensitivity_shape_2, prior_att_specificity_shape_1, prior_att_specificity_shape_2, prior_att_sensitivity_shape_1, prior_att_sensitivity_shape_2, prior_cex_specificity_shape_1, prior_cex_specificity_shape_2, prior_cex_sensitivity_shape_1, prior_cex_sensitivity_shape_2)
    
    pcr_sen_chain[i] <- obs_parms["pcr_sen"]
    pcr_spec_chain[i] <- obs_parms["pcr_spec"]
    att_sen_chain[i] <- obs_parms["att_sen"]
    att_spec_chain[i] <- obs_parms["att_spec"]
    cex_sen_chain[i] <- obs_parms["cex_sen"]
    cex_spec_chain[i] <- obs_parms["cex_spec"]
    
    
    
    berending_history_2 <- rnorm(1, mean = berending_history, sd = sqrt(mh_hist_ber_cov))
    
    if (berending_history_2 < 0) {
      
      berending_history_2 <- berending_history
      
    }
    
    jali_history_2 <- rnorm(1, mean = jali_history, sd = sqrt(mh_hist_jal_cov))
    
    if (jali_history_2 < 0) {
      
      jali_history_2 <- jali_history
      
    }
    
    tr_parms <- c(trb_ber = berending_background, trv_ber = berending_village,
                  trc_ber = berending_compound, trr_ber = berending_room,
                  trb_jal = jali_background, trv_jal = jali_village,
                  trc_jal = jali_compound, trr_jal = jali_room)
    
    tr_parms_2 <- tr_parms + c(mvrnorm(1, mu = rep(0, 4), Sigma = mh_tr_ber_cov), mvrnorm(1, mu = rep(0, 4), Sigma = mh_tr_jal_cov))
    
    if (any(tr_parms_2[1:4] < 0)) {
      
      tr_parms_2[1:4] <- tr_parms[1:4]
      
      ber_ss_2 <- ber_ss
      
      ber_dd_2 <- ber_dd
      
    } else {
      
      ber_ss_2 <- exp(-(tr_parms_2[1] + 
                          inf_base[, 1] * tr_parms_2[2] + 
                          inf_base[, 2] * tr_parms_2[3] + 
                          inf_base[, 3] * tr_parms_2[4]))
      
      
      ber_dd_2 <- exp(-reinfection * (tr_parms_2[1] + 
                                        inf_base[, 1] * tr_parms_2[2] + 
                                        inf_base[, 2] * tr_parms_2[3] + 
                                        inf_base[, 3] * tr_parms_2[4]))
      
    }
    
    if (any(tr_parms_2[5:8] < 0)) {
      
      tr_parms_2[5:8] <- tr_parms[5:8]
      
      jal_ss_2 <- jal_ss
      
      jal_dd_2 <- jal_dd
      
    } else {
      
      jal_ss_2 <- exp(-(tr_parms_2[5] + 
                          inf_base[, 1] * tr_parms_2[6] + 
                          inf_base[, 2] * tr_parms_2[7] + 
                          inf_base[, 3] * tr_parms_2[8]))
      
      
      jal_dd_2 <- exp(-reinfection * (tr_parms_2[5] + 
                                        inf_base[, 1] * tr_parms_2[6] + 
                                        inf_base[, 2] * tr_parms_2[7] + 
                                        inf_base[, 3] * tr_parms_2[8]))
      
    }

    
    i_parms <- mvrnorm(1, mu = c(i_to_id_prob, i_to_id_shape), Sigma = mh_i_cov)

    i_to_id_prob_2 <- i_parms[1]
    i_to_id_shape_2 <- i_parms[2]
    
    i_to_id_mean_2 <- 1 + (1 - i_to_id_prob_2) * i_to_id_shape_2 / i_to_id_prob_2
    
    if (i_to_id_prob_2 < 0 | i_to_id_prob_2 > 1 | i_to_id_shape_2 < 1 | i_to_id_mean_2 > 10000) {
      
      i_to_id_prob_2 <- i_to_id_prob
      i_to_id_shape_2 <- i_to_id_shape
      
      i_to_id_mean_2 <- i_to_id_mean
      
    }
    

    id_parms_1 <- mvrnorm(1, mu = c(id_to_d_prob[1], id_to_d_shape[1]), Sigma = mh_id_cov_1)
    id_parms_2 <- mvrnorm(1, mu = c(id_to_d_prob[2], id_to_d_shape[2]), Sigma = mh_id_cov_2)

    id_to_d_prob_2 <- c(id_parms_1[1], id_parms_2[1])
    id_to_d_shape_2 <- c(id_parms_1[2], id_parms_2[2])
    
    id_to_d_mean_2 <- 1 + (1 - id_to_d_prob_2) * id_to_d_shape_2 / id_to_d_prob_2
    
    for (g in 1:disease_groups) {
      
      if (id_to_d_prob_2[g] < 0 | id_to_d_prob_2[g] > 1 | id_to_d_shape_2[g] < 1 | id_to_d_mean_2[g] > 10000) {
        
        id_to_d_prob_2[g] <- id_to_d_prob[g]
        id_to_d_shape_2[g] <- id_to_d_shape[g]
        
        id_to_d_mean_2[g] <- id_to_d_mean[g]
        
      }
      
    }
    

    d_parms_1 <- mvrnorm(1, mu = c(d_to_s_prob[1], d_to_s_shape[1]), Sigma = mh_d_cov_1)
    d_parms_2 <- mvrnorm(1, mu = c(d_to_s_prob[2], d_to_s_shape[2]), Sigma = mh_d_cov_2)
    
    d_to_s_prob_2 <- c(d_parms_1[1], d_parms_2[1])
    d_to_s_shape_2 <- c(d_parms_1[2], d_parms_2[2])
    
    d_to_s_mean_2 <- 1 + (1 - d_to_s_prob_2) * d_to_s_shape_2 / d_to_s_prob_2
    
    for (g in 1:disease_groups) {
      
      if (d_to_s_prob_2[g] < 0 | d_to_s_prob_2[g] > 1 | d_to_s_shape_2[g] < 1 | d_to_s_mean_2[g] > 10000) {
        
        d_to_s_prob_2[g] <- d_to_s_prob[g]
        d_to_s_shape_2[g] <- d_to_s_shape[g]
        
        d_to_s_mean_2[g] <- d_to_s_mean[g]
        
      }
      
    }
    
    # Length biased distributions
    
    lb_i_pdf_2 <- pnbinom(1:time_count - 2, size = i_to_id_shape_2, prob = i_to_id_prob_2, lower.tail = F) / i_to_id_mean_2
    lb_i_pdf_2[time_count] <- 1 - sum(lb_i_pdf_2[1:(time_count - 1)])
    lb_i_cen_2 <- rev(cumsum(rev(lb_i_pdf_2)))
    
    # Length biased distributions
    
    lb_id_pdf_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_id_cen_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    for (g in 1:disease_groups) {
      
      lb_id_pdf_2[, g] <- pnbinom(1:time_count - 2, size = id_to_d_shape_2[g], prob = id_to_d_prob_2[g], lower.tail = F) / id_to_d_mean_2[g]
      lb_id_pdf_2[time_count, g] <- 1 - sum(lb_id_pdf_2[1:(time_count - 1), g])
      lb_id_cen_2[, g] <- rev(cumsum(rev(lb_id_pdf_2[, g])))
      
    }
    
    # Length biased distributions
    lb_d_pdf_ber_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_cen_ber_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_pdf_jal_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_cen_jal_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    
    for (g in 1:disease_groups) {
      
      # Approximated historic holding times for Berending and Jali
      hist_pdf_ber_2 <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g]) *
        pgeom(1:hist_horiz - 1, prob = reinfec_hist_ber, lower.tail = F) +
        pnbinom(1:hist_horiz - 2, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g], lower.tail = F) *
        dgeom(1:hist_horiz - 1, prob = reinfec_hist_ber)
      
      hist_pdf_jal_2 <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g]) *
        pgeom(1:hist_horiz - 1, prob = reinfec_hist_jal, lower.tail = F) +
        pnbinom(1:hist_horiz - 2, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g], lower.tail = F) *
        dgeom(1:hist_horiz - 1, prob = reinfec_hist_jal)
      
      hist_cdf_rt_ber_2 <- rev(cumsum(rev(hist_pdf_ber_2)))
      hist_mean_ber_2 <- sum((1:hist_horiz) * hist_pdf_ber_2)
      hist_lb_pdf_ber_2 <- hist_cdf_rt_ber_2 / hist_mean_ber_2
      
      hist_cdf_rt_jal_2 <- rev(cumsum(rev(hist_pdf_jal_2)))
      hist_mean_jal_2 <- sum((1:hist_horiz) * hist_pdf_jal_2)
      hist_lb_pdf_jal_2 <- hist_cdf_rt_jal_2 / hist_mean_jal_2
      
      nb_pdf_2 <- dnbinom(1:(time_count + hist_horiz - 1) - 1, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g])
      for (h in 1:hist_horiz) {
        
        if (pnbinom(h - 2, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g], lower.tail = F) == 0) {
          
          break
          
        }
        
        lb_d_pdf_ber_2[, g] <- lb_d_pdf_ber_2[, g] + hist_lb_pdf_ber_2[h] * nb_pdf_2[1:time_count + h - 1] /
          pnbinom(h - 2, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g], lower.tail = F)
        
        lb_d_pdf_jal_2[, g] <- lb_d_pdf_jal_2[, g] + hist_lb_pdf_jal_2[h] * nb_pdf_2[1:time_count + h - 1] /
          pnbinom(h - 2, size = d_to_s_shape_2[g], prob = d_to_s_prob_2[g], lower.tail = F)
        
      }
      
      lb_d_pdf_ber_2[time_count, g] <- 1 - sum(lb_d_pdf_ber_2[1:(time_count - 1), g])
      lb_d_cen_ber_2[, g] <- rev(cumsum(rev(lb_d_pdf_ber_2[, g])))
      
      lb_d_pdf_jal_2[time_count, g] <- 1 - sum(lb_d_pdf_jal_2[1:(time_count - 1), g])
      lb_d_cen_jal_2[, g] <- rev(cumsum(rev(lb_d_pdf_jal_2[, g])))
      
    }
    
    prob_hist <- rep(1, 2)
    prob_tr <- rep(1, 2)
    prob_i <- 1
    prob_id <- rep(1, disease_groups)
    prob_d <- rep(1, disease_groups)
    
    prob_mh <- rep(1, 8)
    
    MH_1(individuals_count_total,
         village_num,
         villages_count,
         compound_num,
         compounds_count,
         room_num,
         rooms_count,
         age,
         time_count,
         disease_groups,
         disease_comps,
         berending_history,
         berending_history_2,
         jali_history,
         jali_history_2,
         ber_ss,
         ber_ss_2,
         ber_dd,
         ber_dd_2,
         jal_ss,
         jal_ss_2,
         jal_dd,
         jal_dd_2,
         i_to_id_prob,
         i_to_id_prob_2,
         i_to_id_shape,
         i_to_id_shape_2,
         id_to_d_prob,
         id_to_d_prob_2,
         id_to_d_shape,
         id_to_d_shape_2,
         d_to_s_prob,
         d_to_s_prob_2,
         d_to_s_shape,
         d_to_s_shape_2,
         lb_i_pdf,
         lb_i_pdf_2,
         lb_id_pdf,
         lb_id_pdf_2,
         lb_d_pdf_ber,
         lb_d_pdf_ber_2,
         lb_d_cen_ber,
         lb_d_cen_ber_2,
         lb_d_pdf_jal,
         lb_d_pdf_jal_2,
         lb_d_cen_jal,
         lb_d_cen_jal_2,
         reinfection,
         states,
         v_ip,
         c_ip,
         r_ip,
         prob_hist,
         prob_tr,
         prob_i,
         prob_id,
         prob_d)
    
    prob_hist_ber <- prob_hist[1]
    prob_hist_jal <- prob_hist[2]
    prob_tr_ber <- prob_tr[1]
    prob_tr_jal <- prob_tr[2]
    
    acc_prob <- prob_hist_ber * (dgamma(berending_history_2, shape = prior_history_shape, rate = prior_history_rate) / dgamma(berending_history, shape = prior_history_shape, rate = prior_history_rate))
    
    if (runif(1) < acc_prob) {
      
      berending_history <- berending_history_2
      acc_hist_ber <- acc_hist_ber + 1
      
    }
    
    berending_history_chain[i] <- berending_history
    
    
    acc_prob <- prob_hist_jal * (dgamma(jali_history_2, shape = prior_history_shape, rate = prior_history_rate) / dgamma(jali_history, shape = prior_history_shape, rate = prior_history_rate))
    
    if (runif(1) < acc_prob) {
      
      jali_history <- jali_history_2
      acc_hist_jal <- acc_hist_jal + 1
      
    }
    
    jali_history_chain[i] <- jali_history
    
    
    acc_prob <- prob_tr_ber *
      (dexp(tr_parms_2["trb_ber"], rate = prior_background_transmission) / dexp(tr_parms["trb_ber"], rate = prior_background_transmission)) *
      (dexp(tr_parms_2["trv_ber"], rate = prior_village_transmission) / dexp(tr_parms["trv_ber"], rate = prior_village_transmission)) *
      (dexp(tr_parms_2["trc_ber"], rate = prior_compound_transmission) / dexp(tr_parms["trc_ber"], rate = prior_compound_transmission)) *
      (dexp(tr_parms_2["trr_ber"], rate = prior_room_transmission) / dexp(tr_parms["trr_ber"], rate = prior_room_transmission))
    
    if (runif(1) < acc_prob) {
      
      berending_background <- unname(tr_parms_2["trb_ber"])
      berending_village <- unname(tr_parms_2["trv_ber"])
      berending_compound <- unname(tr_parms_2["trc_ber"])
      berending_room <- unname(tr_parms_2["trr_ber"])
      
      tr_parms[1:4] <- tr_parms_2[1:4]
      
      ber_ss <- ber_ss_2
      ber_dd <- ber_dd_2
      
      acc_tr_ber <- acc_tr_ber + 1
      
    }
    
    berending_background_chain[i] <- berending_background
    berending_village_chain[i] <- berending_village
    berending_compound_chain[i] <- berending_compound
    berending_room_chain[i] <- berending_room
    
    
    acc_prob <- prob_tr_jal *
      (dexp(tr_parms_2["trb_jal"], rate = prior_background_transmission) / dexp(tr_parms["trb_jal"], rate = prior_background_transmission)) *
      (dexp(tr_parms_2["trv_jal"], rate = prior_village_transmission) / dexp(tr_parms["trv_jal"], rate = prior_village_transmission)) *
      (dexp(tr_parms_2["trc_jal"], rate = prior_compound_transmission) / dexp(tr_parms["trc_jal"], rate = prior_compound_transmission)) *
      (dexp(tr_parms_2["trr_jal"], rate = prior_room_transmission) / dexp(tr_parms["trr_jal"], rate = prior_room_transmission))
    
    if (runif(1) < acc_prob) {
      
      jali_background <- unname(tr_parms_2["trb_jal"])
      jali_village <- unname(tr_parms_2["trv_jal"])
      jali_compound <- unname(tr_parms_2["trc_jal"])
      jali_room <- unname(tr_parms_2["trr_jal"])
      
      tr_parms[5:8] <- tr_parms_2[5:8]
      
      jal_ss <- jal_ss_2
      jal_dd <- jal_dd_2
      
      acc_tr_jal <- acc_tr_jal + 1
      
    }
    
    
    jali_background_chain[i] <- jali_background
    jali_village_chain[i] <- jali_village
    jali_compound_chain[i] <- jali_compound
    jali_room_chain[i] <- jali_room
    
    
    acc_prob <- prob_i *
      (dbeta(i_to_id_prob_2, shape1 = prior_state_shape_1, shape2 = prior_state_shape_2) / dbeta(i_to_id_prob, shape1 = prior_state_shape_1, shape2 = prior_state_shape_2)) *
      (dexp(i_to_id_shape_2, rate = prior_state_rate) / dexp(i_to_id_shape, rate = prior_state_rate))
    
    if (runif(1) < acc_prob) {
      
      i_to_id_prob <- i_to_id_prob_2
      i_to_id_shape <- i_to_id_shape_2
      i_to_id_mean <- i_to_id_mean_2
      
      lb_i_pdf <- lb_i_pdf_2
      
      acc_i <- acc_i + 1
      
    }
    
    i_to_id_prob_chain[i] <- i_to_id_prob
    i_to_id_shape_chain[i] <- i_to_id_shape
    
    
    acc_prob <- prob_id[1] *
      (dbeta(id_to_d_prob_2[1], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2) / dbeta(id_to_d_prob[1], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2)) *
      (dexp(id_to_d_shape_2[1], rate = prior_state_rate) / dexp(id_to_d_shape[1], rate = prior_state_rate))
    
    if (runif(1) < acc_prob & id_to_d_mean_2[1] > id_to_d_mean[2]) {
      
      id_to_d_prob[1] <- id_to_d_prob_2[1]
      id_to_d_shape[1] <- id_to_d_shape_2[1]
      id_to_d_mean[1] <- id_to_d_mean_2[1]
      
      lb_id_pdf[, 1] <- lb_id_pdf_2[, 1]
      
      acc_id_1 <- acc_id_1 + 1
      
    }
    
    acc_prob <- prob_id[2] *
      (dbeta(id_to_d_prob_2[2], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2) / dbeta(id_to_d_prob[2], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2)) *
      (dexp(id_to_d_shape_2[2], rate = prior_state_rate) / dexp(id_to_d_shape[2], rate = prior_state_rate))
    
    if (runif(1) < acc_prob & id_to_d_mean_2[2] < id_to_d_mean[1]) {
      
      id_to_d_prob[2] <- id_to_d_prob_2[2]
      id_to_d_shape[2] <- id_to_d_shape_2[2]
      id_to_d_mean[2] <- id_to_d_mean_2[2]
      
      lb_id_pdf[, 2] <- lb_id_pdf_2[, 2]
      
      acc_id_2 <- acc_id_2 + 1
      
    }
    
    id_to_d_prob_chain[i,] <- id_to_d_prob
    id_to_d_shape_chain[i,] <- id_to_d_shape
    
    
    acc_prob <- prob_d[1] *
      (dbeta(d_to_s_prob_2[1], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2) / dbeta(d_to_s_prob[1], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2)) *
      (dexp(d_to_s_shape_2[1], rate = prior_state_rate) / dexp(d_to_s_shape[1], rate = prior_state_rate))
    
    if (runif(1) < acc_prob & d_to_s_mean_2[1] > d_to_s_mean[2]) {
      
      d_to_s_prob[1] <- d_to_s_prob_2[1]
      d_to_s_shape[1] <- d_to_s_shape_2[1]
      d_to_s_mean[1] <- d_to_s_mean_2[1]
      
      lb_d_pdf_ber[, 1] <- lb_d_pdf_ber_2[, 1]
      lb_d_cen_ber[, 1] <- lb_d_cen_ber_2[, 1]
      lb_d_pdf_jal[, 1] <- lb_d_pdf_jal_2[, 1]
      lb_d_cen_jal[, 1] <- lb_d_cen_jal_2[, 1]
      
      acc_d_1 <- acc_d_1 + 1
      
    }
    
    acc_prob <- prob_d[2] *
      (dbeta(d_to_s_prob_2[2], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2) / dbeta(d_to_s_prob[2], shape1 = prior_state_shape_1, shape2 = prior_state_shape_2)) *
      (dexp(d_to_s_shape_2[2], rate = prior_state_rate) / dexp(d_to_s_shape[2], rate = prior_state_rate))
    
    if (runif(1) < acc_prob & d_to_s_mean_2[2] < d_to_s_mean[1]) {
      
      d_to_s_prob[2] <- d_to_s_prob_2[2]
      d_to_s_shape[2] <- d_to_s_shape_2[2]
      d_to_s_mean[2] <- d_to_s_mean_2[2]
      
      lb_d_pdf_ber[, 2] <- lb_d_pdf_ber_2[, 2]
      lb_d_cen_ber[, 2] <- lb_d_cen_ber_2[, 2]
      lb_d_pdf_jal[, 2] <- lb_d_pdf_jal_2[, 2]
      lb_d_cen_jal[, 2] <- lb_d_cen_jal_2[, 2]
      
      acc_d_2 <- acc_d_2 + 1
      
    }
    
    d_to_s_prob_chain[i,] <- d_to_s_prob
    d_to_s_shape_chain[i,] <- d_to_s_shape
    
    ################################################################################
 
    reinfection_2 <- rnorm(1, mean = reinfection, sd = sqrt(mh_reinf_cov))
    
    if (reinfection_2 < 0) {
      
      reinfection_2 <- reinfection
      
      ber_dd_2 <- ber_dd
      
      jal_dd_2 <- jal_dd
      
    } else {
      
      ber_dd_2 <- exp(-reinfection_2 * (tr_parms[1] + 
                                          inf_base[, 1] * tr_parms[2] + 
                                          inf_base[, 2] * tr_parms[3] + 
                                          inf_base[, 3] * tr_parms[4]))
      
      jal_dd_2 <- exp(-reinfection_2 * (tr_parms[5] + 
                                          inf_base[, 1] * tr_parms[6] + 
                                          inf_base[, 2] * tr_parms[7] + 
                                          inf_base[, 3] * tr_parms[8]))
      
    }
    
    reinfec_hist_ber_2 <- rnorm(1, mean = reinfec_hist_ber, sd = sqrt(mh_re_hist_ber_cov))
    
    if (reinfec_hist_ber_2 < 0) {
      
      reinfec_hist_ber_2 <- reinfec_hist_ber
      
    }
    
    reinfec_hist_jal_2 <- rnorm(1, mean = reinfec_hist_jal, sd = sqrt(mh_re_hist_jal_cov))
    
    if (reinfec_hist_jal_2 < 0) {
      
      reinfec_hist_jal_2 <- reinfec_hist_jal
      
    }
    
    # Length biased distributions
    lb_d_pdf_ber_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_cen_ber_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_pdf_jal_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    lb_d_cen_jal_2 <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
    
    for (g in 1:disease_groups) {
      
      # Approximated historic holding times for Berending and Jali
      hist_pdf_ber_2 <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g]) *
        pgeom(1:hist_horiz - 1, prob = reinfec_hist_ber_2, lower.tail = F) +
        pnbinom(1:hist_horiz - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F) *
        dgeom(1:hist_horiz - 1, prob = reinfec_hist_ber_2)
      
      hist_pdf_jal_2 <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g]) *
        pgeom(1:hist_horiz - 1, prob = reinfec_hist_jal_2, lower.tail = F) +
        pnbinom(1:hist_horiz - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F) *
        dgeom(1:hist_horiz - 1, prob = reinfec_hist_jal_2)
      
      hist_cdf_rt_ber_2 <- rev(cumsum(rev(hist_pdf_ber_2)))
      hist_mean_ber_2 <- sum((1:hist_horiz) * hist_pdf_ber_2)
      hist_lb_pdf_ber_2 <- hist_cdf_rt_ber_2 / hist_mean_ber_2
      
      hist_cdf_rt_jal_2 <- rev(cumsum(rev(hist_pdf_jal_2)))
      hist_mean_jal_2 <- sum((1:hist_horiz) * hist_pdf_jal_2)
      hist_lb_pdf_jal_2 <- hist_cdf_rt_jal_2 / hist_mean_jal_2
      
      nb_pdf <- dnbinom(1:(time_count + hist_horiz - 1) - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g])
      for (h in 1:hist_horiz) {
        
        if (pnbinom(h - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F) == 0) {
          
          break
          
        }
        
        lb_d_pdf_ber_2[, g] <- lb_d_pdf_ber_2[, g] + hist_lb_pdf_ber_2[h] * nb_pdf[1:time_count + h - 1] /
          pnbinom(h - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F)
        
        lb_d_pdf_jal_2[, g] <- lb_d_pdf_jal_2[, g] + hist_lb_pdf_jal_2[h] * nb_pdf[1:time_count + h - 1] /
          pnbinom(h - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F)
        
      }
      
      lb_d_pdf_ber_2[time_count, g] <- 1 - sum(lb_d_pdf_ber_2[1:(time_count - 1), g])
      lb_d_cen_ber_2[, g] <- rev(cumsum(rev(lb_d_pdf_ber_2[, g])))
      
      lb_d_pdf_jal_2[time_count, g] <- 1 - sum(lb_d_pdf_jal_2[1:(time_count - 1), g])
      lb_d_cen_jal_2[, g] <- rev(cumsum(rev(lb_d_pdf_jal_2[, g])))
      
    }
    
    prob_reinf <- 1
    prob_re_hist_ber <- 1
    prob_re_hist_jal <- 1
    
    MH_2(individuals_count_total,
         village_num,
         villages_count,
         compound_num,
         compounds_count,
         room_num,
         rooms_count,
         time_count,
         disease_groups,
         disease_comps,
         berending_history,
         jali_history,
         ber_dd,
         ber_dd_2,
         jal_dd,
         jal_dd_2,
         i_to_id_prob,
         i_to_id_shape,
         id_to_d_prob,
         id_to_d_shape,
         d_to_s_prob,
         d_to_s_shape,
         lb_d_pdf_ber,
         lb_d_pdf_ber_2,
         lb_d_cen_ber,
         lb_d_cen_ber_2,
         lb_d_pdf_jal,
         lb_d_pdf_jal_2,
         lb_d_cen_jal,
         lb_d_cen_jal_2,
         reinfection,
         reinfection_2,
         states,
         v_ip,
         c_ip,
         r_ip,
         prob_reinf,
         prob_re_hist_ber,
         prob_re_hist_jal)
    
    
    acc_prob <- prob_reinf *
      (dgamma(reinfection_2, shape = prior_reinfection_shape, rate = prior_reinfection_rate) / dgamma(reinfection, shape = prior_reinfection_shape, rate = prior_reinfection_rate))
    
    if (runif(1) < acc_prob) {
      
      reinfection <- reinfection_2
      
      ber_dd <- ber_dd_2
      jal_dd <- jal_dd_2
      
      acc_reinf <- acc_reinf + 1
      
    }
    
    reinfection_chain[i] <- reinfection
    
    
    
    acc_prob <- prob_re_hist_ber *
      (dexp(reinfec_hist_ber_2, rate = prior_re_vil_rate) / dexp(reinfec_hist_ber, rate = prior_re_vil_rate))
    
    if (runif(1) < acc_prob) {
      
      reinfec_hist_ber <- reinfec_hist_ber_2
      
      lb_d_pdf_ber[, 1] <- lb_d_pdf_ber_2[, 1]
      lb_d_cen_ber[, 1] <- lb_d_cen_ber_2[, 1]
      lb_d_pdf_ber[, 2] <- lb_d_pdf_ber_2[, 2]
      lb_d_cen_ber[, 2] <- lb_d_cen_ber_2[, 2]
      
      acc_re_hist_ber <- acc_re_hist_ber + 1
      
    }
    
    reinfec_hist_ber_chain[i] <- reinfec_hist_ber
    
    
    acc_prob <- prob_re_hist_jal *
      (dexp(reinfec_hist_jal_2, rate = prior_re_vil_rate) / dexp(reinfec_hist_jal, rate = prior_re_vil_rate))
    
    if (runif(1) < acc_prob) {
      
      reinfec_hist_jal <- reinfec_hist_jal_2
      
      lb_d_pdf_jal[, 1] <- lb_d_pdf_jal_2[, 1]
      lb_d_cen_jal[, 1] <- lb_d_cen_jal_2[, 1]
      lb_d_pdf_jal[, 2] <- lb_d_pdf_jal_2[, 2]
      lb_d_cen_jal[, 2] <- lb_d_cen_jal_2[, 2]
      
      acc_re_hist_jal <- acc_re_hist_jal + 1
      
    }
    
    reinfec_hist_jal_chain[i] <- reinfec_hist_jal
    
    
    
    # Geometric proposals
    i_geom <- 1 / i_to_id_mean
    id_geom <- 1 / id_to_d_mean
    d_geom <- 1 / d_to_s_mean

    urv_s[,] <- runif(individuals_count_total * time_count)
    urv_m <- runif(individuals_count_total)
    
    iffbs(individuals_count_total,
          village_num,
          villages_count,
          compound_num,
          compounds_count,
          room_num,
          rooms_count,
          age,
          time_count,
          disease_groups,
          disease_comps,
          obs_count,
          obs_days,
          obs_pcr,
          obs_att,
          obs_cex,
          berending_history,
          jali_history,
          berending_initial,
          jali_initial,
          ber_ss,
          ber_dd,
          jal_ss,
          jal_dd,
          i_to_id_prob,
          i_to_id_shape,
          id_to_d_prob,
          id_to_d_shape,
          d_to_s_prob,
          d_to_s_shape,
          lb_i_pdf,
          lb_id_pdf,
          lb_d_pdf_ber,
          lb_d_cen_ber,
          lb_d_pdf_jal,
          lb_d_cen_jal,
          i_geom,
          id_geom,
          d_geom,
          obs_parms,
          states,
          v_ip,
          c_ip,
          r_ip,
          acceptances,
          init_probs,
          pred_probs,
          filt_probs,
          prop_states,
          urv_s,
          urv_m)
    
    acc_iffbs_mean[i] <- mean(acceptances)
    acc_iffbs_ind <- acc_iffbs_ind + acceptances
    
    if (!(i %% mcmc_states_thin)) {
      
      states_chain[as.integer(i / mcmc_states_thin), , ] <- states
      
    }
    
    
  }
  
})[3]

save.image(file = "Chain1_Batch2.rdata", sep = "")
