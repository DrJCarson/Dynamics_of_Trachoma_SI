# To do

# Load libraries
library(MASS)
library(Rcpp)
library(MCMCpack)

# Source c++ code
Rcpp::sourceCpp('Trachoma_cpp_rif.cpp')

# Source MCMC R code
source("Trachoma_MCMC_rif.R")

# Set seed for reproducibility
seed <- 1
set.seed(seed)

# Read data
full_data <- read.csv(file = "Gambia_Data.csv")

################################################################################
# Individual info

# Total number of individuals
individuals_count_total <- as.integer(dim(full_data)[1])

# Indices of villagers from each village
individuals_berending <- as.integer(which(full_data[,2]=="Berending"))
individuals_jali <- as.integer(which(full_data[,2]=="Jali"))

# Number of villagers from each village
individuals_count_berending <- as.integer(length(individuals_berending))
individuals_count_jali <- as.integer(length(individuals_jali))

# Village numbers for each individual
village_num <- integer(individuals_count_total)
village_num[individuals_berending] <- as.integer(1)
village_num[individuals_jali] <- as.integer(2)

# Number of villages
villages_count <- as.integer(2)

# Compound numbers for each individual
compound_num <- as.integer(full_data[,3])

# Max number of compounds
compounds_count <- length(unique(compound_num))

# Map room numbers to 0,...,N for each compound
room_num_raw <- full_data[,4]
room_num <- integer(individuals_count_total)
room_num[which(is.na(room_num_raw))] <- as.integer(0)

for (v in 1:villages_count) {

  for (c in unique(compound_num[which(village_num==v)])) {

    rooms_raw <- sort(unique(room_num_raw[which(village_num == v & compound_num == c)]))

    if (length(rooms_raw) > 0) {

      room <- 1

      for (r in rooms_raw) {

        if (r != 99) {

          room_num[which(village_num == v & compound_num == c & room_num_raw == r)] <- as.integer(room)
          room <- room + 1

        } else{

          room_num[which(village_num == v & compound_num == c & room_num_raw == r)] <- as.integer(0)

        }

      }

    }

  }

}

# Max number of rooms
rooms_count <- length(unique(room_num))

# Individual ages
age_raw <- full_data[,5]
age <- integer(individuals_count_total)
age[which(!is.na(age_raw))] <- as.integer(age_raw[which(!is.na(age_raw))])

# Set unknown ages to median value
age[which(is.na(age_raw))] <- as.integer(14)

missing_ages <- which(is.na(age_raw))


################################################################################

# Observation info

# Number of observations
obs_count <- as.integer(13)

# Observations
obs_pcr <- array(integer(individuals_count_total * obs_count),
                 dim = c(individuals_count_total, obs_count))
obs_pcr[,] <- as.integer(-1)
obs_pcr[which(!is.na(full_data[, 8])), 1] <- as.integer(full_data[which(!is.na(full_data[, 8])), 8])

obs_att <- array(integer(individuals_count_total * obs_count),
                 dim = c(individuals_count_total, obs_count))
obs_att[,] <- as.integer(-1)
obs_att[which(!is.na(full_data[, 7])), 1] <- as.integer(full_data[which(!is.na(full_data[, 7])), 7])

obs_cex <- array(integer(individuals_count_total * obs_count),
                 dim = c(individuals_count_total, obs_count))
obs_cex[,] <- as.integer(-1)
obs_cex[which(!is.na(full_data[, 6])), 1]<-as.integer(full_data[which(!is.na(full_data[, 6])), 6])

for(y in 2:obs_count){
  
  obs_att[which(!is.na(full_data[, y + 8])), y] <- as.integer(full_data[which(!is.na(full_data[, y + 8])), y + 8])
  
  obs_cex[which(!is.na(full_data[, y + 21])), y] <- as.integer(full_data[which(!is.na(full_data[, y + 21])), y + 21])
  
}

obs_cex[which(!is.na(full_data[, 35])), 13] <- as.integer(full_data[which(!is.na(full_data[, 35])), 35])

# Remove data for individuals younger than 1
obs_pcr[which(age == 0),] <- as.integer(-1)
obs_att[which(age == 0),] <- as.integer(-1)
obs_cex[which(age == 0),] <- as.integer(-1)

################################################################################
# Model info

# Ladder structure of model
disease_groups <- as.integer(2)
disease_comps <- as.integer(4)
disease_states <- as.integer(disease_groups * disease_comps)

int_steps <- as.integer(2)

time_count <- as.integer((obs_count - 1) * int_steps + 1)

obs_days <- as.integer(which((0:(time_count - 1)) %% int_steps == 0) - 1)


################################################################################
# Prior distribution hyperparameters


# Rate of acquiring previous infections (Poisson process)
prior_history_shape <- 1.5
prior_history_rate <- 0.75

# Concentration hyperparameters for initial probability distributions conditioned on infection history
prior_initial <- rep(2 / disease_comps, disease_comps)

# Exponential hyperparameters for transmission parameters
prior_background_transmission <- 1
prior_village_transmission <- 1
prior_compound_transmission <- 1
prior_room_transmission <- 1

# Beta hyperparameters for Negative-Binomial state transition probabilities
prior_state_shape_1 <- 1
prior_state_shape_2 <- 1

# Exponential hyperparameters for Negative-Binomial state transition shapes
prior_state_rate <- 1 / 50

# Gamma hyperparameters for relative reinfection risk
prior_reinfection_shape <- 2
prior_reinfection_rate <- 2

prior_re_vil_rate <- 20

# Beta hyperparameters for observations specificities and sensitivities
prior_pcr_specificity_shape_1 <- 20
prior_pcr_specificity_shape_2 <- 1

prior_pcr_sensitivity_shape_1 <- 5
prior_pcr_sensitivity_shape_2 <- 2

prior_att_specificity_shape_1 <- 20
prior_att_specificity_shape_2 <- 1

prior_att_sensitivity_shape_1 <- 5
prior_att_sensitivity_shape_2 <- 2

prior_cex_specificity_shape_1 <- 25
prior_cex_specificity_shape_2 <- 1

prior_cex_sensitivity_shape_1 <- 25
prior_cex_sensitivity_shape_2 <- 1


################################################################################
# MCMC intialisation

# MCMC blocks and chain length
mcmc_iterations <- 1
mcmc_length <- 50000
mcmc_states_length <- 500
mcmc_states_thin <- mcmc_length / mcmc_states_length

# Storage for MCMC chains
berending_history_chain <- numeric(mcmc_length)
jali_history_chain <- numeric(mcmc_length)

berending_initial_chain <- array(numeric(mcmc_length * disease_groups * disease_comps),
                                 dim = c(mcmc_length, disease_comps, disease_groups))
jali_initial_chain <- array(numeric(mcmc_length * disease_groups * disease_comps),
                            dim = c(mcmc_length, disease_comps, disease_groups))

berending_background_chain <- numeric(mcmc_length)
berending_village_chain <- numeric(mcmc_length)
berending_compound_chain <- numeric(mcmc_length)
berending_room_chain <- numeric(mcmc_length)
jali_background_chain <- numeric(mcmc_length)
jali_village_chain <- numeric(mcmc_length)
jali_compound_chain <- numeric(mcmc_length)
jali_room_chain <- numeric(mcmc_length)

i_to_id_prob_chain <- numeric(mcmc_length)
i_to_id_shape_chain <- numeric(mcmc_length)

id_to_d_prob_chain <- array(numeric(mcmc_length * disease_groups), dim = c(mcmc_length, disease_groups))
id_to_d_shape_chain <- array(numeric(mcmc_length * disease_groups), dim = c(mcmc_length, disease_groups))

d_to_s_prob_chain <- array(numeric(mcmc_length * disease_groups), dim = c(mcmc_length, disease_groups))
d_to_s_shape_chain <- array(numeric(mcmc_length * disease_groups), dim = c(mcmc_length, disease_groups))

reinfec_hist_ber_chain <- numeric(mcmc_length)
reinfec_hist_jal_chain <- numeric(mcmc_length)

reinfection_chain <- numeric(mcmc_length)

pcr_spec_chain <- numeric(mcmc_length)
pcr_sen_chain <- numeric(mcmc_length)
att_spec_chain <- numeric(mcmc_length)
att_sen_chain <- numeric(mcmc_length)
cex_spec_chain <- numeric(mcmc_length)
cex_sen_chain <- numeric(mcmc_length)

states_chain <- array(integer(mcmc_states_length * individuals_count_total * time_count),
                      dim = c(mcmc_states_length, individuals_count_total, time_count))

missing_ages_chain <- array(integer(mcmc_length * length(missing_ages)), dim = c(mcmc_length, length(missing_ages)))

################################################################################
# Initial parameter values

berending_history <- rgamma(1, shape = prior_history_shape, rate = prior_history_rate)
jali_history <- rgamma(1, shape = prior_history_shape, rate = prior_history_rate)

berending_initial <- array(0, dim = c(disease_comps, disease_groups))
jali_initial <- array(0, dim = c(disease_comps, disease_groups))

ber_exp_1 <- rgamma(4, shape = prior_initial, rate = 1)
ber_exp_2 <- rgamma(4, shape = prior_initial, rate = 1)
jal_exp_1 <- rgamma(4, shape = prior_initial, rate = 1)
jal_exp_2 <- rgamma(4, shape = prior_initial, rate = 1)

berending_initial[, 1] <- ber_exp_1 / sum(ber_exp_1)
berending_initial[, 2] <- ber_exp_2 / sum(ber_exp_2)
jali_initial[, 1] <- jal_exp_1 / sum(jal_exp_1)
jali_initial[, 2] <- jal_exp_2 / sum(jal_exp_2)

berending_background <- 0
berending_village <- rexp(1, 500)
berending_compound <- rexp(1, 500)
berending_room <- rexp(1, 500)
jali_background <- 0
jali_village <- rexp(1, 500)
jali_compound <- rexp(1, 500)
jali_room <- rexp(1, 500)

i_to_id_prob <- runif(1)
i_to_id_shape <- 2

# Length biased distributions
i_to_id_mean <- 1 + (1 - i_to_id_prob) * i_to_id_shape / i_to_id_prob
lb_i_pdf <- pnbinom(1:time_count - 2, size = i_to_id_shape, prob = i_to_id_prob, lower.tail = F) / i_to_id_mean
lb_i_pdf[time_count] <- 1 - sum(lb_i_pdf[1:(time_count - 1)])
lb_i_cen <- rev(cumsum(rev(lb_i_pdf)))

# Geometric proposal
i_geom <- 1 / i_to_id_mean

id_to_d_prob <- sort(runif(2))
id_to_d_shape <- c(2, 2)

# Length biased distributions
id_to_d_mean <- 1 + (1 - id_to_d_prob) * id_to_d_shape / id_to_d_prob
lb_id_pdf <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
lb_id_cen <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
for (g in 1:disease_groups) {
  
  lb_id_pdf[, g] <- pnbinom(1:time_count - 2, size = id_to_d_shape[g], prob = id_to_d_prob[g], lower.tail = F) / id_to_d_mean[g]
  lb_id_pdf[time_count, g] <- 1 - sum(lb_id_pdf[1:(time_count - 1), g])
  lb_id_cen[, g] <- rev(cumsum(rev(lb_id_pdf[, g])))
  
}

# Geometric proposal
id_geom <- 1 / id_to_d_mean


d_to_s_prob <- sort(runif(2))
d_to_s_shape <- c(2, 2)

d_to_s_mean <- 1 + (1 - d_to_s_prob) * d_to_s_shape / d_to_s_prob

# Geometric proposal
d_geom <- 1 / d_to_s_mean


reinfec_hist_ber <- rexp(1, rate = 20)
reinfec_hist_jal <- rexp(1, rate = 20)

# Length biased distributions
lb_d_pdf_ber <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
lb_d_cen_ber <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
lb_d_pdf_jal <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))
lb_d_cen_jal <- array(numeric(time_count * disease_groups), dim = c(time_count, disease_groups))

# Time horizon for numerical approximation
hist_horiz <- 300

for (g in 1:disease_groups) {
  
  # Approximated historic holding times for Berending and Jali
  hist_pdf_ber <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g]) *
    pgeom(1:hist_horiz - 1, prob = reinfec_hist_ber, lower.tail = F) +
    pnbinom(1:hist_horiz - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F) *
    dgeom(1:hist_horiz - 1, prob = reinfec_hist_ber)
  
  hist_pdf_jal <- dnbinom(1:hist_horiz - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g]) *
    pgeom(1:hist_horiz - 1, prob = reinfec_hist_jal, lower.tail = F) +
    pnbinom(1:hist_horiz - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F) *
    dgeom(1:hist_horiz - 1, prob = reinfec_hist_jal)
  
  hist_cdf_rt_ber <- rev(cumsum(rev(hist_pdf_ber)))
  hist_mean_ber <- sum((1:hist_horiz) * hist_pdf_ber)
  hist_lb_pdf_ber <- hist_cdf_rt_ber / hist_mean_ber
  
  hist_cdf_rt_jal <- rev(cumsum(rev(hist_pdf_jal)))
  hist_mean_jal <- sum((1:hist_horiz) * hist_pdf_jal)
  hist_lb_pdf_jal <- hist_cdf_rt_jal / hist_mean_jal
  
  nb_pdf <- dnbinom(1:(time_count + hist_horiz - 1) - 1, size = d_to_s_shape[g], prob = d_to_s_prob[g])
  for (h in 1:hist_horiz) {
    
    lb_d_pdf_ber[, g] <- lb_d_pdf_ber[, g] + hist_lb_pdf_ber[h] * nb_pdf[1:time_count + h - 1] /
      pnbinom(h - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F)
    
    lb_d_pdf_jal[, g] <- lb_d_pdf_jal[, g] + hist_lb_pdf_jal[h] * nb_pdf[1:time_count + h - 1] /
      pnbinom(h - 2, size = d_to_s_shape[g], prob = d_to_s_prob[g], lower.tail = F)
    
  }
  
  lb_d_pdf_ber[time_count, g] <- 1 - sum(lb_d_pdf_ber[1:(time_count - 1), g])
  lb_d_cen_ber[, g] <- rev(cumsum(rev(lb_d_pdf_ber[, g])))
  
  lb_d_pdf_jal[time_count, g] <- 1 - sum(lb_d_pdf_jal[1:(time_count - 1), g])
  lb_d_cen_jal[, g] <- rev(cumsum(rev(lb_d_pdf_jal[, g])))
  
}

reinfection <- rgamma(1, shape = 2, rate = 2)

pcr_spec <- rbeta(1, shape1 = 20, shape2 = 1)
pcr_sen <- rbeta(1, shape1 = 5, shape2 = 2)
att_spec <- rbeta(1, shape1 = 20, shape2 = 1)
att_sen <- rbeta(1, shape1 = 5, shape2 = 2)
cex_spec <- rbeta(1, shape1 = 25, shape2 = 1)
cex_sen <- rbeta(1, shape1 = 25, shape2 = 1)

# Store initial values

berending_history_chain[1] <- berending_history
jali_history_chain[1] <- jali_history

berending_initial_chain[1,,] <- berending_initial
jali_initial_chain[1,,] <- jali_initial

berending_background_chain[1] <- berending_background
berending_village_chain[1] <- berending_village
berending_compound_chain[1] <- berending_compound
berending_room_chain[1] <- berending_room
jali_background_chain[1] <- jali_background
jali_village_chain[1] <- jali_village
jali_compound_chain[1] <- jali_compound
jali_room_chain[1] <- jali_room

i_to_id_prob_chain[1] <- i_to_id_prob
i_to_id_shape_chain[1] <- i_to_id_shape

id_to_d_prob_chain[1,] <- id_to_d_prob
id_to_d_shape_chain[1,] <- id_to_d_shape

d_to_s_prob_chain[1,] <- d_to_s_prob
d_to_s_shape_chain[1,] <- d_to_s_shape

reinfec_hist_ber_chain[1] <- reinfec_hist_ber
reinfec_hist_jal_chain[1] <- reinfec_hist_jal

reinfection_chain[1] <- reinfection

pcr_spec_chain[1] <- pcr_spec
pcr_sen_chain[1] <- pcr_sen
att_spec_chain[1] <- att_spec
att_sen_chain[1] <- att_sen
cex_spec_chain[1] <- cex_spec
cex_sen_chain[1] <- cex_sen

missing_ages_chain[1,] <- age[missing_ages]


################################################################################

# MH Proposals

mh_hist_ber_cov <- 0.0032
mh_hist_jal_cov <- 0.0015

mh_tr_ber_cov <- matrix(c(0, 0, 0, 0, 0, 2.75e-9, -3.62e-8, 3.8e-10, 0, -3.62e-8, 1.44e-6, -2.24e-7, 0, 3.8e-10, -2.24e-7, 1.38e-6), byrow = T, nrow = 4, ncol = 4)
mh_tr_jal_cov <- matrix(c(0, 0, 0, 0, 0, 4.62e-10, -7e-10, -1.79e-8, 0, -7e-10, 3.09e-8, -4.11e-8, 0, -1.79e-8, -4.11e-8, 2.85e-6), byrow = T, nrow = 4, ncol = 4)

mh_i_cov <- diag(c(0.18, 0))

mh_id_cov_1 <- diag(c(0.00077, 0))
mh_id_cov_2 <- diag(c(0.03, 0))

mh_d_cov_1 <- diag(c(0.0015, 0))
mh_d_cov_2 <- diag(c(0.077, 0))

mh_reinf_cov <- 2.11

mh_re_hist_ber_cov <- 0.016
mh_re_hist_jal_cov <- 0.02

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


################################################################################

  
inf_base <- array(dim = c(1430400, 3))

for (i in 0:893) {
  
  for (j in 0:79) {
    
    for (k in 0:19) {
      
      inf_idx <- 1600 * i + 20 * j + k + 1
      
      inf_base[inf_idx, ] <- c(i, j, k)

    }
    
  }
  
}
  

################################################################################
# Simulating states

ber_ss <- exp(-(berending_background + 
                   inf_base[, 1] * berending_village + 
                        inf_base[, 2] * berending_compound + 
                        inf_base[, 3] * berending_room))


ber_dd <- exp(-reinfection * (berending_background + 
                  inf_base[, 1] * berending_village + 
                  inf_base[, 2] * berending_compound + 
                  inf_base[, 3] * berending_room))


jal_ss <- exp(-(jali_background + 
                  inf_base[, 1] * jali_village + 
                  inf_base[, 2] * jali_compound + 
                  inf_base[, 3] * jali_room))


jal_dd <- exp(-reinfection * (jali_background + 
                                inf_base[, 1] * jali_village + 
                                inf_base[, 2] * jali_compound + 
                                inf_base[, 3] * jali_room))


sim <- model_sim(individuals_count_total,
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
                 lb_d_pdf_jal)


states <- array(sim$states, dim = c(individuals_count_total, time_count))

v_ip <- array(sim$village_pressures, dim = c(villages_count, time_count))
c_ip <- array(sim$compound_pressures, dim = c(villages_count, compounds_count, time_count))
r_ip <- array(sim$room_pressures, dim = c(villages_count, compounds_count, rooms_count, time_count))



################################################################################



acceptances <- integer(individuals_count_total)

init_probs <- numeric(disease_states)

pred_probs <- array(numeric(disease_states * time_count), dim = c(disease_states, time_count))
filt_probs <- array(numeric(disease_states * time_count), dim = c(disease_states, time_count))

prop_states <- integer(time_count)

urv_s <- array(runif(individuals_count_total * time_count), dim = c(individuals_count_total, time_count))
urv_m <- runif(individuals_count_total)

obs_parms <- c(pcr_sen = pcr_sen, pcr_spec = pcr_spec,
               att_sen = att_sen, att_spec = att_spec,
               cex_sen = cex_sen, cex_spec = cex_spec)

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

acc_iffbs_mean[1] <- mean(acceptances)
acc_iffbs_ind <- acc_iffbs_ind + acceptances

runtime <- system.time({

  for (i in 2:mcmc_length) {
  
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

save.image(file = "Chain1_Batch1.rdata")

