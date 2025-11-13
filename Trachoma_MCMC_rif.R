# Conditioned initial states (Gibbs)
init_gibbs <- function(init_states, v, v_num, v_count, groups, comps, prior_conc) {
  
  init_prob <- array(numeric(groups * comps), dim = c(comps, groups))
  
  counts <- integer(comps)
  
  for (g in 1:groups) {
    
    g_idx <- which(v_num == v & as.integer(init_states / comps) == (g - 1))
    
    for (c in 1:comps) {
      
      counts[c] <- length(which(init_states[g_idx] %% comps == (c - 1)))
      
    }
    
    post_conc <- prior_conc + counts
    
    init_prob[, g] <- rdirichlet(1, alpha = post_conc)
    
  }
  
  return(init_prob)
  
}


obs_gibbs <- function(pcr_obs, att_obs, cex_obs, obs_states, groups, comps, prior_pcr_spec_s1, prior_pcr_spec_s2, prior_pcr_sen_s1, prior_pcr_sen_s2, prior_att_spec_s1, prior_att_spec_s2, prior_att_sen_s1, prior_att_sen_s2, prior_cex_spec_s1, prior_cex_spec_s2, prior_cex_sen_s1, prior_cex_sen_s2) {
  
  pcr_spec_trials <- length(which((obs_states %% comps == 0 | obs_states %% comps == 3) & (pcr_obs == 0 | pcr_obs == 1)))
  pcr_spec_succs <- length(which((obs_states %% comps == 0 | obs_states %% comps == 3) & (pcr_obs == 0)))
  
  pcr_sen_trials <- length(which((obs_states %% comps == 1 | obs_states %% comps == 2) & (pcr_obs == 0 | pcr_obs == 1)))
  pcr_sen_succs <- length(which((obs_states %% comps == 1 | obs_states %% comps == 2) & (pcr_obs == 1)))
  
  att_spec_trials <- length(which((obs_states %% comps == 0 | obs_states %% comps == 3) & (att_obs == 0 | att_obs == 1)))
  att_spec_succs <- length(which((obs_states %% comps == 0 | obs_states %% comps == 3) & (att_obs == 0)))
  
  att_sen_trials <- length(which((obs_states %% comps == 1 | obs_states %% comps == 2) & (att_obs == 0 | att_obs == 1)))
  att_sen_succs <- length(which((obs_states %% comps == 1 | obs_states %% comps == 2) & (att_obs == 1)))
  
  cex_spec_trials <- length(which((obs_states %% comps == 0 | obs_states %% comps == 1) & (cex_obs == 0 | cex_obs == 1)))
  cex_spec_succs <- length(which((obs_states %% comps == 0 | obs_states %% comps == 1) & (cex_obs == 0)))
  
  cex_sen_trials <- length(which((obs_states %% comps == 2 | obs_states %% comps == 3) & (cex_obs == 0 | cex_obs == 1)))
  cex_sen_succs <- length(which((obs_states %% comps == 2 | obs_states %% comps == 3) & (cex_obs == 1)))
  
  post_pcr_spec_s1 <- prior_pcr_spec_s1 + pcr_spec_succs
  post_pcr_spec_s2 <- prior_pcr_spec_s2 + (pcr_spec_trials - pcr_spec_succs)
  
  post_pcr_sen_s1 <- prior_pcr_sen_s1 + pcr_sen_succs
  post_pcr_sen_s2 <- prior_pcr_sen_s2 + (pcr_sen_trials - pcr_sen_succs)
  
  post_att_spec_s1 <- prior_att_spec_s1 + att_spec_succs
  post_att_spec_s2 <- prior_att_spec_s2 + (att_spec_trials - att_spec_succs)
  
  post_att_sen_s1 <- prior_att_sen_s1 + att_sen_succs
  post_att_sen_s2 <- prior_att_sen_s2 + (att_sen_trials - att_sen_succs)
  
  post_cex_spec_s1 <- prior_cex_spec_s1 + cex_spec_succs
  post_cex_spec_s2 <- prior_cex_spec_s2 + (cex_spec_trials - cex_spec_succs)
  
  post_cex_sen_s1 <- prior_cex_sen_s1 + cex_sen_succs
  post_cex_sen_s2 <- prior_cex_sen_s2 + (cex_sen_trials - cex_sen_succs)
 
#  return(list(pcr_sen = rbeta(1, shape1 = post_pcr_sen_s1, shape2 = post_pcr_sen_s2),
#              pcr_spec = rbeta(1, shape1 = post_pcr_spec_s1, shape2 = post_pcr_spec_s2),
#              att_sen = rbeta(1, shape1 = post_att_sen_s1, shape2 = post_att_sen_s2),
#              att_spec = rbeta(1, shape1 = post_att_spec_s1, shape2 = post_att_spec_s2),
#              cex_sen = rbeta(1, shape1 = post_cex_sen_s1, shape2 = post_cex_sen_s2),
#              cex_spec = rbeta(1, shape1 = post_cex_spec_s1, shape2 = post_cex_spec_s2)))
   
  return(c(pcr_sen = rbeta(1, shape1 = post_pcr_sen_s1, shape2 = post_pcr_sen_s2),
           pcr_spec = rbeta(1, shape1 = post_pcr_spec_s1, shape2 = post_pcr_spec_s2),
           att_sen = rbeta(1, shape1 = post_att_sen_s1, shape2 = post_att_sen_s2),
           att_spec = rbeta(1, shape1 = post_att_spec_s1, shape2 = post_att_spec_s2),
           cex_sen = rbeta(1, shape1 = post_cex_sen_s1, shape2 = post_cex_sen_s2),
           cex_spec =  rbeta(1, shape1 = post_cex_spec_s1, shape2 = post_cex_spec_s2)))
  
}
  






