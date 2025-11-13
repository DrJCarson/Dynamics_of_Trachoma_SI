#include <Rcpp.h>

// Function to convert 2D array indices into a 1D array index
int arr2(int i1, int i2, int d1, int d2) {
  
  int v = i2 * d1 + i1;
  
  return v;
  
}

// Function to convert 3D array indices into a 1D array index
int arr3(int i1, int i2, int i3, int d1, int d2, int d3) {
  
  int v = i3 * d2 * d1 + i2 * d1 + i1;
  
  return v;
  
}

// Function to convert 3D array indices into a 1D array index
int arr4(int i1, int i2, int i3, int i4, int d1, int d2, int d3, int d4) {
  
  int v = i4 * d3 * d2 * d1 + i3 * d2 * d1 + i2 * d1 + i1;
  
  return v;
  
}

// Sample state from vector of probabilities given a uniform random variable
int sample_state(Rcpp::NumericVector p, double u) {
  
  int s = 0;
  double c = p(0);

  while (c < u) {
    
    ++s;
    c += p(s);
    
  }
  
  return s;
    
}

// Sample state from an array of probabilities given a uniform random variable
// and column index
int sample_state_2(Rcpp::NumericVector p, int j, int d1, int d2, double u) {
  
  int s = 0;
  double c = p(arr2(0, j, d1, d2));
  
  while (c < u) {
    
    ++s;
    c += p(arr2(s, j, d1, d2));
    
  }
  
  return s;
  
}

// Return infectious status from state
int is_infec(int s, int comps) {
  
  if (s % comps == 1 || s % comps == 2) {
    
    return 1;
    
  } else {
    
    return 0;
    
  }
  
}

// Return susceptible status from state
int is_sus(int s, int comps) {
  
  if (s % comps == 0 || s % comps == 3) {
    
    return 1;
    
  } else {
    
    return 0;
    
  }
  
}

// Calculate initial state probabilities from parameters
void calc_init(int groups, int comps, int a, double hist_rate,
               Rcpp::NumericVector init_cond, Rcpp::NumericVector init_probs) {
  
  int g, c;
  double g_prob;

  for (g = 0; g < groups; ++g) {
    
    if (g == groups - 1) {
      
      g_prob = R::ppois(double(g - 1), (double(a) + 0.5) * hist_rate, false, false);
      
    } else {
      
      g_prob = R::dpois(double(g), (double(a) + 0.5) * hist_rate, false);
      
    }
    
    for (c = 0; c < comps; ++c) {
      
      init_probs(arr2(c, g, comps, groups)) = g_prob * init_cond(arr2(c, g, comps, groups));
      
    }
    
  }
  
}


int ss_idx(int vi, int ci, int ri) {
  
  return vi * 1600 + ci * 20 + ri; 
  
}


// Simulate from the forward model
// [[Rcpp::export]]
Rcpp::List model_sim(int ind_count,
                     Rcpp::IntegerVector v_num,
                     int v_count,
                     Rcpp::IntegerVector c_num,
                     int c_count,
                     Rcpp::IntegerVector r_num,
                     int r_count,
                     Rcpp::IntegerVector age,
                     int t_count,
                     int groups,
                     int comps,
                     double hist_ber,
                     double hist_jal,
                     Rcpp::NumericVector init_ber,
                     Rcpp::NumericVector init_jal,
                     Rcpp::NumericVector ber_ss,
                     Rcpp::NumericVector ber_dd,
                     Rcpp::NumericVector jal_ss,
                     Rcpp::NumericVector jal_dd,
                     double i_prob,
                     double i_shape,
                     Rcpp::NumericVector id_prob,
                     Rcpp::NumericVector id_shape,
                     Rcpp::NumericVector d_prob,
                     Rcpp::NumericVector d_shape,
                     Rcpp::NumericVector lb_i_pdf,
                     Rcpp::NumericVector lb_id_pdf,
                     Rcpp::NumericVector lb_d_pdf_ber,
                     Rcpp::NumericVector lb_d_pdf_jal) {

  // Individual states
  Rcpp::IntegerVector states(ind_count * t_count);
  
  // Infection pressures for each village, compound, and room
  Rcpp::IntegerVector v_ip(v_count * t_count), 
    c_ip(v_count * c_count * t_count), 
    r_ip(v_count * c_count * r_count * t_count);
  
  // Current holding times for each individual
  Rcpp::IntegerVector hold_time(ind_count);
  
  // Random numbers for simulation
  Rcpp::NumericVector surv(ind_count * t_count);
  surv = Rcpp::runif(ind_count * t_count);

  // Counters
  int i, a, v, c, r, t, g;

  double prob_ss;

  // Current state
  int s;

  // Vector for initial probabilities
  Rcpp::NumericVector init_probs(groups * comps);
  
  // Initialise
  t = 0;
  
  // Loop through individuals
  for (i = 0; i < ind_count; ++i) {
    
    // Read individual information
    a = age(i);
    v = v_num(i);
    c = c_num(i);
    r = r_num(i);
    
    // Calculate initial probabilities for each state
    if (v == 1) {
      
      calc_init(groups, comps, a, hist_ber, init_ber, init_probs);
      
    } else {
      
      calc_init(groups, comps, a, hist_jal, init_jal, init_probs);
      
    }
    
    s = sample_state(init_probs, surv(arr2(i, t, ind_count, t_count)));
    states(arr2(i, t, ind_count, t_count)) = s;
      
    // Update infection pressures
    if (is_infec(s, comps)) {
      
      ++v_ip(arr2(v - 1, t, v_count, t_count));
      ++c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
      ++r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
      
    }
    
    // Sample holding times based on current state
    if (s % comps == 1) {
      
      hold_time[i] = sample_state(lb_i_pdf, R::runif(0.0, 1.0));
      
    } else if (s % comps == 2) {
      
      g = int(s / comps);
      hold_time[i] = sample_state_2(lb_id_pdf, g, t_count, groups, R::runif(0.0, 1.0));
                                    
    } else if (s % comps == 3) {
      
      g = int(s / comps);
      
      if (v == 1) {
        
        hold_time[i] = sample_state_2(lb_d_pdf_ber, g, t_count, groups, R::runif(0.0, 1.0));
        
      } else {
        
        hold_time[i] = sample_state_2(lb_d_pdf_jal, g, t_count, groups, R::runif(0.0, 1.0));
        
      }
      
    }
    
  }
  
  for (t = 1; t < t_count; ++t) {
    
    for (i = 0; i < ind_count; ++i) {
      
      // Read individual information
      a = age(i);
      v = v_num(i);
      c = c_num(i);
      r = r_num(i);
      
      s = states(arr2(i, t - 1, ind_count, t_count));
      
      if (s % comps == 0) {
        
        if (v == 1) { 
          
          prob_ss = ber_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)),
                                   c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                   r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));

        } else {
          
          prob_ss = jal_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)),
                                  c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                  r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));

        }
        
        if (surv(arr2(i, t, ind_count, t_count)) < prob_ss) {
          
          states(arr2(i, t, ind_count, t_count)) = s;

        } else {
          
          states(arr2(i, t, ind_count, t_count)) = s + 1;
          
          hold_time(i) = R::rnbinom(i_shape, i_prob);

        }
        
      } else if (s % comps == 1) {
        
        if (hold_time(i) > 0) {
          
          states(arr2(i, t, ind_count, t_count)) = s;
          
          --hold_time(i);
          
        } else {
          
          states(arr2(i, t, ind_count, t_count)) = s + 1;
          
          g = int(states(arr2(i, t, ind_count, t_count)) / comps);
          
          hold_time(i) = R::rnbinom(id_shape(g), id_prob(g));
          
        }
        
      } else if (s % comps == 2) {
        
        if (hold_time(i) > 0) {
          
          states(arr2(i, t, ind_count, t_count)) = s;
          
          --hold_time(i);
          
        } else {
          
          states(arr2(i, t, ind_count, t_count)) = s + 1;
          
          g = int(states(arr2(i, t, ind_count, t_count)) / comps);
          
          hold_time(i) = R::rnbinom(d_shape(g), d_prob(g));
          
        }
        
      } else {
        
        if (v == 1) { 
          
          prob_ss = ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)),
                                  c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                  r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
          
        } else {
          
          prob_ss = jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)),
                                  c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                  r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
          
        }
        
        if (surv(arr2(i, t, ind_count, t_count)) < prob_ss) {
          
          if (hold_time(i) > 0) {
            
            states(arr2(i, t, ind_count, t_count)) = s;
            
            --hold_time(i);
            
          } else {
            
            g = int(s / comps);
            
            if (g == (groups - 1)) {
              
              states(arr2(i, t, ind_count, t_count)) = s - 3;
              
            } else {
              
              states(arr2(i, t, ind_count, t_count)) = s + 1;
              
            }
            
          }
          
        } else {
          
          g = int(s / comps);
          
          if (g == (groups - 1)) {
            
            states(arr2(i, t, ind_count, t_count)) = s - 1;
            
          } else {
            
            states(arr2(i, t, ind_count, t_count)) = s + 3;
            
          }
          
          g = int(states(arr2(i, t, ind_count, t_count)) / comps);
          
          hold_time(i) = R::rnbinom(id_shape(g), id_prob(g));
          
        }
        
      } 
      
      if (is_infec(states(arr2(i, t, ind_count, t_count)), comps)) {
        
        ++v_ip(arr2(v - 1, t, v_count, t_count));
        ++c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
        ++r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
        
      }
      
    }
    
  }

  // Return infection states and pressures
  Rcpp::List out(4);
  out[0] = states;
  out[1] = v_ip;
  out[2] = c_ip;
  out[3] = r_ip;
  
  out.names() = Rcpp::CharacterVector::create("states", "village_pressures", 
            "compound_pressures", "room_pressures");
  
  return(out);
  
}




// Function for calculating the likelihood of a single observation using a given specificity and sensitivity
double obs_lik(int obs, int act, double spec, double sen){
  
  // Initialise likelihood
  double likelihood = 1.0;
  
  // If observation is negative
  if(obs == 0){
    
    // False negative
    if(act == 1){
      
      likelihood = 1.0 - sen;
      
      // True negative
    } else{
      
      likelihood = spec;
      
    }
    
    // If observation is positive
  } else if(obs == 1){
    
    // True positive
    if(act == 1){
      
      likelihood = sen;
      
      // False positive
    } else{
      
      likelihood = 1.0 - spec;
      
    }
    
  }
  
  return likelihood;
  
}

// Function to classify if state should return a positive PCR test
int pcr_pos(int s, int comps){
  
  if(s % comps == 1 || s % comps == 2){ // If in state I or ID
    
    return 1; // Test should be positive
    
  } else{ // If in any other state
    
    return 0; // Test should be negative
    
  }
  
}

// Function to classify if state should return a positive Ag test
int att_pos(int s, int comps){
  
  if(s % comps == 1 || s % comps == 2){ // If in state I or ID
    
    return 1; // Test should be positive
    
  } else{ // If in any other state
    
    return 0; // Test should be negative
    
  }
  
}

// Function to classify if state should return a positive Dz test
int cex_pos(int s, int comps){
  
  if(s % comps == 2 || s % comps == 3){ // If in state ID or D
    
    return 1; // Test should be positive
    
  } else{ // If in any other state
    
    return 0; // Test should be negative
    
  }
  
}


// Update hidden states using IFFBS
// [[Rcpp::export]]
void iffbs(int ind_count,
           Rcpp::IntegerVector v_num,
           int v_count,
           Rcpp::IntegerVector c_num,
           int c_count,
           Rcpp::IntegerVector r_num,
           int r_count,
           Rcpp::IntegerVector age,
           int t_count,
           int groups,
           int comps,
           int obs_count,
           Rcpp::IntegerVector obs_days,
           Rcpp::NumericVector obs_pcr,
           Rcpp::NumericVector obs_att,
           Rcpp::NumericVector obs_cex,
           double hist_ber,
           double hist_jal,
           Rcpp::NumericVector init_ber,
           Rcpp::NumericVector init_jal,
           Rcpp::NumericVector ber_ss,
           Rcpp::NumericVector ber_dd,
           Rcpp::NumericVector jal_ss,
           Rcpp::NumericVector jal_dd,
           double i_prob,
           double i_shape,
           Rcpp::NumericVector id_prob,
           Rcpp::NumericVector id_shape,
           Rcpp::NumericVector d_prob,
           Rcpp::NumericVector d_shape,
           Rcpp::NumericVector lb_i_pdf,
           Rcpp::NumericVector lb_id_pdf,
           Rcpp::NumericVector lb_d_pdf_ber,
           Rcpp::NumericVector lb_d_cen_ber,
           Rcpp::NumericVector lb_d_pdf_jal,
           Rcpp::NumericVector lb_d_cen_jal,
           double i_geom,
           Rcpp::NumericVector id_geom,
           Rcpp::NumericVector d_geom,
           Rcpp::NumericVector obs_parms,
           Rcpp::IntegerVector states,
           Rcpp::IntegerVector v_ip,
           Rcpp::IntegerVector c_ip,
           Rcpp::IntegerVector r_ip,
           Rcpp::IntegerVector acceptances,
           Rcpp::NumericVector init_probs,
           Rcpp::NumericVector pred_probs,
           Rcpp::NumericVector filt_probs,
           Rcpp::IntegerVector prop_states,
           Rcpp::NumericVector urv_s,
           Rcpp::NumericVector urv_m) {
  
  double pcr_sen = obs_parms["pcr_sen"], pcr_spec = obs_parms["pcr_spec"];
  double att_sen = obs_parms["att_sen"], att_spec = obs_parms["att_spec"];
  double cex_sen = obs_parms["cex_sen"], cex_spec = obs_parms["cex_spec"];
  
  // Probability modification
  double pmodi;
  
  // counters
  int i, i2, t, y, a, v, c, r, v2, c2, r2, g, s1, s2, z1, z2;
  
  // Number of possible states
  int s_count = groups * comps;
 
  // Normalising constant
  double norm;
  
  // Cumulative probability
  double c_prob;
  
  // MH acceptance probability
  double acc_prob;
  
  // Holding time and length bias status
  int curr_ht, prop_ht, curr_lb, prop_lb;
  
  // Transition probabilities
  double curr_nb, curr_geo, prop_nb, prop_geo;
 
  // Loop through individuals
  for (i = 0; i < ind_count; ++i) {
    
    // Read individual information
    a = age(i);
    v = v_num(i);
    c = c_num(i);
    r = r_num(i);
    
    t = 0;
    y = 0;
    
    s1 = states(arr2(i, t, ind_count, t_count));
    s2 = states(arr2(i, t + 1, ind_count, t_count));
    
    if (is_infec(s1, comps)) {
      
      --v_ip(arr2(v - 1, t, v_count, t_count));
      --c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
      --r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
      
    }
    
    // Calculate initial probabilities for each state
    if (v == 1) {
      
      calc_init(groups, comps, a, hist_ber, init_ber, init_probs);
      
    } else {
      
      calc_init(groups, comps, a, hist_jal, init_jal, init_probs);
      
    }
    
    for (s1 = 0; s1 < s_count; ++s1) {
      
      pred_probs(arr2(s1, t, s_count, t_count)) = init_probs(s1);
      filt_probs(arr2(s1, t, s_count, t_count)) = init_probs(s1);
      
    }
    
    pmodi = 1.0;
    
    for (i2 = 0; i2 < ind_count; ++i2) {
      
      v2 = v_num(i2);
      c2 = c_num(i2);
      r2 = r_num(i2);
      
      if (v != v2 || i == i2) {
        
        continue;
        
      }
      
      z1 = states(arr2(i2, t, ind_count, t_count));
      
      if (z1 % comps == 0) {
        
        z2 = states(arr2(i2, t + 1, ind_count, t_count));
        
        if (v == 1) {
          
          if (c != c2) {
            
            if (z1 == z2) {
            
              pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
            
            } else {
              
              pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          } else if (r != r2 || r == 0) {
            
            if (z1 == z2) {
            
              pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));

            }
            
          } else {
            
            if (z1 == z2) {
              
              pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          }
          
        } else {
          
          if (c != c2) {
            
            if (z1 == z2) {
            
              pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
            
          } else if (r != r2 || r == 0) {
            
            if (z1 == z2) {
              
              pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          } else {
            
            if (z1 == z2) {
              
              pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          }
          
        }
        
      } else if (z1 % comps == 3) {
        
        z2 = states(arr2(i2, t + 1, ind_count, t_count));
        
        if (v == 1) {
          
          if (c != c2) {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          } else if (r != r2 || r == 0) {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
            
            
          } else {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          }
          
        } else {
          
          if (c != c2) {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          } else if (r != r2 || r == 0) {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }
            
          } else {
            
            if (z1 == z2 || z2 % comps == 0) {
              
              pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            } else {
              
              pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                     (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
              
              pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                     c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                     r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
              
            }            
            
          }
          
        }
        
      }
      
    }
    
    // Update filtered probabilities
    for (s1 = 0; s1 < s_count; ++s1) {
      
      if (is_infec(s1, comps)) {
        
        filt_probs(arr2(s1, t, s_count, t_count)) *= pmodi;
        
      }
      
    }
    
    // Assimilate observations
    if (t == obs_days[y]) {
      
      for (s1 = 0; s1 < s_count; ++s1) {
        
        filt_probs(arr2(s1, t, s_count, t_count)) *= 
          obs_lik(obs_pcr(arr2(i, y, ind_count, obs_count)), pcr_pos(s1, comps), pcr_spec, pcr_sen);
        
        filt_probs(arr2(s1, t, s_count, t_count)) *= 
          obs_lik(obs_att(arr2(i, y, ind_count, obs_count)), att_pos(s1, comps), att_spec, att_sen);
        
        filt_probs(arr2(s1, t, s_count, t_count)) *= 
          obs_lik(obs_cex(arr2(i, y, ind_count, obs_count)), cex_pos(s1, comps), cex_spec, cex_sen);
        
      }
      
      ++y;
      
    }
    
    // Normalise filtered probabilities
    norm = 0.0;
    
    for (s1 = 0; s1 < s_count; ++s1) {
      
      norm += filt_probs(arr2(s1, t, s_count, t_count));
      
    }
    
    for (s1 = 0; s1 < s_count; ++s1) {
      
      filt_probs(arr2(s1, t, s_count, t_count)) /= norm;
      
    }
    
    for (t = 1; t < t_count; ++t) {
      
      s1 = states(arr2(i, t, ind_count, t_count));
      
      if (is_infec(s1, comps)) {
        
        --v_ip(arr2(v - 1, t, v_count, t_count));
        --c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
        --r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
        
      }
      
      for (s2 = 0; s2 < s_count; ++s2) {
        
        pred_probs(arr2(s2, t, s_count, t_count)) = 0.0;
        
      }
      
      for (s1 = 0; s1 < s_count; ++s1) {
        
        if (s1 % comps == 0) {
          
          if (v == 1) {
            
            pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
              ber_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                            c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
              (1.0 - ber_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                            c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
          } else {
            
            pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
              jal_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                            c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
              (1.0 - jal_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                   c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                   r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
          }
          
        } else if (s1 % comps == 1) {
          
          pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * (1.0 - i_geom);
          
          pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * i_geom;
          
        } else if (s1 % comps == 2) {
          
          g = int(s1 / comps);
          
          pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * (1.0 - id_geom(g));
          
          pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * id_geom(g);
        
        } else {
          
          g = int(s1 / comps);
          
          if (v == 1) {
            
            if (g == (groups - 1)) {
              
              pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                (1.0 - d_geom(g));
              
              pred_probs(arr2(s1 - 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                (1.0 - ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));

              pred_probs(arr2(s1 - 3, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            } else {
              
              pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                (1.0 - d_geom(g));
              
              pred_probs(arr2(s1 + 3, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                (1.0 - ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
              
              pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            }
            
          } else {
            
            if (g == (groups - 1)) {
              
              pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                (1.0 - d_geom(g));
              
              pred_probs(arr2(s1 - 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                (1.0 - jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
              
              pred_probs(arr2(s1 - 3, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            } else {
              
              pred_probs(arr2(s1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                (1.0 - d_geom(g));
              
              pred_probs(arr2(s1 + 3, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                (1.0 - jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
              
              pred_probs(arr2(s1 + 1, t, s_count, t_count)) += filt_probs(arr2(s1, t - 1, s_count, t_count)) * 
                jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            }
            
          }

        } 

      }
      
      for (s2 = 0; s2 < s_count; ++s2) {
        
        filt_probs(arr2(s2, t, s_count, t_count)) = 
          pred_probs(arr2(s2, t, s_count, t_count));
        
      }

      
      // If not on final day
      if ((t + 1) < t_count) {
        
        pmodi = 1.0;
        
        for (i2 = 0; i2 < ind_count; ++i2) {
          
          v2 = v_num(i2);
          c2 = c_num(i2);
          r2 = r_num(i2);
          
          if (v != v2 || i == i2) {
            
            continue;
            
          }
          
          z1 = states(arr2(i2, t, ind_count, t_count));
          
          if (z1 % comps == 0) {
            
            z2 = states(arr2(i2, t + 1, ind_count, t_count));
            
            if (v == 1) {
              
              if (c != c2) {
                
                if (z1 == z2) {
                  
                  pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else if (r != r2 || r == 0) {
                
                if (z1 == z2) {
                  
                  pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else {
                
                if (z1 == z2) {
                  
                  pmodi *= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              }
              
            } else {
              
              if (c != c2) {
                
                if (z1 == z2) {
                  
                  pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
                
              } else if (r != r2 || r == 0) {
                
                if (z1 == z2) {
                  
                  pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else {
                
                if (z1 == z2) {
                  
                  pmodi *= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_ss(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              }
              
            }
            
          } else if (z1 % comps == 3) {
            
            z2 = states(arr2(i2, t + 1, ind_count, t_count));
            
            if (v == 1) {
              
              if (c != c2) {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else if (r != r2 || r == 0) {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
                
                
              } else {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - ber_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              }
              
            } else {
              
              if (c != c2) {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else if (r != r2 || r == 0) {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }
                
              } else {
                
                if (z1 == z2 || z2 % comps == 0) {
                  
                  pmodi *= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                         (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                         c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                         r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                } else {
                  
                  pmodi *= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)) + 1,
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)) + 1,
                                               (r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) + 1) * int(r2 != 0)));
                  
                  pmodi /= 1.0 - jal_dd(ss_idx(v_ip(arr2(v2 - 1, t, v_count, t_count)),
                                               c_ip(arr3(v2 - 1, c2 - 1, t, v_count, c_count, t_count)),
                                               r_ip(arr4(v2 - 1, c2 - 1, r2, t, v_count, c_count, r_count, t_count)) * int(r2 != 0)));
                  
                }            
                
              }
              
            }
            
          }
          
        }
        
        // Update filtered probabilities
        for (s1 = 0; s1 < s_count; ++s1) {
          
          if (is_infec(s1, comps)) {
            
            filt_probs(arr2(s1, t, s_count, t_count)) *= pmodi;
            
          }
          
        }
        
      }
      
      // Assimilate observations
      if (t == obs_days[y]) {
        
        for (s1 = 0; s1 < s_count; ++s1) {
          
          filt_probs(arr2(s1, t, s_count, t_count)) *= 
            obs_lik(obs_pcr(arr2(i, y, ind_count, obs_count)), pcr_pos(s1, comps), pcr_spec, pcr_sen);
          
          filt_probs(arr2(s1, t, s_count, t_count)) *= 
            obs_lik(obs_att(arr2(i, y, ind_count, obs_count)), att_pos(s1, comps), att_spec, att_sen);
          
          filt_probs(arr2(s1, t, s_count, t_count)) *= 
            obs_lik(obs_cex(arr2(i, y, ind_count, obs_count)), cex_pos(s1, comps), cex_spec, cex_sen);
          
        }
        
        ++y;
        
      }

      // Normalise filtered probabilities
      norm = 0.0;
      
      for (s1 = 0; s1 < s_count; ++s1) {
        
        norm += filt_probs(arr2(s1, t, s_count, t_count));
        
      }
      
      for (s1 = 0; s1 < s_count; ++s1) {
        
        filt_probs(arr2(s1, t, s_count, t_count)) /= norm;
        
      }
      
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Initialise backward sampling on final day
    t = t_count - 1;
    
    // Sample state on final day from filtered probabilities
    s1 = 0;
    c_prob = filt_probs(arr2(s1, t, s_count, t_count));
    
    while (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) {
      
      ++s1;
      c_prob += filt_probs(arr2(s1, t, s_count, t_count));
      
    }
    
    prop_states(t) = s1;
    
    // Loop backwards through time
    for (t = (t_count - 2); t >= 0; --t) {
      
      s2 = prop_states(t + 1);
      
      if (s2 % comps == 0) {
        
        g = int(s2 / comps);
        
        if (g == 0) {
          
          s1 = s2;
          
        } else if (g == (groups - 1)) {
          
          s1 = s2 - 1;
          
          if (v == 1) {
            
            c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
              pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
              ber_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                            c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                              d_geom(g - 1);
            
          } else {
            
            c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
              pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
              jal_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                            c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                              d_geom(g - 1);
            
          }
          
          if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
            
            s1 = s2 + 3;
            
            if (v == 1) {
              
              c_prob += (filt_probs(arr2(s1, t, s_count, t_count)) / 
                pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
                ber_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                              c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            } else {
              
              c_prob += (filt_probs(arr2(s1, t, s_count, t_count)) / 
                pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
                jal_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                              c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                              r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                                d_geom(g);
              
            }
            
            if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
              
              s1 = s2;
              
            }
            
          }
          
        } else {
          
          s1 = s2 - 1;
          
          if (v == 1) {
            
            c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
              pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
              ber_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                            c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                              d_geom(g - 1);
            
          } else {
            
            c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
              pred_probs(arr2(s2, t + 1, s_count, t_count))) * 
              jal_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                            c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                            r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))) * 
                              d_geom(g - 1);
            
          }
          
          if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
            
            s1 = s2;
            
          }
            
        }
        
      } else if (s2 % comps == 1) {
        
        s1 = s2;
        
        c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
          pred_probs(arr2(s2, t + 1, s_count, t_count))) * (1.0 - i_geom);
        
        if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
          
          s1 = s2 - 1;
          
        }
        
      } else if (s2 % comps == 2) {
        
        g = int(s2 / comps);
        
        s1 = s2;
        
        c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
          pred_probs(arr2(s2, t + 1, s_count, t_count))) * (1.0 - id_geom(g));
        
        if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
          
          s1 = s2 - 1;
          
          c_prob += (filt_probs(arr2(s1, t, s_count, t_count)) / 
            pred_probs(arr2(s2, t + 1, s_count, t_count))) * (i_geom);
          
          if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
           
           if (g < (groups - 1)) {
             
             s1 = s2 - 3;
             
           } else {
             
             s1 = s2 + 1;
             
             if (v == 1) {
               
               c_prob += (filt_probs(arr2(s1, t, s_count, t_count)) / 
                 pred_probs(arr2(s2, t + 1, s_count, t_count))) * (1.0 - ber_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                                 c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                                 r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))));
               
             } else {
               
               c_prob += (filt_probs(arr2(s1, t, s_count, t_count)) / 
                 pred_probs(arr2(s2, t + 1, s_count, t_count))) * (1.0 - jal_dd(ss_idx(v_ip(arr2(v - 1, t, v_count, t_count)),
                                 c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count)),
                                 r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count)) * int(r != 0))));
               
             }
             
             if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
               
               s1 = s2 - 3;
               
             }
           
           }
           
          }
          
        }
        
      } else {
        
        g = int(s2 / comps);
        
        s1 = s2 - 1;
        
        c_prob = (filt_probs(arr2(s1, t, s_count, t_count)) / 
          pred_probs(arr2(s2, t + 1, s_count, t_count))) * id_geom(g);
        
        if (urv_s(arr2(i, t, ind_count, t_count)) > c_prob) { 
          
          s1 = s2;
          
        }
        
      } 
      
      prop_states(t) = s1;
      
    }
    
    curr_ht = 0;
    prop_ht = 0;
    
    curr_lb = 1;
    prop_lb = 1;
    
    curr_nb = 1.0;
    curr_geo = 1.0;
    prop_nb = 1.0;
    prop_geo = 1.0;
    
    for (t = 1; t < t_count; ++t) {
      
      // Current transition probabilities
      s1 = states(arr2(i, t - 1, ind_count, t_count));
      s2 = states(arr2(i, t, ind_count, t_count));
      
      if (s1 == s2) {
        
        ++curr_ht;
        
      } else {

        if (s1 % comps == 1) {
          
          if (curr_lb == 1) {
            
            curr_nb *= lb_i_pdf(curr_ht);
            
          } else {
            
            curr_nb *= R::dnbinom(curr_ht, i_shape, i_prob, false);
            
          }
          
        } else if (s1 % comps == 2) {
          
          g = int(s1 / comps);
          
          if (curr_lb == 1) {
            
            curr_nb *= lb_id_pdf(arr2(curr_ht, g, t_count, groups));
            
          } else {
            
            curr_nb *= R::dnbinom(curr_ht, id_shape(g), id_prob(g), false);
            
          }
          
        } else if (s1 % comps == 3) {
          
          g = int(s1 / comps);
          
          if (s2 % comps == 2) {
            
            if (curr_lb == 1) {
              
              if (v == 1) {
                
                curr_nb *= lb_d_cen_ber(arr2(curr_ht, g, t_count, groups));
                
              } else {
                
                curr_nb *= lb_d_cen_jal(arr2(curr_ht, g, t_count, groups));
                
              }
              
            } else {
              
              curr_nb *= 
                R::pnbinom(curr_ht - 1, d_shape(g), d_prob(g), false, false);
              
            }
            
          } else {
            
            if (curr_lb == 1) {
              
              if (v == 1) {
                
                curr_nb *= lb_d_pdf_ber(arr2(curr_ht, g, t_count, groups));
                
              } else {
                
                curr_nb *= lb_d_pdf_jal(arr2(curr_ht, g, t_count, groups));
                
              }
              
            } else {
              
              curr_nb *= R::dnbinom(curr_ht, d_shape(g), d_prob(g), false);
              
            }
            
          }
          
        } 
        
        curr_ht = 0;
        curr_lb = 0;
        
      }
      
      if (s1 % comps == 1) {
        
        // Remain infected
        if (s2 == s1) {
          
          curr_geo *= 1.0 - i_geom;
          
          // Transition to infected and diseased
        } else if (s2 == (s1 + 1)) {
          
          curr_geo *= i_geom;
          
        }
        
        // Transition from infectious and diseased state
      } else if (s1 % comps == 2) {
        
        // Determine disease group
        g = int(s1 / comps);
        
        // Remain infectious and diseased
        if (s2 == s1) {
          
          curr_geo *= 1.0 - id_geom(g);
          
          // Transition to diseased
        } else if (s2 == (s1 + 1)) {
          
          curr_geo *= id_geom(g);
          
        }
        
        // If diseased
      } else if (s1 % comps == 3) { 

        // Determine disease group
        g = int(s1 / comps);
        
        if (s2 % comps == 3) {
          
          curr_geo *= 1.0 - d_geom(g);
          
        } else if (s2 % comps == 0) {
          
          curr_geo *= d_geom(g);
          
        }
        
      } 
      
      // Proposed transition probabilities
      s1 = prop_states(t - 1);
      s2 = prop_states(t);
      
      if (s1 == s2) {
        
        ++prop_ht;
        
      } else {
        
        if (s1 % comps == 1) {
          
          if (prop_lb == 1) {
            
            prop_nb *= lb_i_pdf(prop_ht);
            
          } else {
            
            prop_nb *= R::dnbinom(prop_ht, i_shape, i_prob, false);
            
          }
          
        } else if (s1 % comps == 2) {
          
          g = int(s1 / comps);
          
          if (prop_lb == 1) {
            
            prop_nb *= lb_id_pdf(arr2(prop_ht, g, t_count, groups));
            
          } else {
            
            prop_nb *= R::dnbinom(prop_ht, id_shape(g), id_prob(g), false);
            
          }
          
        } else if (s1 % comps == 3) {
          
          g = int(s1 / comps);
          
          if (s2 % comps == 2) {
            
            if (prop_lb == 1) {
              
              if (v == 1) {
                
                prop_nb *= lb_d_cen_ber(arr2(prop_ht, g, t_count, groups));
                
              } else {
                
                prop_nb *= lb_d_cen_jal(arr2(prop_ht, g, t_count, groups));
                
              }
              
            } else {
              
              prop_nb *= 
                R::pnbinom(prop_ht - 1, d_shape(g), d_prob(g), false, false);
              
            }
            
          } else {
            
            if (prop_lb == 1) {
              
              if (v == 1) {
                
                prop_nb *= lb_d_pdf_ber(arr2(prop_ht, g, t_count, groups));
                
              } else {
                
                prop_nb *= lb_d_pdf_jal(arr2(prop_ht, g, t_count, groups));
                
              }
              
            } else {
              
              prop_nb *= R::dnbinom(prop_ht, d_shape(g), d_prob(g), false);
              
            }
            
          }
          
        } 
        
        prop_ht = 0;
        prop_lb = 0;
        
      }
      
      if (s1 % comps == 1) {
        
        // Remain infected
        if (s2 == s1) {
          
          prop_geo *= 1.0 - i_geom;
          
          // Transition to infected and diseased
        } else if (s2 == (s1 + 1)) {
          
          prop_geo *= i_geom;
          
        }
        
        // Transition from infectious and diseased state
      } else if (s1 % comps == 2) {
        
        // Determine disease group
        g = int(s1 / comps);
        
        // Remain infectious and diseased
        if (s2 == s1) {
          
          prop_geo *= 1.0 - id_geom(g);
          
          // Transition to diseased
        } else if (s2 == (s1 + 1)) {
          
          prop_geo *= id_geom(g);
          
        }
        
        // If diseased
      } else if (s1 % comps == 3) { 
        
        // Determine disease group
        g = int(s1 / comps);
        
        if (s2 % comps == 3) {
          
          prop_geo *= 1.0 - d_geom(g);
          
        } else if (s2 % comps == 0) {
          
          prop_geo *= d_geom(g);
          
        }
        
      } 
      
    }
    
    
    // Current transition density from final time point (censored)
    t = t_count - 1;
    s1 = states(arr2(i, t, ind_count, t_count));
    
    if (s1 % comps == 1) {
      
      if (curr_lb == 1) {
        
        curr_nb *= lb_i_pdf(curr_ht);
        
      } else {
        
        curr_nb *= R::pnbinom(curr_ht - 1, i_shape, i_prob, false, false);
        
      }
      
    } else if (s1 % comps == 2) {
      
      g = int(s1 / comps);
      
      if (curr_lb == 1) {
        
        curr_nb *= lb_id_pdf(arr2(curr_ht, g, t_count, groups));
        
      } else {
        
        curr_nb *= R::pnbinom(curr_ht - 1, id_shape(g), id_prob(g), false, false);
        
      }
      
    } else if (s1 % comps == 3) {
      
      g = int(s1 / comps);
  
      if (curr_lb == 1) {
        
        if (v == 1) {
          
          curr_nb *= lb_d_cen_ber(arr2(curr_ht, g, t_count, groups));
          
        } else {
          
          curr_nb *= lb_d_cen_jal(arr2(curr_ht, g, t_count, groups));
          
        }
        
      } else {
        
        curr_nb *= 
          R::pnbinom(curr_ht - 1, d_shape(g), d_prob(g), false, false);
        
      }
      
    } 
    
    // Proposed
    s1 = prop_states(t);
    
    if (s1 % comps == 1) {
      
      if (prop_lb == 1) {
        
        prop_nb *= lb_i_pdf(prop_ht);
        
      } else {
        
        prop_nb *= R::pnbinom(prop_ht - 1, i_shape, i_prob, false, false);
        
      }
      
    } else if (s1 % comps == 2) {
      
      g = int(s1 / comps);
      
      if (prop_lb == 1) {
        
        prop_nb *= lb_id_pdf(arr2(prop_ht, g, t_count, groups));
        
      } else {
        
        prop_nb *= R::pnbinom(prop_ht - 1, id_shape(g), id_prob(g), false, false);
        
      }
      
    } else if (s1 % comps == 3) {
      
      g = int(s1 / comps);
        
      if (prop_lb == 1) {
        
        if (v == 1) {
          
          prop_nb *= lb_d_cen_ber(arr2(prop_ht, g, t_count, groups));
          
        } else {
          
          prop_nb *= lb_d_cen_jal(arr2(prop_ht, g, t_count, groups));
          
        }
        
      } else {
        
        prop_nb *= 
          R::pnbinom(prop_ht - 1, d_shape(g), d_prob(g), false, false);
        
      }
      
    } 
    
    acc_prob = (prop_nb * curr_geo) / (curr_nb * prop_geo);
    
    // MH update
    if (urv_m(i) <= acc_prob) {
      
      // Accept
      acceptances(i) = 1;
      
      // Loop backwards through time
      for(t = (t_count - 1); t >= 0; --t){
        
        // Write proposed states
        states(arr2(i, t, ind_count, t_count)) = prop_states(t);
        
        // Update summaries
        s1 = states(arr2(i, t, ind_count, t_count));

        if (is_infec(s1, comps)) {
          
          ++v_ip(arr2(v - 1, t, v_count, t_count));
          ++c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
          ++r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
          
        }
 
      }
      
      
    } else {
      
      // Reject
      acceptances(i) = 0;
      
      // Loop backwards through time
      for(t = (t_count - 1); t >= 0; --t){
      
        // Update summaries
        s1 = states(arr2(i, t, ind_count, t_count));
        
        if (is_infec(s1, comps)) {
          
          ++v_ip(arr2(v - 1, t, v_count, t_count));
          ++c_ip(arr3(v - 1, c - 1, t, v_count, c_count, t_count));
          ++r_ip(arr4(v - 1, c - 1, r, t, v_count, c_count, r_count, t_count));
          
        }
        
      }
        
    }

  }
  
}



// MH acceptance probabilities
// [[Rcpp::export]]
void MH_1(int ind_count,
           Rcpp::IntegerVector v_num,
           int v_count,
           Rcpp::IntegerVector c_num,
           int c_count,
           Rcpp::IntegerVector r_num,
           int r_count,
           Rcpp::IntegerVector age,
           int t_count,
           int groups,
           int comps,
           double hist_ber,
           double hist_ber_2,
           double hist_jal,
           double hist_jal_2,
           Rcpp::NumericVector ber_ss,
           Rcpp::NumericVector ber_ss_2,
           Rcpp::NumericVector ber_dd,
           Rcpp::NumericVector ber_dd_2,
           Rcpp::NumericVector jal_ss,
           Rcpp::NumericVector jal_ss_2,
           Rcpp::NumericVector jal_dd,
           Rcpp::NumericVector jal_dd_2,
           double i_prob,
           double i_prob_2,
           double i_shape,
           double i_shape_2,
           Rcpp::NumericVector id_prob,
           Rcpp::NumericVector id_prob_2,
           Rcpp::NumericVector id_shape,
           Rcpp::NumericVector id_shape_2,
           Rcpp::NumericVector d_prob,
           Rcpp::NumericVector d_prob_2,
           Rcpp::NumericVector d_shape,
           Rcpp::NumericVector d_shape_2,
           Rcpp::NumericVector lb_i_pdf,
           Rcpp::NumericVector lb_i_pdf_2,
           Rcpp::NumericVector lb_id_pdf,
           Rcpp::NumericVector lb_id_pdf_2,
           Rcpp::NumericVector lb_d_pdf_ber,
           Rcpp::NumericVector lb_d_pdf_ber_2,
           Rcpp::NumericVector lb_d_cen_ber,
           Rcpp::NumericVector lb_d_cen_ber_2,
           Rcpp::NumericVector lb_d_pdf_jal,
           Rcpp::NumericVector lb_d_pdf_jal_2,
           Rcpp::NumericVector lb_d_cen_jal,
           Rcpp::NumericVector lb_d_cen_jal_2,
           double reinf,
           Rcpp::IntegerVector states,
           Rcpp::IntegerVector v_ip,
           Rcpp::IntegerVector c_ip,
           Rcpp::IntegerVector r_ip,
           Rcpp::NumericVector prob_hist,
           Rcpp::NumericVector prob_tr,
           Rcpp::NumericVector prob_i,
           Rcpp::NumericVector prob_id,
           Rcpp::NumericVector prob_d){

  int lb, ht;
  
  int t, i, v, c, a, r, g, s1, s2;
  
  prob_hist(0) = 1.0;
  prob_hist(1) = 1.0;
  prob_tr(0) = 1.0;
  prob_tr(1) = 1.0;
  prob_i(0) = 1.0;
  
  for (g = 0; g < groups; ++g) {
    
    prob_id(g) = 1.0;
    prob_d(g) = 1.0;
    
  }
  
  for (i = 0; i < ind_count; ++i) {
    
    v = v_num(i);
    c = c_num(i);
    r = r_num(i);
    a = age(i);
    
    t = 0;
    
    s1 = states(arr2(i, t, ind_count, t_count));
    
    g = int(s1 / comps);
    
    if (v == 1) {
      
      if (g == groups - 1) {
        
        prob_hist(0) *= (R::ppois(double(g - 1), (double(a) + 0.5) * hist_ber_2, false, false) / 
          R::ppois(double(g - 1), (double(a) + 0.5) * hist_ber, false, false));
        
      } else {

        prob_hist(0) *= (R::dpois(double(g), (double(a) + 0.5) * hist_ber_2, false) / 
          R::dpois(double(g), (double(a) + 0.5) * hist_ber, false));
        
      }
      
    } else {
      
      if (g == groups - 1) {

        prob_hist(1) *= (R::ppois(double(g - 1), (double(a) + 0.5) * hist_jal_2, false, false) / 
          R::ppois(double(g - 1), (double(a) + 0.5) * hist_jal, false, false));
        
      } else {

        prob_hist(1) *= (R::dpois(double(g), (double(a) + 0.5) * hist_jal_2, false) / 
          R::dpois(double(g), (double(a) + 0.5) * hist_jal, false));
        
      }
      
    }
    
    ht = 0;
    lb = 1;
    
    for (t = 1; t < t_count; ++t) {
      
      s1 = states(arr2(i, t - 1, ind_count, t_count));
      s2 = states(arr2(i, t, ind_count, t_count));
      
      if (s1 % comps == 0) {
        
        if (v == 1) {
          
          if (s2 % comps == 0) { 
            
            prob_tr(0) *= ber_ss_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                              ber_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            ++ht;
            
          } else if (s2 % comps == 1) {
            
            prob_tr(0) *= (1.0 - ber_ss_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                              (1.0 - ber_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            ht = 0;
            lb = 0;
            
          }
          
        } else {
          
          if (s2 % comps == 0) { 
            
            prob_tr(1) *= jal_ss_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                                         jal_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            ++ht;
            
          } else if (s2 % comps == 1) {
            
            prob_tr(1) *= (1.0 - jal_ss_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                                (1.0 - jal_ss(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            ht = 0;
            lb = 0;
            
          }
          
        }
        
      } else if (s1 % comps == 1) {
        
        if (s2 % comps == 1) {
          
          ++ht;
          
        }  else {
          
          if (lb == 1) {

            prob_i(0) *= (lb_i_pdf_2(ht) / lb_i_pdf(ht));
            
            lb = 0;
            
          } else {

            prob_i(0) *= (R::dnbinom(ht, i_shape_2, i_prob_2, false) / 
              R::dnbinom(ht, i_shape, i_prob, false));
            
          }
          
          ht = 0;
          
        }
        
      } else if (s1 % comps == 2) {
        
        if (s2 % comps == 2) {
          
          ++ht;
          
        }  else {
          
          g = int(s1 / comps);
          
          if (lb == 1) {
            
            prob_id(g) *= (lb_id_pdf_2(arr2(ht, g, t_count, groups)) / lb_id_pdf(arr2(ht, g, t_count, groups)));
            
            lb = 0;
            
          } else {
            
            prob_id(g) *= (R::dnbinom(ht, id_shape_2(g), id_prob_2(g), false) / 
              R::dnbinom(ht, id_shape(g), id_prob(g), false));
            
          }
          
          ht = 0;
          
        }
        
      } else {
        
        g = int(s1 / comps);
        
        if (v == 1) {
          
          if (s2 % comps == 2) {
            
            prob_tr(0) *= (1.0 - ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                                (1.0 - ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {
              
              prob_d(g) *= (lb_d_cen_ber_2(arr2(ht, g, t_count, groups)) / lb_d_cen_ber(arr2(ht, g, t_count, groups)));
              
            } else {
              
              prob_d(g) *= (R::pnbinom(ht - 1, d_shape_2(g), d_prob_2(g), false, false) / 
                R::pnbinom(ht - 1, d_shape(g), d_prob(g), false, false));
              
            }
            
            ht = 0;
            lb = 0;
            
          } else if (s2 % comps == 0) {
            
            prob_tr(0) *= ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                                         ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            if (lb == 1) {

              prob_d(g) *= (lb_d_pdf_ber_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_ber(arr2(ht, g, t_count, groups)));
              
            } else {

              prob_d(g) *= (R::dnbinom(ht, d_shape_2(g), d_prob_2(g), false) / 
                R::dnbinom(ht, d_shape(g), d_prob(g), false));
              
            }
            
            ht = 0;
            lb = 0;
            
          } else {
            
            prob_tr(0) *= ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                                         ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            ++ht;
            
          }
          
        } else {
          
          if (s2 % comps == 2) {
            
            prob_tr(1) *= (1.0 - jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                              c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                              r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                                (1.0 - jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                     c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                     r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {
              
              prob_d(g) *= (lb_d_cen_jal_2(arr2(ht, g, t_count, groups)) / lb_d_cen_jal(arr2(ht, g, t_count, groups)));
              
            } else {
              
              prob_d(g) *= (R::pnbinom(ht - 1, d_shape_2(g), d_prob_2(g), false, false) / 
                R::pnbinom(ht - 1, d_shape(g), d_prob(g), false, false));
              
            }
            
            ht = 0;
            lb = 0;
            
          } else if (s2 % comps == 0) {
            
            prob_tr(1) *= jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                                         jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            if (lb == 1) {
              
              prob_d(g) *= (lb_d_pdf_jal_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_jal(arr2(ht, g, t_count, groups)));
              
            } else {
              
              prob_d(g) *= (R::dnbinom(ht, d_shape_2(g), d_prob_2(g), false) / 
                R::dnbinom(ht, d_shape(g), d_prob(g), false));
              
            }
            
            ht = 0;
            lb = 0;
            
          } else {
            
            prob_tr(1) *= jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))) / 
                                         jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)));
            
            ++ht;
            
          }
          
        }
        
      } 
      
    }
    
    t = t_count - 1;
    s1 = states(arr2(i, t, ind_count, t_count));
    
    if (s1 % comps == 1) {
      
      if (lb == 1) {

        prob_i(0) *= (lb_i_pdf_2(ht) / lb_i_pdf(ht));
        
      } else {

        prob_i(0) *= (R::pnbinom(ht - 1, i_shape_2, i_prob_2, false, false) / 
          R::pnbinom(ht - 1, i_shape, i_prob, false, false));
        
      }
      
    } else if (s1 % comps == 2) {
      
      g = int(s1 / comps);
      
      if (lb == 1) {
        
        prob_id(g) *= (lb_id_pdf_2(arr2(ht, g, t_count, groups)) / lb_id_pdf(arr2(ht, g, t_count, groups)));
        
      } else {

        prob_id(g) *= (R::pnbinom(ht - 1, id_shape_2(g), id_prob_2(g), false, false) / 
          R::pnbinom(ht - 1, id_shape(g), id_prob(g), false, false));
        
      }
      
    } else if (s1 % comps == 3) {
      
      g = int(s1 / comps);
      
      if (v == 1) {
        
        if (lb == 1) {
          
          prob_d(g) *= (lb_d_pdf_ber_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_ber(arr2(ht, g, t_count, groups)));
          
        } else {
          
          prob_d(g) *= (R::pnbinom(ht - 1, d_shape_2(g), d_prob_2(g), false, false) / 
            R::pnbinom(ht - 1, d_shape(g), d_prob(g), false, false));
          
        }
        
      } else {

        if (lb == 1) {
          
          prob_d(g) *= (lb_d_pdf_jal_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_jal(arr2(ht, g, t_count, groups)));
          
        } else {
          
          prob_d(g) *= (R::pnbinom(ht - 1, d_shape_2(g), d_prob_2(g), false, false) / 
            R::pnbinom(ht - 1, d_shape(g), d_prob(g), false, false));
          
        }
        
      }
      
    } 
      
  }
  
}


// MH acceptance probabilities
// [[Rcpp::export]]
void MH_2(int ind_count,
          Rcpp::IntegerVector v_num,
          int v_count,
          Rcpp::IntegerVector c_num,
          int c_count,
          Rcpp::IntegerVector r_num,
          int r_count,
          int t_count,
          int groups,
          int comps,
          double hist_ber,
          double hist_jal,
          Rcpp::NumericVector ber_dd,
          Rcpp::NumericVector ber_dd_2,
          Rcpp::NumericVector jal_dd,
          Rcpp::NumericVector jal_dd_2,
          double i_prob,
          double i_shape,
          Rcpp::NumericVector id_prob,
          Rcpp::NumericVector id_shape,
          Rcpp::NumericVector d_prob,
          Rcpp::NumericVector d_shape,
          Rcpp::NumericVector lb_d_pdf_ber,
          Rcpp::NumericVector lb_d_pdf_ber_2,
          Rcpp::NumericVector lb_d_cen_ber,
          Rcpp::NumericVector lb_d_cen_ber_2,
          Rcpp::NumericVector lb_d_pdf_jal,
          Rcpp::NumericVector lb_d_pdf_jal_2,
          Rcpp::NumericVector lb_d_cen_jal,
          Rcpp::NumericVector lb_d_cen_jal_2,
          double reinf,
          double reinf_2,
          Rcpp::IntegerVector states,
          Rcpp::IntegerVector v_ip,
          Rcpp::IntegerVector c_ip,
          Rcpp::IntegerVector r_ip,
          Rcpp::NumericVector prob_reinf,
          Rcpp::NumericVector prob_re_hist_ber,
          Rcpp::NumericVector prob_re_hist_jal) {
  
  int lb, ht;
  
  int t, i, v, c, r, g, s1, s2;

  prob_reinf(0) = 1.0;
  prob_re_hist_ber(0) = 1.0;
  prob_re_hist_jal(0) = 1.0;
  
  for (i = 0; i < ind_count; ++i) {
    
    v = v_num(i);
    c = c_num(i);
    r = r_num(i);
    
    ht = 0;
    lb = 1;
    
    for (t = 1; t < t_count; ++t) {
      
      s1 = states(arr2(i, t - 1, ind_count, t_count));
      s2 = states(arr2(i, t, ind_count, t_count));
      
      if (s1 % comps == 0) {
        
        if (s2 % comps == 0) {
          
          ++ht;
          
        } else {
          
          ht = 0;
          lb = 0;
          
        }
        
      } else if (s1 % comps == 1) {
        
        if (s2 % comps == 1) {
          
          ++ht;
          
        } else {
          
          ht = 0;
          lb = 0;
          
        }
        
      } else if (s1 % comps == 2) {
        
        if (s2 % comps == 2) {
          
          ++ht;
          
        } else {
          
          ht = 0;
          lb = 0;
          
        }
        
      } else {

        g = int(s1 / comps);
        
        if (v == 1) {
          
          if (s2 % comps == 2) {
            
            prob_reinf(0) *= (1.0 - ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                          c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                          r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                            (1.0 - ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                 c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                 r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {
              
              prob_re_hist_ber(0) *= (lb_d_cen_ber_2(arr2(ht, g, t_count, groups)) / lb_d_cen_ber(arr2(ht, g, t_count, groups)));
              
            } 
            
            ht = 0;
            lb = 0;
            
          } else if (s2 % comps == 0) {
            
            prob_reinf(0) *= (ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                             c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                             r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                               (ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                    c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                    r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {

              prob_re_hist_ber(0) *= (lb_d_pdf_ber_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_ber(arr2(ht, g, t_count, groups)));
              
            } 
            
            ht = 0;
            lb = 0;
            
          } else {
            
            prob_reinf(0) *= (ber_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                         (ber_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                        c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                        r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            ++ht;
            
          }
          
        } else {
          
          if (s2 % comps == 2) {
            
            prob_reinf(0) *= (1.0 - jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                             c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                             r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                               (1.0 - jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                                    c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                                    r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {
              
              prob_re_hist_jal(0) *= (lb_d_cen_jal_2(arr2(ht, g, t_count, groups)) / lb_d_cen_jal(arr2(ht, g, t_count, groups)));
              
            } 
            
            ht = 0;
            lb = 0;
            
          } else if (s2 % comps == 0) {
            
            prob_reinf(0) *= (jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                         (jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                        c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                        r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            if (lb == 1) {
              
              prob_re_hist_jal(0) *= (lb_d_pdf_jal_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_jal(arr2(ht, g, t_count, groups)));
              
            } 
            
            ht = 0;
            lb = 0;
            
          } else {
            
            prob_reinf(0) *= (jal_dd_2(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                       c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                       r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0)))) / 
                                         (jal_dd(ss_idx(v_ip(arr2(v - 1, t - 1, v_count, t_count)), 
                                                        c_ip(arr3(v - 1, c - 1, t - 1, v_count, c_count, t_count)),
                                                        r_ip(arr4(v - 1, c - 1, r, t - 1, v_count, c_count, r_count, t_count)) * int(r != 0))));
            
            ++ht;
            
          }
          
        }
          
      } 
      
    }
      
    t = t_count - 1;
    s1 = states(arr2(i, t, ind_count, t_count));
    
    if (s1 % comps == 3) {
      
      g = int(s1 / comps);
      
      if (v == 1) {
        
        if (lb == 1) {

          prob_re_hist_ber(0) *= (lb_d_pdf_ber_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_ber(arr2(ht, g, t_count, groups)));
          
        } 
        
      } else {
        
        if (lb == 1) {

          prob_re_hist_jal(0) *= (lb_d_pdf_jal_2(arr2(ht, g, t_count, groups)) / lb_d_pdf_jal(arr2(ht, g, t_count, groups)));
          
        } 
        
      }
        
    }
    
  }
  
}


