### Generate OTU table
gen_OTUtab = function(n = 100,
                       OTU_prop,
                       theta,
                       mu = 1000,
                       size = 25) 
{
  OTU_prop_over = OTU_prop * (1-theta) / theta
  n_OTUs = length(OTU_prop_over)
  OTU_tab = matrix(NA, nrow = n, ncol = length(OTU_prop_over))
  colnames(OTU_tab) = names(OTU_prop_over)
  n_seq = rnbinom(n, mu = mu, size = size)
  for (i in 1:n) {
    temp_p = rgamma(n_OTUs, OTU_prop_over)
    temp_p = temp_p / sum(temp_p)
    OTU_tab[i,] = rmultinom(1, n_seq[i], prob = temp_p)
  }
  return(OTU_tab)
}

### Generate binary outcome
gen_data = function(OTUtab, # OTU count table
                    info_OTUs_ls, # a list of informative OTUs; each list is a set of OTU names
                    beta_vec, # a vector of betas for the groups of OTUs in info_OTUs_ls; its length = length(info_OTU_ls)
                    model_vec)
{
  if (1 %in% model_vec) { # if there is abundance scenario
    abun_OTUtab = OTUtab / rowSums(OTUtab) # OTU table of abundances
  }
  if (2 %in% model_vec) { # if there is presence/absence scenario
    pres_OTUtab = ifelse(OTUtab > 0, 1, 0) # note that this require OTUtab is the count table; otherwise, even 0 may become 1
  }
  
  ### Generate eta
  eta_mat = matrix(NA, nrow = nrow(OTUtab), ncol = length(info_OTUs_ls))
  for (k in 1:length(info_OTUs_ls)) { # for each cluster of informative OTUs
    info_OTUs = info_OTUs_ls[[k]] # current set of info OTUs
    if (model_vec[k] == 1) { # if using abundance
      Z = abun_OTUtab[, info_OTUs, drop = F]
    } else if (model_vec[k] == 2) { # if using presence
      Z = pres_OTUtab[, info_OTUs, drop = F]
    }
    eta = rowSums(Z)
    eta_mat[, k] = (eta-mean(eta))/sd(eta)
  }
  
  ### Generate outcomes
  prob = 1/(1 + exp(- eta_mat %*% beta_vec))
  Y = rbinom(n, 1, prob) # 0 negative; 1 positive
  data_mat = cbind(OTUtab, Y)
  colnames(data_mat) = c(colnames(OTUtab), 'Y')
  return(list(data_mat = data_mat, 
              prob = prob))
}

### Bray-Curtis distance
UniFrac_BC = function(OTUtab = OTUtab, w = 2) { return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/w) }

