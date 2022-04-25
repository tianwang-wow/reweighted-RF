##########################################################################################################
### Generate OTU clusters and info OTUs 
##########################################################################################################
### Get tree and OTU abundance from smoking data
library(MiSPU)
data(dd)
data(throat.tree)
tree = throat.tree
OTU_prop = dd$pi
tree$tip.label = paste0('OTU', tree$tip.label)
names(OTU_prop) = paste0('OTU', names(OTU_prop))
OTU_prop = OTU_prop[tree$tip.label]
# OTU_prop_over = OTU_prop * (1-dd$theta) / dd$theta

### Do clustering
n_clus = 20
clus_results = pam(as.dist(cophenetic(tree)), n_clus, diss = TRUE)$clustering
clus_prop = sort(tapply(OTU_prop, clus_results, sum), decreasing = T)

### Get informative OTUs
all_info_OTUs_ls = vector('list', length = 2)
names(all_info_OTUs_ls) = c('Phy_Related', 'Phy_Unrelated')
### Phylogentically-related OTUs
info_cluster_id = c(2, 6, 15, 18, 19)
for (j in 1:5) {
  info_OTUs = names(which(clus_results == names(clus_prop)[info_cluster_id[j]]))
  all_info_OTUs_ls$Phy_Related[[j]] = info_OTUs
}
### Phylogentically-Unrelated OTUs
un_OTU_id_vec = c(11, 31, 71, 131, 201)
ranked_OTU_prop = sort(OTU_prop, decreasing = T)
for (j in 1:5) {
  temp_id = un_OTU_id_vec[j]
  info_OTUs = names(ranked_OTU_prop)[temp_id]
  current_clus = clus_results[names(ranked_OTU_prop[temp_id])]
  while (length(info_OTUs) <= 9) {
    temp_id = temp_id + 1
    temp_OTU = names(ranked_OTU_prop[temp_id])
    temp_clus = clus_results[temp_OTU]
    if (!temp_clus %in% current_clus) {
      info_OTUs = c(info_OTUs, temp_OTU)
      current_clus = c(current_clus, temp_clus)
    }
  }
  all_info_OTUs_ls$Phy_Unrelated[[j]] = info_OTUs
}

### Check selected OTUs
if (F) {
  for (ll in 1:2) { # Related or Unrelated OTUs
    print('-----')
    for (j in 1:5) { # Abundance High -> Low
      temp_OTUs = all_info_OTUs_ls[[ll]][[j]] # selected info OTUs
      print(sum(OTU_prop[temp_OTUs])) # show abundance sum
      print(length(unique(clus_results[temp_OTUs]))) # show number of unique clusters
    }
  }
}

### Set simulation settings
# j1: set of Phy-Related OTUs
# j2: set of Phy-Unrelated OTUs
# beta1: effect size for j1
# beta2: effect size for j2
# model1: for j1 set: 1 (abundance), 2 (presence)
# model2: for j2 set: 1 (abundance), 2 (presence)
pm_mat = matrix(NA, nrow = 18, ncol = 6)
colnames(pm_mat) =  c('j1', 'j2', 'beta1', 'beta2', 'model1', 'model2')
pm_mat[, 'j1'] = pm_mat[, 'j2'] = c(1, 3, 5)
pm_mat[, 'beta1'] = c(rep(3, 6), rep(0, 6), rep(3, 6))
pm_mat[, 'beta2'] = c(rep(0, 6), rep(3, 6), rep(3, 6))
pm_mat[, 'model1'] = rep(c(rep(1, 3), rep(2, 3)), 3)
pm_mat[, 'model2'] = pm_mat[, 'model1']
print(pm_mat) # check settings, looks OK

save(list = c('tree', 'dd', 'OTU_prop', 'clus_results', 'all_info_OTUs_ls'), file = './R_files/simulation_settings.RData')

