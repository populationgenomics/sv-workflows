library(mashr)
library(tidyverse)

#load in the estrs
estr = read.csv('~/Downloads/str_associatr_final_freeze_tob_n950_and_bioheart_n975_strs_meta_results_meta_with_fixed_mashr_beta_se_all_chrom_all_celltypes_beta_se.tsv', sep = '\t')
estr <- na.omit(estr)
estr_beta = estr %>% select(-contains('_se'))
estr_se = estr %>% select(-contains('_beta'))

#load in all TRs tested on chr2 (random set of tests)
null_data = read.csv('~/Downloads/str_associatr_final_freeze_tob_n950_and_bioheart_n975_strs_meta_results_meta_with_fixed_mashr_chr2_null_beta_se_all_celltypes_beta_se.tsv', sep = '\t')
null_data_beta = null_data %>% select(-contains('_se'))
null_data_se = null_data %>% select(-contains('_beta'))

#subset 150,000 variants from chr2 
set.seed(123)
null_data <- null_data[sample(nrow(null_data), 150000), ]
null_data_beta = null_data %>% select(-contains('_se'))
null_data_se = null_data %>% select(-contains('_beta'))

#wrangle the data into a format amenable for mashR
df_subset<- subset(estr_beta, select = -locus)
beta_x = as.matrix(df_subset)
rownames(beta_x) <- estr_beta$locus

df_subset_se<- subset(estr_se, select = -locus)
se_x = as.matrix(df_subset_se)
rownames(se_x) <- estr_se$locus

df_subset<- subset(null_data_beta, select = -locus)
beta_null = as.matrix(df_subset)
rownames(beta_null) <- null_data_beta$locus

df_subset_se<- subset(null_data_se, select = -locus)
se_null = as.matrix(df_subset_se)
rownames(se_null) <- null_data_se$locus

#prepare inputs for mashR
strong.subset= mash_set_data(beta_x, se_x,zero_Bhat_Shat_reset =1e3)
random.subset = mash_set_data(beta_null, se_null,zero_Bhat_Shat_reset = 1e3)
Vhat = estimate_null_correlation_simple(random.subset)

data.strong= mash_set_data(beta_x, se_x, V= Vhat,zero_Bhat_Shat_reset = 1e3)
data.random = mash_set_data(beta_null, se_null, V =Vhat,zero_Bhat_Shat_reset = 1e3)

U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

#save outputs 
as.data.frame(m2$result$PosteriorMean)%>% write.csv('~/Desktop/etrs_fixed_1e3_chr2_posterior_mean_mash.csv')
as.data.frame(m2$result$PosteriorSD)%>% write.csv('~/Desktop/etrs_fixed_1e3_chr2_posterior_sd_mash.csv')
as.data.frame(m2$result$lfsr)%>% write.csv('~/Desktop/etrs_fixed_1e3_chr2_lfsr_mash.csv')
as.data.frame(m2$result$NegativeProb)%>% write.csv('~/Desktop/etrs_fixed_1e3_chr2_negativeprob_mash.csv')
