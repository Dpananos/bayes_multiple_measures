data{
  
  int y[4]; // Data containing the table counts.  Entries are (PP, NP, PN, NN)
  int N; // Total number of obserbations.  Equal to sum(y)
  
  real p_mu; // Prior mean for prevalence.  We place a beta prior on the prevalence.  
  real p_kappa; // Prior precision parameter
  
  real mu_se_1; //prior mean for admin senitivity
  real sd_se_1; //prior standard deviation admin sensitivity
  
  real mu_se_2; //prior mean for survey sensitivity
  real sd_se_2; //prior standard deviation for survey sensitivity
  
  real mu_sp_1; //admin specificity
  real sd_sp_1;
  
  real mu_sp_2; //survey specificity
  real sd_sp_2;
}
parameters{
  
  real<lower=0, upper=1> p; //population prevalence
  real<lower=0, upper=1> se_1; //admin sensitivity
  real<lower=0, upper=1> se_2; //survey sensitivity
  
  real<lower=0, upper=1> sp_1; //admin specifitiy
  real<lower=0, upper=1> sp_2; //survey specificity
  
}
transformed parameters{
  // parameter for proportion of sample falling in each cell of the 2x2 table
  // The two by two table is flattened so the cells are c(PP, NP, PN, NN).
  // Here, the first letter indicates the result of the admin test and the second the result of the survey.
  // The algebra to determine how the sensitivities and specificities are combined for the multinomial likelihood
  // is as follows.  Using the result NP as an example...
  //
  // P(NP) = P(NP|D+)P(D+) + P(NP|D-)P(D-)
  //
  // Note that P(D+) = pi is the prevalence.  Assuming the test results are independent
  //
  // = P(N|D+)P(P|D+)P(D+) + P(N|D-)P(P|D-)P(D-) = (1-se_1) x se_2 x pi + sp_1 x (1-sp_2) x (1-pi)
  //
  // Similar computations can be performed for the remaining three cells.
  //
  // The code is vectorized in order to do all the computations simultanously.
  // Thus, the mulitnomial probability for the second cell is the SE1[2]*SE2[2]*pi + SP1[2]*SP2[2]*(1-pi)
  vector[4] SE1 = to_vector({se_1, 1-se_1, se_1, 1-se_1});
  vector[4] SE2 = to_vector({se_2, se_2, 1-se_2, 1-se_2});
  vector[4] SP1 = to_vector({1-sp_1, sp_1, 1-sp_1, sp_1});
  vector[4] SP2 = to_vector({1-sp_2, 1-sp_2, sp_2, sp_2});
  
  simplex[4] p_sample = p * SE1 .* SE2 + (1-p) * SP1 .* SP2;
}
model{
  // Model likelihood
  y ~ multinomial(p_sample);
  
  // Model priors
  p ~ beta_proportion(p_mu,p_kappa);
  se_1 ~ normal(mu_se_1, sd_se_1);
  se_2 ~ normal(mu_se_2, sd_se_2);
  
  sp_1 ~ normal(mu_sp_1, sd_sp_1);
  sp_2 ~ normal(mu_sp_2, sd_sp_2);
  
}
generated quantities{
  // Posterior predictive distribution
  // Take the learned multinomial sample probabilities
  // Draw from the multinomial distribution
  // 2x2 tables which are probable under the model.
  int y_ppc[4] = multinomial_rng(p_sample, N);
  
  real log_lik = multinomial_lpmf(y|p_sample);
}
