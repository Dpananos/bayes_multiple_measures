data{
  int y[4];
  int N;
  
  real p_mu;
  real p_kappa;
  
  real mu_se_1;
  real sd_se_1;
  
  real mu_se_2;
  real sd_se_2;
  
  real mu_sp_1;
  real sd_sp_1;
  
  real mu_sp_2;
  real sd_sp_2;
}
parameters{
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> se_1;
  real<lower=0, upper=1> se_2;
  
  real<lower=0, upper=1> sp_1;
  real<lower=0, upper=1> sp_2;
  
}
transformed parameters{
  simplex[4] p_sample;
  p_sample[1] = p*se_1*se_2 + (1-p)*(1-sp_1)*(1-sp_2);
  p_sample[2] = p*se_1*(1-se_2) + (1-p)*(1-sp_1)*sp_2;
  p_sample[3] = p*(1-se_1)*se_2 + (1-p)*sp_1*(1-sp_2);
  p_sample[4] = p*(1-se_1)*(1-se_2) + (1-p)*sp_1*sp_2;
}
model{
  
  y ~ multinomial(p_sample);
  
  p ~ beta_proportion(p_mu,p_kappa);
  se_1 ~ normal(mu_se_1, sd_se_1);
  se_2 ~ normal(mu_se_2, sd_se_2);
  
  sp_1 ~ normal(mu_sp_1, sd_sp_1);
  sp_2 ~ normal(mu_sp_2, sd_sp_2);
  
}
generated quantities{
  int y_ppc[4] = multinomial_rng(p_sample, N);
}
