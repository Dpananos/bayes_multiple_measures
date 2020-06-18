data{
  int y[4];
  int N;
  
  real p_mu;
  real p_kappa;
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
  se_1 ~ normal(0.692, 0.02 );
  se_2 ~ normal(0.553, 0.07 );
  
  sp_1 ~ normal(0.938, 0.005);
  sp_2 ~ normal(0.937, 0.02);
  
}
generated quantities{
  int y_ppc[4] = multinomial_rng(p_sample, N);
}
