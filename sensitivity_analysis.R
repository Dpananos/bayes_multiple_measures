library(rstan)
library(tidyverse)
library(tidybayes)

param_grid <- crossing(y = list(c(164, 415, 268, 3310)), 
           N = 4157,
           p_mu = 0.5, 
           p_kappa = 2,
           mu_se_1 = seq(0.65, 0.75, 0.025),
           sd_se_1 = seq(0.01, 0.05, 0.01),
           mu_se_2 = 0.553, 
           sd_se_2 = 0.07,
           mu_sp_1 = 0.938,
           sd_sp_1 = 0.005,
           mu_sp_2 = 0.937,
           sd_sp_2 = 0.02) 
  

dl = transpose(param_grid)
  
model <- stan_model(file = 'stan/combined_tests.stan')

fit_model<-function(dl){
  fit <- sampling(model, data = dl, chains = 4)
  
  fit_summ<- fit %>% 
             spread_draws(p) %>% 
             mean_qi()
  
  rhats<- bayesplot::rhat(fit)
  fit_summ$rhat <- all(rhats<1.01)
  
  fit_summ
}

sa = map_dfr(dl, fit_model)


sa %>% 
  bind_cols(param_grid) %>% 
  ggplot(aes(sd_se_1, p, ymin = .lower, ymax = .upper))+
  geom_pointrange()+
  facet_wrap(~mu_se_1, nrow = 1)

